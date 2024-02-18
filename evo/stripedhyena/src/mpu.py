import torch
import deepspeed
from deepspeed import comm as dist
from deepspeed.utils import logger
from deepspeed.runtime.utils import bwc_tensor_model_parallel_rank
import contextlib
from contextlib import contextmanager

from .utils import divide
from copy import copy

# Model parallel group that the current rank belongs to.
_MODEL_PARALLEL_GROUP = None
# Data parallel group that the current rank belongs to.
_DATA_PARALLEL_GROUP = None
# Pipeline parallel group that the current rank belongs to.
_PIPE_PARALLEL_GROUP = None

# A group used to sync during the IO process. Usually this is data_parallel_group(),
# but with pipeline parallelism it must also involve the last stage (which is not in the
# DP group of rank 0)
_IO_PARALLEL_GROUP = None

# These values enable us to change the mpu sizes on the fly.
_MPU_WORLD_SIZE = None
_MPU_RANK = None

# Used to query 3D topology
_MPU_TOPOLOGY = None

# Get fp32_allreduce flag
_FP32_ALLREDUCE = None


# TODO: extract these
get_cuda_rng_tracker = deepspeed.checkpointing.get_cuda_rng_tracker
_set_cuda_rng_state = deepspeed.checkpointing._set_cuda_rng_state
model_parallel_cuda_manual_seed = deepspeed.checkpointing.model_parallel_cuda_manual_seed
rng_str = "rng_state"


def print_rank_0(message, debug=False, end="\n"):
    """Print from rank 0 only."""
    if torch.distributed.get_rank() == 0:
        print(message, flush=True, end=end)


def _get_extra_te_kwargs(config):
    extra_transformer_engine_kwargs = {}
    from importlib.metadata import version

    from pkg_resources import packaging

    te_version = packaging.version.Version(version("transformer-engine"))
    if te_version >= packaging.version.Version("0.12.0"):
        if config.use_cpu_initialization:
            extra_transformer_engine_kwargs["device"] = "cpu"
        else:
            extra_transformer_engine_kwargs["device"] = torch.cuda.current_device()
    return extra_transformer_engine_kwargs


def get_model_parallel_group(check_initialized=True):
    """Get the model parallel group the caller rank belongs to."""
    if check_initialized:
        assert _MODEL_PARALLEL_GROUP is not None, "model parallel group is not initialized"
    return _MODEL_PARALLEL_GROUP


def get_model_parallel_rank():
    """Return my rank for the model parallel group."""
    global _MPU_RANK
    if _MPU_RANK is not None:
        return _MPU_RANK
    return torch.distributed.get_rank(group=get_model_parallel_group())


def ensure_divisibility(numerator, denominator):
    """Ensure that numerator is divisible by the denominator."""
    assert numerator % denominator == 0, "{} is not divisible by {}".format(numerator, denominator)


def get_data_parallel_group():
    """Get the data parallel group the caller rank belongs to."""
    assert _DATA_PARALLEL_GROUP is not None, "data parallel group is not initialized"
    return _DATA_PARALLEL_GROUP


def initialize_model_parallel(model_parallel_size, topology=None, fp32_allreduce=False):
    """
    Initialize model data parallel groups.
    """
    if torch.distributed.get_rank() == 0:
        print("> initializing model parallel with size {}".format(model_parallel_size))
    # Get world size and rank. Ensure some consistencies.
    assert torch.distributed.is_initialized()
    world_size = torch.distributed.get_world_size()
    if world_size < model_parallel_size:
        raise ValueError("world size cannot be smaller than model parallel size")
    ensure_divisibility(world_size, model_parallel_size)
    rank = torch.distributed.get_rank()

    global _MPU_TOPOLOGY
    if topology:
        _MPU_TOPOLOGY = topology

    # Build the data parallel groups.
    global _DATA_PARALLEL_GROUP
    assert _DATA_PARALLEL_GROUP is None, "data parallel group is already initialized"
    if topology:
        for dp_group in topology.get_axis_comm_lists("data"):
            group = torch.distributed.new_group(ranks=dp_group)
            if rank == 0:
                print(f"MPU DP:", dp_group)
            if rank in dp_group:
                _DATA_PARALLEL_GROUP = group
    else:
        for i in range(model_parallel_size):
            ranks = range(i, world_size, model_parallel_size)
            group = torch.distributed.new_group(ranks)
            if i == (rank % model_parallel_size):
                _DATA_PARALLEL_GROUP = group

    # Build pipeline parallel group
    if topology is not None:
        global _PIPE_PARALLEL_GROUP
        for pp_group in topology.get_axis_comm_lists("pipe"):
            group = torch.distributed.new_group(ranks=pp_group)
            if rank == 0:
                print(f"MPU PP:", pp_group)
            if rank in pp_group:
                _PIPE_PARALLEL_GROUP = group

    # Build IO group
    global _IO_PARALLEL_GROUP
    if topology and topology.get_dim("pipe") > 1:
        io_stages = [0, topology.get_dim("pipe") - 1]
        io_group = []
        for stage in io_stages:
            io_group.extend(topology.filter_match(pipe=stage, model=0))
        if rank == 0:
            print(f"MPU IO:", io_group)
        group = torch.distributed.new_group(ranks=io_group)
        if rank in io_group:
            _IO_PARALLEL_GROUP = group
    else:
        _IO_PARALLEL_GROUP = get_data_parallel_group()

    # Build the model parallel groups.
    global _MODEL_PARALLEL_GROUP
    assert _MODEL_PARALLEL_GROUP is None, "model parallel group is already initialized"
    if topology:
        # Short circuit case without model parallelism.
        # TODO: it would be nice  to avoid this branching case?
        if model_parallel_size == 1:
            for group_rank in range(world_size):
                group = torch.distributed.new_group(ranks=[group_rank])
                if rank == 0:
                    print(f"MPU MP:", [group_rank])
                if rank == group_rank:
                    _MODEL_PARALLEL_GROUP = group
            return

        for mp_group in topology.get_axis_comm_lists("model"):
            group = torch.distributed.new_group(ranks=mp_group)
            if rank == 0:
                print(f"MPU MP:", mp_group)
            if rank in mp_group:
                _MODEL_PARALLEL_GROUP = group

    else:
        for i in range(world_size // model_parallel_size):
            ranks = range(i * model_parallel_size, (i + 1) * model_parallel_size)
            group = torch.distributed.new_group(ranks)
            if i == (rank // model_parallel_size):
                _MODEL_PARALLEL_GROUP = group

    global _FP32_ALLREDUCE
    assert _FP32_ALLREDUCE is None, "fp32_allreduce is already initialized"
    _FP32_ALLREDUCE = fp32_allreduce


_MODEL_PARALLEL_GROUP = None
_FP32_ALLREDUCE = None
_MODEL_PARALLEL_RNG_TRACKER_NAME = deepspeed.checkpointing._MODEL_PARALLEL_RNG_TRACKER_NAME


class CudaRNGStatesTracker:
    """Tracker for the cuda RNG states.

    Using the `add` method, a cuda rng state is initialized based on
    the input `seed` and is assigned to `name`. Later, by forking the
    rng state, we can perform operations and return to our starting
    cuda state.
    """

    def __init__(self):
        # Map from a string name to the cuda rng state.
        self.states_ = {}
        # Seeds are just for book keeping and ensure no seed is set twice.
        self.seeds_ = set()

    def reset(self):
        """Set to the initial state (no tracker)."""
        self.states_ = {}
        self.seeds_ = set()

    def get_states(self):
        """Get rng states. Copy the dictionary so we have direct
        pointers to the states, not just a pointer to the dictionary."""
        states = {}
        for name in self.states_:
            states[name] = self.states_[name]
        return states

    def set_states(self, states):
        """Set the rng states. For efficiency purposes, we do not check
        the size of seed for compatibility."""
        self.states_ = states

    def add(self, name, seed):
        """Track the rng state."""
        # Check seed is not already used.
        if seed in self.seeds_:
            raise Exception("seed {} already exists".format(seed))
        self.seeds_.add(seed)
        # Check that state is not already defined.
        if name in self.states_:
            raise Exception("cuda rng state {} already exists".format(name))
        # Get the current rng state.
        orig_rng_state = torch.cuda.get_rng_state()
        # Set the new state and store it.
        torch.cuda.manual_seed(seed)
        self.states_[name] = torch.cuda.get_rng_state()
        # Reset rng state to what it was.
        _set_cuda_rng_state(orig_rng_state)

    @contextlib.contextmanager
    def fork(self, name=rng_str):
        """Fork the cuda rng state, perform operations, and exit with
        the original state."""
        # Check if we have added the state
        if name not in self.states_:
            raise Exception("cuda rng state {} is not added".format(name))
        # Store current rng state.
        orig_cuda_rng_state = torch.cuda.get_rng_state()
        # Set rng state to the desired one
        _set_cuda_rng_state(self.states_[name])
        # Do the stuff we wanted to do.
        try:
            yield
        finally:
            # Update the current rng state for later use.
            self.states_[name] = torch.cuda.get_rng_state()
            # And set the state to the original state we started with.
            _set_cuda_rng_state(orig_cuda_rng_state)


_CUDA_RNG_STATE_TRACKER = CudaRNGStatesTracker()


def get_model_parallel_group(check_initialized=True):
    """Get the model parallel group the caller rank belongs to."""
    if check_initialized:
        assert _MODEL_PARALLEL_GROUP is not None, "model parallel group is not initialized"
    return _MODEL_PARALLEL_GROUP


def get_model_parallel_rank():
    """Return my rank for the model parallel group."""
    global _MPU_RANK
    if _MPU_RANK is not None:
        return _MPU_RANK
    return torch.distributed.get_rank(group=get_model_parallel_group())


def get_model_parallel_world_size():
    """Return world size for the model parallel group."""
    global _MPU_WORLD_SIZE
    if _MPU_WORLD_SIZE is not None:
        return _MPU_WORLD_SIZE
    return torch.distributed.get_world_size(group=get_model_parallel_group())


def get_fp32_allreduce():
    """Get the fp32 allreduce flag"""
    assert _FP32_ALLREDUCE is not None, "fp32_allreduce is not Initialized"
    return _FP32_ALLREDUCE


def _reduce(input_):
    """All-reduce the the input tensor across model parallel group."""

    # Bypass the function if we are using only 1 GPU.
    if get_model_parallel_world_size() == 1:
        return input_

    # Bf16 convert
    dt = input_.dtype
    if dt == torch.bfloat16 and get_fp32_allreduce():
        input_ = input_.float()

    # All-reduce.
    torch.distributed.all_reduce(input_, group=get_model_parallel_group())

    # Bf16 convert
    if dt == torch.bfloat16 and get_fp32_allreduce():
        input_ = input_.bfloat16()

    return input_


def split_tensor_along_last_dim(tensor, num_partitions, contiguous_split_chunks=False):
    """Split a tensor along its last dimension.
    Arguments:
        tensor: input tensor.
        num_partitions: number of partitions to split the tensor
        contiguous_split_chunks: If True, make each chunk contiguous
                                 in memory.
    """
    # Get the size and dimension.
    last_dim = tensor.dim() - 1
    last_dim_size = divide(tensor.size()[last_dim], num_partitions)
    # Split.
    tensor_list = torch.split(tensor, last_dim_size, dim=last_dim)
    # Note: torch.split does not create contiguous tensors by default.
    if contiguous_split_chunks:
        return tuple(chunk.contiguous() for chunk in tensor_list)

    return tensor_list


def _split(input_):
    """Split the tensor along its last dimension and keep the
    corresponding slice."""

    world_size = get_model_parallel_world_size()
    # Bypass the function if we are using only 1 GPU.
    if world_size == 1:
        return input_

    # Bf16 convert
    dt = input_.dtype
    if dt == torch.bfloat16 and get_fp32_allreduce():
        input_ = input_.float()

    # Split along last dimension.
    input_list = split_tensor_along_last_dim(input_, world_size)

    # Note: torch.split does not create contiguous tensors by default.
    rank = get_model_parallel_rank()
    output = input_list[rank].contiguous()

    # Bf16 convert
    if dt == torch.bfloat16 and get_fp32_allreduce():
        output = output.bfloat16()

    return output


def _gather(input_):
    """Gather tensors and concatinate along the last dimension."""

    world_size = get_model_parallel_world_size()
    # Bypass the function if we are using only 1 GPU.
    if world_size == 1:
        return input_

    # Bf16 convert
    dt = input_.dtype
    if dt == torch.bfloat16 and get_fp32_allreduce():
        input_ = input_.float()

    # Size and dimension.
    last_dim = input_.dim() - 1
    rank = get_model_parallel_rank()

    tensor_list = [torch.empty_like(input_) for _ in range(world_size)]
    tensor_list[rank] = input_
    torch.distributed.all_gather(tensor_list, input_, group=get_model_parallel_group())

    # Note: torch.cat already creates a contiguous tensor.
    output = torch.cat(tensor_list, dim=last_dim).contiguous()

    # Bf16 convert
    if dt == torch.bfloat16 and get_fp32_allreduce():
        output = output.bfloat16()

    return output


class _CopyToModelParallelRegion(torch.autograd.Function):
    """Pass the input to the model parallel region."""

    @staticmethod
    def symbolic(graph, input_):
        return input_

    @staticmethod
    def forward(ctx, input_):
        return input_

    @staticmethod
    def backward(ctx, grad_output):
        return _reduce(grad_output)


class _ReduceFromModelParallelRegion(torch.autograd.Function):
    """All-reduce the input from the model parallel region."""

    @staticmethod
    def symbolic(graph, input_):
        return _reduce(input_)

    @staticmethod
    def forward(ctx, input_):
        return _reduce(input_)

    @staticmethod
    def backward(ctx, grad_output):
        return grad_output


class _ScatterToModelParallelRegion(torch.autograd.Function):
    """Split the input and keep only the corresponding chuck to the rank."""

    @staticmethod
    def symbolic(graph, input_):
        return _split(input_)

    @staticmethod
    def forward(ctx, input_):
        return _split(input_)

    @staticmethod
    def backward(ctx, grad_output):
        return _gather(grad_output)


class _GatherFromModelParallelRegion(torch.autograd.Function):
    """Gather the input from model parallel region and concatinate."""

    @staticmethod
    def symbolic(graph, input_):
        return _gather(input_)

    @staticmethod
    def forward(ctx, input_):
        return _gather(input_)

    @staticmethod
    def backward(ctx, grad_output):
        return _split(grad_output)


def copy_to_model_parallel_region(input_):
    return _CopyToModelParallelRegion.apply(input_)


def reduce_from_model_parallel_region(input_):
    return _ReduceFromModelParallelRegion.apply(input_)


def scatter_to_model_parallel_region(input_):
    return _ScatterToModelParallelRegion.apply(input_)


def gather_from_model_parallel_region(input_):
    return _GatherFromModelParallelRegion.apply(input_)


def get_dtype_from_string(dtype_str):
    if type(dtype_str) == str:
        if dtype_str == "float32" or dtype_str == "fp32":
            return torch.float32
        elif dtype_str == "float16" or dtype_str == "fp16":
            return torch.float16
        elif dtype_str == "bfloat16" or dtype_str == "bf16":
            return torch.bfloat16
        else:
            raise ValueError(f"Unrecognized dtype {dtype_str}")
    else:
        return dtype_str


def get_init_from_string(init_str):
    if type(init_str) == str:
        if init_str == "torch.nn.init.zeros_":
            return torch.nn.init.zeros_


def _initialize_affine_weight_gpu(weight, init_method, partition_dim, stride=1, rng_fork=True):
    """Initialize affine weight for model parallel on GPU."""

    weight.model_parallel = True
    weight.partition_dim = partition_dim
    weight.partition_stride = stride

    if rng_fork:
        with deepspeed.checkpointing.get_cuda_rng_tracker().fork():
            init_method(weight)
    else:
        init_method(weight)


def _initialize_affine_weight_cpu(
    global_config,
    weight,
    output_size,
    input_size,
    per_partition_size,
    partition_dim,
    init_method,
    stride=1,
    return_master_weight=False,
):
    """Initialize affine weight for model parallel.

    Build the master weight on all processes and scatter
    the relevant chunk."""

    weight.model_parallel = True
    weight.partition_dim = partition_dim
    weight.partition_stride = stride

    # Initialize master weight
    master_weight = torch.empty(output_size, input_size, dtype=torch.float, requires_grad=False)
    init_method(master_weight)
    master_weight = master_weight.to(dtype=global_config.params_dtype)

    # Split and copy
    per_partition_per_stride_size = divide(per_partition_size, stride)
    weight_list = torch.split(master_weight, per_partition_per_stride_size, dim=partition_dim)
    rank = get_model_parallel_rank()
    world_size = get_model_parallel_world_size()
    my_weight_list = weight_list[rank::world_size]

    with torch.no_grad():
        torch.cat(my_weight_list, dim=partition_dim, out=weight)
    if return_master_weight:
        return master_weight
    return None


def _vocab_size_with_padding(orig_vocab_size, args):
    """Pad vocab size so it is divisible by model parallel size and
    still having GPU friendly size."""

    after = orig_vocab_size
    multiple = args.make_vocab_size_divisible_by * args.model_parallel_size
    while (after % multiple) != 0:
        after += 1
    if args.rank == 0:
        print(
            " > padded vocab (size: {}) with {} dummy tokens "
            "(new size: {})".format(orig_vocab_size, after - orig_vocab_size, after),
            flush=True,
        )
    return after
