# Copyright (c) Together
# This software is distributed under the terms of the Apache License, Version 2.0
# Author: Michael Poli
import gc

import torch
import torch.nn as nn
import torch.nn.functional as F

try:
    import conv1d_cpp
except:
    pass
from .utils import column_split

IIR_PREFILL_MODES = [
    "recurrence",
    "modal-fft",
    "hybrid-modal-recurrence",
    "modal-scan",
    "canonical-fft",
    "iir-fir-caching",
]


def canonicalize_modal_system(poles, residues):
    """Canonicalize a modal system.

    Args:
        poles (Tensor): The poles of the system.
        residues (Tensor): The residues of the system.

    Returns:
        Tuple[Tensor, Tensor]: The canonicalized poles and residues.
    """
    raise NotImplementedError


def list_tensors(idx):
    for obj in gc.get_objects():
        try:
            if torch.is_tensor(obj) and isinstance(obj, torch.Tensor):
                # dump to log
                print(type(obj), obj.size())
                el = obj[0]
                with open(f"tensors_{idx}.txt", "a") as f:
                    f.write(f"{type(obj)} {obj.size()} {el}\n")
        except Exception as e:
            pass


class HyenaInferenceEngine:
    def __init__(
        self,
        fir_fn=None,
        iir_prefill_style="modal-fft",
        layer_idx=None,
    ) -> None:
        self.fir_fn = fir_fn
        assert iir_prefill_style in IIR_PREFILL_MODES, f"iir_prefill_style must be one of {IIR_PREFILL_MODES}"
        self.iir_prefill_style = iir_prefill_style
        self.layer_idx = layer_idx
        self.low_mem_mode = False

    def parallel_fir(
        self,
        fir_fn,
        u,
        weight,
        bias,
        L,
        fir_length=3,
        inference_params=None,
        prefill_mode=None,
        padding_mask=None,
    ):
        """Compute the output state of the long convolutional filter."""
        # prepare input layout, dimensions and dispatch to fir kernel
        if fir_fn != torch.nn.functional.conv1d:
            z_pre = fir_fn(u)[:, :L]  # B, L, D
            z_pre = z_pre.permute(0, 2, 1)
        else:
            u = u.permute(0, 2, 1)  # B, D, L
            z_pre = fir_fn(
                u,
                weight,
                bias=None,  # don't pass it here, add manually instead!  source of small error
                stride=1,
                padding=fir_length - 1,
                groups=u.shape[1],
            )[..., :L]

            # add manually instead!  source of small error
            z_pre = z_pre + bias[None, :, None]

        # handle padding post fir, the only place with biases
        if type(padding_mask) == torch.Tensor:
            z_pre = z_pre * padding_mask[:, None]

        if inference_params is not None:
            # handle seqlen last and dim last cases for `u`
            if fir_fn != torch.nn.functional.conv1d:
                fir_state = u[:, -fir_length + 1 :].permute(0, 2, 1)
            else:
                fir_state = u[..., -fir_length + 1 :]
        else:
            fir_state = None

        return z_pre, fir_state

    def parallel_iir(
        self,
        z_pre,
        h,
        D,
        L,
        poles,
        residues,
        t,
        dims,
        layer_idx,
        inference_params=None,
        prefill_style="fft",
        fftconv_fn=None,
        padding_mask=None,
        use_flashfft=False,
        column_split_hyena=False,
        long_fir_threshold=None,
    ):
        """Compute the output state of the short convolutional filter."""
        fft_size = 2 * L
        hidden_size, num_attention_heads, hidden_size_per_attention_head, _, _ = dims
        # Compatibility with training infra that column splits the projections
        if column_split_hyena:
            z = z_pre.reshape(
                z_pre.shape[0],
                num_attention_heads,
                3 * hidden_size_per_attention_head,
                z_pre.shape[2],
            )
            x2, x1, v = (
                z[:, :, :hidden_size_per_attention_head],
                z[
                    :,
                    :,
                    hidden_size_per_attention_head : 2 * hidden_size_per_attention_head,
                ],
                z[:, :, 2 * hidden_size_per_attention_head :],
            )
            x2, x1, v = (
                x2.reshape(x2.shape[0], -1, x2.shape[-1]),
                x1.reshape(x1.shape[0], -1, x1.shape[-1]),
                v.reshape(v.shape[0], -1, v.shape[-1]),
            )
        else:
            x2, x1, v = z_pre.split([hidden_size, hidden_size, hidden_size], dim=1)

        x1v = x1 * v

        if inference_params is not None and prefill_style == "recurrence":
            y = self.prefill_via_direct_recurrence(
                inference_params=inference_params,
                x1v=x1v,
                L=L,
                poles=poles,
                residues=residues,
            )

        else:
            if use_flashfft and (L % 2) == 0:  # only works with even L
                y = fftconv_fn(
                    x1v.to(dtype=torch.bfloat16).contiguous(),
                    h.to(dtype=torch.float32),
                )
                X_s = None

            elif long_fir_threshold is None:
                H = torch.fft.rfft(h.to(dtype=torch.float32), n=fft_size) / fft_size
                X_s = torch.fft.fft(x1v.to(dtype=torch.float32), n=fft_size)
                X = X_s[..., : H.shape[-1]]
                if len(z_pre.shape) > 3:
                    H = H.unsqueeze(1)
                y = torch.fft.irfft(X * H, n=fft_size, norm="forward")[..., :L]

            else:
                assert h.shape[0] == 1, "batch size must be 1 for long_fir_threshold"
                h = h[0][:, None]  # rearrange to d, 1, l for depthwise conv1d
                h = h[..., :long_fir_threshold]
                y = F.conv1d(
                    x1v,
                    h.to(dtype=x1v.dtype),
                    stride=1,
                    groups=x1v.shape[1],
                    padding=h.shape[-1] - 1,
                )[..., :L]

        y = y.to(dtype=x1v.dtype)
        y = (y + x1v * D.unsqueeze(-1)) * x2

        if inference_params is not None:
            if prefill_style == "fft":
                self.prefill_via_modal_fft(
                    inference_params=inference_params,
                    x1v=x1v,
                    X_s=X_s,
                    L=L,
                    t=t,
                    poles=poles,
                    dims=dims,
                    layer_idx=layer_idx,
                    use_flashfft=use_flashfft,
                    fftconv_fn=fftconv_fn,
                )

            elif prefill_style == "recurrence":
                # recurrent prefill is done before
                pass
            else:
                raise NotImplementedError
            if self.low_mem_mode:
                # TODO: smarter gc
                del z_pre, x2, x1, v, x1v, h, poles, residues
                torch.cuda.empty_cache()

        return y.permute(0, 2, 1)

    def step_fir(self, u, fir_state, weight, bias=None):
        """Step the FIR filter.

        Note:
        `fir_state` contains the last `short_filter_length - 1` elements of `u`: `u_(L-2), u_{L-1), ...`
        We assume dimensions of `short_filter_weight` to be `[d, 1, short_filter_len]` (SISO / multi SISO layout).
        """
        h0, h = weight[..., 0, -1], weight[..., 0, :-1]
        h0, h = h0[None], h[None]
        y = h0 * u + torch.sum(fir_state * h, dim=-1) + bias

        # update
        fir_state = torch.roll(fir_state, -1, dims=2)
        fir_state[..., -1] = u
        return y, fir_state

    def step_iir(self, x2, x1, v, D, residues, poles, iir_state, iir_groups=1):
        x1v = x1 * v

        residues, poles = (
            torch.view_as_complex(residues.to(torch.float32)),
            torch.view_as_complex(poles.to(torch.float32)),
        )
        # squeeze the dummy seqlen dimension
        # D, state_dim, 1 -> 1, D, state_dim
        residues, poles = residues[..., 0][None], poles[..., 0][None]
        iir_state = poles * iir_state + x1v[..., None]

        res_state = torch.sum(residues * iir_state, dim=-1).real

        if iir_groups > 1:
            raise NotImplementedError
        y = x2 * (res_state + D * x1v)

        return y, iir_state

    def prefill_via_fir_caching(self, u, inference_params, L, *args, **kwargs):
        """Turns the IIR filter into a FIR and uses a cache for decoding."""
        raise NotImplementedError(":)")

    def prefill_via_direct_recurrence(
        self, inference_params, x1v, L, residues, poles, *args, **kwargs
    ) -> torch.Tensor:
        """
        Compute the IIR state via explicit SSM recurrence (modal form)

        This is the most memory efficient prefilling method for Hyena filters.

        Note:
            dtypes: [state: float32, poles: float32, x1v: bfloat16, output: bfloat16]
        """
        state_dim = poles.shape[1]
        x1v_ = x1v[..., None, None]  # b, d, l, sdim, reim
        x1v_ = x1v_.repeat(1, 1, 1, state_dim, 2)  # b, d, l, sdim, reim
        x1v_[..., 1] = 0

        state = 0 * x1v_[:, :, 0]
        output = 0 * x1v_[:, :, :, 0, 0]  # b, d, l

        # suppress dummy seqlen dimension
        poles = poles[:, :, 0][None]
        residues = residues[:, :, 0][None].repeat(x1v_.shape[0], 1, 1, 1)  # b, d, sdim, reim

        # state: b, d, sdim, reim
        # poles: 1, d, sdim, reim
        # x1v_: b, d, l, sdim, reim
        for i in range(L):
            state[..., 0] = poles[..., 0] * state[..., 0] - poles[..., 1] * state[..., 1] + x1v_[:, :, i, :, 0]
            state[..., 1] = poles[..., 0] * state[..., 1] + poles[..., 1] * state[..., 0] + x1v_[:, :, i, :, 1] 
            output[:, :, i] = torch.sum(residues * state, dim=-2)[..., 0]  # .real
            
        inference_params.state_dict[self.layer_idx] = torch.view_as_complex(state.to(dtype=torch.float32))

        return output

    def prefill_via_hybrid_recurrence(self, inference_params, u, log_poles, x1v_f_a, L, *args, **kwargs):
        """
        Compute the IIR state via hybrid recurrence-convolution over blocks
        """
        raise NotImplementedError(":)")

    def prefill_via_scan(self, u, inference_params=None, *args, **kwargs):
        raise NotImplementedError

    def prefill_via_canonical_fft(self, u, inference_params=None, *args, **kwargs):
        """
        Compute the IIR state via a single FFT with the denominator of the SSM in companion form.

        This is the most memory efficient "parallelized" prefilling method for Hyena.

        From: https://arxiv.org/abs/2310.18780
        """
        raise NotImplementedError(":)")

    def prefill_via_modal_fft(
        self,
        inference_params,
        x1v,
        L,
        poles,
        t,
        dims,
        layer_idx,
        X_s=None,
        use_flashfft=False,
        fftconv_fn=None,
        state_dtype=torch.complex64,
        *args,
        **kwargs,
    ):
        """
        Compute the IIR state via a single FFT, using the poles of the SSM in modal form.
        """
        # When the model has a long convolution derived from a SSM in modal form and prefill_style is "fft",
        # we split the filter into poles and residues and reuse FFT computation on the input.
        # This optimization is currently not supported when using flashfftconv.
        hidden_size, _, _, state_size, hyena_filter_groups = dims

        if use_flashfft:
            # using real states
            poles = poles.squeeze().reshape(poles.shape[0], -1)[..., None]

            state_s = poles**t
            if hyena_filter_groups > 1:
                raise NotImplementedError

            x1v = x1v[:, :, None].repeat(1, 1, 2 * state_size, 1)
            x1v = x1v.reshape(x1v.shape[0], -1, x1v.shape[-1])
            state_s = state_s[None]

            state = fftconv_fn(
                x1v.contiguous(),
                state_s.to(dtype=torch.float32),
            )
            state = state[..., L - 1].reshape(x1v.shape[0], hidden_size, state_size, 2)
            state = torch.view_as_complex(state.contiguous().to(dtype=torch.float32))
            inference_params.state_dict[self.layer_idx] = state
        else:
            assert X_s is not None
            bs = x1v.shape[0]
            fft_size = 2 * L
            poles = torch.view_as_complex(poles.to(torch.float32))
            state_s = poles**t
            state_S = torch.fft.fft(state_s, n=fft_size).repeat(bs, 1, 1, 1)  # B, D, state_dim, 2 * L
            if hyena_filter_groups > 1:
                state_S = state_S.repeat_interleave(hidden_size // hyena_filter_groups, 1)
            state = torch.fft.ifft(X_s[..., None, :] * state_S, n=fft_size)
            inference_params.state_dict[layer_idx] = state[..., L - 1].to(dtype=state_dtype)

    def _compute_state(self, log_poles, u, t, L, *args, **kwargs):
        """
        Compute the IIR state given an input `u` and log_poles of the modal system.
        """
        bs = u.shape[0]
        fft_size = 2 * L
        U = torch.fft.rfft(u.to(torch.float32), n=fft_size)
        fft_size = 2 * L
        x = (log_poles * t).exp()
        # [batch, hidden_size, state_dim, 2 * seqlen]
        X = torch.fft.fft(x, n=fft_size).repeat(bs, 1, 1, 1)
        state = torch.fft.ifft(U[..., None, :] * X, n=fft_size)[..., :L]
        return state
