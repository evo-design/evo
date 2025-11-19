import argparse
import ast
import multiprocessing as mp
import time
from functools import partial
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pandas as pd
from sklearn.feature_extraction.text import CountVectorizer


def parse_hairpins(hairpin_str: str) -> List[Tuple[int, int, int, int]]:
    """Return parsed hairpin annotations from their string representation."""
    if pd.isna(hairpin_str) or hairpin_str == "[]":
        return []
    try:
        return ast.literal_eval(hairpin_str)
    except Exception:
        return []


def extract_structure_features(structure: str, mfe: float, hairpins: str) -> np.ndarray:
    """Convert dot-bracket structures into a numeric feature vector for filtering."""
    if not structure or pd.isna(structure):
        return np.zeros(10)

    # Basic structure features
    length = len(structure)
    num_pairs = structure.count("(")
    num_unpaired = structure.count(".")
    pairing_ratio = num_pairs / length if length > 0 else 0

    # Stem and loop regions
    stems = []
    loops = []
    current_stem = 0
    current_loop = 0

    for char in structure:
        if char == "(":
            if current_loop > 0:
                loops.append(current_loop)
                current_loop = 0
            current_stem += 1
        elif char == ")":
            if current_stem > 0:
                stems.append(current_stem)
                current_stem = 0
        else:  # '.'
            if current_stem > 0:
                stems.append(current_stem)
                current_stem = 0
            current_loop += 1

    # Add final regions
    if current_stem > 0:
        stems.append(current_stem)
    if current_loop > 0:
        loops.append(current_loop)

    # Statistical features
    avg_stem_length = np.mean(stems) if stems else 0
    avg_loop_length = np.mean(loops) if loops else 0
    max_stem_length = max(stems) if stems else 0
    num_stems = len(stems)

    # MFE features
    mfe_val = float(mfe) if not pd.isna(mfe) else 0
    mfe_per_nucleotide = mfe_val / length if length > 0 else 0

    # Hairpin features
    hairpin_list = parse_hairpins(hairpins)
    num_hairpins = len(hairpin_list)

    return np.array(
        [
            length,
            num_pairs,
            pairing_ratio,
            avg_stem_length,
            avg_loop_length,
            max_stem_length,
            num_stems,
            mfe_val,
            mfe_per_nucleotide,
            num_hairpins,
        ]
    )


def create_structure_index(
    df: pd.DataFrame, structure_type: str = "both"
) -> Dict[str, np.ndarray]:
    """Prepare feature matrices for the requested structure types."""
    print("Creating structure index...")

    indices = {}

    if structure_type in ["consensus", "both"]:
        print("  Indexing consensus structures...")
        consensus_features = []
        for _, row in df.iterrows():
            features = extract_structure_features(
                row["Consensus Structure"],
                row["Consensus MFE"],
                row["Consensus Hairpins"],
            )
            consensus_features.append(features)
        indices["consensus"] = np.array(consensus_features)

    if structure_type in ["full", "both"]:
        print("  Indexing full structures...")
        full_features = []
        for _, row in df.iterrows():
            features = extract_structure_features(
                row["Full Structure"], row["Full MFE"], row["Full Hairpins"]
            )
            full_features.append(features)
        indices["full"] = np.array(full_features)

    return indices


def cosine_similarity_filter(
    query_features: np.ndarray, target_features: np.ndarray, threshold: float = 0.8
) -> np.ndarray:
    """Apply cosine similarity between feature vectors and return candidates above a threshold."""
    # Normalize features to handle different scales
    query_norm = query_features / (np.linalg.norm(query_features) + 1e-8)
    target_norm = target_features / (
        np.linalg.norm(target_features, axis=1, keepdims=True) + 1e-8
    )

    # Cosine similarity
    similarities = np.dot(target_norm, query_norm)

    # Return indices of candidates above threshold
    return np.where(similarities >= threshold)[0]


def structure_to_kmer_features(structure: str, k: int = 3) -> str:
    """Represent a dot-bracket structure as a bag of structural k-mers."""
    if not structure:
        return ""

    # Create k-mers from STRUCTURE (dot-bracket notation)
    # This captures structural motifs like "(((" or ").))"
    kmers = []
    for i in range(len(structure) - k + 1):
        kmer = structure[i : i + k]
        # Only include k-mers with structural information
        if any(c in kmer for c in "()"):
            kmers.append(kmer)

    return " ".join(kmers)


def create_structural_motif_index(
    df: pd.DataFrame, structure_type: str = "both"
) -> Dict[str, Dict[str, object]]:
    """Vectorize structural motifs for later similarity computations."""
    print("Creating structural motif index...")

    indices = {}

    if structure_type in ["consensus", "both"]:
        print("  Indexing consensus structural motifs...")
        consensus_motifs = []
        for _, row in df.iterrows():
            motifs = extract_structural_motifs(row["Consensus Structure"])
            consensus_motifs.append(" ".join(motifs))

        # Create binary matrix of structural motifs
        vectorizer = CountVectorizer(analyzer="word", token_pattern=r"\S+", binary=True)
        consensus_matrix = vectorizer.fit_transform(consensus_motifs)
        indices["consensus"] = {
            "matrix": consensus_matrix,
            "vectorizer": vectorizer,
            "motifs": consensus_motifs,
        }

    if structure_type in ["full", "both"]:
        print("  Indexing full structural motifs...")
        full_motifs = []
        for _, row in df.iterrows():
            motifs = extract_structural_motifs(row["Full Structure"])
            full_motifs.append(" ".join(motifs))

        vectorizer = CountVectorizer(analyzer="word", token_pattern=r"\S+", binary=True)
        full_matrix = vectorizer.fit_transform(full_motifs)
        indices["full"] = {
            "matrix": full_matrix,
            "vectorizer": vectorizer,
            "motifs": full_motifs,
        }

    return indices


def extract_structural_motifs(structure: str) -> List[str]:
    """Return a list of structural motifs (stems, hairpins, bulges, etc.)."""
    if not structure:
        return []

    motifs = []

    # 1. Stem-loop patterns
    i = 0
    while i < len(structure):
        if structure[i] == "(":
            # Found start of stem
            stem_start = i
            stem_count = 1
            i += 1

            # Count consecutive opening brackets
            while i < len(structure) and structure[i] == "(":
                stem_count += 1
                i += 1

            # Count loop region
            loop_count = 0
            while i < len(structure) and structure[i] == ".":
                loop_count += 1
                i += 1

            # Count closing brackets
            close_count = 0
            while i < len(structure) and structure[i] == ")":
                close_count += 1
                i += 1

            # Create motif descriptor
            if close_count > 0:
                if loop_count == 0:
                    # Perfect stem
                    motifs.append(f"STEM_{min(stem_count, close_count)}")
                else:
                    # Hairpin: stem_length:loop_length
                    motifs.append(
                        f"HAIRPIN_{min(stem_count, close_count)}:{loop_count}"
                    )
        else:
            i += 1

    # 2. Bulge patterns (internal loops)
    # Look for patterns like (((.)))
    bulge_pattern = ""
    in_stem = False
    for char in structure:
        if char == "(":
            if not in_stem:
                bulge_pattern = "("
                in_stem = True
            else:
                bulge_pattern += char
        elif char == ")":
            bulge_pattern += char
            if bulge_pattern.count("(") == bulge_pattern.count(")"):
                # Complete bulge pattern
                dots = bulge_pattern.count(".")
                stems = bulge_pattern.count("(")
                if dots > 0 and stems > 1:
                    motifs.append(f"BULGE_{stems}:{dots}")
                bulge_pattern = ""
                in_stem = False
        elif char == "." and in_stem:
            bulge_pattern += char

    # 3. Multi-loop patterns
    # Count branching complexity
    max_depth = 0
    current_depth = 0
    for char in structure:
        if char == "(":
            current_depth += 1
            max_depth = max(max_depth, current_depth)
        elif char == ")":
            current_depth -= 1

    if max_depth > 0:
        motifs.append(f"DEPTH_{max_depth}")

    # 4. Consecutive unpaired regions
    unpaired_lengths = []
    current_unpaired = 0
    for char in structure:
        if char == ".":
            current_unpaired += 1
        else:
            if current_unpaired > 0:
                unpaired_lengths.append(current_unpaired)
                current_unpaired = 0
    if current_unpaired > 0:
        unpaired_lengths.append(current_unpaired)

    # Add unpaired region motifs
    for length in unpaired_lengths:
        if length >= 3:  # Only significant unpaired regions
            motifs.append(f"UNPAIRED_{min(length, 10)}")  # Cap at 10

    return motifs


def motif_similarity_filter(
    query_motifs: str, target_motif_matrix, vectorizer, threshold: float = 0.3
) -> np.ndarray:
    """Calculate motif overlap (Jaccard-like) and return indices above the threshold."""
    if not query_motifs.strip():
        return np.array([])

    # Transform query motifs to vector
    query_vector = vectorizer.transform([query_motifs])

    # Calculate Jaccard-like similarity with all targets
    # |intersection| / |union| approximated by dot product / (sum - dot product)
    intersection = query_vector.dot(target_motif_matrix.T).toarray().flatten()
    query_sum = query_vector.sum()
    target_sums = np.array(target_motif_matrix.sum(axis=1)).flatten()

    # Avoid division by zero
    union = query_sum + target_sums - intersection
    similarities = np.divide(
        intersection,
        union,
        out=np.zeros_like(intersection, dtype=float),
        where=union != 0,
    )

    # Return indices above threshold
    return np.where(similarities >= threshold)[0]


def pure_structural_pre_filter(
    query_structure: str, target_structures: List[str], threshold: float = 0.4
) -> np.ndarray:
    """Perform a structural property check (length/pairing ratios) before heavier comparisons."""
    if not query_structure:
        return np.array([])

    # Query structural properties
    q_len = len(query_structure)
    q_pairs = query_structure.count("(")
    q_unpaired = query_structure.count(".")
    q_complexity = len(set(query_structure))  # Number of different symbols

    candidates = []

    for i, target_struct in enumerate(target_structures):
        if not target_struct:
            continue

        # Target structural properties
        t_len = len(target_struct)
        t_pairs = target_struct.count("(")
        t_unpaired = target_struct.count(".")

        # Quick similarity checks
        # 1. Length similarity
        len_sim = 1.0 - abs(q_len - t_len) / max(q_len, t_len)
        if len_sim < 0.5:  # Very different lengths
            continue

        # 2. Pairing ratio similarity
        q_pair_ratio = q_pairs / q_len if q_len > 0 else 0
        t_pair_ratio = t_pairs / t_len if t_len > 0 else 0
        pair_sim = 1.0 - abs(q_pair_ratio - t_pair_ratio)

        # 3. Combined quick score
        quick_score = 0.6 * len_sim + 0.4 * pair_sim

        if quick_score >= threshold:
            candidates.append(i)

    return np.array(candidates)


def base_pair_distance(struct1: str, struct2: str) -> float:
    """Compute one minus the normalized overlap between base-pair sets of two structures."""
    if not struct1 or not struct2:
        return 1.0

    # Quick length check
    len_diff = abs(len(struct1) - len(struct2))
    max_len = max(len(struct1), len(struct2))
    if len_diff / max_len > 0.5:  # Very different lengths
        return 1.0

    def get_pairs(structure: str) -> set[Tuple[int, int]]:
        chars = np.array(list(structure))
        opens = np.where(chars == "(")[0]
        closes = np.where(chars == ")")[0]

        if len(opens) != len(closes):
            return set()

        pairs = set()
        stack = []
        for i, char in enumerate(structure):
            if char == "(":
                stack.append(i)
            elif char == ")" and stack:
                pairs.add((stack.pop(), i))
        return pairs

    pairs1 = get_pairs(struct1)
    pairs2 = get_pairs(struct2)

    if not pairs1 and not pairs2:
        return 0.0

    intersection = len(pairs1.intersection(pairs2))
    union = len(pairs1.union(pairs2))

    return 1.0 - (intersection / union) if union > 0 else 1.0


def optimized_similarity_search(
    query_row_data,
    target_df: pd.DataFrame,
    target_indices: Dict[str, np.ndarray],
    motif_indices: Dict[str, Dict[str, object]],
    structure_type: str = "both",
    min_similarity: float = 0.7,
    pre_filter_threshold: float = 0.4,
) -> List[Dict[str, object]]:
    """Run the multi-stage similarity search for a single query row."""
    query_idx, query_row = query_row_data
    results = []

    # Process each structure type requested
    structure_types = []
    if structure_type == "both":
        structure_types = ["consensus", "full"]
    else:
        structure_types = [structure_type]

    for struct_type in structure_types:
        query_structure = query_row[f"{struct_type.title()} Structure"]

        if not query_structure or pd.isna(query_structure):
            continue

        # Stage 1: quick structural property filter
        # Extract ALL target structures of this type
        target_structures = []
        for _, target_row in target_df.iterrows():
            target_struct = target_row[f"{struct_type.title()} Structure"]
            target_structures.append(
                target_struct if not pd.isna(target_struct) else ""
            )

        quick_candidates = pure_structural_pre_filter(
            query_structure, target_structures, threshold=0.3
        )

        if len(quick_candidates) == 0:
            continue

        print(
            f"  {struct_type.title()} Stage 1: {len(quick_candidates)}/{len(target_structures)} candidates after structural pre-filter"
        )

        # Stage 2: Structural motif filtering
        query_motifs = " ".join(extract_structural_motifs(query_structure))
        motif_index = motif_indices[struct_type]

        if len(quick_candidates) > 0:
            # Get motif matrix for candidates only
            candidate_motif_matrix = motif_index["matrix"][quick_candidates]
            motif_candidates_rel = motif_similarity_filter(
                query_motifs,
                candidate_motif_matrix,
                motif_index["vectorizer"],
                threshold=0.2,
            )
            motif_candidates = (
                quick_candidates[motif_candidates_rel]
                if len(motif_candidates_rel) > 0
                else np.array([])
            )
        else:
            motif_candidates = np.array([])

        if len(motif_candidates) == 0:
            continue

        print(
            f"  {struct_type.title()} Stage 2: {len(motif_candidates)} candidates after motif filter"
        )

        # Stage 3: Feature-based filtering
        query_features = extract_structure_features(
            query_structure,
            query_row[f"{struct_type.title()} MFE"],
            query_row[f"{struct_type.title()} Hairpins"],
        )

        if len(motif_candidates) > 0:
            # Get feature vectors for motif candidates only
            target_features = target_indices[struct_type][motif_candidates]
            feature_candidates_rel = cosine_similarity_filter(
                query_features, target_features, pre_filter_threshold
            )
            final_candidates = (
                motif_candidates[feature_candidates_rel]
                if len(feature_candidates_rel) > 0
                else np.array([])
            )
        else:
            final_candidates = np.array([])

        if len(final_candidates) == 0:
            continue

        print(
            f"  {struct_type.title()} Stage 3: {len(final_candidates)} candidates after feature filter"
        )

        # Stage 4: Detailed similarity calculation
        detailed_matches = 0
        for target_idx in final_candidates:
            target_row = target_df.iloc[target_idx]

            # Skip self-comparison
            if query_row["Root ID"] == target_row["Root ID"]:
                continue

            # Get target structure for this specific type
            target_structure = target_row[f"{struct_type.title()} Structure"]
            if not target_structure or pd.isna(target_structure):
                continue

            # Detailed base pair distance calculation
            struct_sim = 1.0 - base_pair_distance(query_structure, target_structure)

            if struct_sim >= min_similarity:
                # Calculate other similarities
                query_mfe = query_row[f"{struct_type.title()} MFE"]
                target_mfe = target_row[f"{struct_type.title()} MFE"]

                if not pd.isna(query_mfe) and not pd.isna(target_mfe):
                    mfe_diff = abs(float(query_mfe) - float(target_mfe))
                    mfe_sim = max(0.0, 1.0 - mfe_diff / 10.0)
                else:
                    mfe_sim = 0.0

                # Combined score
                combined_score = 0.7 * struct_sim + 0.3 * mfe_sim

                if combined_score >= min_similarity:
                    # Calculate complementary structure similarity for comparison
                    complementary_type = (
                        "full" if struct_type == "consensus" else "consensus"
                    )
                    comp_query_structure = query_row.get(
                        f"{complementary_type.title()} Structure", ""
                    )
                    comp_target_structure = target_row.get(
                        f"{complementary_type.title()} Structure", ""
                    )

                    # Calculate complementary structure similarity if both structures exist
                    comp_struct_sim = 0.0
                    comp_mfe_sim = 0.0
                    comp_combined_score = 0.0

                    if (
                        comp_query_structure
                        and comp_target_structure
                        and not pd.isna(comp_query_structure)
                        and not pd.isna(comp_target_structure)
                    ):
                        comp_struct_sim = 1.0 - base_pair_distance(
                            comp_query_structure, comp_target_structure
                        )

                        comp_query_mfe = query_row.get(
                            f"{complementary_type.title()} MFE", None
                        )
                        comp_target_mfe = target_row.get(
                            f"{complementary_type.title()} MFE", None
                        )

                        if not pd.isna(comp_query_mfe) and not pd.isna(comp_target_mfe):
                            comp_mfe_diff = abs(
                                float(comp_query_mfe) - float(comp_target_mfe)
                            )
                            comp_mfe_sim = max(0.0, 1.0 - comp_mfe_diff / 10.0)

                        comp_combined_score = 0.7 * comp_struct_sim + 0.3 * comp_mfe_sim

                    result = {
                        "Query_ID": query_row["Root ID"],
                        "Target_ID": target_row["Root ID"],
                        "Structure_Type": struct_type,
                        "Similarity_Score": combined_score,
                        "Structure_Similarity": struct_sim,
                        "MFE_Similarity": mfe_sim,
                        "Query_MFE": query_mfe,
                        "Target_MFE": target_mfe,
                        "Query_Structure": query_structure,
                        "Target_Structure": target_structure,
                        "Query_Hairpins": query_row[f"{struct_type.title()} Hairpins"],
                        "Target_Hairpins": target_row[
                            f"{struct_type.title()} Hairpins"
                        ],
                        # Add complementary structure information
                        "Complementary_Type": complementary_type,
                        "Complementary_Similarity_Score": comp_combined_score,
                        "Complementary_Structure_Similarity": comp_struct_sim,
                        "Complementary_MFE_Similarity": comp_mfe_sim,
                        "Query_Complementary_Structure": comp_query_structure,
                        "Target_Complementary_Structure": comp_target_structure,
                        "Query_Complementary_MFE": query_row.get(
                            f"{complementary_type.title()} MFE", None
                        ),
                        "Target_Complementary_MFE": target_row.get(
                            f"{complementary_type.title()} MFE", None
                        ),
                        "Query_Complementary_Hairpins": query_row.get(
                            f"{complementary_type.title()} Hairpins", []
                        ),
                        "Target_Complementary_Hairpins": target_row.get(
                            f"{complementary_type.title()} Hairpins", []
                        ),
                        # Add both consensus and full sequences for query
                        "Query_Consensus_RNA": query_row.get("Consensus RNA", ""),
                        "Query_Consensus_DNA": query_row.get("Sequence", ""),
                        "Query_Full_RNA": query_row.get("Full RNA", ""),
                        "Query_Full_DNA": query_row.get("Full TRF Region", ""),
                        # Add both consensus and full sequences for target
                        "Target_Consensus_RNA": target_row.get("Consensus RNA", ""),
                        "Target_Consensus_DNA": target_row.get("Sequence", ""),
                        "Target_Full_RNA": target_row.get("Full RNA", ""),
                        "Target_Full_DNA": target_row.get("Full TRF Region", ""),
                    }
                    results.append(result)
                    detailed_matches += 1

        print(f"  {struct_type.title()} Stage 4: {detailed_matches} final matches")

    return results


def optimized_structural_search(
    query_csv: Path,
    target_csv: Path,
    output_csv: Path,
    structure_type: str = "both",
    min_similarity: float = 0.7,
    max_results: Optional[int] = None,
    n_cpus: Optional[int] = None,
    pre_filter_threshold: float = 0.6,
    batch_size: int = 100,
) -> pd.DataFrame:
    """Coordinate indexing, batching, and filtering across the full query set."""
    if n_cpus is None:
        n_cpus = min(mp.cpu_count(), 16)  # Don't use all CPUs for indexing

    print(f"Loading data...")
    start_time = time.time()

    query_df = pd.read_csv(query_csv)
    target_df = pd.read_csv(target_csv)

    print(f"Query sequences: {len(query_df)}")
    print(f"Target sequences: {len(target_df)}")
    print(f"Total potential comparisons: {len(query_df) * len(target_df):,}")

    # Create indices
    print(f"Creating search indices...")
    target_feature_indices = create_structure_index(target_df, structure_type)
    target_motif_indices = create_structural_motif_index(target_df, structure_type)

    index_time = time.time() - start_time
    print(f"Index creation completed in {index_time:.1f} seconds")

    # Process in batches to manage memory
    print(f"Starting optimized search with {n_cpus} CPUs...")
    search_start = time.time()

    all_results = []
    total_batches = (len(query_df) + batch_size - 1) // batch_size

    for batch_num in range(total_batches):
        batch_start = batch_num * batch_size
        batch_end = min((batch_num + 1) * batch_size, len(query_df))

        print(
            f"Processing batch {batch_num + 1}/{total_batches} (sequences {batch_start}-{batch_end})"
        )

        batch_queries = query_df.iloc[batch_start:batch_end]

        # Parallel processing of batch
        search_func = partial(
            optimized_similarity_search,
            target_df=target_df,
            target_indices=target_feature_indices,
            motif_indices=target_motif_indices,
            structure_type=structure_type,
            min_similarity=min_similarity,
            pre_filter_threshold=pre_filter_threshold,
        )

        with mp.Pool(n_cpus) as pool:
            batch_results = pool.map(search_func, batch_queries.iterrows())

        # Flatten batch results
        for result_list in batch_results:
            all_results.extend(result_list)

        print(
            f"Batch {batch_num + 1} completed. Running total: {len(all_results)} matches"
        )

    search_time = time.time() - search_start
    total_time = time.time() - start_time

    print(f"Search completed in {search_time:.1f} seconds")
    print(f"Total runtime: {total_time:.1f} seconds")

    if not all_results:
        print("No matches found above similarity threshold")
        return pd.DataFrame()

    # Convert to DataFrame and process results
    results_df = pd.DataFrame(all_results)
    results_df = results_df.sort_values("Similarity_Score", ascending=False)

    if max_results:
        results_df = results_df.head(max_results)

    # Save results
    results_df.to_csv(output_csv, index=False)

    print(f"\nPerformance Summary:")
    print(f"Matches found: {len(results_df)}")
    print(f"Results saved to: {output_csv}")

    return results_df


def main() -> Optional[pd.DataFrame]:
    parser = argparse.ArgumentParser(
        description="Optimized RNA Structure Similarity Search"
    )
    parser.add_argument("--query", required=True, help="Query CSV file")
    parser.add_argument("--target", required=True, help="Target CSV file")
    parser.add_argument("--output", required=True, help="Output CSV file")
    parser.add_argument(
        "--structure-type",
        choices=["consensus", "full", "both"],
        default="both",
        help="Structure types to compare",
    )
    parser.add_argument(
        "--min-similarity", type=float, default=0.7, help="Minimum similarity threshold"
    )
    parser.add_argument(
        "--pre-filter-threshold",
        type=float,
        default=0.6,
        help="Pre-filter threshold (lower = more candidates)",
    )
    parser.add_argument("--max-results", type=int, help="Maximum results to return")
    parser.add_argument("--cpus", type=int, help="Number of CPU cores")
    parser.add_argument(
        "--batch-size", type=int, default=100, help="Batch size for processing"
    )

    args = parser.parse_args()

    results_df = optimized_structural_search(
        query_csv=args.query,
        target_csv=args.target,
        output_csv=args.output,
        structure_type=args.structure_type,
        min_similarity=args.min_similarity,
        max_results=args.max_results,
        n_cpus=args.cpus,
        pre_filter_threshold=args.pre_filter_threshold,
        batch_size=args.batch_size,
    )

    return results_df


if __name__ == "__main__":
    main()
