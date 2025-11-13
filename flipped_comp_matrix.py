import os
import numpy as np
from scipy.stats import pearsonr
import pandas as pd
from scipy.cluster.hierarchy import dendrogram, linkage
import matplotlib
matplotlib.use('Agg')  # Using a non-GUI backend at PARAM himalaya
import matplotlib.pyplot as plt
from scipy.spatial.distance import squareform
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations


# we calculated the best Pearson correlation for a pair of sequences with original and flipped shorter sequence
def calculate_best_pearson(seq1, seq2):
    def get_best_corr(short, long):
        shifts = range(0, len(long) - len(short) + 1, 2)
        return max((pearsonr(short, long[s:s + len(short)])[0] for s in shifts), default=-1)

    if len(seq1) <= len(seq2):
        short, long = seq1, seq2
    else:
        short, long = seq2, seq1

    flipped_short = short[::-1]  # flipping the sequence

    best_original = get_best_corr(short, long)
    best_flipped = get_best_corr(flipped_short, long)

    return max(best_original, best_flipped)


# Parallelize Pearson coefficient calculations using chunked parallelism
def parallel_calculate(sequences, sequence_names, chunk_size=500):
    results = []
    sequence_pairs = list(combinations(range(len(sequence_names)), 2))  # we generated all sequence pairs
    chunks = [sequence_pairs[i:i + chunk_size] for i in range(0, len(sequence_pairs), chunk_size)]

    with ProcessPoolExecutor() as executor:
        for chunk in chunks:
            futures = [
                executor.submit(
                    calculate_best_pearson,
                    sequences[sequence_names[i]],
                    sequences[sequence_names[j]]
                )
                for i, j in chunk
            ]
            for future, (i, j) in zip(futures, chunk):
                results.append((sequence_names[i], sequence_names[j], future.result()))
    return results


def process_file(input_file, output_dir):
    print(f"Processing {input_file}...")

    # Read sequences from Excel
    df = pd.read_excel(input_file, index_col=0)
    sequences = {row: df.loc[row].dropna().tolist() for row in df.index}
    sequence_names = list(sequences.keys())
    print(f"Sequences read from {input_file}...")

    # Calculate Pearson coefficients
    pearson_results = parallel_calculate(sequences, sequence_names, chunk_size=500)
    pearson_df = pd.DataFrame(pearson_results,
                              columns=["Comparison Sequence", "Target Sequence", "Best Pearson Coefficient"])

    # Prepare distance matrix
    num_sequences = len(sequence_names)
    distance_matrix = np.zeros((num_sequences, num_sequences))  # Initialized with 0
    sequence_to_index = {seq: idx for idx, seq in enumerate(sequence_names)}

    for _, row in pearson_df.iterrows():
        comp_idx = sequence_to_index[row["Comparison Sequence"]]
        targ_idx = sequence_to_index[row["Target Sequence"]]
        coeff = row["Best Pearson Coefficient"]
        distance = 1 - coeff  # Converted Pearson to distance
        distance_matrix[comp_idx, targ_idx] = distance
        distance_matrix[targ_idx, comp_idx] = distance  # Symmetric matrix

    # Generating output file names based on input file
    base_name = os.path.splitext(os.path.basename(input_file))[0]
    output_dendrogram_file = os.path.join(output_dir, f"{base_name}_dendrogram.png")
    output_linkage_csv_file = os.path.join(output_dir, f"{base_name}_linkage_matrix.csv")
    output_full_distance_csv_file = os.path.join(output_dir, f"{base_name}_full_distance_matrix.csv")

    # Saving the full distance matrix
    distance_df = pd.DataFrame(distance_matrix, index=sequence_names, columns=sequence_names)
    distance_df.to_csv(output_full_distance_csv_file)
    print(f"Full distance matrix saved to {output_full_distance_csv_file}")

    # Performing hierarchical clustering
    condensed_matrix = squareform(distance_matrix)
    linkage_matrix = linkage(condensed_matrix, method="complete")
    print("Linkage matrix generated...")

    # Save the linkage matrix to a CSV file
    linkage_df = pd.DataFrame(
        linkage_matrix,
        columns=["Cluster 1", "Cluster 2", "Distance", "Number of Original Observations"]
    )
    linkage_df.to_csv(output_linkage_csv_file, index=False)
    print(f"Linkage matrix saved to {output_linkage_csv_file}")

    # Plot the dendrogram
    plt.figure(figsize=(20, 10))
    dendrogram(linkage_matrix, labels=sequence_names)
    plt.title("Hierarchical Clustering Dendrogram")
    plt.xlabel("Sequences")
    plt.ylabel("1 - Pearson Coefficient (Distance)")
    plt.tight_layout()
    plt.savefig(output_dendrogram_file)
    print(f"Dendrogram saved to {output_dendrogram_file}")


if __name__ == "__main__":
    # Define input and output directories
    input_dir = r"/scratch/shipras.sbb.iitmandi/flipped_90chrcomb/"
    output_dir = r"/scratch/shipras.sbb.iitmandi/flipped_90o/"

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Process each file in the input directory
    for file_name in os.listdir(input_dir):
        if file_name.endswith(".xlsx"):
            input_file_path = os.path.join(input_dir, file_name)
            process_file(input_file_path, output_dir)
