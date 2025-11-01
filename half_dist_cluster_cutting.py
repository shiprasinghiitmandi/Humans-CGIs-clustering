import pandas as pd
import numpy as np
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster, dendrogram
import matplotlib
matplotlib.use('Agg')  # GUI backend
import matplotlib.pyplot as plt
import os

def cut_clusters_from_upper_triangular_matrix(input_distance_file, height_threshold, output_file):
    """
    Cuts clusters directly from an upper triangular distance matrix at a specified height threshold.

    Parameters:
        input_distance_file (str): Path to the upper triangular distance matrix file (CSV format).
        height_threshold (float): The height at which to cut the dendrogram.
        output_file (str): Path to save the cluster information and dendrogram.
    """
    # Load the upper triangular distance matrix and fill empty cells with zeros
    upper_triangular_matrix = pd.read_csv(input_distance_file, index_col=0).fillna(0)
    sequence_names = upper_triangular_matrix.columns.tolist()

    # Extract the condensed form of the distance matrix
    condensed_matrix = squareform(upper_triangular_matrix.values, checks=False)

    # Generate a linkage matrix using the condensed distance matrix
    linkage_matrix = linkage(condensed_matrix, method='complete')

    # Perform clustering by cutting the dendrogram at the specified height
    cluster_labels = fcluster(linkage_matrix, t=height_threshold, criterion='distance')

    # Group sequences into clusters
    clusters = {}
    for instance, cluster_id in enumerate(cluster_labels):
        if cluster_id not in clusters:
            clusters[cluster_id] = []
        clusters[cluster_id].append(sequence_names[instance])

    # Save cluster information to the output file
    with open(output_file, 'w') as f:
        f.write(f"Height Threshold: {height_threshold}\n")
        f.write(f"Number of Clusters Found: {len(clusters)}\n\n")
        for cluster_id, sequences in clusters.items():
            f.write(f"Cluster {cluster_id} ({len(sequences)} instances): {sequences}\n")

    print(f"Cluster information saved to {output_file}")

    # Plot and save the dendrogram with single cluster labels
    plt.figure(figsize=(18, 10), dpi=150)

    # Generate dendrogram
    dendrogram(linkage_matrix, color_threshold=height_threshold)
    plt.axhline(y=height_threshold, color='r', linestyle='--', label=f'Cut Height: {height_threshold}')
    plt.legend()
    plt.title('flipped_90thperc_clustercut_1.15')
    plt.xlabel('Cluster_IDs')
    plt.ylabel('Distance')
    dendrogram_output_file = os.path.splitext(output_file)[0] + '_dendrogram.png'
    plt.savefig(dendrogram_output_file)
    plt.close()
    plt.show()
    print(f"Dendrogram with cut height saved to {dendrogram_output_file}")

# Example usage
input_distance_file = r"/scratch/shipras.sbb.iitmandi/1tilly_90th_flipped.csv" # Replace with your distance matrix file path
height_threshold = 1.15
output_file = r"/scratch/shipras.sbb.iitmandi/flipped_90th_clustcut/fl_90thperc1.15.txt"  # Replace with desired output file path

cut_clusters_from_upper_triangular_matrix(input_distance_file, height_threshold, output_file)

