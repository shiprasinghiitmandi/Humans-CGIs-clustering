import os
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score



input_distance_file = (
    r"1tilly_90th_flipped.csv" # distance matrix of 90th percentile as input file
)

output_directory = (
    r"KMeans_90d" # output file containing silhouette scores obtained after k-means clustering
)

min_clusters = 2
max_clusters = 20


current_cluster_number = 5



# K-means settings
random_state = 42
n_init = 50


os.makedirs(output_directory, exist_ok=True)


# ============================================================
# LOAD DISTANCE MATRIX
# ============================================================

print("Loading distance matrix...")

distance_df = pd.read_csv(
    input_distance_file,
    index_col=0
)

sequence_names = distance_df.columns.tolist()

raw_matrix = distance_df.to_numpy(dtype=float)


# ============================================================
# RECONSTRUCT SYMMETRIC DISTANCE MATRIX
# ============================================================

upper_triangle = np.triu(raw_matrix, k=1)

distance_matrix = upper_triangle + upper_triangle.T

np.fill_diagonal(distance_matrix, 0)


# ============================================================
# CHECK DISTANCE MATRIX
# ============================================================

print("\nDistance matrix information:")
print("--------------------------------")

print("Number of sequences:", len(sequence_names))
print("Matrix shape:", distance_matrix.shape)

print(
    "Symmetric:",
    np.allclose(distance_matrix, distance_matrix.T)
)

print(
    "Diagonal zero:",
    np.allclose(np.diag(distance_matrix), 0)
)

print(
    "Minimum distance:",
    np.min(distance_matrix)
)

print(
    "Maximum distance:",
    np.max(distance_matrix)
)


# ============================================================
# PCoA / CLASSICAL MDS
# ============================================================

print("\nPerforming classical PCoA / MDS...")

n = distance_matrix.shape[0]

# Squared distance matrix
D_squared = distance_matrix ** 2

# Centering matrix
J = np.eye(n) - np.ones((n, n)) / n

# Double-centered matrix
B = -0.5 * J @ D_squared @ J


# ============================================================
# EIGENVALUE DECOMPOSITION
# ============================================================

eigenvalues, eigenvectors = np.linalg.eigh(B)

# Sort from largest to smallest
order = np.argsort(eigenvalues)[::-1]

eigenvalues = eigenvalues[order]
eigenvectors = eigenvectors[:, order]


# ============================================================
# POSITIVE EIGENVALUES
# ============================================================

positive_mask = eigenvalues > 1e-10

positive_eigenvalues = eigenvalues[positive_mask]

positive_eigenvectors = eigenvectors[:, positive_mask]

# Automatically use all positive PCoA axes
number_of_pcoa_axes = len(positive_eigenvalues)

print("\nEigenvalue information:")
print("--------------------------------")

print(
    "Positive eigenvalues:",
    len(positive_eigenvalues)
)

print(
    "Negative eigenvalues:",
    np.sum(eigenvalues < -1e-10)
)

print(
    "Approximately zero:",
    np.sum(np.abs(eigenvalues) <= 1e-10)
)


print("\nEigenvalues:")

for i, value in enumerate(eigenvalues, start=1):

    print(
        f"PCoA{i:2d}: {value:.10f}"
    )


# ============================================================
# CHECK REQUESTED NUMBER OF PCoA AXES
# ============================================================

if number_of_pcoa_axes > len(positive_eigenvalues):

    raise ValueError(
        f"You requested {number_of_pcoa_axes} PCoA axes, "
        f"but only {len(positive_eigenvalues)} positive "
        f"eigenvalues are available."
    )


# ============================================================
# CALCULATE PCoA COORDINATES
# ============================================================

coordinates = (
    positive_eigenvectors[:, :number_of_pcoa_axes]
    *
    np.sqrt(
        positive_eigenvalues[:number_of_pcoa_axes]
    )
)


# ============================================================
# SAVE PCoA COORDINATES
# ============================================================

coordinate_columns = [
    f"PCoA{i + 1}"
    for i in range(number_of_pcoa_axes)
]

coordinates_df = pd.DataFrame(
    coordinates,
    columns=coordinate_columns
)

coordinates_df.insert(
    0,
    "Sequence_ID",
    sequence_names
)


coordinates_file = os.path.join(
    output_directory,
    "KMeans_PCoA_coordinates.csv"
)

coordinates_df.to_csv(
    coordinates_file,
    index=False
)

print("\nPCoA coordinates saved to:")
print(coordinates_file)


# ============================================================
# K-MEANS CLUSTERING
# ============================================================

print("\n========================================")
print("K-MEANS CLUSTERING")
print("========================================")

results = []


for k in range(
    min_clusters,
    max_clusters + 1
):

    print(
        f"\nPerforming K-means for k = {k}..."
    )

    kmeans = KMeans(
        n_clusters=k,
        random_state=random_state,
        n_init=n_init
    )

    # K-means is performed on PCoA coordinates
    kmeans_labels = kmeans.fit_predict(
        coordinates
    )

    # Convert labels from 0-based to 1-based
    # for easier interpretation
    cluster_labels = kmeans_labels + 1


    # ========================================================
    # SILHOUETTE USING NEW DISTANCE MATRIX
    # ========================================================

    silhouette = silhouette_score(
        coordinates,
        cluster_labels,
        metric="euclidean"
    )


    # ========================================================
    # WITHIN-CLUSTER SUM OF SQUARES
    # ========================================================

    inertia = kmeans.inertia_


    results.append(
        {
            "Number_of_Clusters": k,
            "Silhouette_Score": silhouette,
            "KMeans_Inertia": inertia
        }
    )


    print(
        f"k = {k:2d}   "
        f"Silhouette = {silhouette:.6f}   "
        f"Inertia = {inertia:.6f}"
    )


# ============================================================
# SAVE K-MEANS RESULTS
# ============================================================

results_df = pd.DataFrame(results)


results_file = os.path.join(
    output_directory,
    "KMeans_silhouette_scores.csv"
)

results_df.to_csv(
    results_file,
    index=False
)


print("\nK-means results saved to:")
print(results_file)


# ============================================================
# BEST K ACCORDING TO SILHOUETTE
# ============================================================

best_row = results_df.loc[
    results_df["Silhouette_Score"].idxmax()
]

best_k = int(
    best_row["Number_of_Clusters"]
)

best_score = float(
    best_row["Silhouette_Score"]
)


print("\n========================================")
print("BEST K ACCORDING TO SILHOUETTE")
print("========================================")

print(
    "Best k:",
    best_k
)

print(
    "Best silhouette score:",
    best_score
)


# ============================================================
# CURRENT K SILHOUETTE
# ============================================================

print(
    f"\nSilhouette score for your current "
    f"k = {current_cluster_number}:"
)


current_score = results_df.loc[
    results_df["Number_of_Clusters"]
    == current_cluster_number,
    "Silhouette_Score"
]


if not current_score.empty:

    print(
        f"{float(current_score.iloc[0]):.6f}"
    )


# ============================================================
# GET FINAL K-MEANS CLUSTERS FOR CURRENT K
# ============================================================

print(
    f"\nGenerating cluster assignments "
    f"for current k = {current_cluster_number}..."
)


final_kmeans = KMeans(
    n_clusters=current_cluster_number,
    random_state=random_state,
    n_init=n_init
)


final_labels = final_kmeans.fit_predict(
    coordinates
)

final_labels = final_labels + 1


## ============================================================
# SAVE CLUSTER ASSIGNMENTS
# ============================================================

cluster_assignment_df = pd.DataFrame({
    "Sequence_ID": sequence_names,
    "KMeans_Cluster": final_labels
})


cluster_file = os.path.join(
    output_directory,
    f"KMeans_cluster_assignments_k{current_cluster_number}.csv"
)

cluster_assignment_df.to_csv(
    cluster_file,
    index=False
)


print("\nCluster assignments saved to:")
print(cluster_file)


# ============================================================
# PRINT CLUSTER SIZES
# ============================================================

print("\nCluster sizes:")

cluster_counts = (
    pd.Series(final_labels)
    .value_counts()
    .sort_index()
)


for cluster, count in cluster_counts.items():

    print(
        f"Cluster {cluster}: {count} sequences"
    )


# ============================================================
# SILHOUETTE PLOT
# ============================================================

plt.figure(
    figsize=(10, 7),
    dpi=150
)


plt.plot(
    results_df["Number_of_Clusters"],
    results_df["Silhouette_Score"],
    marker="o",
    linewidth=2
)


# Best k
plt.scatter(
    best_k,
    best_score,
    s=100,
    label=f"Best k = {best_k}"
)


# Current k
if (
    current_cluster_number
    in results_df["Number_of_Clusters"].values
):

    current_score_value = float(
        results_df.loc[
            results_df["Number_of_Clusters"]
            == current_cluster_number,
            "Silhouette_Score"
        ].iloc[0]
    )

    plt.scatter(
        current_cluster_number,
        current_score_value,
        s=100,
        label=(
            f"Current k = "
            f"{current_cluster_number}"
        )
    )


plt.xlabel(
    "Number of Clusters (k)"
)

plt.ylabel(
    "Mean Silhouette Score"
)

plt.title(
    "K-means Clustering after PCoA"
)

plt.xticks(
    range(
        min_clusters,
        max_clusters + 1
    )
)

plt.grid(
    True,
    alpha=0.3
)

plt.legend()

plt.tight_layout()


plot_file = os.path.join(
    output_directory,
    "KMeans_silhouette_score_vs_k.png"
)


plt.savefig(plot_file)

plt.close()


print("\nSilhouette plot saved to:")
print(plot_file)


# ============================================================
# 2D PCoA PLOT WITH K-MEANS CLUSTERS
# ============================================================

if number_of_pcoa_axes >= 2:

    plt.figure(
        figsize=(8, 7),
        dpi=150
    )


    scatter = plt.scatter(
        coordinates[:, 0],
        coordinates[:, 1],
        c=final_labels,
        s=70,
        cmap="tab10"
    )


    for i, name in enumerate(sequence_names):

        plt.annotate(
            name,
            (
                coordinates[i, 0],
                coordinates[i, 1]
            ),
            xytext=(5, 5),
            textcoords="offset points",
            fontsize=8
        )


    plt.xlabel("PCoA1")

    plt.ylabel("PCoA2")

    plt.title(
        f"K-means Clusters on PCoA Space "
        f"(k = {current_cluster_number})"
    )


    plt.axhline(
        0,
        linewidth=0.8
    )

    plt.axvline(
        0,
        linewidth=0.8
    )


    plt.tight_layout()


    pcoa_cluster_plot = os.path.join(
        output_directory,
        f"KMeans_PCoA_clusters_k{current_cluster_number}.png"
    )


    plt.savefig(
        pcoa_cluster_plot
    )

    plt.close()


    print(
        "\nPCoA K-means cluster plot saved to:"
    )

    print(
        pcoa_cluster_plot
    )


# ============================================================
# COMPLETED
# ============================================================

print("\n========================================")
print("K-MEANS ANALYSIS COMPLETED")
print("========================================")