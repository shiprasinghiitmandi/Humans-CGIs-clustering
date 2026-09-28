import os
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, fcluster

from sklearn.metrics import silhouette_score


input_distance_file = r"1tilly_90th_flipped.csv" # distance matrix of 90th percentile as input file

output_directory = r"silhouette_analysis90_HC" #hierarchical clustering silhouette scores output folder

min_clusters = 2
max_clusters = 20

# Your manuscript's current cluster number
current_cluster_number = 5 # cluster number selected for 90th percentile according to optimum height threshold found using previous method of cluster conservation mentioned in the main text


os.makedirs(output_directory, exist_ok=True)



# LOAD UPPER-TRIANGULAR DISTANCE MATRIX


print("Loading distance matrix...")

distance_df = pd.read_csv(
    input_distance_file,
    index_col=0
)

sequence_names = distance_df.columns.tolist()

raw_matrix = distance_df.to_numpy(dtype=float)



# RECONSTRUCT COMPLETE SYMMETRIC DISTANCE MATRIX


# Keep only the upper triangle
upper_triangle = np.triu(raw_matrix, k=1)

# Mirror the upper triangle
distance_matrix = upper_triangle + upper_triangle.T

# Diagonal must be zero
np.fill_diagonal(distance_matrix, 0)



# CHECK DISTANCE MATRIX


print("\nDistance matrix information:")
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



# CONVERT TO CONDENSED FORM


condensed_matrix = squareform(distance_matrix)



# COMPLETE-LINKAGE HIERARCHICAL CLUSTERING


print("\nPerforming complete-linkage clustering...")

linkage_matrix = linkage(
    condensed_matrix,
    method="complete"
)



# SILHOUETTE ANALYSIS


results = []

print("\nCalculating silhouette scores...\n")

for k in range(min_clusters, max_clusters + 1):

    # Cut dendrogram to obtain exactly k clusters
    cluster_labels = fcluster(
        linkage_matrix,
        t=k,
        criterion="maxclust"
    )

    # Calculate silhouette score using the ORIGINAL
    # pairwise distance matrix
    score = silhouette_score(
        distance_matrix,
        cluster_labels,
        metric="precomputed"
    )

    results.append(
        {
            "Number_of_Clusters": k,
            "Silhouette_Score": score
        }
    )

    print(
        f"Clusters = {k:2d}   "
        f"Silhouette Score = {score:.6f}"
    )



# SAVE RESULTS


results_df = pd.DataFrame(results)

results_file = os.path.join(
    output_directory,
    "silhouette_scores.csv"
)

results_df.to_csv(
    results_file,
    index=False
)

print("\nSilhouette results saved to:")
print(results_file)



# IDENTIFY BEST CLUSTER NUMBER


best_row = results_df.loc[
    results_df["Silhouette_Score"].idxmax()
]

best_k = int(best_row["Number_of_Clusters"])
best_score = float(best_row["Silhouette_Score"])

print("\n========================================")
print("BEST CLUSTER NUMBER")
print("========================================")
print("Best k:", best_k)
print("Best silhouette score:", best_score)

print(
    f"\nSilhouette score for your current "
    f"k = {current_cluster_number}:"
)

current_score = results_df.loc[
    results_df["Number_of_Clusters"] == current_cluster_number,
    "Silhouette_Score"
]

if not current_score.empty:
    print(
        f"{float(current_score.iloc[0]):.6f}"
    )



# PLOT SILHOUETTE SCORE VS NUMBER OF CLUSTERS


plt.figure(figsize=(10, 7), dpi=150)

plt.plot(
    results_df["Number_of_Clusters"],
    results_df["Silhouette_Score"],
    marker="o",
    linewidth=2
)

# Mark the best k
plt.scatter(
    best_k,
    best_score,
    s=100,
    label=f"Best k = {best_k}"
)

# Mark your current k = 5
if current_cluster_number in results_df["Number_of_Clusters"].values:

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
        label=f"Current k = {current_cluster_number}"
    )


plt.xlabel("Number of Clusters (k)")
plt.ylabel("Mean Silhouette Score")

plt.title(
    "Silhouette Analysis for Complete-Linkage Clustering"
)

plt.xticks(
    range(min_clusters, max_clusters + 1)
)

plt.grid(True, alpha=0.3)

plt.legend()

plt.tight_layout()




plot_file = os.path.join(
    output_directory,
    "silhouette_score_vs_cluster_number.png"
)

plt.savefig(plot_file)

plt.close()

print("\nPlot saved to:")
print(plot_file)

print("\nAnalysis completed successfully.")