# CpG clustering and methylation maps formation
# Introduction
CpG island (CGI) sequences from each chromosome were converted into strings of numerical dy–dx values. All possible permutations and combinations of these CGI sequences were then compared to evaluate their similarity.
Since each CGI was represented by only two variables (dy and dx), we computed Pearson correlation coefficients for every possible pair of CGI sequences, establishing a measure of linear correlation between them.

Pairwise correlation matrices were generated for each chromosome pair, followed by distance matrix computation, linkage matrix generation, and hierarchical clustering.
Because processing all chromosomes simultaneously was memory-intensive, we instead analyzed all possible chromosome pairs individually. This resulted in 300 chromosome pairs (from 24 chromosomes) and 300 corresponding distance matrices, which were subsequently merged to produce a single comprehensive distance matrix representing a given percentile of CGI sequences.
# Pearson calculation, distance matrix formation and hierarchical clustering
Flipped_comp_matrix.py processes files containing CpG island (CGI) sequences from any given chromosome pair.
For each CGI pair, the script calculates the Pearson correlation coefficient by shifting the shorter sequence over the longer one in 2 bp increments (since the minimum dy value is 2 bp and each CGI begins with dy, ensuring that equivalent elements from both CGIs overlap at a time).

To ensure no optimal alignment is missed, the flipped orientation of each sequence is also evaluated, and Pearson correlations are computed for all possible flipped alignments as well.
Among all correlations calculated, the highest Pearson coefficient is considered the best score (representing the strongest similarity between two CGIs), and the corresponding position is recorded as the best shift.

This process is repeated for every possible CGI pair, after which:

1) A distance matrix is generated based on the best Pearson scores.

2) The linkage matrix is then derived from this distance matrix.

3) Agglomerative hierarchical clustering is performed, where each CGI initially forms an individual cluster, and clusters iteratively merge with their nearest neighbors until a single cluster remains.

The same workflow was applied to all 300 chromosome pairs, producing 300 distance matrices in total.
These were then combined using 1tillY_distance_combined.py to generate a single unified distance matrix covering all chromosomes (1 through Y).

Finally, cluster cutting was performed using half_dist_cluster_cutting.py, which used the upper half of this final matrix to form clusters at a defined cut threshold.
# Cluster-wise sequence alignment
flipped_align2.py is used to align all CpG islands (CGIs) within each cluster based on their best Pearson correlation scores and corresponding best shift positions.
The algorithm for computing the best Pearson correlation between all possible CGI sequence pairs follows the same approach as used previously during distance matrix generation.
However, in this step, the correlations are calculated only among CGIs belonging to the same cluster.

For each pair, the best Pearson score, shift position, and flip orientation (normal or reversed) are recorded to enable accurate alignment of sequences within the cluster.
flipped_align3.py uses the recorded alignment data from flipped_align2.py to perform the actual sequence alignment for all CGIs within each cluster.
Multiple alignment scenarios are handled in this script, each of which is explicitly addressed in the code.

The CGI pairs are sorted from highest to lowest best Pearson scores, ensuring that the most similar pairs are aligned first.
The best pair is aligned initially using the previously recorded flip orientation and shift information, followed by alignment of the remaining sequences according to their respective cases.

Example:
Different alignment cases account for specific scenarios or exceptions encountered during pairwise alignment (details available in the code comments).
<img width="1118" height="511" alt="image" src="https://github.com/user-attachments/assets/4f2c26bf-878c-416c-8dab-5f6fc589a3e0" />


# Obtaining cluster-wise methylation patterns

