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
<img width="1118" height="511" alt="image" src="https://github.com/user-attachments/assets/d03df890-63e9-4326-a6ee-abadd85c735a" />



# Obtaining cluster-wise methylation patterns
Step1: Filtering out only 100 methylation fraction containing sites from bedgraph file of each tissue using the 1_methylation_data_mining.py code and preparing a master file.

Step2: Extracting each methylated chromosome individually from the master file using 2_extracting_methylated_chr.py code.

Step3: Comparing the obtained methylated positions with CGI coordinates to obtain methylated positions only within CGIs using 3_CGI_extracted_methylated_pos.py code.

Step4: Preparing a CGI docx with highlighted methylated positions for better visibility and further analyses using 4_highlighting_CGIdocx.py

Step5: Renaming the headers of the highlighted docx so as to assign the same CGI sequence Id as we used during the clustering analyses using 5_renaming_methhltdocx.py code.

Step6: Recording absolute chromosome-wise dy dx methylation levels using 6_methylated_dy_dx_table.py code.

Step7: Recording the flipped or unflipped versions for every aligned sequences within each cluster using 11_alignedmeth.py which is further used during methylation mapping 
process.

Step8: Mapping the methylated sites with absolute methylation levels recorded on different CGIs over the aligned CGIS of each cluster using the code 12_alignedmeth.py, considering the flipped status recorded earlier.

Step9: Calculating column-wise methyaltion levels in each cluster using 13_alignedmethplot.py code.

Step10: Calculation of cluster-wise methylation levels for each tissue using 14_methylationlevels_histogram.py code.

Step11: Random methylation simulation code: 15_random methylation simulation.py

