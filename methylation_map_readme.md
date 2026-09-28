# Methylation Map Formation, Simulation, and Spatial Profiling Pipeline

This repository contains an end-to-end Python pipeline for **genomic CpG methylation data mining, CGI-sequence annotation, orientation-aware spatial mapping, stochastic baseline simulation, and multi-tissue profiling**.

The pipeline bridges raw bedGraph methylation scores with structural cluster alignments to generate spatial methylation maps across CpG-rich ($dy$) and non-CpG ($dx$) sequence domains.

##  Environment Requirements & Dependencies

* **Python:** 3.8 or higher

* **Required Libraries:**

  * `numpy`

  * `pandas`

  * `matplotlib`

  * `python-docx` (for Word document sequence markup handling)

  * `openpyxl` (for Excel sheet read/write operations)

Install all python dependencies via:

```
pip install numpy pandas matplotlib python-docx openpyxl

```

## Execution Order & Script Details

### Module 1: Genomic Data Mining & CGI Filtering

#### 1. `1_ file-methylation data mining.py`

* **Purpose:** Filters raw whole-genome bedGraph methylation files to retain sites with 100% methylation fraction.

* **Input:** `CpG_methylation_levels_Lung.bedgraph`

* **Output:** `methylated chr all_Lung.txt`

#### 2. `2_ extracting each methylated chromosome.py`

* **Purpose:** Extracts chromosome-specific methylation records (e.g., `chrY` or `chr20`) from the filtered genome-wide file.

* **Input:** `methylated chr all.txt`

* **Output:** `methylated chr Y.txt`

#### 3. `3_CGI extracted methylated positions.py`

* **Purpose:** Intersects chromosome-specific methylation coordinates with genomic CpG Island (CGI) coordinate boundaries.

* **Input:**

  * CGI coordinate range file: `chrY.txt` (BED format: `chrom`, `start`, `end`)

  * Chromosome methylated positions file: `methylated chr Y.txt`

* **Output:** Filtered CGI-specific bedGraph positions: `chrY_methcgi_strt_end.txt`

### Module 2: Sequence Annotation & Feature Extraction

#### 4. `4_a_highlightingCGIdocx.py`

* **Purpose:** Maps single-base methylation coordinates onto raw FASTA sequence strings and exports Microsoft Word (`.docx`) documents with yellow highlight formatting applied to methylated positions.

* **Input:**

  * Sequence text file: `chr20_replaced_ATCG.txt`

  * Methylation coordinates Excel file: `chr20_cg_meth.xlsx`

* **Output:** Formatted DOCX file: `chr20highlighted_sequences.docx`

#### 5. `5_renaming_methhltdocx.py`

* **Purpose:** Standardizes sequence headers within highlighted `.docx` documents from genomic coordinates to unified global `sequence <ID>` identifiers based on chromosome range lookup tables.

* **Input:**

  * Highlighted sequence document: `chr20highlighted_sequences.docx`

  * Lookup ID table: `sequence_id_ranges.xlsx`

* **Output:** `chr20highlighted_sequences_renamed.docx`

#### 6. `6_methylated_dx_dy_table_9april.py`

* **Purpose:** Parses highlighted Word documents to quantify methylated site counts across CpG dinucleotide runs ($dy$) and non-CpG/spacer segments ($dx$), outputting structural feature tables per sequence.

* **Input Directory:** Directory containing renamed highlighted `.docx` files (`cancerlung_renamed_docx/`)

* **Output Directory:** Excel directory containing absolute $dy$/$dx$ counts (`cancerlung_abs_dydx_excel/*_dydx.xlsx`)

### Module 3: Alignment Matrix & Spatial Mapping

#### 7. `7_alignedmeth.py` (Orientation Detection)

* **Purpose:** Compares unaligned raw cluster sequences against Multiple Sequence Alignment (MSA) matrices to determine whether individual sequences were maintained in forward orientation (`flipped_status = 0`) or reverse-flipped (`flipped_status = 1`).

* **Input Files:**

  * Unaligned sequence CSV: `Cluster 3_90th.csv`

  * Aligned sequence CSV: `Cluster 3_90thAlign.csv`

* **Output File:** `Cluster 3_90th_flip_status.csv`

#### 8. `8_alignedmeth.py` (Methylation Map Construction)

* **Purpose:** Projects $dy$ and $dx$ methylation counts onto the spatial grid of aligned sequence matrices, respecting orientation adjustments (reversing features for `flipped_status = 1`).

* **Input Files:**

  * Folder containing $dy$/$dx$ Excel tables (`lung_abs_dydx_excel/`)

  * Alignment CSV with flip status: `Cluster 3_90th_flip_status.csv`

* **Output File:** `Cluster 3_90th_methmap.csv`

### Module 4: Null Baseline Stochastic Simulation

#### `random_methylation_simulation.py`

* **Purpose:** Simulates a null baseline methylation map via binomial sampling. For each $dy$ feature containing $n$ nucleotides (equivalent to $n/2$ CpG sites), methylation events are sampled from $B(n/2, p)$, where $p$ is the empirical global methylation probability. Results are averaged over $M$ single-cell replicates.

* **Key Parameters:**

  * `p_methylation`:  methylation probability fixed according to experimental value obtained for the same tissue cluster under the specific pathological condition(the real methylation normalised sum value) $p$ (e.g., `0.046`)

  * `num_replicates`: Simulated cell replicates $M$ (e.g., `100`)

* **Input File:** Aligned sequence matrix `Cluster 5_90thAlign.csv`

* **Output File:** Simulated average methylation map `Cluster 5_90th_methylation_r.csv`

### Module 5: Spatial Profiling & Cross-Tissue Visualization

#### 9. `9_alignedmethplot.py`

* **Purpose:** Computes index-wise total sums and per-sequence normalized averages for $dy$ and $dx$ positions, exports profile summary spreadsheets, and generates spatial profile line plots.

* **Input File:** Experimental or simulated methylation map CSV (`Cluster 1_90th_methylation_r3.csv`)

* **Output Files:**

  * Summary matrix CSV: `Cl1_random_output_with_sums.csv`

  * Spatial profile plot PNG: `Cl1_random_output_with_sums.png`

#### 10. `10_methylationlevels_histogram.py`

* **Purpose:** Aggregates total $dy$ and $dx$ normalized methylation levels across multiple tissues and clusters, generating comparative horizontal stacked bar charts and consolidated CSV summaries.

* **Input Directory:** Root folder containing tissue-specific methylation map folders (`methmaps_tissues_folder/`)

* **Output Files:**

  * Multi-tissue comparison plot: `all_tissues_normalisedsum_methlevels_correct_r.png`

  * Consolidated summary table: `r_dy_dx_methylationsummary.csv`

## Required Repository Input Files & Example Structure

To ensure complete reproducibility, upload the following example data files into an `example_data/` folder in your repository:

| Input File Name / Type | Description | Relevant Pipeline Step | 
 | ----- | ----- | ----- | 
| `CpG_methylation_levels_sample.bedgraph` | Raw bedGraph containing chromosome coordinates and methylation scores. | Script 1 | 
| `chrY_CGI_ranges.txt` | Tab-delimited file containing CpG Island genomic coordinates (`chr`, `start`, `end`). | Script 3 | 
| `chr_sequences.txt` | Text file containing chromosome header tags and nucleotide sequence strings. | Script 4 | 
| `sequence_id_ranges.xlsx` | Mapping table connecting chromosome identifiers to sequence ID ranges. | Script 5 | 
| `Cluster_90th.csv` | Raw (unaligned) cluster sequence matrix. | Script 7 | 
| `Cluster_90thAlign.csv` | Aligned cluster sequence matrix produced by MSA. | Script 7 & Simulation | 
| `methmaps_tissues_folder/` | Directory containing tissue subfolders (e.g., `Lung_methmap_files/`) with cluster methylation maps. | Script 10 | 

## 📄 License & Citation

When using this pipeline or its output maps in academic publications, please cite:

> **\[S. Singh, K. Hungyo\]**, *"A genome-wide sequence analysis of human CpG islands reveals distinct patterns of CpG organization shaping DNA methylation and gene regulatory networks"*, \[Journal Name / Year\]. DOI: `[10.5281/zenodo.22929492]`