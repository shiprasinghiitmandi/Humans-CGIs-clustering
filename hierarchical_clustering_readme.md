# Hierarchical Clustering and Orientation-Aware Sequence Alignment Pipeline

This repository contains a modular Python pipeline designed for **distance matrix calculation, chromosome-wide matrix aggregation, threshold-based hierarchical clustering, and orientation-aware sequence alignment**.

The workflow uses Pearson correlation distance metrics ($1 - r$) across sliding shifts—checking both forward and reverse-flipped orientations—to robustly group and align sequential profile data (e.g., genomic signal distributions, epigenetic tracks, or windowed biological vectors).

## Requirements \& Dependencies

The pipeline requires **Python 3.8+** with the following scientific libraries:

* `numpy`
* `pandas`
* `scipy`
* `matplotlib`
* `openpyxl` (for `.xlsx` input handling)

Install all dependencies via `pip`:

```
pip install numpy pandas scipy matplotlib openpyxl

```

##  Execution Pipeline \& Script Descriptions

### Step 1: Distance Matrix \& Local Clustering Calculation

* **Script:** `flipped\_comp\_matrix.py`
* **Description:** Iterates through individual chromosome profile files (Excel `.xlsx`), calculates the optimal Pearson correlation ($r$) between all pair combinations across sliding shifts and forward/reverse orientations, converts correlation into distance ($1 - r$), generates full symmetric distance matrices, and exports individual dendrogram plots.
* **Input File:**

  * Per-chromosome Excel file (`.xlsx`) where **Row index = Sequence/Locus ID** and **Columns = Ordered numerical vector values**.
  * *Example input file included in repository:* `example\_data/chr2\_50th.xlsx`
* **Output Files:**

  * Full distance matrix CSV (`<basename>\_full\_distance\_matrix.csv`)
  * Linkage matrix CSV (`<basename>\_linkage\_matrix.csv`)
  * Dendrogram plot (`<basename>\_dendrogram.png`)

### Step 2: Chromosome-Wide Distance Matrix Aggregation

* **Script:** `Chr1tillY\_distance\_combined.py`
* **Description:** Reads all individual per-chromosome distance matrices from a designated directory (`flipped\_90dist/`), orders them sequentially (e.g., Chr1 through ChrY), and combines them into a single genome-wide distance matrix while removing redundant pairwise entries.
* **Input Directory \& File:**

  * Directory containing per-chromosome distance CSV files (`flipped\_50dist/`).
  * *Example input files included in repository:* `example\_data/distance\_matrices/chr1\_chr2\_dist.csv`
* **Output File:**

  * Combined upper-triangular distance matrix CSV: `1tilly\_50th\_flipped.csv`

### Step 3: Global Hierarchical Clustering \& Cut-Height Segmentation

* **Script:** `half\_dist\_cluster\_cutting.py`
* **Description:** Takes the combined genome-wide distance matrix, performs complete-linkage hierarchical clustering, cuts the resulting dendrogram at a user-defined distance height threshold ($t$), and groups sequences into discrete cluster IDs.
* **Input File:**

  * Combined upper-triangular distance matrix CSV: `1tilly\_50th\_flipped.csv`
  * *Example input file included in repository:* `example\_data/1tilly\_50th\_flipped.csv`
* **Parameters to Set:**

  * `height\_threshold`: Distance threshold value for cutting clusters (e.g., `1.11`).
* **Output Files:**

  * Text file containing cluster assignments and sequence IDs (`fl\_50thperc1.11.txt`).
  * Genome-wide dendrogram plot with cut line visualization (`fl\_50thperc1.11\_dendrogram.png`).

### Step 4: Cluster-Specific Pairwise Alignment Parameters

* **Script:** `flipped\_align2.py`
* **Description:** Computes parallelized, orientation-aware pairwise shift offsets and Pearson correlation scores for all sequence pairs within an isolated cluster. Tracks whether a sequence needs to be reverse-flipped (`Flipped = 1`) to achieve optimal correlation.
* **Input File:**

  * CSV file containing sequence vectors for a specific cluster (`Cluster 2\_50th.csv`), where column 1 is `SEQUENCE\_ID` and remaining columns represent numerical profile values.
  * *Example input file included in repository:* `example\_data/Cluster2\_50th.csv`
* **Output File:**

  * Pairwise alignment parameter table CSV (`Cluster 1\_50thp.csv`) containing headers: `\[Short\_Sequence, Long\_Sequence, Best\_Shift, Best\_Correlation, Flipped]`.

### Step 5: Iterative Multiple Sequence Alignment Assembly

* **Script:** `flipped\_align3.py`
* **Description:** Takes the raw cluster sequence vectors and the pairwise alignment parameter output from Step 4. Iteratively positions, shifts, and flips individual sequences relative to a reference coordinate frame to build a unified Multiple Sequence Alignment (MSA) matrix padded with zeros (`0`).
* **Input Files:**

  * Raw cluster profile CSV (`Cluster 2\_50th.csv`)
  * Alignment parameter CSV from Step 4 (`Cluster 2\_50thp.csv`)
* **Output File:**

  * Aligned sequence matrix CSV (`Cluster 1\_50thAlign.csv`) where all sequences share an aligned spatial index.

##  

## 📄 License \& Citation

If you use this pipeline or adapt these code scripts in your research, please cite this repository and associated manuscript:

> \*\*\\\[S. Singh, K. Hungyo\\]\*\*, \*"\*\*\*A genome-wide sequence analysis of human CpG islands reveals distinct patterns of CpG organization shaping DNA methylation and gene regulatory networks\*\*\*"\*, \\\[Journal Name / Year\\]. DOI: `\[`\*\*10.5281/zenodo.22929492\*\*`]`

