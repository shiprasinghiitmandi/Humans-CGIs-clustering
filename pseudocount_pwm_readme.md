# Pseudocount Position Weight Matrix (PWM) Construction, Sequence Scanning, and Score Calculation Pipeline

This repository provides an end-to-end computational pipeline for constructing **Pseudocount-adjusted Position Weight Matrices (PWMs)** from structural alignments of CpG Islands (CGIs), scanning these PWM models across target sequence sets (**real-included, real-excluded, and synthetic random CGIs**), and determining optimal alignment scores under phase-preserving orientation search.

The framework processes structural sequence matrices represented as alternating CpG dinucleotide runs ($dy$) and non-CpG spacer lengths ($dx$).

##  Environment Requirements & Dependencies

* **Python:** 3.8 or higher

* **Dependencies:**

  * `numpy`

  * `pandas`

  * `openpyxl` (for reading and writing Excel spreadsheets)

Install all required Python packages:

```
pip install numpy pandas openpyxl

```

##  Mathematical Framework

### 1. Pseudocount Position Frequency Matrix (PFM)

To handle small sample sizes and prevent zero-probability artifacts without arbitrary pseudocount constants, background probabilities $f_0(i)$ from genome-wide distributions are incorporated as prior counts:

$$
f_{ij}^o=PFM_{ij}=\frac{\left(PCM+f_i^e\right)}{n+1}
$$

Where:

* $PCM_{i,j}$ is the observed frequency count of feature value $i$ at column position $j$.

* fe is the normalized background expectation for value $i$ in the corresponding domain ($dx$ or $dy$).

* $N$ is the total number of aligned sequences.

### 2. Log-Odds Position Weight Matrix (PWM)

Weights $PWM_{i,j}$ are computed as log-odds scores relative to background frequencies:

$$
PWM_{ij}\ =\log \left(\frac{f_{ij}^o}{f_i^e}\right)
$$

*Note:* For $dy$ columns, odd-valued lengths (e.g., 1, 3, 5) are physically unallowable due to dinucleotide symmetry and are assigned `NaN`.

### 3. Phase-Preserving Sequence Scanning

When scanning target sequences, a sliding window steps by **2 positions at a time** ($\Delta s = 2$) to enforce structural phase alignment between alternating $dy$ and $dx$ positions. Both original and reverse-flipped PWM configurations are evaluated to derive the maximum alignment score:

$$
\text{Score}(S_{\text{window}}) = \sum_{j=1}^{L_{\text{PWM}}} W_{S[j], j}
$$

## Execution Order & Workflow Script Details

```
Filter Effective Range
       │
       ▼
Compute PCM & Pseudocount PFM
       │
       ▼
Calculate Log-Odds PWM
       │
  ┌────┴──────────────────────────┐
  ▼                            ▼
PWM scanning           Target Sets
                 (Real-Included / Excluded/ Random CGIs)
  │                            │
  └────────────┬─────────────────┘
              ▼
Scan Sequences & Compute Best Scores

```

### Step 1: Align Feature Trimming & Range Selection

#### Script: `PseudocountA_PWMrevisedCBCa.py`

* **Purpose:** Trims raw sequence alignment matrices to predefined cluster-specific core structural boundaries ($dy_{\text{start}}$ to $dy_{\text{end}}$) while removing statistical summary rows.

* **Input Directory:** `/50aln_plots_input/`

* **Output Directory:** `/50aln_effectivesize/`

* **Key Operations:** Extracts designated $dy$ ranges per cluster (e.g., `dy25` to `dy154` for Cluster 1) and preserves the `SEQUENCE_ID` column.

### Step 2: Position Count & Pseudocount PFM Generation

#### Script: `PseudocountB_PWM(PCM_PFM)revisedCBC.py`

* **Purpose:** Counts feature values across matrix positions to construct Position Count Matrices (PCMs) and transforms counts into pseudocount-adjusted Position Frequency Matrices (PFMs).

* **Inputs:**

  * Filtered CSVs from `/50aln_effectivesize/`

  * Background frequency reference: `combined_cgiandnoncgi_normalized.csv`

* **Outputs:**

  * PCM files saved to `/50aln_pcm_outputs/`

  * PFM files saved to `/50aln_pfm_outputs/`

### Step 3: Log-Odds Weight Matrix Calculation

#### Script: `PseudocountC_PWM_finalrevisedCBC.py`

* **Purpose:** Converts PFM frequencies into log-odds Position Weight Matrices (PWMs).

* **Inputs:**

  * PFM matrices from `/50aln_pfm_outputs/`

  * Background reference file: `combined_cgiandnoncgi_normalized.csv`

* **Output Directory:** `/50aln_PWM_outputs/`

* **Key Handling:** Automatically handles invalid odd $dy$ entries and zero background probabilities.

### Step 4: Random CGI Generation

#### Script: `PWMscanning_random_cgis_generation.py`

* **Purpose:** Generates synthetic, randomly constructed CpG island sequences to serve as negative/random baseline controls during scanning evaluation.

* **Mechanism:** Samples alternating $dy$ (lengths 2, 4, 6, 8) and $dx$ (lengths 1–15) values from empirical probability distributions ($p_{dy}$ and $p_{dx}$) until target genomic length bounds (1500–3000 bp) are reached.

* **Output File:** `randomCGIs_matrix.csv`

### Step 5: Sequence Scanning & Score Computation

#### Script: `PseudocountD_PWM_scanning_sequencesRevisedCBC.py`

* **Purpose:** Scans target sequence sets (Real-Included CGIs, Real-Excluded CGIs, or Random CGIs) across forward and flipped PWM orientations using step size = 2 to retain phase.

* **Inputs:**

  * Target PWM file (e.g., `Cluster 1_90thAlign_PWM.xlsx`)

  * Candidate sequences matrix (e.g., `excludedseqs_matrix_padded2.xlsx` or `randomCGIs_matrix.csv`)

  * Expected background probabilities: `combined_cgiandnoncgi_normalized.csv`

* **Output:** CSV spreadsheet containing best match scores, optimal shift indices, and flip orientation status (`cluster1_score.csv`).

  ## 

### Sample Format for `combined_cgiandnoncgi_normalized.csv`:

```
Value,dx Normalized Frequency,dy Normalized Frequency
1,0.1234,0.0000
2,0.2150,0.7800
3,0.0890,0.0000
4,0.1120,0.1500

```

## Summary Output Format (`cluster1_score.csv`)

The final scanning step outputs a standardized table detailing optimal alignment parameters per sequence:

| Sequence_ID | PWM_Name | Best_PWM_Score | Best_Shift_Position | Flipped | 
 | ----- | ----- | ----- | ----- | ----- | 
| `CGI_1002` | `Cluster 1_90thAlign_PWM` | `42.8512` | `4` | `0` | 
| `CGI_1003` | `Cluster 1_90thAlign_PWM` | `38.1204` | `12` | `1` | 

* **`Best_Shift_Position`**: Starting column index of the highest-scoring alignment window.

* **`Flipped`**: `0` for forward orientation; `1` if the reverse-flipped PWM yielded a higher score.

## 📄 License & Citation

If you use this pseudocount PWM scanning workflow in your research, please cite:

> $$
> S. Singh, K. Hungyo
> $$
>
> , *"A genome-wide sequence analysis of human CpG islands reveals distinct patterns of CpG organization shaping DNA methylation and gene regulatory networks"*, 
>
> $$
> Journal Name / Year
> $$
>
> . DOI: `[10.5281/zenodo.22929492]`