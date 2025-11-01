from itertools import combinations
from scipy.stats import pearsonr
from concurrent.futures import ProcessPoolExecutor
import pandas as pd
import csv

# Function to calculate the best Pearson correlation for a single pair (with flip check)
def compute_best_pearson_for_pair(pair, sequences):
    seq1_id, seq2_id = pair
    seq1 = sequences[seq1_id]
    seq2 = sequences[seq2_id]

    if len(seq1) >= len(seq2):
        long_seq, short_seq = seq1, seq2
        long_id, short_id = seq1_id, seq2_id
    else:
        long_seq, short_seq = seq2, seq1
        long_id, short_id = seq2_id, seq1_id

    best_pearson = -1
    best_shift = 0
    flipped_flag = 0

    # Original short sequence
    for shift in range(0, len(long_seq) - len(short_seq) + 1, 2):
        window = long_seq[shift:shift + len(short_seq)]
        if len(window) == len(short_seq):
            try:
                correlation, _ = pearsonr(window, short_seq)
                if correlation > best_pearson:
                    best_pearson = correlation
                    best_shift = shift
                    flipped_flag = 0
            except Exception:
                continue

    # Flipped short sequence
    flipped_short = short_seq[::-1]
    for shift in range(0, len(long_seq) - len(short_seq) + 1, 2):
        window = long_seq[shift:shift + len(short_seq)]
        if len(window) == len(short_seq):
            try:
                correlation, _ = pearsonr(window, flipped_short)
                if correlation > best_pearson:
                    best_pearson = correlation
                    best_shift = shift
                    flipped_flag = 1
            except Exception:
                continue

    if best_pearson > 0:
        return (short_id, long_id, best_shift, best_pearson, flipped_flag)
    else:
        return None

# Function to calculate all best Pearson results in parallel
def calculate_all_pearsons(sequences, max_workers=4):
    pairs = list(combinations(sequences.keys(), 2))
    results = []

    with ProcessPoolExecutor(max_workers=max_workers) as executor:
        futures = [executor.submit(compute_best_pearson_for_pair, pair, sequences) for pair in pairs]
        for future in futures:
            result = future.result()
            if result:
                results.append(result)

    return results

# Function to read input sequences from CSV
def read_sequences(csv_path):
    df = pd.read_csv(csv_path)
    df.set_index(df.columns[0], inplace=True)
    sequences = {}
    for index, row in df.iterrows():
        sequences[str(index)] = row.dropna().tolist()
    return sequences

# Function to save Pearson results
def save_pearson_results(results, output_csv):
    with open(output_csv, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(["Short_Sequence", "Long_Sequence", "Best_Shift", "Best_Correlation", "Flipped"])
        for row in results:
            writer.writerow(row)

# Main function
def run_pearson_calculation(input_csv, output_csv):
    print(f"Reading sequences from {input_csv}")
    sequences = read_sequences(input_csv)

    print("Calculating best Pearson correlations...")
    results = calculate_all_pearsons(sequences)

    print(f"Saving results to {output_csv}")
    save_pearson_results(results, output_csv)


if __name__ == "__main__":
    input_csv = r"/scratch/shipras.sbb.iitmandi/flipped_1.11/90/Cluster 3_90th.csv"
    output_csv = r"/scratch/shipras.sbb.iitmandi/flipped_1.11/90/Cluster 3_90thp.csv"
    run_pearson_calculation(input_csv, output_csv)
