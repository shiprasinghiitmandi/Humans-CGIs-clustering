import pandas as pd
import numpy as np

def clean_sequence_elements(row):
    """Extract and clean sequence elements from a DataFrame row (skip ID)."""
    # Convert to list of ints, skipping NaNs
    elements = [int(e) for e in row[1:] if not pd.isna(e)]
    return elements

def get_trimmed_elements(elements):
    """Trim leading and trailing zeros from elements."""
    # Remove leading/trailing zeros
    while elements and elements[0] == 0:
        elements.pop(0)
    while elements and elements[-1] == 0:
        elements.pop()
    return elements

def compare_sequences(file1, file2, output_file):
    df1 = pd.read_csv(file1)
    df2 = pd.read_csv(file2)

    # Create dictionary of sequences from file1
    seq1_dict = {}
    for _, row in df1.iterrows():
        seq_id = row.iloc[0].strip()
        elements = clean_sequence_elements(row)
        seq1_dict[seq_id] = elements

    flipped_status = []

    for _, row in df2.iterrows():
        seq_id = row.iloc[0].strip()
        elements = clean_sequence_elements(row)
        trimmed = get_trimmed_elements(elements.copy())

        if seq_id in seq1_dict:
            original = seq1_dict[seq_id]
            if original == trimmed:
                flipped_status.append(0)
            elif original == trimmed[::-1]:
                flipped_status.append(1)
            else:
                flipped_status.append(None)  # No match
        else:
            flipped_status.append(None)  # Sequence not in file1

    # Add flipped_status to second file
    df2['flipped_status'] = flipped_status

    # Save result
    df2.to_csv(output_file, index=False)



compare_sequences(r"F:\flipped_matrices\methmaps and fractions_liver\90\Cluster 3_90th.csv", r"F:\flipped_matrices\methmaps and fractions_liver\90\Cluster 3_90thAlign.csv", r"F:\flipped_matrices\methmaps and fractions_liver\90\Cluster 3_90th_flip_status.csv")
