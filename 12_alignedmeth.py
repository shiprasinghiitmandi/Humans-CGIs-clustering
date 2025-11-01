import pandas as pd
import os
import glob

def process_multiple_methylation_files(methylation_folder, aligned_csv, output_csv):
    # Read and combine all .xlsx methylation files
    methylation_files = glob.glob(os.path.join(methylation_folder, "*.xlsx"))
    all_methylation_data = []

    for file in methylation_files:
        df = pd.read_excel(file)
        all_methylation_data.append(df)

    combined_df = pd.concat(all_methylation_data, ignore_index=True)

    # Read aligned sequence file with flipped status
    df2 = pd.read_csv(aligned_csv)
    output_df = df2.copy()
    df2['SEQUENCE_ID'] = df2['SEQUENCE_ID'].astype(str)

    # Track which sequences were processed
    processed_sequences = set()

    # Group methylation data by SEQUENCE_ID
    grouped = combined_df.groupby('Sequence_ID')

    for seq_id, group in grouped:
        sequence_label = f"sequence {int(seq_id)}"
        seq_row = df2[df2['SEQUENCE_ID'] == sequence_label]

        if seq_row.empty:
            continue

        processed_sequences.add(sequence_label)

        flipped_status = seq_row.iloc[0, -1]
        sequence_values = seq_row.iloc[0, 1:-1].tolist()

        if flipped_status == 1:
            sequence_values = sequence_values[::-1]

        non_zero_indices = [i for i, val in enumerate(sequence_values) if val != 0]

        logical_map = {}
        labels = ['dy', 'dx']
        for i, idx in enumerate(non_zero_indices):
            label = labels[i % 2] + str((i // 2) + 1)
            logical_map[label] = idx

        updated_sequence = [0] * len(sequence_values)

        for _, row in group.iterrows():
            dy_pos = str(row['dy_position']) if pd.notna(row['dy_position']) else None
            dy_val = row['dy_values'] if pd.notna(row['dy_values']) else None
            dx_pos = str(row['dx_position']) if pd.notna(row['dx_position']) else None
            dx_val = row['dx_values'] if pd.notna(row['dx_values']) else None

            if dy_pos in logical_map and dy_val is not None:
                updated_sequence[logical_map[dy_pos]] = dy_val

            if dx_pos in logical_map and dx_val is not None:
                updated_sequence[logical_map[dx_pos]] = dx_val

        if flipped_status == 1:
            updated_sequence = updated_sequence[::-1]

        output_df.loc[output_df['SEQUENCE_ID'] == sequence_label, output_df.columns[1:-1]] = updated_sequence

    # Zero out sequences not present in methylation files
    all_sequence_labels = df2['SEQUENCE_ID'].tolist()
    missing_sequences = set(all_sequence_labels) - processed_sequences

    for missing_seq in missing_sequences:
        output_df.loc[output_df['SEQUENCE_ID'] == missing_seq, output_df.columns[1:-1]] = [0] * (output_df.shape[1] - 2)

    # Save final output
    output_df.to_csv(output_csv, index=False)



process_multiple_methylation_files(
    methylation_folder=r"F:\flipped_matrices\methmaps and fractions_lung\90\lung_abs_dydx_excel",
    aligned_csv= r"F:\flipped_matrices\methmaps and fractions_lung\90\Cluster 3_90th_flip_status.csv",
    output_csv=  r"F:\flipped_matrices\methmaps and fractions_lung\90\Cluster 3_90th_methmap.csv"
)



