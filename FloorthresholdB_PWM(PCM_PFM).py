import os
import pandas as pd


input_directory =  r"\90aln_effectivesize"
pcm_output_directory = r"\90aln_pcm_outputs"
pfm_output_directory = r"\90aln_pfm_outputs"

os.makedirs(pcm_output_directory, exist_ok=True)
os.makedirs(pfm_output_directory, exist_ok=True)

for filename in os.listdir(input_directory):
    if filename.endswith(".csv"):
        input_path = os.path.join(input_directory, filename)
        print(f"Processing {filename}...")

        # Load the data
        df = pd.read_csv(input_path)

        # Remove first column (sequence IDs)
        data = df.iloc[:, 1:].copy()
        data = data.apply(pd.to_numeric, errors='coerce')

        # Remove last 3 rows if they are summary stats
        data_no_footer = data.copy()
        if any(df.iloc[-1, 0].strip().lower().startswith(stat) for stat in ['sum', 'sum/count', 'std']):
            data_no_footer = data_no_footer.iloc[:-4, :]

        # Total number of sequences (used for PFM)
        total_sequences = data_no_footer.shape[0]
        print(f"Total sequences (rows): {total_sequences}")

        # Determine max value for PCM row count
        max_value = int(data_no_footer.max().max())
        print(max_value)

        # Create PCM
        pcm_data = []
        for i in range(1, max_value + 1):
            row_counts = (data_no_footer == i).sum(axis=0)
            pcm_data.append(row_counts.values)

        pcm_df = pd.DataFrame(pcm_data, columns=data.columns)
        pcm_df.insert(0, 'Value', list(range(1, max_value + 1)))

        # Save PCM
        pcm_path = os.path.join(pcm_output_directory, filename)
        pcm_df.to_csv(pcm_path, index=False)
        print(f"Saved PCM: {pcm_path}")

        # Create PFM (normalize by number of sequences)
        pfm_df = pcm_df.copy()
        pfm_df.iloc[:, 1:] = pfm_df.iloc[:, 1:] / total_sequences

        # Apply minimum threshold: set all values < 0.000000000001 to 0.000000000001
        pfm_df.iloc[:, 1:] = pfm_df.iloc[:, 1:].clip(lower=0.000000000001)

        # Save PFM
        pfm_path = os.path.join(pfm_output_directory, filename)
        pfm_df.to_csv(pfm_path, index=False)
        print(f"Saved PFM: {pfm_path}")
