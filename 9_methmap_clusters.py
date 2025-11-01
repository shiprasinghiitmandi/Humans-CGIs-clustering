import pandas as pd
import re
import os

def extract_numeric_id(text):
    """Extract numeric part from a string like 'sequence 3991'."""
    match = re.search(r'\d+', str(text))
    return int(match.group()) if match else None

def filter_from_directory(csv_file, xlsx_dir, output_file):
    # Read CSV and extract numeric IDs
    csv_df = pd.read_csv(csv_file)
    csv_df['Numeric_ID'] = csv_df['SEQUENCE_ID'].apply(extract_numeric_id)
    csv_seq_ids = set(csv_df['Numeric_ID'].dropna().astype(int))

    # Prepare a DataFrame to collect all matched rows
    all_matches = []

    # Loop over all .xlsx files in the directory
    for filename in os.listdir(xlsx_dir):
        if filename.endswith('.xlsx'):
            file_path = os.path.join(xlsx_dir, filename)
            try:
                xlsx_df = pd.read_excel(file_path)
                xlsx_df['Sequence_ID'] = xlsx_df['Sequence_ID'].fillna(0).astype(int)
                filtered_df = xlsx_df[xlsx_df['Sequence_ID'].isin(csv_seq_ids)]
                if not filtered_df.empty:
                    all_matches.append(filtered_df)
            except Exception as e:
                print(f"?? Error reading {filename}: {e}")

    # Concatenate all matched DataFrames
    if all_matches:
        result_df = pd.concat(all_matches, ignore_index=True)
        result_df.to_csv(output_file, index=False)
        print(f"? Combined results saved to: {output_file}")
    else:
        print("?? No matching sequences found in any .xlsx file.")

# Example usage:
csv_input = r"/scratch/shipras.sbb.iitmandi/flipped_1.11/90/Cluster 5_90th.csv"
xlsx_directory = r"/scratch/shipras.sbb.iitmandi/meth_norm_allchr/"
csv_output = r"/scratch/shipras.sbb.iitmandi/flipped_norm_methmap/90/cluster5_methmap.csv"

filter_from_directory(csv_input, xlsx_directory, csv_output)
