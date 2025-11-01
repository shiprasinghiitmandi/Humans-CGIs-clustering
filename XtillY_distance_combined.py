import os
import pandas as pd

# Directory containing the distance matrix files
input_directory = r"/scratch/shipras.sbb.iitmandi/flipped_90dist/"
output_file = r"/scratch/shipras.sbb.iitmandi/1tilly_90th_flipped.csv"


def extract_chrom_number(filename):
    """
    We extracted chromosome numbers from filenames safely.

    """
    parts = filename.split('_')
    try:
        chrom1 = int(parts[0][3:]) if parts[0].startswith('chr') and parts[0][3:].isdigit() else float('inf')
        chrom2 = int(parts[1][3:]) if len(parts) > 1 and parts[1].startswith('chr') and parts[1][
                                                                                        3:].isdigit() else float('inf')
    except (IndexError, ValueError):
        return (float('inf'), float('inf'))
    return (chrom1, chrom2)


def combine_distance_matrices(input_directory, output_file):


    files = [f for f in os.listdir(input_directory) if f.endswith(".csv")]


    files.sort(key=extract_chrom_number)

    # Initializing an empty DataFrame for the combined matrix
    combined_matrix = pd.DataFrame()

    # Keeping track of which sequence pairs have already been written
    processed_sequences = set()

    # Processing each file in order
    for file in files:
        file_path = os.path.join(input_directory, file)
        print(f"Processing {file}...")  # Debugging print

        try:
            # Reading the distance matrix
            df = pd.read_csv(file_path, index_col=0)

            # Checking for empty files
            if df.empty:
                print(f"Skipping {file}, as it is empty.")
                continue


            df.columns = df.columns.astype(str)
            df.index = df.index.astype(str)


            for row in df.index:
                for col in df.columns:
                    seq_pair = tuple(sorted([row, col]))
                    if seq_pair not in processed_sequences:
                        combined_matrix.loc[row, col] = df.at[row, col]
                        processed_sequences.add(seq_pair)

        except Exception as e:
            print(f"Error processing {file}: {e}")


    combined_matrix = combined_matrix.fillna(0)

    # Saving the final combined matrix
    combined_matrix.to_csv(output_file)
    print(f"Combined distance matrix saved to {output_file}")



combine_distance_matrices(input_directory, output_file)
