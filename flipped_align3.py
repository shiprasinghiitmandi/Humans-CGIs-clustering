import pandas as pd
import csv

# Function to flip a sequence
def flip_sequence(sequence):
    return sequence[::-1]

# Function to align sequences based on Pearson results
def align_sequences(sequences, best_pearson_results):
    sorted_pairs = sorted(best_pearson_results.items(), key=lambda x: -x[1][1])
    aligned_sequences = {}
    sequence_shifts = {}
    aligned_sequences_flipped = {}  # Tracking flipped status of aligned sequences
    placed_sequences = set()

    if not sorted_pairs:
        return {}

    # Process the first pair
    (short_id, long_id), (best_shift, _, flip) = sorted_pairs[0]
    
    # Handling flipping if necessary
    if flip == 1:  # Flipping short sequence if needed
        sequences[short_id] = flip_sequence(sequences[short_id])

    sequence_shifts[long_id] = 0
    sequence_shifts[short_id] = best_shift

    max_length = max(len(sequences[long_id]), len(sequences[short_id]) + best_shift)
    aligned_sequences[long_id] = [0] * (max_length - len(sequences[long_id])) + sequences[long_id]
    aligned_sequences[short_id] = [0] * best_shift + sequences[short_id] + [0] * (max_length - best_shift - len(sequences[short_id]))

    placed_sequences.update([short_id, long_id])
    aligned_sequences_flipped[short_id] = flip  # Tracking flip status
    aligned_sequences_flipped[long_id] = 0  # Long sequence is not flipped
    remaining_sequences = set(sequences.keys()) - placed_sequences

    while remaining_sequences:
        best_candidate = None
        best_match = None
        best_shift = 0
        best_value = -1
        best_flip = 0

        # Finding best candidate to align with already placed sequences
        for candidate in remaining_sequences:
            for aligned in placed_sequences:
                pair = (candidate, aligned) if (candidate, aligned) in best_pearson_results else (aligned, candidate)
                if pair in best_pearson_results:
                    shift, value, flip = best_pearson_results[pair]
                    if value > best_value:
                        best_candidate = candidate
                        best_match = aligned
                        best_shift = shift
                        best_value = value
                        best_flip = flip

        if best_candidate is None:
            break

        new_seq = sequences[best_candidate]
        ref_seq = sequences[best_match]
        ref_shift = sequence_shifts[best_match]
        new_len = len(new_seq)
        ref_len = len(ref_seq)
#mentioned all the cases
        # CASE 1: ONE SEQUENCE FROM THE PAIR IS ALREADY PRESENT IN UNFLIP FORM IN THE ALIGNED MATRIX

        # SITUATION 1: shorter seq already present in unflip form
        if aligned_sequences_flipped[best_match] == 0 and best_flip == 1: 
            new_seq = flip_sequence(new_seq)
            new_shift = len(ref_seq) - (len(new_seq) + best_shift)
            start_pos = ref_shift + new_shift
        # SITUATION 2: longer sequence already present in unflip form
        elif aligned_sequences_flipped[best_match] == 0 and best_flip == 0: 
            start_pos = ref_shift + best_shift
        # SITUATION 3: short sequence already present unflipped, align unflipped pair
        elif aligned_sequences_flipped[best_match] == 0 and best_flip == 0:
            start_pos = ref_shift + best_shift
        # SITUATION 4: longer seq already present in unflip form, align unflipped pair
        elif aligned_sequences_flipped[best_match] == 0 and best_flip == 1:
            new_seq = flip_sequence(new_seq)
            start_pos = ref_shift + best_shift

        # CASE 2: ONE OF THE SEQUENCES FROM THE PAIR IS ALREADY PRESENT IN THE ALIGNED MATRIX IN FLIPPED FORM
        # SITUATION 1: short seq already flipped, align in flipped version
        if aligned_sequences_flipped[best_match] == 1 and best_flip == 1:
            start_pos = ref_shift + best_shift
        # SITUATION 2: long seq already flipped, pair align unflipped
        elif aligned_sequences_flipped[best_match] == 1 and best_flip == 0:
            new_seq = flip_sequence(new_seq)
            new_shift = len(ref_seq) - (len(new_seq) + best_shift)
            start_pos = ref_shift + new_shift
        # SITUATION 3: long seq already flipped, align in flipped version
        elif aligned_sequences_flipped[best_match] == 1 and best_flip == 1:
            start_pos = ref_shift + best_shift
        # SITUATION 4: short seq already flipped, align unflipped pair
        elif aligned_sequences_flipped[best_match] == 1 and best_flip == 0:
            new_seq = flip_sequence(new_seq)
            start_pos = ref_shift + new_shift


        # Handle the alignment with the current matrix
        if len(new_seq) >= len(ref_seq):
            start_pos = ref_shift - best_shift
        else:
            start_pos = ref_shift + best_shift

        min_start = min(0, start_pos)
        shift_all = -min_start
        if shift_all > 0:
            for key in aligned_sequences:
                aligned_sequences[key] = [0] * shift_all + aligned_sequences[key]
            for k in sequence_shifts:
                sequence_shifts[k] += shift_all
            start_pos += shift_all

        current_width = max(len(v) for v in aligned_sequences.values())
        end_pos = start_pos + len(new_seq)
        final_width = max(current_width, end_pos)

        for k in aligned_sequences:
            aligned_sequences[k] += [0] * (final_width - len(aligned_sequences[k]))

        new_aligned = [0] * start_pos + new_seq + [0] * (final_width - start_pos - len(new_seq))
        aligned_sequences[best_candidate] = new_aligned
        sequence_shifts[best_candidate] = start_pos
        placed_sequences.add(best_candidate)
        remaining_sequences.remove(best_candidate)
        aligned_sequences_flipped[best_candidate] = best_flip  # Track flip status

    return aligned_sequences


# Function to read sequences from CSV file using .iterrows()
def read_sequences(csv_path):
    df = pd.read_csv(csv_path)
    df.set_index(df.columns[0], inplace=True)
    sequences = {}
    for index, row in df.iterrows():
        sequences[str(index)] = row.dropna().tolist()
    return sequences

# Function to read best Pearson results from CSV (now with Flipped column)
def read_best_pearsons(csv_path):
    best_pearsons = {}
    with open(csv_path, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            long_seq = row['Long_Sequence']
            short_seq = row['Short_Sequence']
            best_shift = int(row['Best_Shift'])
            best_corr = float(row['Best_Correlation'])
            flip = int(row['Flipped'])  # Read flip status
            best_pearsons[(short_seq, long_seq)] = (best_shift, best_corr, flip)
    return best_pearsons

# Function to save aligned sequences
def save_aligned_sequences(aligned_sequences, output_path):
    df = pd.DataFrame.from_dict(aligned_sequences, orient='index')
    df.index.name = 'SEQUENCE_ID'
    df.reset_index(inplace=True)
    df.to_csv(output_path, index=False)

# Main function
def run_alignment(input_csv, pearson_csv, output_csv):
    print(f"Reading sequences from: {input_csv}")
    sequences = read_sequences(input_csv)

    print(f"Reading Pearson results from: {pearson_csv}")
    best_pearsons = read_best_pearsons(pearson_csv)

    print("Aligning sequences...")
    aligned = align_sequences(sequences, best_pearsons)

    print(f"Saving aligned sequences to: {output_csv}")
    save_aligned_sequences(aligned, output_csv)

# Entry point
if __name__ == "__main__":
    input_csv = r"/scratch/shipras.sbb.iitmandi/flipped_1.11/90/Cluster 3_90th.csv"
    pearson_csv = r"/scratch/shipras.sbb.iitmandi/flipped_1.11/90/Cluster 3_90thp.csv"
    output_csv = r"/scratch/shipras.sbb.iitmandi/flipped_1.11/90/Cluster 3_90thAlign.csv"
    run_alignment(input_csv, pearson_csv, output_csv)
