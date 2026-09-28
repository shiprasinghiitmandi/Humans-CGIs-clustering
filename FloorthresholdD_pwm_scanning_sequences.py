import pandas as pd
import numpy as np
import math
import os


def load_pwm(file_path):
    pwm = pd.read_excel(file_path, index_col=0)
    return pwm


def load_sequences(file_path):
    sequences = pd.read_excel(file_path)
    return sequences


def load_expected_probs(file_path):
    df = pd.read_excel(file_path)
    expected_probs = {
        'dx': df.set_index('Value')['dx Normalized Frequency'].to_dict(),
        'dy': df.set_index('Value')['dy Normalized Frequency'].to_dict()
    }
    return expected_probs


def compute_weight(value, channel, pwm, col_name, expected_probs):
    try:
        return pwm.loc[value, col_name]
    except KeyError:
        obs_prob = 0.000000000001
        exp_prob = expected_probs[channel].get(value, 0.000000000001)
        if exp_prob < 0.000000000001:
            exp_prob = 0.000000000001
        return math.log(obs_prob / exp_prob)


def scan_sequence(seq_row, pwm, expected_probs, mode='original'):
    seq_id = seq_row.iloc[0]
    seq_values = seq_row.iloc[1:]
    seq_cols = seq_values.index.tolist()

    pwm_cols = pwm.columns.tolist()
    pwm_len = len(pwm_cols)

    if len(seq_cols) < pwm_len:
        return seq_id, None, None, 'Too short'

    best_score = None
    best_shift = None

    for shift in range(0, len(seq_cols) - pwm_len + 1, 2):
        window = seq_cols[shift:shift + pwm_len]
        score = 0.0
        valid = True
        for col, seq_col in zip(pwm_cols, window):
            channel = 'dy' if 'dy' in col else 'dx'
            val = seq_values[seq_col]
            try:
                val = int(val)
            except:
                valid = False
                break
            weight = compute_weight(val, channel, pwm, col, expected_probs)
            print(f"Seq: {seq_id}, Shift start: {shift}, Col: {col}, Value: {val}, Weight: {weight:.4f}")
            score += weight
        if valid:
            if best_score is None or score > best_score:
                best_score = score
                best_shift = shift

    return seq_id, best_score, best_shift, None


def flip_pwm(pwm):
    cols = pwm.columns.tolist()
    flipped_cols = cols[::-1]
    return pwm[flipped_cols]


def main(pwm_path, seq_path, prob_path, output_path):
    pwm = load_pwm(pwm_path)
    sequences = load_sequences(seq_path)
    expected_probs = load_expected_probs(prob_path)
    pwm_name = os.path.splitext(os.path.basename(pwm_path))[0]  # Extract PWM name

    results = []

    for i, row in sequences.iterrows():
        seq_id, best_score_1, shift1, reason1 = scan_sequence(row, pwm, expected_probs, 'original')
        _, best_score_2, shift2, reason2 = scan_sequence(row, flip_pwm(pwm), expected_probs, 'flipped')

        if reason1 == 'Too short' and reason2 == 'Too short':
            print(f"Skipping {seq_id} — too short")
            continue

        # Determine best score and where it came from
        if best_score_1 is not None and (best_score_2 is None or best_score_1 >= best_score_2):
            best_score = best_score_1
            best_shift = shift1
            flipped = 0
        else:
            best_score = best_score_2
            best_shift = shift2
            flipped = 1

        results.append([seq_id, pwm_name, best_score, best_shift, flipped])

    # Save results
    result_df = pd.DataFrame(results, columns=["Sequence_ID", "PWM_Name", "Best_PWM_Score", "Best_Shift_Position", "Flipped"])
    result_df.to_csv(output_path, index=False)
    print(f"Scores saved to {output_path}")


main(
    pwm_path=r"\Cluster 2_75thAlign_pwm.xlsx",
    seq_path=r"\excludedseqs_matrix_padded2.xlsx",
    prob_path=r"\combined_cgiandnoncgi_normalized2.xlsx",
    output_path=r"\cluster2_score.csv"
)
