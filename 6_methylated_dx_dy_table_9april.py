import pandas as pd
from docx import Document
import re
import os
from glob import glob


def read_highlighted_positions(file_path):
    doc = Document(file_path)
    sequences = {}
    seq_id = None

    for para in doc.paragraphs:
        text = para.text.strip()

        match = re.match(r"sequence\s+(\d+)", text, re.IGNORECASE)
        if match:
            seq_id = int(match.group(1))
            continue

        if seq_id is None or not text:
            continue

        sequence = ""
        highlights = set()
        seq_start = 0

        for run in para.runs:
            run_text = run.text
            sequence += run_text
            if run.font.highlight_color == 6:
                highlights.update(range(seq_start, seq_start + len(run_text)))
            seq_start += len(run_text)

        sequences[seq_id] = (sequence, highlights)
        seq_id = None

    return sequences


def extract_dy_dx_highlights(sequences):
    rows = []

    for seq_id, (seq, highlights) in sequences.items():
        i = 0
        dy_count = 1
        dx_count = 1

        while i < len(seq):
            if seq[i:i + 2] == 'CG':
                start = i
                is_highlighted = False
                length = 0
                while i < len(seq) - 1 and seq[i:i + 2] == 'CG':
                    if i in highlights or (i + 1) in highlights:
                        is_highlighted = True
                    i += 2
                    length += 2

                if is_highlighted:
                    highlight_len = sum(1 for j in range(start, start + length) if j in highlights)
                    rows.append({
                        "Sequence_ID": seq_id,
                        "dy_values": highlight_len,
                        "dy_position": f"dy{dy_count}",
                        "dx_values": None,
                        "dx_position": None
                    })
                dy_count += 1

            elif seq[i] == 'N':
                start = i
                is_highlighted = False
                length = 0
                while i < len(seq) and seq[i] == 'N':
                    if i in highlights:
                        is_highlighted = True
                    i += 1
                    length += 1

                if is_highlighted:
                    highlight_len = sum(1 for j in range(start, start + length) if j in highlights)
                    rows.append({
                        "Sequence_ID": seq_id,
                        "dy_values": None,
                        "dy_position": None,
                        "dx_values": highlight_len,
                        "dx_position": f"dx{dx_count}"
                    })
                dx_count += 1
            else:
                i += 1

    return pd.DataFrame(rows)


def save_to_excel(df, output_file):
    df.to_excel(output_file, index=False)
    print(f" Saved to: {output_file}")


def main():
    input_dir = r"/scratch/shipras.sbb.iitmandi/meth_map_codes/cancerlung_renamed_docx/"
    output_dir = r"/scratch/shipras.sbb.iitmandi/meth_map_codes/cancerlung_abs_dydx_excel/"

    os.makedirs(output_dir, exist_ok=True)

    docx_files = glob(os.path.join(input_dir, "*.docx"))

    if not docx_files:
        print(" No .docx files found in input directory.")
        return

    for docx_file in docx_files:
        base_name = os.path.splitext(os.path.basename(docx_file))[0]
        output_excel = os.path.join(output_dir, f"{base_name}_dydx.xlsx")

        print(f"Processing: {docx_file}")
        sequences = read_highlighted_positions(docx_file)
        highlight_df = extract_dy_dx_highlights(sequences)
        save_to_excel(highlight_df, output_excel)


if __name__ == "__main__":
    main()
