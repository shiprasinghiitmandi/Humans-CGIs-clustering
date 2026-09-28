import os
import pandas as pd
import numpy as np

pfm_input_directory = r"90aln_pfm_outputs" # Directory with PFM CSVs
expected_probabilities_path = r"combined_cgiandnoncgi_normalized.csv"
pwm_output_directory = r"90aln_PWM_outputs"

os.makedirs(pwm_output_directory, exist_ok=True)

#  Load expected probabilities
expected_df = pd.read_csv(expected_probabilities_path)
expected_df.columns = ['Value', 'dx Normalized Frequency', 'dy Normalized Frequency']
expected_dict_dx = expected_df.set_index('Value')['dx Normalized Frequency'].to_dict()
expected_dict_dy = expected_df.set_index('Value')['dy Normalized Frequency'].to_dict()

#  Process Each PFM File
for filename in os.listdir(pfm_input_directory):
    if filename.endswith(".csv"):
        pfm_path = os.path.join(pfm_input_directory, filename)
        print(f"Processing PFM: {filename}")

        pfm_df = pd.read_csv(pfm_path)

        # Extract values and headers
        values = pfm_df['Value'].tolist()
        data_only = pfm_df.drop(columns=['Value'])
        pwm_data = []

        for i, row in data_only.iterrows():
            pwm_row = []
            val = values[i]  # This is the row header (1, 2, 3,...)
            for col_name in data_only.columns:
                obs = row[col_name]

                # Determine if it's a dx or dy column
                if col_name.lower().startswith('dy'):
                    expected = expected_dict_dy.get(val,0.000000000001)
                elif col_name.lower().startswith('dx'):
                    expected = expected_dict_dx.get(val, 0.000000000001)
                else:
                    expected =0.000000000001 # Default in case of unknown column type

                # Enforce minimum threshold
                expected = max(expected, 0.000000000001)

                # Calculate log(observed/expected)
                if obs == 0:
                    pwm_row.append(np.log(0.000000000001/ expected))
                else:
                    pwm_row.append(np.log(obs / expected))
            pwm_data.append(pwm_row)

        # Create PWM DataFrame
        pwm_df = pd.DataFrame(pwm_data, columns=data_only.columns)
        pwm_df.insert(0, 'Value', values)

        # Save to output
        pwm_filename = os.path.join(pwm_output_directory, filename)
        pwm_df.to_csv(pwm_filename, index=False)
        print(f"Saved PWM: {pwm_filename}")
