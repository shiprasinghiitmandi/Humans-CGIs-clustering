import pandas as pd
import numpy as np
import os


INPUT_PWM = r"D:\CGI_patterns_manuscript\Revision_CBC_codesandfiles\shuffledPWMs_pseudocountBonferoni\Cluster 2_90thAlign_PWM.xlsx"

OUTPUT_DIRECTORY = r"D:\CGI_patterns_manuscript\Revision_CBC_codesandfiles\shuffledPWMs_pseudocountBonferoni\shuffled_PWMs_dy_dx_random2"

NUMBER_OF_PERMUTATIONS = 5



os.makedirs(OUTPUT_DIRECTORY, exist_ok=True)



pwm = pd.read_excel(INPUT_PWM)

if "Value" not in pwm.columns:
    raise ValueError("The PWM file must contain a 'Value' column.")



pwm_columns = [col for col in pwm.columns if col != "Value"]

dy_columns = [col for col in pwm_columns if "dy" in str(col).lower()]
dx_columns = [col for col in pwm_columns if "dx" in str(col).lower()]

print(f"Number of dy columns: {len(dy_columns)}")
print(f"Number of dx columns: {len(dx_columns)}")

print("\nOriginal dy columns:")
print(dy_columns)

print("\nOriginal dx columns:")
print(dx_columns)


# No seed = different random permutations each time
rng = np.random.default_rng()



permutation_records = []

for i in range(1, NUMBER_OF_PERMUTATIONS + 1):


    # Shuffle dy columns ONLY among dy columns


    shuffled_dy = rng.permutation(dy_columns).tolist()


    # Shuffle dx columns ONLY among dx columns


    shuffled_dx = rng.permutation(dx_columns).tolist()


    # Reconstruct PWM in ORIGINAL dy/dx positional pattern


    shuffled_columns = []

    dy_index = 0
    dx_index = 0

    for column in pwm_columns:

        if "dy" in str(column).lower():

            shuffled_columns.append(
                shuffled_dy[dy_index]
            )

            dy_index += 1

        elif "dx" in str(column).lower():

            shuffled_columns.append(
                shuffled_dx[dx_index]
            )

            dx_index += 1

        else:

            raise ValueError(
                f"Unexpected PWM column: {column}"
            )


    # Create shuffled PWM


    shuffled_pwm = pwm[["Value"] + shuffled_columns]


    # Save shuffled PWM


    output_file = os.path.join(
        OUTPUT_DIRECTORY,
        f"shuffled_PWM_{i}.xlsx"
    )

    shuffled_pwm.to_excel(
        output_file,
        index=False
    )


    # Record permutation


    permutation_records.append(
        [i] + shuffled_columns
    )

    print(f"\nPermutation {i}")
    print("Shuffled column order:")
    print(shuffled_columns)

    print(f"Saved to:")
    print(output_file)

    print("-" * 70)


# ==========================================================
# SAVE PERMUTATION INFORMATION
# ==========================================================

permutation_df = pd.DataFrame(
    permutation_records,
    columns=[
        "Permutation"
    ] + [
        f"Position_{j+1}"
        for j in range(len(pwm_columns))
    ]
)

permutation_file = os.path.join(
    OUTPUT_DIRECTORY,
    "permutation_column_orders_dy_dx.xlsx"
)

permutation_df.to_excel(
    permutation_file,
    index=False
)


# ==========================================================
# FINISHED
# ==========================================================

print("\nAll 5 dy/dx-restricted shuffled PWMs generated successfully.")
print(f"Permutation information saved to:")
print(permutation_file)