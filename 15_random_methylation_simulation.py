import pandas as pd
import numpy as np

#step 1
alignment_file = r"F:\flipped_matrices\random_methylation_simulation\Cancerous_lung\Cluster 5_90thAlign.csv"  # <-- change per cluster
output_file = r"F:\flipped_matrices\random_methylation_simulation\Cancerous_lung\Cluster 5_90th_methylation_r.csv"

df_align = pd.read_csv(alignment_file)

# Step 2: Simulation parameters
p_methylation = 0.025  # probability of CG methylation (change as needed)
num_replicates = 100  # number of simulated cells (change as needed)

print(f"Simulating with p={p_methylation}, replicates={num_replicates}")


# Step 3: Simulation function
def simulate_cell(col_name, value):
    """ We Simulated methylation counts for a single dy/dx cell value."""
    if value == 0:
        return 0

    if col_name.startswith("dy"):
        trials = value // 2  # number of CGs
        return np.random.binomial(trials, p_methylation)

    elif col_name.startswith("dx"):
        return 0  # dx methylation fixed at 0

    else:
        return value  # keep non-dx/dy columns unchanged


# Step 4: Running replicates(increasing the cell counts)
all_replicates = []

for r in range(num_replicates):
    df_sim = df_align.copy()
    for col in df_sim.columns:
        if col.startswith("dy") or col.startswith("dx"):
            df_sim[col] = df_sim[col].apply(lambda v: simulate_cell(col, v))
    all_replicates.append(df_sim)

# Step 5: Averaging across cells (replicates)
# Keep only dy/dx numeric columns averaged, keep SEQUENCE_ID as first column
df_avg = df_align.copy()

for col in df_avg.columns:
    if col.startswith("dy") or col.startswith("dx"):
        df_avg[col] = np.mean([rep[col].values for rep in all_replicates], axis=0)

#Step 6: Saving output
df_avg.to_csv(output_file, index=False)
print(f"Averaged simulated methylation levels written to {output_file}")
