import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def compute_normalized_sums(file_path):
    df = pd.read_csv(file_path)
    df_trimmed = df.iloc[:, :-1]  # Drop flipped_status (last column)
    num_sequences = df_trimmed.shape[0]

    dy_cols = [col for col in df_trimmed.columns if col.startswith("dy")]
    dx_cols = [col for col in df_trimmed.columns if col.startswith("dx")]

    dy_sum = df_trimmed[dy_cols].replace(0, np.nan).sum().sum(skipna=True)
    dx_sum = df_trimmed[dx_cols].replace(0, np.nan).sum().sum(skipna=True)

    dy_norm = dy_sum / num_sequences
    dx_norm = dx_sum / num_sequences

    return dy_norm, dx_norm

def process_tissue_folders(base_input_folder):
    result_data = []

    for folder in sorted(os.listdir(base_input_folder)):
        tissue_path = os.path.join(base_input_folder, folder)
        if not os.path.isdir(tissue_path):
            continue

        tissue_name = folder.replace("_methmap_files", "").strip()

        for file in sorted(os.listdir(tissue_path)):
            if file.endswith(".csv") and file.startswith("Cluster"):
                file_path = os.path.join(tissue_path, file)

                cluster_part = file.split('_')[0]  # 'Cluster 1'
                cluster_num = int(cluster_part.split()[1])  # extract number

                dy_norm, dx_norm = compute_normalized_sums(file_path)

                result_data.append({
                    "Tissue": tissue_name,
                    "Cluster": cluster_part,
                    "ClusterNum": cluster_num,
                    "dy": dy_norm,
                    "dx": dx_norm
                })

    return pd.DataFrame(result_data)

def plot_custom_horizontal_bar_chart(df_summary, output_path):
    df_summary = df_summary.sort_values(by=["Tissue", "ClusterNum"], ascending=[True, True]).reset_index(drop=True)

    # Prepare bar positions
    y_labels = []
    tissue_positions = {}
    dy_values = []
    dx_values = []

    y_pos = []
    current_y = 0

    for tissue, group in df_summary.groupby("Tissue"):
        group = group.sort_values("ClusterNum", ascending=True)
        n = len(group)

        midpoint = current_y + (n - 1) / 2.0
        tissue_positions[tissue] = midpoint

        for _, row in group.iterrows():
            label = f"{row['Cluster']}"
            y_labels.append(label)
            y_pos.append(current_y)
            dy_values.append(row['dy'])
            dx_values.append(row['dx'])
            current_y += 1

    # Plot
    fig, ax = plt.subplots(figsize=(10, 12), dpi=500)
    bar_height = 0.35

    ax.barh([y - bar_height/2 for y in y_pos], dy_values, height=bar_height, color='darkkhaki', label='dy')
    ax.barh([y + bar_height/2 for y in y_pos], dx_values, height=bar_height, color='teal', label='dx')

    ax.set_yticks(y_pos)
    ax.set_yticklabels(y_labels, fontsize=18)
    ax.set_xlabel('Normalized total Sum', fontsize=18)
    #ax.set_title('Normalized dy and dx Sums Across Clusters and Tissues', fontsize=14)
    ax.legend(fontsize=16)
    ax.grid(axis='x', linestyle='--', alpha=0.5)


    for tissue, y_center in tissue_positions.items():
        ax.text(-0.15, y_center, tissue, fontsize=18, rotation=90, va='center', ha='right', transform=ax.get_yaxis_transform())
    ax.tick_params(axis='x', labelsize=18)  # You can adjust the number (e.g., 16, 18)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.show()


base_input_folder = r"F:\flipped_matrices\random_methylation_simulation\methmaps_tissues_folder"  # contains folders like "Healthy Liver_methmap_files"
output_plot = r"F:\flipped_matrices\random_methylation_simulation\all_tissues_normalisedsum_methlevels_correct_r.png"

summary_df = process_tissue_folders(base_input_folder)
plot_custom_horizontal_bar_chart(summary_df, output_plot)
# === SAVE EXTENDED DATAFRAME TO CSV ===
def save_detailed_csv(df_summary, base_input_folder, output_csv_path):
    detailed_rows = []

    for folder in sorted(os.listdir(base_input_folder)):
        tissue_path = os.path.join(base_input_folder, folder)
        if not os.path.isdir(tissue_path):
            continue

        tissue_name = folder.replace("_methmap_files", "").strip()

        for file in sorted(os.listdir(tissue_path)):
            if file.endswith(".csv") and file.startswith("Cluster"):
                file_path = os.path.join(tissue_path, file)

                cluster_part = file.split('_')[0]  # 'Cluster 1'
                cluster_num = int(cluster_part.split()[1])  # extract number

                df = pd.read_csv(file_path)
                df_trimmed = df.iloc[:, :-1]
                num_sequences = df_trimmed.shape[0]

                dy_cols = [col for col in df_trimmed.columns if col.startswith("dy")]
                dx_cols = [col for col in df_trimmed.columns if col.startswith("dx")]

                dy_sum = df_trimmed[dy_cols].replace(0, np.nan).sum().sum(skipna=True)
                dx_sum = df_trimmed[dx_cols].replace(0, np.nan).sum().sum(skipna=True)

                dy_norm = dy_sum / num_sequences
                dx_norm = dx_sum / num_sequences

                detailed_rows.append({
                    "Tissue": tissue_name,
                    "Cluster": cluster_part,
                    "ClusterNum": cluster_num,
                    "NumSequences": num_sequences,
                    "dy_sum": dy_sum,
                    "dx_sum": dx_sum,
                    "dy_normalized": dy_norm,
                    "dx_normalized": dx_norm
                })

    detailed_df = pd.DataFrame(detailed_rows)
    detailed_df = detailed_df.sort_values(by=["Tissue", "ClusterNum"])
    detailed_df.to_csv(output_csv_path, index=False)
    print(f"Detailed CSV saved to: {output_csv_path}")


output_csv_path = r"F:\flipped_matrices\random_methylation_simulation\r_dy_dx_methylationsummary.csv"
save_detailed_csv(summary_df, base_input_folder, output_csv_path)
