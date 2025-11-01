import pandas as pd
import matplotlib
import matplotlib.pyplot as plt

def calculate_and_plot(file1, file2, excel_output, plot_output):
    # Read input CSVs
    df1 = pd.read_csv(file1, sep=None, engine='python')
    df2 = pd.read_csv(file2, sep=None, engine='python')

    # -------- Step 1: Extract and Sum from File 1 --------
    a1 = df1['dy_values'].fillna(0).sum()
    a2 = df1['dx_values'].fillna(0).sum()
    a = a1 + a2

    # -------- Step 2: Extract and Sum from File 2 --------
    dy_columns = [col for col in df2.columns if col.lower().startswith('dy')]
    dx_columns = [col for col in df2.columns if col.lower().startswith('dx')]

    n1 = df2[dy_columns].fillna(0).values.sum()
    n2 = df2[dx_columns].fillna(0).values.sum()
    n = n1 + n2

    # -------- Step 1: Save All Values to Excel --------
    results = {
        'Metric': ['dy+dx (a)', 'dy+dx (n)', 'dy only (a1)', 'dy only (n1)', 'dx only (a2)', 'dx only (n2)'],
        'Value': [a, n, a1, n1, a2, n2]
    }
    result_df = pd.DataFrame(results)
    result_df.to_excel(excel_output, index=False)

    # -------- Step 4: Plotting dy and dx only --------
    plt.figure(figsize=(6, 5))

    plt.bar(0, a1 / n1, color='darkblue', width=0.4, label='p_dy')
    plt.bar(0.5, a2 / n2, color='darkred', width=0.4, label='p_dx')

    plt.xticks([0, 0.5], ['p_dy', 'p_dx'])
    plt.ylabel("Normalized Methylation fraction")
    plt.title("90th_perc cluster_6 methylation(healthy)")
    plt.legend()
    plt.tight_layout()
    plt.savefig(plot_output)
    plt.show()



csv_file1 = r"F:\1_normal_methdata_n_mapping\unflipped_norm_methmap_1.17\90\cluster6_methmap.csv"
csv_file2 = r"F:\1_normal_methdata_n_mapping\unflipped_norm_methmap_1.17\1.17o\90\Cluster 6_90th.csv"
excel_output = r"F:\1_normal_methdata_n_mapping\unflipped_norm_methmap_1.17\90_methclust_plots\90clust6_calculations.xlsx"
plot_output = r"F:\1_normal_methdata_n_mapping\unflipped_norm_methmap_1.17\90_methclust_plots\90clust6_meth_histogram.png"

calculate_and_plot(csv_file1, csv_file2, excel_output, plot_output)
