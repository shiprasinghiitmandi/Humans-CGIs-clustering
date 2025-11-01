import pandas as pd
import matplotlib.pyplot as plt
import os

def process_and_plot(input_csv, output_csv):
    # Read the input CSV
    df = pd.read_csv(input_csv)

    # Drop the flipped_status column(as it was not required while plotting curves)
    df_trimmed = df.iloc[:, :-1]

    # Extract only numeric columns (exclude SEQUENCE_ID)
    data_only = df_trimmed.iloc[:, 1:]

    # Calculate the number of sequences (rows)
    num_sequences = df_trimmed.shape[0]
    print(f"Total number of sequences: {num_sequences}")

    # Sum each column
    column_sums = data_only.sum()

    # Normalize the column sums
    normalized_sums = column_sums / num_sequences

    # Append SUM and NORMALIZED rows to the DataFrame
    new_rows = pd.DataFrame([["SUM"] + column_sums.tolist(), ["NORMALIZED"] + normalized_sums.tolist()],
                            columns=df_trimmed.columns)
    df_output = pd.concat([df_trimmed, new_rows], ignore_index=True)

    # Save to a new CSV
    df_output.to_csv(output_csv, index=False)
    print(f"Output saved to: {output_csv}")

    # Plotting
    col_names = data_only.columns.tolist()
    y_values = normalized_sums.tolist()

    # Separate dy and dx values
    dy_labels = [col for col in col_names if col.startswith("dy")]
    dx_labels = [col for col in col_names if col.startswith("dx")]

    dy_y = [normalized_sums[col] for col in dy_labels]
    dx_y = [normalized_sums[col] for col in dx_labels]

    # Insert origin point
    dy_x_plot = [0] + list(range(1, len(dy_y)+1))
    dx_x_plot = [0] + list(range(1, len(dx_y)+1))
    dy_y_plot = [0] + dy_y
    dx_y_plot = [0] + dx_y

    # Plot
    plt.figure(figsize=(14, 6))

    # Plot dy line and dots
    plt.plot(dy_x_plot, dy_y_plot, color='darkkhaki', linewidth=1.5, label='dy')
    plt.scatter(dy_x_plot[1:], dy_y_plot[1:], color='darkkhaki')  # skip origin for dot

    # Plot dx line and dots
    plt.plot(dx_x_plot, dx_y_plot, color='teal', linewidth=1.5, label='dx')
    plt.scatter(dx_x_plot[1:], dx_y_plot[1:], color='teal')  # skip origin for dot

    # Custom legend
    from matplotlib.lines import Line2D
    legend_elements = [
        Line2D([0], [0], color='darkkhaki', marker='o', label='dy', markersize=8),
        Line2D([0], [0], color='teal', marker='o', label='dx', markersize=8)
    ]
    plt.legend(handles=legend_elements)

    plt.xlabel("Indices")
    plt.ylabel("Normalized Sum")
    plt.title("90th cluster 1_random methylation (Healthy Brain)")
    plt.grid(False)
    plt.tight_layout()

    # Save the figure
    plot_path = os.path.splitext(output_csv)[0] + ".png"
    plt.savefig(plot_path, dpi=300)
    print(f"Plot saved to: {plot_path}")

    # Show the plot
    plt.show()


process_and_plot(
    input_csv=r"F:\flipped_matrices\random_methylation_simulation\Cluster 1_90th_methylation_r3.csv",
    output_csv= r"F:\flipped_matrices\random_methylation_simulation\Cl1_random_output_with_sums.csv"
)
