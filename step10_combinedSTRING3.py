import os
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt

from matplotlib.lines import Line2D
from scipy.stats import norm
from scipy.stats import normaltest



# CONFIGURATION


BASE_DIR = (
    r"path to the directory containing STRING network files"
)


CLUSTERS = {

    "Cluster 1": {
        "node_file": os.path.join(
            BASE_DIR,
            "cl1_string",
            "STRING_cl1Node.csv"
        ),
        "edge_file": os.path.join(
            BASE_DIR,
            "cl1_string",
            "STRING_cl1edge.csv"
        )
    },

    "Cluster 2": {
        "node_file": os.path.join(
            BASE_DIR,
            "cl2_string",
            "STRING_cl2Node.csv"
        ),
        "edge_file": os.path.join(
            BASE_DIR,
            "cl2_string",
            "STRING_cl2edge.csv"
        )
    },

    "Cluster 3": {
        "node_file": os.path.join(
            BASE_DIR,
            "cl3_string",
            "STRING_cl3Node.csv"
        ),
        "edge_file": os.path.join(
            BASE_DIR,
            "cl3_string",
            "STRING_cl3edge.csv"
        )
    },

    "Cluster 4": {
        "node_file": os.path.join(
            BASE_DIR,
            "cl4_string",
            "STRING_cl4Node.csv"
        ),
        "edge_file": os.path.join(
            BASE_DIR,
            "cl4_string",
            "STRING_cl4edge.csv"
        )
    },

    "Cluster 5": {
        "node_file": os.path.join(
            BASE_DIR,
            "cl5_string",
            "STRING_cl5Node.csv"
        ),
        "edge_file": os.path.join(
            BASE_DIR,
            "cl5_string",
            "STRING_cl5edge.csv"
        )
    }
}



# OUTPUT DIRECTORY


OUTPUT_DIR = os.path.join(
    BASE_DIR,
    "Combined_STRING_DPM_Zscore"
)

os.makedirs(
    OUTPUT_DIR,
    exist_ok=True
)



# RANDOMIZATION SETTINGS


N_RANDOM = 1000

SWAP_FACTOR = 10

RANDOM_SEED = 42



# PLOT COLORS


RANDOM_COLOR = "rosybrown"
REAL_COLOR = "darksalmon"



# GLOBAL PLOT FORMATTING


plt.rcParams.update({

    "font.family": "Times New Roman",

    "font.size": 14,

    "text.color": "black",

    "axes.labelcolor": "black",

    "axes.titlecolor": "black",

    "xtick.color": "black",

    "ytick.color": "black",

    "axes.edgecolor": "black",

    "pdf.fonttype": 42,

    "ps.fonttype": 42
})



# HELPER — Z-SCORE STATISTICS


def z_to_normal_statistics(z):
    """
    Convert a Z-score to:

    1. Absolute Z-score
    2. Standard normal CDF: Phi(|Z|)
    3. Standard normal survival function: SF(|Z|)
    4. Two-sided normal-theory p-value

    p = 2 * SF(|Z|)
    """

    if not np.isfinite(z):

        return (
            np.nan,
            np.nan,
            np.nan,
            np.nan
        )


    abs_z_score = abs(z)


    normal_cdf = norm.cdf(
        abs_z_score
    )


    normal_sf = norm.sf(
        abs_z_score
    )


    p_value = (
        2.0 *
        normal_sf
    )


    return (
        abs_z_score,
        normal_cdf,
        normal_sf,
        p_value
    )



# COMMON PLOT FORMATTING


def format_axis(
    ax,
    bottom_row=False
):

    """
    Apply common formatting to every subplot.

    The x-axis category label is shown only on
    the bottom row of the combined figures.
    """


    # Spines


    ax.spines["top"].set_visible(False)

    ax.spines["right"].set_visible(False)

    ax.spines["left"].set_color("black")

    ax.spines["bottom"].set_color("black")

    ax.spines["left"].set_linewidth(1.0)

    ax.spines["bottom"].set_linewidth(1.0)



    # Tick formatting


    ax.tick_params(
        axis="both",
        which="major",
        labelsize=9,
        colors="black",
        width=1.0
    )



    # Make tick labels black and bold


    for label in ax.get_xticklabels():

        label.set_fontname(
            "Times New Roman"
        )

        label.set_fontweight(
            "bold"
        )

        label.set_color(
            "black"
        )


    for label in ax.get_yticklabels():

        label.set_fontname(
            "Times New Roman"
        )

        label.set_fontweight(
            "bold"
        )

        label.set_color(
            "black"
        )



    # X-axis label


    ax.set_xticks(
        [1]
    )


    if bottom_row:

        ax.set_xticklabels(
            ["Random networks"],
            fontname="Times New Roman",
            fontsize=9,
            fontweight="bold",
            color="black"
        )

        ax.tick_params(
            axis="x",
            which="major",
            length=4
        )

    else:

        ax.set_xticklabels([])

        ax.tick_params(
            axis="x",
            which="major",
            length=0
        )



    # Grid


    ax.grid(
        axis="y",
        linestyle="--",
        linewidth=0.6,
        alpha=0.3,
        color="black"
    )



# HELPER — AXIS LABEL


def set_ylabel(
    ax,
    text
):

    ax.set_ylabel(
        text,
        fontsize=10,
        fontname="Times New Roman",
        fontweight="bold",
        color="black"
    )



# HELPER — SUBPLOT TITLE


def set_subplot_title(
    ax,
    cluster_name,
    centrality_name
):

    ax.set_title(
        f"{cluster_name}",
        fontsize=11,
        fontname="Times New Roman",
        fontweight="bold",
        color="black",
        pad=7
    )



# HELPER — STATISTICAL ANNOTATION


def add_statistics(
    ax,
    difference,
    z_score,
    p_value
):

    ax.text(
        0.05,
        0.90,
        (
            f"Δ = {difference:.3g}\n"
            f"Z = {z_score:.3g}\n"
            f"p = {p_value:.3g}"
        ),
        transform=ax.transAxes,
        fontsize=8,
        fontname="Times New Roman",
        fontweight="bold",
        color="black",
        verticalalignment="top"
    )



# LEGEND HANDLE


real_network_legend = Line2D(
    [0],
    [0],
    marker="o",
    color="none",
    markerfacecolor=REAL_COLOR,
    markeredgecolor="black",
    markersize=5,
    label="Real network"
)



# STORAGE


all_results = {}

summary_rows = []



# PROCESS EACH CLUSTER


for cluster_name, cluster_info in CLUSTERS.items():

    print("\n")
    print("=" * 75)
    print(f"PROCESSING {cluster_name}")
    print("=" * 75)


    NODE_FILE = cluster_info["node_file"]

    EDGE_FILE = cluster_info["edge_file"]



    # 1. READ NODE TABLE


    print("\nReading node table...")

    nodes = pd.read_csv(
        NODE_FILE
    )


    required_node_columns = [
        "@id",
        "Degree",
        "BetweennessCentrality",
        "ClosenessCentrality"
    ]


    for col in required_node_columns:

        if col not in nodes.columns:

            raise ValueError(
                f"'{col}' not found in {NODE_FILE}"
            )


    print(
        f"Number of nodes in node table: "
        f"{len(nodes)}"
    )



    # 2. READ EDGE TABLE


    print("\nReading edge table...")

    edges = pd.read_csv(
        EDGE_FILE
    )


    if "name" not in edges.columns:

        raise ValueError(
            f"'name' column not found in {EDGE_FILE}"
        )



    # 3. EXTRACT NODE IDS FROM EDGE NAMES


    def extract_nodes(edge_name):

        if pd.isna(edge_name):

            return None, None


        parts = str(
            edge_name
        ).split(" (pp) ")


        if len(parts) != 2:

            return None, None


        node_a = parts[0].strip()

        node_b = parts[1].strip()


        return node_a, node_b


    edges[
        ["Node_A", "Node_B"]
    ] = edges["name"].apply(
        lambda x: pd.Series(
            extract_nodes(x)
        )
    )



    # 4. CLEAN EDGES


    edges = edges.dropna(
        subset=[
            "Node_A",
            "Node_B"
        ]
    )



    # Remove self-loops


    edges = edges[
        edges["Node_A"] != edges["Node_B"]
    ]



    # Remove duplicate undirected edges


    edges["Edge_key"] = edges.apply(
        lambda row: tuple(
            sorted(
                [
                    row["Node_A"],
                    row["Node_B"]
                ]
            )
        ),
        axis=1
    )


    edges = edges.drop_duplicates(
        subset="Edge_key"
    )


    print(
        f"Number of edges after cleaning: "
        f"{len(edges)}"
    )



    # 5. CONSTRUCT REAL NETWORK


    G_real = nx.Graph()


    G_real.add_edges_from(
        edges[
            ["Node_A", "Node_B"]
        ].itertuples(
            index=False,
            name=None
        )
    )


    print("\nReal network:")

    print(
        f"Nodes : "
        f"{G_real.number_of_nodes()}"
    )

    print(
        f"Edges : "
        f"{G_real.number_of_edges()}"
    )



    # 6. MATCH NODE TABLE WITH NETWORK


    nodes["@id"] = (
        nodes["@id"]
        .astype(str)
        .str.strip()
    )


    nodes["STRING_ID"] = (
        nodes["@id"]
        .str.replace(
            "stringdb:",
            "",
            regex=False
        )
    )


    network_nodes = set(
        G_real.nodes()
    )


    table_nodes = set(
        nodes["STRING_ID"]
    )


    common_nodes = (
        network_nodes.intersection(
            table_nodes
        )
    )


    print("\nNode matching:")

    print(
        f"Network nodes  : "
        f"{len(network_nodes)}"
    )

    print(
        f"Node-table IDs : "
        f"{len(table_nodes)}"
    )

    print(
        f"Common IDs     : "
        f"{len(common_nodes)}"
    )


    if len(common_nodes) == 0:

        raise ValueError(
            f"No common node IDs found "
            f"for {cluster_name}."
        )



    # Keep only matched nodes


    nodes = nodes[
        nodes["STRING_ID"].isin(
            common_nodes
        )
    ].copy()



    # Keep only matched network nodes


    G_real = G_real.subgraph(
        common_nodes
    ).copy()


    print(
        f"Matched network: "
        f"{G_real.number_of_nodes()} nodes, "
        f"{G_real.number_of_edges()} edges"
    )



    # 7. REAL NETWORK CENTRALITY


    real_betweenness = (
        pd.to_numeric(
            nodes[
                "BetweennessCentrality"
            ],
            errors="coerce"
        )
        .mean()
    )


    real_closeness = (
        pd.to_numeric(
            nodes[
                "ClosenessCentrality"
            ],
            errors="coerce"
        )
        .mean()
    )


    print("\nReal-network centrality:")

    print(
        f"Mean Betweenness : "
        f"{real_betweenness:.10g}"
    )

    print(
        f"Mean Closeness   : "
        f"{real_closeness:.10g}"
    )



    # 8. DEGREE SEQUENCE


    real_degree_sequence = sorted(
        dict(
            G_real.degree()
        ).values()
    )



    # 9. RANDOM NUMBER GENERATOR


    cluster_seed = (
        RANDOM_SEED
        +
        int(
            cluster_name.split()[-1]
        )
    )


    rng = np.random.default_rng(
        cluster_seed
    )



    # 10. RANDOM NETWORK STORAGE


    random_betweenness = []

    random_closeness = []

    successful_networks = 0


    number_of_edges = (
        G_real.number_of_edges()
    )


    n_swaps = max(
        number_of_edges * SWAP_FACTOR,
        10
    )


    max_tries = max(
        n_swaps * 20,
        100
    )



    # 11. GENERATE DEGREE-PRESERVING RANDOM NETWORKS


    print(
        "\nGenerating degree-preserving "
        "random networks..."
    )


    for i in range(
        N_RANDOM
    ):

        G_random = G_real.copy()


        try:

            nx.double_edge_swap(
                G_random,
                nswap=n_swaps,
                max_tries=max_tries,
                seed=int(
                    rng.integers(
                        0,
                        2**32 - 1
                    )
                )
            )


        except nx.NetworkXAlgorithmError:

            print(
                f"Random network "
                f"{i + 1} failed; skipping."
            )

            continue


        #
        # Verify exact degree sequence
        #

        random_degree_sequence = sorted(
            dict(
                G_random.degree()
            ).values()
        )


        if (
            random_degree_sequence
            != real_degree_sequence
        ):

            print(
                f"Degree sequence mismatch "
                f"in random network "
                f"{i + 1}; skipping."
            )

            continue


        #
        # Calculate betweenness
        #

        betweenness = (
            nx.betweenness_centrality(
                G_random
            )
        )


        #
        # Calculate closeness
        #

        closeness = (
            nx.closeness_centrality(
                G_random
            )
        )


        #
        # Mean centrality across nodes
        #

        random_betweenness.append(
            np.mean(
                list(
                    betweenness.values()
                )
            )
        )


        random_closeness.append(
            np.mean(
                list(
                    closeness.values()
                )
            )
        )


        successful_networks += 1


        if (
            successful_networks % 100
            == 0
        ):

            print(
                f"Successful random networks: "
                f"{successful_networks}/"
                f"{N_RANDOM}"
            )



    # 12. CONVERT RANDOM RESULTS TO ARRAYS


    random_betweenness = np.asarray(
        random_betweenness,
        dtype=float
    )


    random_closeness = np.asarray(
        random_closeness,
        dtype=float
    )


    if successful_networks == 0:

        raise RuntimeError(
            f"No randomized networks "
            f"were successfully generated "
            f"for {cluster_name}."
        )


    print(
        f"\nSuccessfully generated "
        f"{successful_networks} "
        f"randomized networks."
    )



    # 13. RANDOM DISTRIBUTION SUMMARY


    random_betweenness_mean = (
        np.mean(
            random_betweenness
        )
    )


    random_betweenness_sd = (
        np.std(
            random_betweenness,
            ddof=1
        )
    )


    random_closeness_mean = (
        np.mean(
            random_closeness
        )
    )


    random_closeness_sd = (
        np.std(
            random_closeness,
            ddof=1
        )
    )



    # 14. DIFFERENCE BETWEEN REAL AND RANDOM MEAN


    betweenness_difference = (
        real_betweenness
        -
        random_betweenness_mean
    )


    closeness_difference = (
        real_closeness
        -
        random_closeness_mean
    )



    # 15. Z-SCORE


    if random_betweenness_sd == 0:

        raise RuntimeError(
            f"Random betweenness SD is zero "
            f"for {cluster_name}."
        )


    if random_closeness_sd == 0:

        raise RuntimeError(
            f"Random closeness SD is zero "
            f"for {cluster_name}."
        )


    betweenness_z = (
        betweenness_difference
        /
        random_betweenness_sd
    )


    closeness_z = (
        closeness_difference
        /
        random_closeness_sd
    )



    # 16. Z-SCORE → NORMAL-DISTRIBUTION STATISTICS


    (
        bet_abs_z,
        bet_normal_cdf,
        bet_normal_sf,
        p_betweenness
    ) = z_to_normal_statistics(
        betweenness_z
    )


    (
        close_abs_z,
        close_normal_cdf,
        close_normal_sf,
        p_closeness
    ) = z_to_normal_statistics(
        closeness_z
    )



    # 17. NORMALITY DIAGNOSTIC


    if len(random_betweenness) >= 20:

        bet_normality_stat, bet_normality_p = (
            normaltest(
                random_betweenness
            )
        )

    else:

        bet_normality_stat = np.nan

        bet_normality_p = np.nan


    if len(random_closeness) >= 20:

        close_normality_stat, close_normality_p = (
            normaltest(
                random_closeness
            )
        )

    else:

        close_normality_stat = np.nan

        close_normality_p = np.nan



    # 18. PRINT STATISTICAL SUMMARY


    print("\n" + "-" * 70)

    print(
        f"{cluster_name}"
    )

    print("-" * 70)



    # BETWEENNESS


    print("\nBETWEENNESS")

    print(
        f"Real mean          : "
        f"{real_betweenness:.10g}"
    )

    print(
        f"Random mean        : "
        f"{random_betweenness_mean:.10g}"
    )

    print(
        f"Random SD          : "
        f"{random_betweenness_sd:.10g}"
    )

    print(
        f"Difference         : "
        f"{betweenness_difference:.10g}"
    )

    print(
        f"Z-score            : "
        f"{betweenness_z:.10g}"
    )

    print(
        f"|Z|-score          : "
        f"{bet_abs_z:.10g}"
    )

    print(
        f"Normal CDF         : "
        f"{bet_normal_cdf:.10g}"
    )

    print(
        f"Normal SF          : "
        f"{bet_normal_sf:.10g}"
    )

    print(
        f"Two-sided normal p : "
        f"{p_betweenness:.10g}"
    )

    print(
        f"Normality test p   : "
        f"{bet_normality_p:.10g}"
    )



    # CLOSENESS


    print("\nCLOSENESS")

    print(
        f"Real mean          : "
        f"{real_closeness:.10g}"
    )

    print(
        f"Random mean        : "
        f"{random_closeness_mean:.10g}"
    )

    print(
        f"Random SD          : "
        f"{random_closeness_sd:.10g}"
    )

    print(
        f"Difference         : "
        f"{closeness_difference:.10g}"
    )

    print(
        f"Z-score            : "
        f"{closeness_z:.10g}"
    )

    print(
        f"|Z|-score          : "
        f"{close_abs_z:.10g}"
    )

    print(
        f"Normal CDF         : "
        f"{close_normal_cdf:.10g}"
    )

    print(
        f"Normal SF          : "
        f"{close_normal_sf:.10g}"
    )

    print(
        f"Two-sided normal p : "
        f"{p_closeness:.10g}"
    )

    print(
        f"Normality test p   : "
        f"{close_normality_p:.10g}"
    )



    # 19. STORE RESULTS


    all_results[cluster_name] = {

        "real_betweenness":
            real_betweenness,

        "random_betweenness":
            random_betweenness,

        "random_betweenness_mean":
            random_betweenness_mean,

        "random_betweenness_sd":
            random_betweenness_sd,

        "betweenness_difference":
            betweenness_difference,

        "betweenness_z":
            betweenness_z,

        "bet_abs_z":
            bet_abs_z,

        "bet_normal_cdf":
            bet_normal_cdf,

        "bet_normal_sf":
            bet_normal_sf,

        "p_betweenness":
            p_betweenness,

        "bet_normality_p":
            bet_normality_p,

        "real_closeness":
            real_closeness,

        "random_closeness":
            random_closeness,

        "random_closeness_mean":
            random_closeness_mean,

        "random_closeness_sd":
            random_closeness_sd,

        "closeness_difference":
            closeness_difference,

        "closeness_z":
            closeness_z,

        "close_abs_z":
            close_abs_z,

        "close_normal_cdf":
            close_normal_cdf,

        "close_normal_sf":
            close_normal_sf,

        "p_closeness":
            p_closeness,

        "close_normality_p":
            close_normality_p,

        "successful_networks":
            successful_networks
    }



    # 20. SUMMARY ROW — BETWEENNESS


    summary_rows.append({

        "Cluster":
            cluster_name,

        "Centrality":
            "Betweenness",

        "Real_Network_Mean":
            real_betweenness,

        "Random_Network_Mean":
            random_betweenness_mean,

        "Random_Network_SD":
            random_betweenness_sd,

        "Difference_in_Means":
            betweenness_difference,

        "Z_Score":
            betweenness_z,

        "Absolute_Z_Score":
            bet_abs_z,

        "Standard_Normal_CDF":
            bet_normal_cdf,

        "Standard_Normal_SF":
            bet_normal_sf,

        "Normal_Theory_P_Value":
            p_betweenness,

        "Null_Distribution_Normality_P":
            bet_normality_p,

        "N_Random_Networks":
            successful_networks
    })



    # 21. SUMMARY ROW — CLOSENESS


    summary_rows.append({

        "Cluster":
            cluster_name,

        "Centrality":
            "Closeness",

        "Real_Network_Mean":
            real_closeness,

        "Random_Network_Mean":
            random_closeness_mean,

        "Random_Network_SD":
            random_closeness_sd,

        "Difference_in_Means":
            closeness_difference,

        "Z_Score":
            closeness_z,

        "Absolute_Z_Score":
            close_abs_z,

        "Standard_Normal_CDF":
            close_normal_cdf,

        "Standard_Normal_SF":
            close_normal_sf,

        "Normal_Theory_P_Value":
            p_closeness,

        "Null_Distribution_Normality_P":
            close_normality_p,

        "N_Random_Networks":
            successful_networks
    })



# GLOBAL Y-AXIS LIMIT FOR ALL PLOTS


all_y_values = []


for cluster_name in all_results:

    result = all_results[
        cluster_name
    ]



    # Betweenness random values


    all_y_values.extend(
        result["random_betweenness"]
    )



    # Betweenness real value


    all_y_values.append(
        result["real_betweenness"]
    )



    # Closeness random values


    all_y_values.extend(
        result["random_closeness"]
    )



    # Closeness real value


    all_y_values.append(
        result["real_closeness"]
    )



# Remove non-finite values


all_y_values = np.asarray(
    all_y_values,
    dtype=float
)

all_y_values = all_y_values[
    np.isfinite(all_y_values)
]


if len(all_y_values) == 0:

    raise RuntimeError(
        "No valid y-axis values were found "
        "for global y-axis calculation."
    )



# Maximum value across all clusters
# and both centralities


maximum_y_value = np.max(
    all_y_values
)



# Add 5% headroom


global_y_max = (
    maximum_y_value * 1.05
)


global_y_min = 0


print("\n" + "=" * 75)

print("GLOBAL Y-AXIS LIMIT")

print("=" * 75)

print(
    f"Maximum recorded y-value : "
    f"{maximum_y_value:.10g}"
)

print(
    f"Global y-axis minimum     : "
    f"{global_y_min:.10g}"
)

print(
    f"Global y-axis maximum     : "
    f"{global_y_max:.10g}"
)



# 22. SAVE SUMMARY CSV


summary = pd.DataFrame(
    summary_rows
)


summary_file = os.path.join(
    OUTPUT_DIR,
    "all_clusters_degree_preserving_Zscore_summary.csv"
)


summary.to_csv(
    summary_file,
    index=False
)


print("\n" + "=" * 75)

print("COMBINED Z-SCORE SUMMARY CSV")

print("=" * 75)

print(summary_file)



# 23. COMBINED DOT / STRIP PLOT


fig, axes = plt.subplots(
    5,
    2,
    figsize=(8, 12)
)


plot_rng = np.random.default_rng(
    12345
)


for row, cluster_name in enumerate(
    CLUSTERS.keys()
):

    result = all_results[
        cluster_name
    ]


    bottom_row = (
        row == len(CLUSTERS) - 1
    )



    # BETWEENNESS


    ax = axes[row, 0]


    random_values = (
        result["random_betweenness"]
    )


    real_value = (
        result["real_betweenness"]
    )


    x_jitter = plot_rng.uniform(
        -0.08,
        0.08,
        size=len(random_values)
    )


    ax.scatter(
        np.ones(
            len(random_values)
        ) + x_jitter,
        random_values,
        s=10,
        alpha=0.40,
        color=RANDOM_COLOR,
        edgecolors="none"
    )


    ax.scatter(
        1,
        real_value,
        s=40,
        color=REAL_COLOR,
        edgecolor="black",
        linewidth=1.0,
        zorder=5
    )


    set_ylabel(
        ax,
        "Mean Betweenness"
    )


    set_subplot_title(
        ax,
        cluster_name,
        "Mean Betweenness Centrality"
    )


    add_statistics(
        ax,
        result["betweenness_difference"],
        result["betweenness_z"],
        result["p_betweenness"]
    )


    format_axis(
        ax,
        bottom_row=bottom_row
    )



    # GLOBAL Y-AXIS


    ax.set_ylim(
        global_y_min,
        global_y_max
    )



    # CLOSENESS


    ax = axes[row, 1]


    random_values = (
        result["random_closeness"]
    )


    real_value = (
        result["real_closeness"]
    )


    x_jitter = plot_rng.uniform(
        -0.08,
        0.08,
        size=len(random_values)
    )


    ax.scatter(
        np.ones(
            len(random_values)
        ) + x_jitter,
        random_values,
        s=10,
        alpha=0.40,
        color=RANDOM_COLOR,
        edgecolors="none"
    )


    ax.scatter(
        1,
        real_value,
        s=40,
        color=REAL_COLOR,
        edgecolor="black",
        linewidth=1.0,
        zorder=5
    )


    set_ylabel(
        ax,
        "Mean Closeness"
    )


    set_subplot_title(
        ax,
        cluster_name,
        "Mean Closeness Centrality"
    )


    add_statistics(
        ax,
        result["closeness_difference"],
        result["closeness_z"],
        result["p_closeness"]
    )


    format_axis(
        ax,
        bottom_row=bottom_row
    )



    # GLOBAL Y-AXIS


    ax.set_ylim(
        global_y_min,
        global_y_max
    )



# FIGURE TITLE


fig.suptitle(
    "Degree-Preserving Randomization",
    fontsize=14,
    fontname="Times New Roman",
    fontweight="bold",
    color="black",
    y=0.995
)



# LEGEND


legend = fig.legend(
    handles=[
        real_network_legend
    ],
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        0.978
    ),
    frameon=True,
    fontsize=8
)


for text in legend.get_texts():

    text.set_fontname(
        "Times New Roman"
    )

    text.set_fontweight(
        "bold"
    )

    text.set_color(
        "black"
    )



# LAYOUT


fig.tight_layout(
    rect=[
        0,
        0,
        1,
        0.965
    ]
)



# SAVE DOT PLOT


dot_png = os.path.join(
    OUTPUT_DIR,
    "all_clusters_degree_preserving_Zscore_dot_plot.png"
)


dot_pdf = os.path.join(
    OUTPUT_DIR,
    "all_clusters_degree_preserving_Zscore_dot_plot.pdf"
)


fig.savefig(
    dot_png,
    dpi=600,
    bbox_inches="tight"
)


fig.savefig(
    dot_pdf,
    bbox_inches="tight"
)


print("\nCombined dot plot saved:")

print(dot_png)

print(dot_pdf)


plt.show()

plt.close(fig)



# 24. COMBINED BOXPLOT


fig, axes = plt.subplots(
    5,
    2,
    figsize=(8, 12)
)


for row, cluster_name in enumerate(
    CLUSTERS.keys()
):

    result = all_results[
        cluster_name
    ]


    bottom_row = (
        row == len(CLUSTERS) - 1
    )



    # BETWEENNESS


    ax = axes[row, 0]


    ax.boxplot(
        result[
            "random_betweenness"
        ],
        positions=[1],
        widths=0.45,
        patch_artist=True,
        boxprops=dict(
            facecolor=RANDOM_COLOR,
            edgecolor="black",
            linewidth=1.0,
            alpha=0.65
        ),
        medianprops=dict(
            color="black",
            linewidth=2
        ),
        whiskerprops=dict(
            color="black",
            linewidth=1.2
        ),
        capprops=dict(
            color="black",
            linewidth=1.2
        ),
        flierprops=dict(
            marker="o",
            markerfacecolor=RANDOM_COLOR,
            markeredgecolor="black",
            markersize=3,
            alpha=0.3
        )
    )


    ax.scatter(
        1,
        result[
            "real_betweenness"
        ],
        s=40,
        color=REAL_COLOR,
        edgecolor="black",
        linewidth=1.0,
        zorder=5
    )


    set_ylabel(
        ax,
        "Mean Betweenness"
    )


    set_subplot_title(
        ax,
        cluster_name,
        "Mean Betweenness Centrality"
    )


    add_statistics(
        ax,
        result["betweenness_difference"],
        result["betweenness_z"],
        result["p_betweenness"]
    )


    format_axis(
        ax,
        bottom_row=bottom_row
    )



    # GLOBAL Y-AXIS


    ax.set_ylim(
        global_y_min,
        global_y_max
    )



    # CLOSENESS


    ax = axes[row, 1]


    ax.boxplot(
        result[
            "random_closeness"
        ],
        positions=[1],
        widths=0.45,
        patch_artist=True,
        boxprops=dict(
            facecolor=RANDOM_COLOR,
            edgecolor="black",
            linewidth=1.0,
            alpha=0.65
        ),
        medianprops=dict(
            color="black",
            linewidth=2
        ),
        whiskerprops=dict(
            color="black",
            linewidth=1.2
        ),
        capprops=dict(
            color="black",
            linewidth=1.2
        ),
        flierprops=dict(
            marker="o",
            markerfacecolor=RANDOM_COLOR,
            markeredgecolor="black",
            markersize=3,
            alpha=0.3
        )
    )


    ax.scatter(
        1,
        result[
            "real_closeness"
        ],
        s=40,
        color=REAL_COLOR,
        edgecolor="black",
        linewidth=1.0,
        zorder=5
    )


    set_ylabel(
        ax,
        "Mean Closeness"
    )


    set_subplot_title(
        ax,
        cluster_name,
        "Mean Closeness Centrality"
    )


    add_statistics(
        ax,
        result["closeness_difference"],
        result["closeness_z"],
        result["p_closeness"]
    )


    format_axis(
        ax,
        bottom_row=bottom_row
    )



    # GLOBAL Y-AXIS


    ax.set_ylim(
        global_y_min,
        global_y_max
    )



# FIGURE TITLE


fig.suptitle(
    "Degree-Preserving Randomization",
    fontsize=14,
    fontname="Times New Roman",
    fontweight="bold",
    color="black",
    y=0.995
)



# LEGEND


legend = fig.legend(
    handles=[
        real_network_legend
    ],
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        0.978
    ),
    frameon=True,
    fontsize=8
)


for text in legend.get_texts():

    text.set_fontname(
        "Times New Roman"
    )

    text.set_fontweight(
        "bold"
    )

    text.set_color(
        "black"
    )



# LAYOUT


fig.tight_layout(
    rect=[
        0,
        0,
        1,
        0.965
    ]
)



# SAVE BOXPLOT


box_png = os.path.join(
    OUTPUT_DIR,
    "all_clusters_degree_preserving_Zscore_boxplot.png"
)


box_pdf = os.path.join(
    OUTPUT_DIR,
    "all_clusters_degree_preserving_Zscore_boxplot.pdf"
)


fig.savefig(
    box_png,
    dpi=600,
    bbox_inches="tight"
)


fig.savefig(
    box_pdf,
    bbox_inches="tight"
)


print("\nCombined boxplot saved:")

print(box_png)

print(box_pdf)


plt.show()

plt.close(fig)



# 25. COMBINED VIOLIN PLOT


fig, axes = plt.subplots(
    5,
    2,
    figsize=(8, 12)
)


for row, cluster_name in enumerate(
    CLUSTERS.keys()
):

    result = all_results[
        cluster_name
    ]


    bottom_row = (
        row == len(CLUSTERS) - 1
    )



    # BETWEENNESS


    ax = axes[row, 0]


    parts = ax.violinplot(
        result[
            "random_betweenness"
        ],
        positions=[1],
        widths=0.65,
        showmeans=False,
        showmedians=True,
        showextrema=True
    )


    for body in parts["bodies"]:

        body.set_facecolor(
            RANDOM_COLOR
        )

        body.set_edgecolor(
            "black"
        )

        body.set_linewidth(
            1.0
        )

        body.set_alpha(
            0.65
        )


    if "cmedians" in parts:

        parts[
            "cmedians"
        ].set_color(
            "black"
        )

        parts[
            "cmedians"
        ].set_linewidth(
            2
        )


    if "cmins" in parts:

        parts[
            "cmins"
        ].set_color(
            "black"
        )


    if "cmaxes" in parts:

        parts[
            "cmaxes"
        ].set_color(
            "black"
        )


    if "cbars" in parts:

        parts[
            "cbars"
        ].set_color(
            "black"
        )


    ax.scatter(
        1,
        result[
            "real_betweenness"
        ],
        s=40,
        color=REAL_COLOR,
        edgecolor="black",
        linewidth=1.0,
        zorder=5
    )


    set_ylabel(
        ax,
        "Mean Betweenness"
    )


    set_subplot_title(
        ax,
        cluster_name,
        "Mean Betweenness Centrality"
    )


    add_statistics(
        ax,
        result["betweenness_difference"],
        result["betweenness_z"],
        result["p_betweenness"]
    )


    format_axis(
        ax,
        bottom_row=bottom_row
    )



    # GLOBAL Y-AXIS


    ax.set_ylim(
        global_y_min,
        global_y_max
    )



    # CLOSENESS


    ax = axes[row, 1]


    parts = ax.violinplot(
        result[
            "random_closeness"
        ],
        positions=[1],
        widths=0.65,
        showmeans=False,
        showmedians=True,
        showextrema=True
    )


    for body in parts["bodies"]:

        body.set_facecolor(
            RANDOM_COLOR
        )

        body.set_edgecolor(
            "black"
        )

        body.set_linewidth(
            1.0
        )

        body.set_alpha(
            0.65
        )


    if "cmedians" in parts:

        parts[
            "cmedians"
        ].set_color(
            "black"
        )

        parts[
            "cmedians"
        ].set_linewidth(
            2
        )


    if "cmins" in parts:

        parts[
            "cmins"
        ].set_color(
            "black"
        )


    if "cmaxes" in parts:

        parts[
            "cmaxes"
        ].set_color(
            "black"
        )


    if "cbars" in parts:

        parts[
            "cbars"
        ].set_color(
            "black"
        )


    ax.scatter(
        1,
        result[
            "real_closeness"
        ],
        s=40,
        color=REAL_COLOR,
        edgecolor="black",
        linewidth=1.0,
        zorder=5
    )


    set_ylabel(
        ax,
        "Mean Closeness"
    )


    set_subplot_title(
        ax,
        cluster_name,
        "Mean Closeness Centrality"
    )


    add_statistics(
        ax,
        result["closeness_difference"],
        result["closeness_z"],
        result["p_closeness"]
    )


    format_axis(
        ax,
        bottom_row=bottom_row
    )



    # GLOBAL Y-AXIS


    ax.set_ylim(
        global_y_min,
        global_y_max
    )



# FIGURE TITLE


fig.suptitle(
    "Degree-Preserving Randomization",
    fontsize=14,
    fontname="Times New Roman",
    fontweight="bold",
    color="black",
    y=0.995
)



# LEGEND


legend = fig.legend(
    handles=[
        real_network_legend
    ],
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        0.978
    ),
    frameon=True,
    fontsize=8
)


for text in legend.get_texts():

    text.set_fontname(
        "Times New Roman"
    )

    text.set_fontweight(
        "bold"
    )

    text.set_color(
        "black"
    )



# LAYOUT


fig.tight_layout(
    rect=[
        0,
        0,
        1,
        0.965
    ]
)



# SAVE VIOLIN PLOT


violin_png = os.path.join(
    OUTPUT_DIR,
    "all_clusters_degree_preserving_Zscore_violin_plot.png"
)


violin_pdf = os.path.join(
    OUTPUT_DIR,
    "all_clusters_degree_preserving_Zscore_violin_plot.pdf"
)


fig.savefig(
    violin_png,
    dpi=600,
    bbox_inches="tight"
)


fig.savefig(
    violin_pdf,
    bbox_inches="tight"
)


print("\nCombined violin plot saved:")

print(violin_png)

print(violin_pdf)


plt.show()

plt.close(fig)



# 26. FINAL MESSAGE


print("\n" + "=" * 75)

print("ALL CLUSTERS PROCESSED")

print("=" * 75)


print("\nSummary CSV:")

print(summary_file)


print("\nCombined dot plot:")

print(dot_png)

print(dot_pdf)


print("\nCombined boxplot:")

print(box_png)

print(box_pdf)


print("\nCombined violin plot:")

print(violin_png)

print(violin_pdf)


print("\nAnalysis completed successfully.")