import os
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA

# -----------------------------
# 1. Read data
# -----------------------------

expression = pd.read_csv(
    "data/raw/GSE193258/GSE193258_RNAseq_log2TPM_abundance.tsv.gz",
    sep="\t",
    compression="gzip"
)

metadata = pd.read_csv(
    "data/metadata/GSE193258_metadata.csv"
)

expression = expression.set_index("gene")

cell_lines = ["PC9", "H1975", "HCC827", "HCC2935"]

os.makedirs("results/pca_distances", exist_ok=True)

all_results = []


# -----------------------------
# 2. Analyze each cell line
# -----------------------------

for cell_line in cell_lines:

    print("\n" + "=" * 60)
    print("Cell line:", cell_line)
    print("=" * 60)

    # Metadata for this cell line
    meta_sub = metadata[
        metadata["cell_line"] == cell_line
    ].copy()

    sample_names = meta_sub["expression_sample_name"].tolist()

    # Expression data for this cell line
    expr_sub = expression[sample_names]

    # -----------------------------
    # 3. Select top 2000 variable genes
    # -----------------------------

    gene_variance = expr_sub.var(axis=1)

    top_genes = gene_variance.nlargest(2000).index

    expr_top = expr_sub.loc[top_genes]

    # rows = samples
    # columns = genes
    X = expr_top.T

    # -----------------------------
    # 4. Run PCA with all possible components
    # -----------------------------

    pca = PCA()

    pca_result = pca.fit_transform(X)

    cumulative_variance = np.cumsum(
        pca.explained_variance_ratio_
    )

    # Find minimum number of PCs needed for >= 80% variance
    n_pcs = np.argmax(cumulative_variance >= 0.80) + 1

    explained = cumulative_variance[n_pcs - 1] * 100

    print("PCs needed for >=80% variance:", n_pcs)
    print("Variance explained:", round(explained, 2), "%")

    # Keep only selected PCs
    pca_selected = pca_result[:, :n_pcs]

    pc_columns = [
        f"PC{i+1}" for i in range(n_pcs)
    ]

    pca_df = pd.DataFrame(
        pca_selected,
        columns=pc_columns,
        index=X.index
    )

    pca_df["expression_sample_name"] = pca_df.index

    # Add treatment information
    pca_df = pca_df.merge(
        meta_sub[
            ["expression_sample_name", "treatment_state"]
        ],
        on="expression_sample_name",
        how="left"
    )

    # -----------------------------
    # 5. Calculate centroids
    # -----------------------------

    centroids = (
        pca_df
        .groupby("treatment_state")[pc_columns]
        .mean()
    )

    # Euclidean distance using all selected PCs
    def distance(state1, state2):

        point1 = centroids.loc[state1].values
        point2 = centroids.loc[state2].values

        return np.linalg.norm(point1 - point2)

    # -----------------------------
    # 6. Calculate distances
    # -----------------------------

    short_to_control = distance(
        "short_washout",
        "control"
    )

    short_to_dtp = distance(
        "short_washout",
        "DTP"
    )

    long_to_control = distance(
        "long_washout",
        "control"
    )

    long_to_dtp = distance(
        "long_washout",
        "DTP"
    )

    short_closer_to = (
        "Control"
        if short_to_control < short_to_dtp
        else "DTP"
    )

    long_closer_to = (
        "Control"
        if long_to_control < long_to_dtp
        else "DTP"
    )

    print("\nShort washout:")
    print("Distance to Control:", round(short_to_control, 2))
    print("Distance to DTP:", round(short_to_dtp, 2))
    print("Closer to:", short_closer_to)

    print("\nLong washout:")
    print("Distance to Control:", round(long_to_control, 2))
    print("Distance to DTP:", round(long_to_dtp, 2))
    print("Closer to:", long_closer_to)

    # -----------------------------
    # 7. Save results
    # -----------------------------

    all_results.append({
        "cell_line": cell_line,
        "n_pcs_for_80_percent": n_pcs,
        "variance_explained_percent": explained,
        "short_to_control": short_to_control,
        "short_to_DTP": short_to_dtp,
        "short_closer_to": short_closer_to,
        "long_to_control": long_to_control,
        "long_to_DTP": long_to_dtp,
        "long_closer_to": long_closer_to
    })


# -----------------------------
# 8. Save summary
# -----------------------------

results_df = pd.DataFrame(all_results)

results_df.to_csv(
    "results/pca_distances/GSE193258_multiPC_distances.csv",
    index=False
)

print("\n" + "=" * 60)
print("FINAL SUMMARY")
print("=" * 60)

print(results_df)