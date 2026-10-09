import os
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA

# 1. Read expression data
expression = pd.read_csv(
    "data/raw/GSE193258/GSE193258_RNAseq_log2TPM_abundance.tsv.gz",
    sep="\t",
    compression="gzip"
)

# 2. Read metadata
metadata = pd.read_csv(
    "data/metadata/GSE193258_metadata.csv"
)

# 3. Set gene names as index
expression = expression.set_index("gene")

# Output folders
os.makedirs("figures/pca_by_cell_line", exist_ok=True)
os.makedirs("results/pca_by_cell_line", exist_ok=True)

# Cell lines to analyze
cell_lines = ["PC9", "H1975", "HCC827", "HCC2935"]

for cell_line in cell_lines:

    print("\n" + "=" * 50)
    print("Cell line:", cell_line)
    print("=" * 50)

    # 4. Select metadata for this cell line
    meta_sub = metadata[
        metadata["cell_line"] == cell_line
    ].copy()

    sample_names = meta_sub["expression_sample_name"].tolist()

    print("Number of samples:", len(sample_names))

    # 5. Select expression samples for this cell line
    expr_sub = expression[sample_names]

    # 6. Find the 2000 most variable genes within this cell line
    gene_variance = expr_sub.var(axis=1)

    top_genes = gene_variance.nlargest(2000).index

    expr_top = expr_sub.loc[top_genes]

    print("Genes used for PCA:", len(top_genes))

    # 7. Transpose:
    # rows = samples
    # columns = genes
    X = expr_top.T

    # 8. Run PCA
    pca = PCA(n_components=2)

    pca_result = pca.fit_transform(X)

    # 9. Create PCA dataframe
    pca_df = pd.DataFrame(
        pca_result,
        columns=["PC1", "PC2"],
        index=X.index
    )

    pca_df["expression_sample_name"] = pca_df.index

    # Add metadata
    pca_df = pca_df.merge(
        meta_sub,
        on="expression_sample_name",
        how="left"
    )

    # 10. Print explained variance
    pc1_var = pca.explained_variance_ratio_[0] * 100
    pc2_var = pca.explained_variance_ratio_[1] * 100

    print("PC1 explained variance:", round(pc1_var, 2), "%")
    print("PC2 explained variance:", round(pc2_var, 2), "%")

    # 11. Save PCA coordinates
    pca_df.to_csv(
        f"results/pca_by_cell_line/{cell_line}_PCA_scores.csv",
        index=False
    )

    # 12. Plot by treatment state
    plt.figure(figsize=(8, 6))

    for state, df in pca_df.groupby("treatment_state"):
        plt.scatter(
            df["PC1"],
            df["PC2"],
            label=state,
            s=70
        )

    plt.xlabel(f"PC1 ({pc1_var:.1f}%)")
    plt.ylabel(f"PC2 ({pc2_var:.1f}%)")

    plt.title(f"GSE193258 PCA - {cell_line}")

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        f"figures/pca_by_cell_line/{cell_line}_PCA.png",
        dpi=300
    )

    plt.close()

print("\nFinished.")