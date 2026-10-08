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

# 4. Calculate variance of each gene
gene_variance = expression.var(axis=1)

# 5. Select the 2000 most variable genes
top_genes = gene_variance.nlargest(2000).index

expression_top = expression.loc[top_genes]

print("Genes used for PCA:", expression_top.shape[0])

# 6. Transpose:
# rows = samples
# columns = genes
X = expression_top.T

print("Samples used for PCA:", X.shape[0])

# 7. Run PCA
pca = PCA(n_components=2)

pca_result = pca.fit_transform(X)

# 8. Create PCA dataframe
pca_df = pd.DataFrame(
    pca_result,
    columns=["PC1", "PC2"],
    index=X.index
)

pca_df["expression_sample_name"] = pca_df.index

# 9. Add metadata
pca_df = pca_df.merge(
    metadata,
    on="expression_sample_name",
    how="left"
)

# 10. Print explained variance
print(
    "PC1 explained variance:",
    round(pca.explained_variance_ratio_[0] * 100, 2),
    "%"
)

print(
    "PC2 explained variance:",
    round(pca.explained_variance_ratio_[1] * 100, 2),
    "%"
)

# Create output folders
os.makedirs("results/pca", exist_ok=True)
os.makedirs("figures", exist_ok=True)

# Save PCA coordinates
pca_df.to_csv(
    "results/pca/GSE193258_pca_scores.csv",
    index=False
)

# 11. PCA colored/grouped by cell line
plt.figure(figsize=(8, 6))

for group, df in pca_df.groupby("cell_line"):
    plt.scatter(
        df["PC1"],
        df["PC2"],
        label=group,
        s=60
    )

plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)"
)
plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)"
)

plt.title("GSE193258 PCA by Cell Line")
plt.legend()
plt.tight_layout()

plt.savefig(
    "figures/GSE193258_PCA_cell_line.png",
    dpi=300
)

plt.show()


# 12. PCA grouped by treatment state
plt.figure(figsize=(8, 6))

for group, df in pca_df.groupby("treatment_state"):
    plt.scatter(
        df["PC1"],
        df["PC2"],
        label=group,
        s=60
    )

plt.xlabel(
    f"PC1 ({pca.explained_variance_ratio_[0] * 100:.1f}%)"
)
plt.ylabel(
    f"PC2 ({pca.explained_variance_ratio_[1] * 100:.1f}%)"
)

plt.title("GSE193258 PCA by Treatment State")
plt.legend()
plt.tight_layout()

plt.savefig(
    "figures/GSE193258_PCA_treatment_state.png",
    dpi=300
)

plt.show()