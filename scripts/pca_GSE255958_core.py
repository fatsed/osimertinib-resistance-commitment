import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from sklearn.decomposition import PCA


# ============================================================
# Paths
# ============================================================

tpm_file = "data/processed/GSE255958_core_TPM.tsv.gz"
metadata_file = "data/metadata/GSE255958_core_metadata.csv"

results_dir = "results/pca"
figures_dir = "figures"

os.makedirs(results_dir, exist_ok=True)
os.makedirs(figures_dir, exist_ok=True)


# ============================================================
# Read data
# ============================================================

tpm = pd.read_csv(
    tpm_file,
    sep="\t",
    index_col=0
)

metadata = pd.read_csv(metadata_file)

print("TPM shape:", tpm.shape)
print("Metadata samples:", len(metadata))


# ============================================================
# Check sample matching
# ============================================================

sample_ids = metadata["sample_id"].tolist()

missing_samples = set(sample_ids) - set(tpm.columns)
extra_samples = set(tpm.columns) - set(sample_ids)

print("\nMissing samples:", missing_samples)
print("Extra samples:", extra_samples)

assert len(missing_samples) == 0
assert len(extra_samples) == 0


# Reorder expression columns to match metadata
tpm = tpm[sample_ids]


# ============================================================
# Remove genes with zero TPM in all samples
# ============================================================

tpm = tpm.loc[
    tpm.sum(axis=1) > 0
]

print("\nGenes after removing all-zero genes:")
print(tpm.shape[0])


# ============================================================
# Log transformation
# ============================================================

log_tpm = np.log2(
    tpm + 1
)


# ============================================================
# Select top 2000 variable genes
# ============================================================

gene_variance = log_tpm.var(
    axis=1
)

top_genes = gene_variance.nlargest(
    2000
).index

pca_input = log_tpm.loc[
    top_genes
].T


# ============================================================
# PCA
# ============================================================

pca = PCA(
    n_components=2
)

scores = pca.fit_transform(
    pca_input
)

pc1_var = (
    pca.explained_variance_ratio_[0] * 100
)

pc2_var = (
    pca.explained_variance_ratio_[1] * 100
)

print(
    f"\nPC1 explained variance: {pc1_var:.2f}%"
)

print(
    f"PC2 explained variance: {pc2_var:.2f}%"
)


# ============================================================
# Build PCA result table
# ============================================================

pca_df = metadata.copy()

pca_df["PC1"] = scores[:, 0]
pca_df["PC2"] = scores[:, 1]

pca_df.to_csv(
    os.path.join(
        results_dir,
        "GSE255958_core_pca_scores.csv"
    ),
    index=False
)


# ============================================================
# PCA plot by cell line
# ============================================================

plt.figure(
    figsize=(8, 6)
)

for cell_line in metadata[
    "cell_line"
].unique():

    subset = pca_df[
        pca_df["cell_line"] == cell_line
    ]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=cell_line,
        s=70
    )

plt.xlabel(
    f"PC1 ({pc1_var:.2f}%)"
)

plt.ylabel(
    f"PC2 ({pc2_var:.2f}%)"
)

plt.title(
    "GSE255958 Core PCA - Cell Line"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    os.path.join(
        figures_dir,
        "GSE255958_core_PCA_cell_line.png"
    ),
    dpi=300
)

plt.close()


# ============================================================
# PCA plot by treatment state
# ============================================================

plt.figure(
    figsize=(8, 6)
)

for state in [
    "control",
    "acute",
    "DTP",
    "resistant"
]:

    subset = pca_df[
        pca_df["treatment_state"] == state
    ]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=state,
        s=70
    )

plt.xlabel(
    f"PC1 ({pc1_var:.2f}%)"
)

plt.ylabel(
    f"PC2 ({pc2_var:.2f}%)"
)

plt.title(
    "GSE255958 Core PCA - Treatment State"
)

plt.legend()

plt.tight_layout()

plt.savefig(
    os.path.join(
        figures_dir,
        "GSE255958_core_PCA_treatment_state.png"
    ),
    dpi=300
)

plt.close()


print(
    "\nPCA analysis completed successfully."
)