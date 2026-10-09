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

results_dir = "results/pca_by_cell_line/GSE255958"
figures_dir = "figures/pca_by_cell_line/GSE255958"

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
# Treatment-state order
# ============================================================

state_order = [
    "control",
    "acute",
    "DTP",
    "resistant"
]


# ============================================================
# Run PCA separately for each cell line
# ============================================================

for cell_line in [
    "PC9",
    "H3122",
    "H358"
]:

    print("\n" + "=" * 50)
    print("Cell line:", cell_line)
    print("=" * 50)

    # --------------------------------------------------------
    # Select metadata
    # --------------------------------------------------------

    meta_sub = metadata[
        metadata["cell_line"] == cell_line
    ].copy()

    sample_ids = meta_sub[
        "sample_id"
    ].tolist()

    print(
        "Number of samples:",
        len(sample_ids)
    )

    print(
        "\nSamples per treatment state:"
    )

    print(
        meta_sub[
            "treatment_state"
        ].value_counts()
    )


    # --------------------------------------------------------
    # Select expression data
    # --------------------------------------------------------

    tpm_sub = tpm[
        sample_ids
    ].copy()


    # --------------------------------------------------------
    # Remove genes with zero TPM
    # in all samples of this cell line
    # --------------------------------------------------------

    tpm_sub = tpm_sub.loc[
        tpm_sub.sum(axis=1) > 0
    ]

    print(
        "\nGenes after zero filtering:",
        tpm_sub.shape[0]
    )


    # --------------------------------------------------------
    # log2(TPM + 1)
    # --------------------------------------------------------

    log_tpm = np.log2(
        tpm_sub + 1
    )


    # --------------------------------------------------------
    # Select top 2000 variable genes
    # --------------------------------------------------------

    gene_variance = log_tpm.var(
        axis=1
    )

    n_top = min(
        2000,
        len(gene_variance)
    )

    top_genes = gene_variance.nlargest(
        n_top
    ).index

    pca_input = log_tpm.loc[
        top_genes
    ].T


    # --------------------------------------------------------
    # PCA
    # --------------------------------------------------------

    pca = PCA(
        n_components=2
    )

    scores = pca.fit_transform(
        pca_input
    )

    pc1_var = (
        pca.explained_variance_ratio_[0]
        * 100
    )

    pc2_var = (
        pca.explained_variance_ratio_[1]
        * 100
    )

    print(
        f"PC1 explained variance: {pc1_var:.2f}%"
    )

    print(
        f"PC2 explained variance: {pc2_var:.2f}%"
    )


    # --------------------------------------------------------
    # Save PCA scores
    # --------------------------------------------------------

    pca_df = meta_sub.copy()

    pca_df["PC1"] = scores[:, 0]
    pca_df["PC2"] = scores[:, 1]

    pca_df.to_csv(
        os.path.join(
            results_dir,
            f"{cell_line}_PCA_scores.csv"
        ),
        index=False
    )


    # --------------------------------------------------------
    # Plot by treatment state
    # --------------------------------------------------------

    plt.figure(
        figsize=(8, 6)
    )

    for state in state_order:

        subset = pca_df[
            pca_df["treatment_state"] == state
        ]

        plt.scatter(
            subset["PC1"],
            subset["PC2"],
            label=state,
            s=80
        )

    plt.xlabel(
        f"PC1 ({pc1_var:.2f}%)"
    )

    plt.ylabel(
        f"PC2 ({pc2_var:.2f}%)"
    )

    plt.title(
        f"GSE255958 {cell_line} PCA"
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        os.path.join(
            figures_dir,
            f"{cell_line}_PCA.png"
        ),
        dpi=300
    )

    plt.close()


print(
    "\nAll cell-line PCA analyses completed successfully."
)