import os
import pandas as pd


# ============================================================
# Paths
# ============================================================

metadata_file = "data/metadata/GSE255958_core_metadata.csv"

expression_dir = "data/raw/GSE255958/extracted"

output_dir = "data/processed"

os.makedirs(
    output_dir,
    exist_ok=True
)


# ============================================================
# Read metadata
# ============================================================

metadata = pd.read_csv(metadata_file)

print("Number of metadata samples:", len(metadata))


# ============================================================
# Containers for expression data
# ============================================================

tpm_series = []
count_series = []

reference_genes = None


# ============================================================
# Read all 36 samples
# ============================================================

for _, row in metadata.iterrows():

    sample_id = row["sample_id"]
    filename = row["expression_file"]

    file_path = os.path.join(
        expression_dir,
        filename
    )

    df = pd.read_csv(
        file_path,
        sep="\t",
        compression="gzip"
    )

    # ----------------------------------------
    # Basic checks
    # ----------------------------------------

    if df["gene_id"].duplicated().any():
        raise ValueError(
            f"Duplicated gene IDs found in {sample_id}"
        )

    current_genes = df["gene_id"].tolist()

    # Use first sample as reference gene list
    if reference_genes is None:

        reference_genes = current_genes

        print(
            "Genes in reference sample:",
            len(reference_genes)
        )

    else:

        if current_genes != reference_genes:

            raise ValueError(
                f"Gene list/order differs in {sample_id}"
            )

    # ----------------------------------------
    # TPM
    # ----------------------------------------

    tpm = df.set_index("gene_id")["TPM"]

    tpm.name = sample_id

    tpm_series.append(tpm)

    # ----------------------------------------
    # Expected counts
    # ----------------------------------------

    counts = df.set_index("gene_id")[
        "expected_count"
    ]

    counts.name = sample_id

    count_series.append(counts)


# ============================================================
# Combine samples into matrices
# ============================================================

tpm_matrix = pd.concat(
    tpm_series,
    axis=1
)

count_matrix = pd.concat(
    count_series,
    axis=1
)


# ============================================================
# Quality checks
# ============================================================

print("\nTPM matrix shape:")
print(tpm_matrix.shape)

print("\nExpected-count matrix shape:")
print(count_matrix.shape)

print(
    "\nMissing values in TPM:",
    tpm_matrix.isna().sum().sum()
)

print(
    "Missing values in expected counts:",
    count_matrix.isna().sum().sum()
)

print(
    "\nDuplicate genes:",
    tpm_matrix.index.duplicated().sum()
)


# ============================================================
# All-zero genes
# ============================================================

tpm_all_zero = (
    tpm_matrix.sum(axis=1) == 0
).sum()

count_all_zero = (
    count_matrix.sum(axis=1) == 0
).sum()

print(
    "\nGenes with TPM = 0 in all 36 samples:",
    tpm_all_zero
)

print(
    "Genes with expected_count = 0 in all 36 samples:",
    count_all_zero
)


# ============================================================
# Save matrices
# ============================================================

tpm_file = os.path.join(
    output_dir,
    "GSE255958_core_TPM.tsv.gz"
)

count_file = os.path.join(
    output_dir,
    "GSE255958_core_expected_counts.tsv.gz"
)

tpm_matrix.to_csv(
    tpm_file,
    sep="\t",
    compression="gzip"
)

count_matrix.to_csv(
    count_file,
    sep="\t",
    compression="gzip"
)


# ============================================================
# Final message
# ============================================================

print("\nMatrices created successfully.")

print("\nTPM saved to:")
print(tpm_file)

print("\nExpected counts saved to:")
print(count_file)