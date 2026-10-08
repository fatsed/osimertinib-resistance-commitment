import pandas as pd

# 1. Read metadata
metadata = pd.read_csv(
    "data/metadata/GSE193258_metadata.csv"
)

metadata_samples = set(metadata["expression_sample_name"])


def check_file(file_path, label):
    # Read only column names
    data = pd.read_csv(
        file_path,
        sep="\t",
        compression="gzip",
        nrows=0
    )

    expression_columns = set(data.columns)

    matched_samples = metadata_samples.intersection(expression_columns)
    missing_in_expression = metadata_samples - expression_columns
    extra_in_expression = expression_columns - metadata_samples

    print("=" * 50)
    print(label)
    print("=" * 50)

    print("Metadata samples:", len(metadata_samples))
    print("Expression columns:", len(expression_columns))
    print("Matched samples:", len(matched_samples))

    print("\nMissing in expression:")
    print(missing_in_expression)

    print("\nExtra columns in expression:")
    print(extra_in_expression)

    print()


# 2. Check estimated counts
check_file(
    "data/raw/GSE193258/GSE193258_RNAseq_estimated_counts.tsv.gz",
    "Estimated counts"
)

# 3. Check log2TPM
check_file(
    "data/raw/GSE193258/GSE193258_RNAseq_log2TPM_abundance.tsv.gz",
    "log2TPM"
)