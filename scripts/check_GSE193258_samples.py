import pandas as pd

# 1. Read metadata
metadata = pd.read_csv(
    "data/metadata/GSE193258_metadata.csv"
)

# 2. Read only the header of the expression file
expression = pd.read_csv(
    "data/raw/GSE193258/GSE193258_RNAseq_estimated_counts.tsv.gz",
    sep="\t",
    compression="gzip",
    nrows=0
)

# 3. Get sample names
metadata_samples = set(metadata["expression_sample_name"])
expression_columns = set(expression.columns)

# 4. Find matches
matched_samples = metadata_samples.intersection(expression_columns)

# Samples present in metadata but not expression
missing_in_expression = metadata_samples - expression_columns

# Columns present in expression but not metadata
extra_in_expression = expression_columns - metadata_samples

# 5. Print results
print("Metadata samples:", len(metadata_samples))
print("Expression columns:", len(expression_columns))
print("Matched samples:", len(matched_samples))

print("\nMissing in expression:")
print(missing_in_expression)

print("\nExtra columns in expression:")
print(extra_in_expression)