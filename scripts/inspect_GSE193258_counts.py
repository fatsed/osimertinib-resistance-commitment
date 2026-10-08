import pandas as pd

# 1. Read the count matrix
counts = pd.read_csv(
    "data/raw/GSE193258/GSE193258_RNAseq_estimated_counts.tsv.gz",
    sep="\t",
    compression="gzip"
)

# 2. Basic information
print("Shape:", counts.shape)
print()

print("First 5 rows:")
print(counts.head())
print()

print("First 10 column names:")
print(counts.columns[:10].tolist())
print()

# 3. Missing values
print("Total missing values:")
print(counts.isna().sum().sum())
print()

# 4. Duplicate genes
print("Duplicated gene IDs:")
print(counts["gene"].duplicated().sum())
print()

# 5. Check expression values
sample_columns = counts.columns[1:]

print("Minimum count:")
print(counts[sample_columns].min().min())

print("Maximum count:")
print(counts[sample_columns].max().max())

print()

# 6. How many genes have zero counts in all samples?
all_zero_genes = (counts[sample_columns].sum(axis=1) == 0).sum()

print("Genes with zero counts in all samples:")
print(all_zero_genes)

# 7. Remove genes with zero expression in all samples
filtered_counts = counts[
    counts[sample_columns].sum(axis=1) > 0
].copy()

print("Shape after removing all-zero genes:")
print(filtered_counts.shape)

# 8. Save filtered matrix
filtered_counts.to_csv(
    "data/processed/GSE193258_estimated_counts_no_allzero.tsv",
    sep="\t",
    index=False
)