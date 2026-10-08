import pandas as pd

tpm = pd.read_csv(
    "data/raw/GSE193258/GSE193258_RNAseq_log2TPM_abundance.tsv.gz",
    sep="\t",
    compression="gzip"
)

print("Shape:", tpm.shape)
print()

print("First 5 rows:")
print(tpm.head())
print()

print("First 10 columns:")
print(tpm.columns[:10].tolist())
print()

print("Total missing values:")
print(tpm.isna().sum().sum())
print()

print("Duplicated gene IDs:")
print(tpm["gene"].duplicated().sum())
print()

sample_columns = tpm.columns[1:]

print("Minimum value:")
print(tpm[sample_columns].min().min())

print("Maximum value:")
print(tpm[sample_columns].max().max())

print()

print("Genes with zero in all samples:")
print((tpm[sample_columns].sum(axis=1) == 0).sum())