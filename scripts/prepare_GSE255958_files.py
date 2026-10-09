import os
import tarfile
import gzip
import pandas as pd

archive_path = "data/raw/GSE255958/GSE255958_RAW.tar"
extract_dir = "data/raw/GSE255958/extracted"

os.makedirs(extract_dir, exist_ok=True)

# --------------------------------------------------
# 1. Extract all processed expression files
# --------------------------------------------------

with tarfile.open(archive_path, "r") as tar:
    tar.extractall(extract_dir)

print("Archive extracted successfully.")


# --------------------------------------------------
# 2. Inspect one representative file
# --------------------------------------------------

example_file = os.path.join(
    extract_dir,
    "GSM8083374_PC9_D_1.genes.results.gz"
)

df = pd.read_csv(
    example_file,
    sep="\t",
    compression="gzip"
)

print("\nExample file:")
print(example_file)

print("\nShape:")
print(df.shape)

print("\nColumns:")
print(df.columns.tolist())

print("\nFirst 5 rows:")
print(df.head())