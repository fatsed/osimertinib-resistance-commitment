import pandas as pd

file_path = (
    "results/generic_DTP/GSE255958/"
    "GSE255958_generic_DTP_profile.csv"
)

df = pd.read_csv(file_path)

genes = df["gene"].astype(str)

print("Total genes:", len(genes))

# Example IDs
print("\nFirst 30 gene IDs:")
print(genes.head(30).to_string(index=False))


# --------------------------------------------------
# Check duplicated-symbol pattern:
# ABC_ABC
# --------------------------------------------------

def repeated_symbol(gene):
    parts = gene.split("_")

    return (
        len(parts) == 2
        and parts[0] == parts[1]
    )


repeated = genes.apply(repeated_symbol)

print(
    "\nIDs with SYMBOL_SYMBOL pattern:",
    repeated.sum()
)

print(
    "Percentage:",
    round(
        100 * repeated.mean(),
        2
    ),
    "%"
)


# --------------------------------------------------
# Show IDs that contain underscore
# but are NOT simple SYMBOL_SYMBOL duplicates
# --------------------------------------------------

underscore = genes.str.contains("_", regex=False)

unusual = genes[
    underscore & ~repeated
]

print(
    "\nIDs containing underscore but not SYMBOL_SYMBOL:",
    len(unusual)
)

print("\nExamples:")
print(
    unusual.head(50).to_string(index=False)
)