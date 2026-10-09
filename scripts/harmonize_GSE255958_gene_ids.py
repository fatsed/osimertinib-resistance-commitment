import os
import pandas as pd


# ============================================================
# Paths
# ============================================================

input_file = (
    "results/generic_DTP/GSE255958/"
    "GSE255958_generic_DTP_profile.csv"
)

output_dir = "results/generic_DTP/GSE255958"

output_file = os.path.join(
    output_dir,
    "GSE255958_generic_DTP_profile_harmonized.csv"
)


# ============================================================
# Read data
# ============================================================

df = pd.read_csv(input_file)

print("Genes before harmonization:", len(df))


# ============================================================
# Harmonize gene IDs
#
# ABC_ABC -> ABC
# Anything else remains unchanged
# ============================================================

def harmonize_gene_id(gene):

    gene = str(gene)

    parts = gene.split("_")

    # Check whether the whole ID consists of
    # two identical halves.
    #
    # Examples:
    # PADI3_PADI3 -> PADI3
    # C4B_2_C4B_2 -> C4B_2

    if len(parts) % 2 == 0:

        half = len(parts) // 2

        first_half = parts[:half]
        second_half = parts[half:]

        if first_half == second_half:
            return "_".join(first_half)

    return gene

# Preserve original ID

df["gene_original"] = df["gene"]

df["gene"] = df["gene"].apply(
    harmonize_gene_id
)


# ============================================================
# Check changes
# ============================================================

changed = (
    df["gene"] != df["gene_original"]
)

print(
    "IDs harmonized:",
    changed.sum()
)

print(
    "IDs unchanged:",
    (~changed).sum()
)


# ============================================================
# Show unusual IDs
# ============================================================

unusual = df.loc[
    ~changed &
    df["gene_original"].str.contains(
        "_",
        regex=False
    ),
    [
        "gene_original",
        "gene"
    ]
]

print("\nUnusual underscore-containing IDs kept unchanged:")

print(
    unusual.to_string(
        index=False
    )
)


# ============================================================
# Check whether harmonization created duplicate gene symbols
# ============================================================

duplicate_mask = df["gene"].duplicated(
    keep=False
)

n_duplicate_rows = duplicate_mask.sum()

print(
    "\nRows involved in duplicated harmonized gene symbols:",
    n_duplicate_rows
)

if n_duplicate_rows > 0:

    print(
        "\nExamples of duplicated harmonized IDs:"
    )

    print(
        df.loc[
            duplicate_mask,
            [
                "gene_original",
                "gene"
            ]
        ].head(30).to_string(
            index=False
        )
    )


# ============================================================
# Save harmonized profile
# ============================================================

df.to_csv(
    output_file,
    index=False
)

print(
    "\nHarmonized profile saved to:"
)

print(output_file)


# ============================================================
# Re-create harmonized candidate sets
# ============================================================

core = df[
    df["generic_DTP_core"]
].copy()

recurrent = df[
    df["generic_DTP_2plus"]
].copy()


core.to_csv(
    os.path.join(
        output_dir,
        "GSE255958_generic_DTP_core_3of3_harmonized.csv"
    ),
    index=False
)

recurrent.to_csv(
    os.path.join(
        output_dir,
        "GSE255958_generic_DTP_recurrent_2plus_harmonized.csv"
    ),
    index=False
)


print(
    "\nCore 3/3 genes:",
    len(core)
)

print(
    "Recurrent >=2/3 genes:",
    len(recurrent)
)