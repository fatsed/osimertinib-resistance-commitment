import os
import numpy as np
import pandas as pd


# ============================================================
# Settings
# ============================================================

cell_lines = [
    "PC9",
    "H3122",
    "H358"
]

base_dir = "results/differential_expression/GSE255958"

output_dir = "results/generic_DTP/GSE255958"

os.makedirs(
    output_dir,
    exist_ok=True
)


# ============================================================
# Read and prepare one cell line
# ============================================================

def prepare_cell_line(cell_line):

    file_path = os.path.join(
        base_dir,
        cell_line,
        "DTP_vs_Control.csv"
    )

    df = pd.read_csv(file_path)

    df = df[
        [
            "gene",
            "logFC",
            "adj.P.Val"
        ]
    ].copy()

    df = df.rename(
        columns={
            "logFC": f"logFC_{cell_line}",
            "adj.P.Val": f"FDR_{cell_line}"
        }
    )

    return df


# ============================================================
# Read all three models
# ============================================================

tables = [
    prepare_cell_line(cell_line)
    for cell_line in cell_lines
]


# ============================================================
# Outer merge
#
# Some genes may have been filtered out in one model,
# so absence is treated as "not tested", not as "no change".
# ============================================================

combined = tables[0]

for table in tables[1:]:

    combined = combined.merge(
        table,
        on="gene",
        how="outer"
    )


# ============================================================
# Model-level evidence
# ============================================================

for cell_line in cell_lines:

    logfc_col = f"logFC_{cell_line}"
    fdr_col = f"FDR_{cell_line}"

    # Was this gene tested after expression filtering?
    combined[
        f"tested_{cell_line}"
    ] = (
        combined[logfc_col].notna() &
        combined[fdr_col].notna()
    )

    # DTP-responsive definition
    combined[
        f"DTP_responsive_{cell_line}"
    ] = (
        combined[f"tested_{cell_line}"] &
        (combined[fdr_col] < 0.05) &
        (combined[logfc_col].abs() >= 0.5)
    )


# ============================================================
# Count how many models tested/responded
# ============================================================

tested_cols = [
    f"tested_{cell_line}"
    for cell_line in cell_lines
]

responsive_cols = [
    f"DTP_responsive_{cell_line}"
    for cell_line in cell_lines
]

combined["n_models_tested"] = combined[
    tested_cols
].sum(axis=1)

combined["n_DTP_responsive"] = combined[
    responsive_cols
].sum(axis=1)


# ============================================================
# Direction consistency
# ============================================================

def direction_summary(row):

    directions = []

    for cell_line in cell_lines:

        responsive = row[
            f"DTP_responsive_{cell_line}"
        ]

        if responsive:

            logfc = row[
                f"logFC_{cell_line}"
            ]

            directions.append(
                np.sign(logfc)
            )

    if len(directions) == 0:
        return "none"

    if all(x > 0 for x in directions):
        return "up"

    if all(x < 0 for x in directions):
        return "down"

    return "mixed"


combined["responsive_direction"] = combined.apply(
    direction_summary,
    axis=1
)


# ============================================================
# Strongest generic DTP evidence
#
# Core:
# - tested in all 3 models
# - DTP-responsive in all 3
# - same direction in all 3
# ============================================================

combined["generic_DTP_core"] = (
    (combined["n_models_tested"] == 3) &
    (combined["n_DTP_responsive"] == 3) &
    combined["responsive_direction"].isin(
        ["up", "down"]
    )
)


# ============================================================
# Broader recurrent evidence
#
# At least 2/3 models responsive
# in the same direction
# ============================================================

combined["generic_DTP_2plus"] = (
    (combined["n_DTP_responsive"] >= 2) &
    combined["responsive_direction"].isin(
        ["up", "down"]
    )
)


# ============================================================
# Effect-size summaries
# ============================================================

logfc_cols = [
    f"logFC_{cell_line}"
    for cell_line in cell_lines
]

combined["median_logFC"] = combined[
    logfc_cols
].median(
    axis=1,
    skipna=True
)

combined["median_abs_logFC"] = combined[
    logfc_cols
].abs().median(
    axis=1,
    skipna=True
)

combined["min_abs_logFC"] = combined[
    logfc_cols
].abs().min(
    axis=1,
    skipna=True
)


# ============================================================
# Sort strongest evidence first
# ============================================================

combined = combined.sort_values(
    by=[
        "generic_DTP_core",
        "n_DTP_responsive",
        "min_abs_logFC"
    ],
    ascending=[
        False,
        False,
        False
    ]
)


# ============================================================
# Save full profile
# ============================================================

combined.to_csv(
    os.path.join(
        output_dir,
        "GSE255958_generic_DTP_profile.csv"
    ),
    index=False
)


# ============================================================
# Save core 3/3 generic DTP genes
# ============================================================

core = combined[
    combined["generic_DTP_core"]
].copy()

core.to_csv(
    os.path.join(
        output_dir,
        "GSE255958_generic_DTP_core_3of3.csv"
    ),
    index=False
)


# ============================================================
# Save recurrent 2+/3 genes
# ============================================================

recurrent = combined[
    combined["generic_DTP_2plus"]
].copy()

recurrent.to_csv(
    os.path.join(
        output_dir,
        "GSE255958_generic_DTP_recurrent_2plus.csv"
    ),
    index=False
)


# ============================================================
# Summary
# ============================================================

print("\nGeneric DTP profile created.")

print(
    "\nTotal genes in combined profile:",
    len(combined)
)

print(
    "Genes tested in all 3 models:",
    (combined["n_models_tested"] == 3).sum()
)

print(
    "\nDTP-responsive in 3/3 models:",
    (combined["n_DTP_responsive"] == 3).sum()
)

print(
    "DTP-responsive in 2/3 models:",
    (combined["n_DTP_responsive"] == 2).sum()
)

print(
    "DTP-responsive in 1/3 models:",
    (combined["n_DTP_responsive"] == 1).sum()
)

print(
    "\nCore generic DTP genes (3/3 + same direction):",
    len(core)
)

print(
    "Recurrent generic DTP genes (>=2/3 + same direction):",
    len(recurrent)
)

print("\nCore direction counts:")

print(
    core["responsive_direction"].value_counts()
)

print("\nTop 20 core generic DTP genes:")

print(
    core[
        [
            "gene",
            "responsive_direction",
            "logFC_PC9",
            "logFC_H3122",
            "logFC_H358",
            "min_abs_logFC"
        ]
    ].head(20).to_string(index=False)
)