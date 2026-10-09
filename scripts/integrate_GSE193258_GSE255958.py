import os
import numpy as np
import pandas as pd


# ============================================================
# Paths
# ============================================================

persistence_file = (
    "results/persistence/"
    "GSE193258_cross_cellline_persistence_profile.csv"
)

generic_file = (
    "results/generic_DTP/GSE255958/"
    "GSE255958_generic_DTP_profile_harmonized.csv"
)

output_dir = "results/cross_dataset_integration"

os.makedirs(
    output_dir,
    exist_ok=True
)


# ============================================================
# Read data
# ============================================================

persistence = pd.read_csv(persistence_file)
generic = pd.read_csv(generic_file)

print(
    "GSE193258 genes:",
    len(persistence)
)

print(
    "GSE255958 genes:",
    len(generic)
)


# ============================================================
# GSE193258 direction among persistence-supporting cell lines
# ============================================================

cell_lines_193258 = [
    "H1975",
    "PC9",
    "HCC827",
    "HCC2935"
]


def persistence_direction(row):

    directions = []

    for cell_line in cell_lines_193258:

        support = row[
            f"Persistence_support_{cell_line}"
        ]

        logfc = row[
            f"DTP_vs_Control_logFC_{cell_line}"
        ]

        if support and pd.notna(logfc):

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


persistence[
    "persistence_direction"
] = persistence.apply(
    persistence_direction,
    axis=1
)


# ============================================================
# Direction-consistent persistence sets
# ============================================================

persistence[
    "persistent_core_4of4"
] = (
    (persistence["n_persistence_support"] == 4) &
    persistence["persistence_direction"].isin(
        ["up", "down"]
    )
)


persistence[
    "persistent_robust_3plus"
] = (
    (persistence["n_persistence_support"] >= 3) &
    persistence["persistence_direction"].isin(
        ["up", "down"]
    )
)


# ============================================================
# Keep useful generic-DTP columns
# ============================================================

generic_small = generic[
    [
        "gene",
        "generic_DTP_core",
        "generic_DTP_2plus",
        "responsive_direction",
        "n_DTP_responsive",
        "median_logFC",
        "median_abs_logFC"
    ]
].copy()


generic_small = generic_small.rename(
    columns={
        "responsive_direction":
            "generic_DTP_direction",

        "n_DTP_responsive":
            "generic_n_models_responsive",

        "median_logFC":
            "generic_median_logFC",

        "median_abs_logFC":
            "generic_median_abs_logFC"
    }
)


# ============================================================
# Merge datasets
# ============================================================

integrated = persistence.merge(
    generic_small,
    on="gene",
    how="left"
)


# Genes absent from GSE255958 profile
integrated["generic_DTP_core"] = (
    integrated["generic_DTP_core"]
    .fillna(False)
    .astype(bool)
)

integrated["generic_DTP_2plus"] = (
    integrated["generic_DTP_2plus"]
    .fillna(False)
    .astype(bool)
)
print("\ngeneric_DTP_core dtype:", integrated["generic_DTP_core"].dtype)
print("generic_DTP_2plus dtype:", integrated["generic_DTP_2plus"].dtype)

# ============================================================
# Direction agreement between datasets
# ============================================================

integrated[
    "cross_dataset_same_direction"
] = (
    integrated["persistence_direction"] ==
    integrated["generic_DTP_direction"]
)


# ============================================================
# Strongest overlap:
#
# GSE193258:
# persistence support in 4/4
# same persistence direction
#
# GSE255958:
# generic DTP in 3/3
# same direction
#
# Across datasets:
# directions agree
# ============================================================

integrated[
    "generic_persistent_core"
] = (
    integrated["persistent_core_4of4"] &
    integrated["generic_DTP_core"] &
    integrated["cross_dataset_same_direction"]
)


# ============================================================
# Broader overlap:
#
# persistence >=3/4
# generic DTP >=2/3
# matching direction
# ============================================================

integrated[
    "generic_persistent_robust"
] = (
    integrated["persistent_robust_3plus"] &
    integrated["generic_DTP_2plus"] &
    integrated["cross_dataset_same_direction"]
)


# ============================================================
# Persistent signals without generic support
#
# Important:
# This does NOT mean osimertinib-specific.
# It only means generic DTP support was not observed here.
# ============================================================

integrated[
    "persistent_without_generic_core_support"
] = (
    integrated["persistent_core_4of4"] &
    ~integrated["generic_DTP_core"]
)


integrated[
    "persistent_without_generic_recurrent_support"
] = (
    integrated["persistent_robust_3plus"] &
    ~integrated["generic_DTP_2plus"]
)


# ============================================================
# Save full integrated table
# ============================================================

integrated.to_csv(
    os.path.join(
        output_dir,
        "GSE193258_GSE255958_integrated_profile.csv"
    ),
    index=False
)


# ============================================================
# Save strongest generic-persistent genes
# ============================================================

generic_core = integrated[
    integrated["generic_persistent_core"]
].copy()

generic_core.to_csv(
    os.path.join(
        output_dir,
        "generic_persistent_core.csv"
    ),
    index=False
)


# ============================================================
# Save broader generic-persistent genes
# ============================================================

generic_robust = integrated[
    integrated["generic_persistent_robust"]
].copy()

generic_robust.to_csv(
    os.path.join(
        output_dir,
        "generic_persistent_robust.csv"
    ),
    index=False
)


# ============================================================
# Save persistent genes without generic support
# ============================================================

nongeneric_core = integrated[
    integrated[
        "persistent_without_generic_core_support"
    ]
].copy()

nongeneric_core.to_csv(
    os.path.join(
        output_dir,
        "persistent_core_without_generic_support.csv"
    ),
    index=False
)


# ============================================================
# Summary
# ============================================================

print("\n======================================")
print("GSE193258 persistence")
print("======================================")

print(
    "4/4 persistence before direction check:",
    (
        integrated["n_persistence_support"] == 4
    ).sum()
)

print(
    "4/4 persistence + consistent direction:",
    integrated[
        "persistent_core_4of4"
    ].sum()
)

print(
    ">=3/4 persistence + consistent direction:",
    integrated[
        "persistent_robust_3plus"
    ].sum()
)


print("\n======================================")
print("Cross-dataset integration")
print("======================================")

print(
    "Core generic-persistent genes:",
    integrated[
        "generic_persistent_core"
    ].sum()
)

print(
    "Robust generic-persistent genes:",
    integrated[
        "generic_persistent_robust"
    ].sum()
)

print(
    "Core persistent genes without "
    "3/3 generic DTP support:",
    integrated[
        "persistent_without_generic_core_support"
    ].sum()
)

print(
    "Robust persistent genes without "
    ">=2/3 generic DTP support:",
    integrated[
        "persistent_without_generic_recurrent_support"
    ].sum()
)


# ============================================================
# Direction counts
# ============================================================

print("\nCore generic-persistent direction:")

print(
    generic_core[
        "persistence_direction"
    ].value_counts()
)


# ============================================================
# Top genes
# ============================================================

print(
    "\nTop 20 core generic-persistent genes:"
)

cols_to_show = [
    "gene",
    "persistence_direction",
    "n_persistence_support",
    "median_long_retention",
    "generic_n_models_responsive",
    "generic_median_logFC"
]

print(
    generic_core[
        cols_to_show
    ]
    .sort_values(
        by="median_long_retention",
        ascending=False
    )
    .head(20)
    .to_string(index=False)
)