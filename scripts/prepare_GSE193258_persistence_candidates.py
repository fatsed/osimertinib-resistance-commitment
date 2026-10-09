import os
import pandas as pd

# --------------------------------------------------
# 1. Read cross-cell-line persistence profile
# --------------------------------------------------

df = pd.read_csv(
    "results/persistence/"
    "GSE193258_cross_cellline_persistence_profile.csv"
)

print("Total genes:", len(df))


# --------------------------------------------------
# 2. Strongest candidates: persistent in 4/4
# --------------------------------------------------

core_4of4 = df[
    df["n_persistence_support"] == 4
].copy()

print(
    "Persistent in 4/4 cell lines:",
    len(core_4of4)
)


# --------------------------------------------------
# 3. Robust candidates: persistent in at least 3/4
# --------------------------------------------------

robust_3plus = df[
    df["n_persistence_support"] >= 3
].copy()

print(
    "Persistent in at least 3/4 cell lines:",
    len(robust_3plus)
)


# --------------------------------------------------
# 4. Moderate evidence: at least 2/4
# --------------------------------------------------

support_2plus = df[
    df["n_persistence_support"] >= 2
].copy()

print(
    "Persistent in at least 2/4 cell lines:",
    len(support_2plus)
)


# --------------------------------------------------
# 5. Sort candidates
# --------------------------------------------------

core_4of4 = core_4of4.sort_values(
    by="median_long_retention",
    ascending=False
)

robust_3plus = robust_3plus.sort_values(
    by=[
        "n_persistence_support",
        "median_long_retention"
    ],
    ascending=[
        False,
        False
    ]
)


# --------------------------------------------------
# 6. Save candidate sets
# --------------------------------------------------

output_dir = "results/persistence/candidates"

os.makedirs(
    output_dir,
    exist_ok=True
)

core_4of4.to_csv(
    f"{output_dir}/GSE193258_core_4of4.csv",
    index=False
)

robust_3plus.to_csv(
    f"{output_dir}/GSE193258_robust_3plus.csv",
    index=False
)

support_2plus.to_csv(
    f"{output_dir}/GSE193258_support_2plus.csv",
    index=False
)


# --------------------------------------------------
# 7. Print top candidates
# --------------------------------------------------

print("\nTop 20 core persistence candidates:")

print(
    core_4of4[
        [
            "gene",
            "n_persistence_support",
            "median_long_retention"
        ]
    ].head(20)
)
