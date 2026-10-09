import os
import pandas as pd


# ============================================================
# Paths
# ============================================================

extract_dir = "data/raw/GSE255958/extracted"
output_file = "data/metadata/GSE255958_core_metadata.csv"


# ============================================================
# Helper function
# ============================================================

rows = []


def add_group(
    gsm_start,
    filenames,
    cell_line,
    driver,
    state,
    drug,
    duration,
    model_status
):
    for i, filename in enumerate(filenames, start=1):

        gsm_number = gsm_start + i - 1
        gsm = f"GSM{gsm_number}"

        rows.append(
            {
                "dataset": "GSE255958",
                "sample_id": gsm,
                "expression_file": filename,
                "cell_line": cell_line,
                "oncogenic_driver": driver,
                "treatment_state": state,
                "drug": drug,
                "treatment_duration": duration,
                "model_status": model_status,
                "replicate": i,
            }
        )


# ============================================================
# PC9
# EGFR-mutant
# ============================================================

add_group(
    8083374,
    [
        "GSM8083374_PC9_D_1.genes.results.gz",
        "GSM8083375_PC9_D_2.genes.results.gz",
        "GSM8083376_PC9_D_3.genes.results.gz",
    ],
    "PC9",
    "EGFR-mutant",
    "control",
    "DMSO",
    "48h",
    "parental",
)

add_group(
    8083377,
    [
        "GSM8083377_PC9_O2_1.genes.results.gz",
        "GSM8083378_PC9_O2_2.genes.results.gz",
        "GSM8083379_PC9_O2_3.genes.results.gz",
    ],
    "PC9",
    "EGFR-mutant",
    "acute",
    "osimertinib",
    "48h",
    "parental",
)

add_group(
    8083380,
    [
        "GSM8083380_PC9_O9_1.genes.results.gz",
        "GSM8083381_PC9_O9_2.genes.results.gz",
        "GSM8083382_PC9_O9_3.genes.results.gz",
    ],
    "PC9",
    "EGFR-mutant",
    "DTP",
    "osimertinib",
    "9d",
    "persister",
)

add_group(
    8083383,
    [
        "GSM8083383_PC9_d47_1.genes.results.gz",
        "GSM8083384_PC9_d47_2.genes.results.gz",
        "GSM8083385_PC9_d47_3.genes.results.gz",
    ],
    "PC9",
    "EGFR-mutant",
    "resistant",
    "osimertinib",
    "48h",
    "acquired_resistant",
)


# ============================================================
# H3122
# ALK fusion-positive
# ============================================================

add_group(
    8083386,
    [
        "GSM8083386_H3122_1.genes.results.gz",
        "GSM8083387_H3122_2.genes.results.gz",
        "GSM8083388_H3122_3.genes.results.gz",
    ],
    "H3122",
    "ALK-fusion",
    "control",
    "DMSO",
    "48h",
    "parental",
)

add_group(
    8083389,
    [
        "GSM8083389_H3122_4.genes.results.gz",
        "GSM8083390_H3122_5.genes.results.gz",
        "GSM8083391_H3122_6.genes.results.gz",
    ],
    "H3122",
    "ALK-fusion",
    "acute",
    "alectinib",
    "48h",
    "parental",
)

add_group(
    8083392,
    [
        "GSM8083392_H3122_7.genes.results.gz",
        "GSM8083393_H3122_8.genes.results.gz",
        "GSM8083394_H3122_9.genes.results.gz",
    ],
    "H3122",
    "ALK-fusion",
    "DTP",
    "alectinib",
    "9d",
    "persister",
)

add_group(
    8083395,
    [
        "GSM8083395_H3122_10.genes.results.gz",
        "GSM8083396_H3122_11.genes.results.gz",
        "GSM8083397_H3122_12.genes.results.gz",
    ],
    "H3122",
    "ALK-fusion",
    "resistant",
    "alectinib",
    "48h",
    "acquired_resistant",
)


# ============================================================
# H358
# KRAS-mutant
# ============================================================

add_group(
    8083398,
    [
        "GSM8083398_H358_1.genes.results.gz",
        "GSM8083399_H358_2.genes.results.gz",
        "GSM8083400_H358_3.genes.results.gz",
    ],
    "H358",
    "KRAS-mutant",
    "control",
    "DMSO",
    "48h",
    "parental",
)

add_group(
    8083401,
    [
        "GSM8083401_H358_4.genes.results.gz",
        "GSM8083402_H358_5.genes.results.gz",
        "GSM8083403_H358_6.genes.results.gz",
    ],
    "H358",
    "KRAS-mutant",
    "acute",
    "RMC-4550",
    "48h",
    "parental",
)

add_group(
    8083404,
    [
        "GSM8083404_H358_7.genes.results.gz",
        "GSM8083405_H358_8.genes.results.gz",
        "GSM8083406_H358_9.genes.results.gz",
    ],
    "H358",
    "KRAS-mutant",
    "DTP",
    "RMC-4550",
    "9d",
    "persister",
)

add_group(
    8083407,
    [
        "GSM8083407_H358_10.genes.results.gz",
        "GSM8083408_H358_11.genes.results.gz",
        "GSM8083409_H358_12.genes.results.gz",
    ],
    "H358",
    "KRAS-mutant",
    "resistant",
    "RMC-4550",
    "48h",
    "acquired_resistant",
)


# ============================================================
# Create metadata table
# ============================================================

metadata = pd.DataFrame(rows)


# ============================================================
# Basic checks
# ============================================================

assert len(metadata) == 36

assert metadata["sample_id"].is_unique

assert metadata["expression_file"].is_unique


# Check that all 36 expression files really exist

missing_files = []

for filename in metadata["expression_file"]:
    path = os.path.join(extract_dir, filename)

    if not os.path.exists(path):
        missing_files.append(filename)


if missing_files:
    raise FileNotFoundError(
        f"Missing expression files: {missing_files}"
    )


# ============================================================
# Save metadata
# ============================================================

metadata.to_csv(
    output_file,
    index=False
)


# ============================================================
# Print summary
# ============================================================

print("\nMetadata created successfully.")

print("\nNumber of samples:")
print(len(metadata))

print("\nSamples by cell line:")
print(metadata["cell_line"].value_counts())

print("\nSamples by treatment state:")
print(metadata["treatment_state"].value_counts())

print("\nCell line x treatment state:")
print(
    pd.crosstab(
        metadata["cell_line"],
        metadata["treatment_state"]
    )
)

print("\nSaved to:")
print(output_file)