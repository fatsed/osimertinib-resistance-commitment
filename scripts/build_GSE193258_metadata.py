import GEOparse
import pandas as pd

# 1. Download GEO metadata
gse = GEOparse.get_GEO(
    geo="GSE193258",
    destdir="data/metadata/"
)

# 2. Extract information for every sample
rows = []

for gsm_id, gsm in gse.gsms.items():

    title = gsm.metadata.get("title", [""])[0]

    row = {
        "dataset": "GSE193258",
        "sample_id": gsm_id,
        "sample_name": title
    }

    rows.append(row)

# 3. Create dataframe
metadata = pd.DataFrame(rows)

# 4. Save as CSV
metadata.to_csv(
    "data/metadata/GSE193258_metadata_basic.csv",
    index=False
)

print(metadata.head())
print("Number of samples:", len(metadata))

# 5. Clean sample names
metadata["sample_name_clean"] = (
    metadata["sample_name"]
    .str.replace(r"\s*\[RNA-seq\]\s*$", "", regex=True)
)

# 6. Extract cell line
def get_cell_line(name):
    for cell_line in ["H1975", "HCC2935", "HCC827", "PC9"]:
        if name.startswith(cell_line):
            return cell_line
    return "unknown"

metadata["cell_line"] = metadata["sample_name_clean"].apply(get_cell_line)

# 7. Extract treatment state
def get_state(name):
    name = name.lower()

    if "long_washout" in name:
        return "long_washout"
    elif "short_washout" in name:
        return "short_washout"
    elif "dtp" in name:
        return "DTP"
    elif "acute" in name:
        return "acute"
    elif "dmso" in name:
        return "control"
    else:
        return "unknown"

metadata["treatment_state"] = metadata["sample_name_clean"].apply(get_state)

# 8. Extract replicate number
metadata["replicate"] = (
    metadata["sample_name_clean"]
    .str.extract(r"_(\d+)$")[0]
)

# 9. Save complete metadata
metadata.to_csv(
    "data/metadata/GSE193258_metadata.csv",
    index=False
)

print(metadata.head())
print(metadata["treatment_state"].value_counts())

# 9. Add drug information
metadata["drug"] = metadata["treatment_state"].apply(
    lambda x: "DMSO" if x == "control" else "osimertinib"
)

# 10. Add treatment time
def get_treatment_time(state):
    if state == "control":
        return "24h"
    elif state == "acute":
        return "24h"
    elif state in ["DTP", "short_washout", "long_washout"]:
        return "21d"
    return ""

metadata["treatment_time"] = metadata["treatment_state"].apply(get_treatment_time)

# 11. Add washout information
metadata["washout"] = metadata["treatment_state"].apply(
    lambda x: "yes" if x in ["short_washout", "long_washout"] else "no"
)

# 12. Add washout time
def get_washout_time(name, state):
    if state == "short_washout":
        return "24h"

    elif state == "long_washout":
        # Examples:
        # H1975_long_washout_72h_1
        # HCC827_long_washout_96h_1
        # PC9_long_washout_7d_1
        # HCC2935_long_washout_10d_1
        parts = name.split("_")

        for part in parts:
            if part.endswith("h") or part.endswith("d"):
                if part not in ["21d"]:
                    return part

    return ""

metadata["washout_time"] = metadata.apply(
    lambda row: get_washout_time(
        row["sample_name_clean"],
        row["treatment_state"]
    ),
    axis=1
)

# Create sample names that match the expression matrix columns
def get_expression_sample_name(name):

    name = name.replace("_short_washout_", "_short_wash_")

    if "_long_washout_" in name:
        parts = name.split("_")

        # Example:
        # PC9_long_washout_7d_1
        # becomes:
        # PC9_long_wash_1

        cell_line = parts[0]
        replicate = parts[-1]

        return f"{cell_line}_long_wash_{replicate}"

    return name


metadata["expression_sample_name"] = (
    metadata["sample_name_clean"]
    .apply(get_expression_sample_name)
)

metadata.to_csv(
    "data/metadata/GSE193258_metadata.csv",
    index=False
)

print(metadata.head(15))