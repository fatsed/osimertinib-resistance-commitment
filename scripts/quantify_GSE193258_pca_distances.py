import os
import pandas as pd
import numpy as np

cell_lines = ["PC9", "H1975", "HCC827", "HCC2935"]

os.makedirs("results/pca_distances", exist_ok=True)

all_results = []

for cell_line in cell_lines:

    print("\n" + "=" * 50)
    print("Cell line:", cell_line)
    print("=" * 50)

    # 1. Read PCA scores
    pca_df = pd.read_csv(
        f"results/pca_by_cell_line/{cell_line}_PCA_scores.csv"
    )

    # 2. Calculate centroid of each treatment state
    centroids = (
        pca_df
        .groupby("treatment_state")[["PC1", "PC2"]]
        .mean()
    )

    print("\nCentroids:")
    print(centroids)

    # 3. Helper function for Euclidean distance
    def distance(state1, state2):

        point1 = centroids.loc[state1].values
        point2 = centroids.loc[state2].values

        return np.linalg.norm(point1 - point2)

    # 4. Distances for short washout
    short_to_control = distance(
        "short_washout",
        "control"
    )

    short_to_dtp = distance(
        "short_washout",
        "DTP"
    )

    # 5. Distances for long washout
    long_to_control = distance(
        "long_washout",
        "control"
    )

    long_to_dtp = distance(
        "long_washout",
        "DTP"
    )

    print("\nShort washout:")
    print("Distance to Control:", round(short_to_control, 2))
    print("Distance to DTP:", round(short_to_dtp, 2))

    print("\nLong washout:")
    print("Distance to Control:", round(long_to_control, 2))
    print("Distance to DTP:", round(long_to_dtp, 2))

    # 6. Which state is closer?
    short_closer_to = (
        "Control"
        if short_to_control < short_to_dtp
        else "DTP"
    )

    long_closer_to = (
        "Control"
        if long_to_control < long_to_dtp
        else "DTP"
    )

    print("\nShort washout is closer to:", short_closer_to)
    print("Long washout is closer to:", long_closer_to)

    # 7. Save results
    all_results.append({
        "cell_line": cell_line,
        "short_to_control": short_to_control,
        "short_to_DTP": short_to_dtp,
        "short_closer_to": short_closer_to,
        "long_to_control": long_to_control,
        "long_to_DTP": long_to_dtp,
        "long_closer_to": long_closer_to
    })


# 8. Save summary table
results_df = pd.DataFrame(all_results)

results_df.to_csv(
    "results/pca_distances/GSE193258_PCA_distances.csv",
    index=False
)

print("\n" + "=" * 50)
print("SUMMARY")
print("=" * 50)
print(results_df)