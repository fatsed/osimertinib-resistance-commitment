# ============================================================
# GSE193258
# Gene trajectory and persistence summary for one cell line
# ============================================================

# Choose cell line
#cell_line <- "PC9"
#cell_line <- "HCC827"
cell_line <- "HCC2935"
# ============================================================
# Paths
# ============================================================

base_dir <- file.path(
  "results/differential_expression",
  cell_line
)

figure_dir <- file.path(
  "figures/differential_expression"
)

dir.create(
  figure_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Read differential expression results
# ============================================================

acute_control <- read.csv(
  file.path(base_dir, "Acute_vs_Control.csv")
)

dtp_control <- read.csv(
  file.path(base_dir, "DTP_vs_Control.csv")
)

dtp_acute <- read.csv(
  file.path(base_dir, "DTP_vs_Acute.csv")
)

short_control <- read.csv(
  file.path(base_dir, "ShortWashout_vs_Control.csv")
)

long_control <- read.csv(
  file.path(base_dir, "LongWashout_vs_Control.csv")
)

long_dtp <- read.csv(
  file.path(base_dir, "LongWashout_vs_DTP.csv")
)


# ============================================================
# Function for keeping important columns
# ============================================================

prepare_result <- function(
    df,
    logfc_name,
    fdr_name
) {
  
  out <- df[
    ,
    c(
      "gene",
      "logFC",
      "adj.P.Val"
    )
  ]
  
  names(out) <- c(
    "gene",
    logfc_name,
    fdr_name
  )
  
  return(out)
}


# ============================================================
# Prepare each comparison
# ============================================================

acute_control_small <- prepare_result(
  acute_control,
  "Acute_vs_Control_logFC",
  "Acute_vs_Control_FDR"
)

dtp_control_small <- prepare_result(
  dtp_control,
  "DTP_vs_Control_logFC",
  "DTP_vs_Control_FDR"
)

dtp_acute_small <- prepare_result(
  dtp_acute,
  "DTP_vs_Acute_logFC",
  "DTP_vs_Acute_FDR"
)

short_control_small <- prepare_result(
  short_control,
  "Short_vs_Control_logFC",
  "Short_vs_Control_FDR"
)

long_control_small <- prepare_result(
  long_control,
  "Long_vs_Control_logFC",
  "Long_vs_Control_FDR"
)

long_dtp_small <- prepare_result(
  long_dtp,
  "Long_vs_DTP_logFC",
  "Long_vs_DTP_FDR"
)


# ============================================================
# Merge all comparisons into one trajectory table
# ============================================================

trajectory <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "gene"
    )
  },
  
  list(
    acute_control_small,
    dtp_control_small,
    dtp_acute_small,
    short_control_small,
    long_control_small,
    long_dtp_small
  )
)


# ============================================================
# Direction consistency
# ============================================================

trajectory$short_same_direction <- (
  sign(
    trajectory$DTP_vs_Control_logFC
  ) ==
    sign(
      trajectory$Short_vs_Control_logFC
    )
)

trajectory$long_same_direction <- (
  sign(
    trajectory$DTP_vs_Control_logFC
  ) ==
    sign(
      trajectory$Long_vs_Control_logFC
    )
)


# ============================================================
# Retention ratios
# ============================================================

trajectory$short_retention_ratio <- (
  abs(
    trajectory$Short_vs_Control_logFC
  ) /
    (
      abs(
        trajectory$DTP_vs_Control_logFC
      ) + 0.001
    )
)

trajectory$long_retention_ratio <- (
  abs(
    trajectory$Long_vs_Control_logFC
  ) /
    (
      abs(
        trajectory$DTP_vs_Control_logFC
      ) + 0.001
    )
)


# ============================================================
# Distance between Long Washout and DTP
# ============================================================

trajectory$distance_from_DTP <- abs(
  trajectory$Long_vs_DTP_logFC
)


# ============================================================
# Save complete trajectory table
# ============================================================

trajectory_file <- file.path(
  base_dir,
  paste0(
    cell_line,
    "_gene_trajectory_table.csv"
  )
)

write.csv(
  trajectory,
  trajectory_file,
  row.names = FALSE
)

cat(
  "\nTrajectory table saved with",
  nrow(trajectory),
  "genes.\n"
)


# ============================================================
# Focus on real DTP-responsive genes
# ============================================================

dtp_candidates_05 <- trajectory[
  trajectory$DTP_vs_Control_FDR < 0.05 &
    abs(
      trajectory$DTP_vs_Control_logFC
    ) >= 0.5,
]

dtp_candidates_1 <- trajectory[
  trajectory$DTP_vs_Control_FDR < 0.05 &
    abs(
      trajectory$DTP_vs_Control_logFC
    ) >= 1,
]


cat(
  "\nDTP candidates with |logFC| >= 0.5:",
  nrow(dtp_candidates_05),
  "\n"
)

cat(
  "DTP candidates with |logFC| >= 1:",
  nrow(dtp_candidates_1),
  "\n"
)


# ============================================================
# Retention summaries
# ============================================================

cat(
  "\nMedian short retention:",
  median(
    dtp_candidates_05$short_retention_ratio,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Median long retention:",
  median(
    dtp_candidates_05$long_retention_ratio,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nShort retention quartiles:\n"
)

print(
  quantile(
    dtp_candidates_05$short_retention_ratio,
    probs = c(
      0.25,
      0.50,
      0.75
    ),
    na.rm = TRUE
  )
)

cat(
  "\nLong retention quartiles:\n"
)

print(
  quantile(
    dtp_candidates_05$long_retention_ratio,
    probs = c(
      0.25,
      0.50,
      0.75
    ),
    na.rm = TRUE
  )
)


# ============================================================
# Scatter plot: DTP vs Short Washout
# ============================================================

png(
  file.path(
    figure_dir,
    paste0(
      cell_line,
      "_DTP_vs_ShortWashout.png"
    )
  ),
  width = 1000,
  height = 900
)

plot(
  dtp_candidates_05$DTP_vs_Control_logFC,
  dtp_candidates_05$Short_vs_Control_logFC,
  pch = 16,
  cex = 0.5,
  xlab = "DTP vs Control logFC",
  ylab = "Short Washout vs Control logFC",
  main = paste0(
    cell_line,
    ": DTP vs Short Washout"
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

abline(
  h = 0,
  v = 0,
  lty = 3
)

dev.off()


# ============================================================
# Scatter plot: DTP vs Long Washout
# ============================================================

png(
  file.path(
    figure_dir,
    paste0(
      cell_line,
      "_DTP_vs_LongWashout.png"
    )
  ),
  width = 1000,
  height = 900
)

plot(
  dtp_candidates_05$DTP_vs_Control_logFC,
  dtp_candidates_05$Long_vs_Control_logFC,
  pch = 16,
  cex = 0.5,
  xlab = "DTP vs Control logFC",
  ylab = "Long Washout vs Control logFC",
  main = paste0(
    cell_line,
    ": DTP vs Long Washout"
  )
)

abline(
  a = 0,
  b = 1,
  lty = 2
)

abline(
  h = 0,
  v = 0,
  lty = 3
)

dev.off()


# ============================================================
# Spearman correlations
# ============================================================

short_cor <- cor(
  dtp_candidates_05$DTP_vs_Control_logFC,
  dtp_candidates_05$Short_vs_Control_logFC,
  method = "spearman",
  use = "complete.obs"
)

long_cor <- cor(
  dtp_candidates_05$DTP_vs_Control_logFC,
  dtp_candidates_05$Long_vs_Control_logFC,
  method = "spearman",
  use = "complete.obs"
)


cat(
  "\nSpearman correlation - DTP vs Short Washout:",
  round(
    short_cor,
    3
  ),
  "\n"
)

cat(
  "Spearman correlation - DTP vs Long Washout:",
  round(
    long_cor,
    3
  ),
  "\n"
)


cat(
  "\n",
  cell_line,
  "trajectory analysis finished.\n"
)