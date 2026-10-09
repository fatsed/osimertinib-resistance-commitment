# ============================================================
# GSE193258 - H1975
# Differential Expression Summary + Gene Trajectory Table
# ============================================================


# ------------------------------------------------------------
# 1. Folder containing differential expression results
# ------------------------------------------------------------

base_dir <- "results/differential_expression/H1975"


# ------------------------------------------------------------
# 2. Read differential expression result files
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 3. Function for summarizing significant genes
# ------------------------------------------------------------

summarize_de <- function(df, name) {
  
  fdr_only <- sum(
    df$adj.P.Val < 0.05,
    na.rm = TRUE
  )
  
  fdr_fc05 <- sum(
    df$adj.P.Val < 0.05 &
      abs(df$logFC) >= 0.5,
    na.rm = TRUE
  )
  
  fdr_fc1 <- sum(
    df$adj.P.Val < 0.05 &
      abs(df$logFC) >= 1,
    na.rm = TRUE
  )
  
  cat("\n", name, "\n")
  
  cat(
    "FDR < 0.05:",
    fdr_only,
    "\n"
  )
  
  cat(
    "FDR < 0.05 and |logFC| >= 0.5:",
    fdr_fc05,
    "\n"
  )
  
  cat(
    "FDR < 0.05 and |logFC| >= 1:",
    fdr_fc1,
    "\n"
  )
}


# ------------------------------------------------------------
# 4. Print summaries
# ------------------------------------------------------------

summarize_de(
  acute_control,
  "Acute vs Control"
)

summarize_de(
  dtp_control,
  "DTP vs Control"
)

summarize_de(
  dtp_acute,
  "DTP vs Acute"
)

summarize_de(
  short_control,
  "Short Washout vs Control"
)

summarize_de(
  long_control,
  "Long Washout vs Control"
)

summarize_de(
  long_dtp,
  "Long Washout vs DTP"
)


# ------------------------------------------------------------
# 5. Keep gene, logFC and FDR from every comparison
# ------------------------------------------------------------


# Acute vs Control

acute_control_small <- acute_control[
  ,
  c(
    "gene",
    "logFC",
    "adj.P.Val"
  )
]

names(acute_control_small) <- c(
  "gene",
  "Acute_vs_Control_logFC",
  "Acute_vs_Control_FDR"
)


# DTP vs Control

dtp_control_small <- dtp_control[
  ,
  c(
    "gene",
    "logFC",
    "adj.P.Val"
  )
]

names(dtp_control_small) <- c(
  "gene",
  "DTP_vs_Control_logFC",
  "DTP_vs_Control_FDR"
)


# DTP vs Acute

dtp_acute_small <- dtp_acute[
  ,
  c(
    "gene",
    "logFC",
    "adj.P.Val"
  )
]

names(dtp_acute_small) <- c(
  "gene",
  "DTP_vs_Acute_logFC",
  "DTP_vs_Acute_FDR"
)


# Short Washout vs Control

short_control_small <- short_control[
  ,
  c(
    "gene",
    "logFC",
    "adj.P.Val"
  )
]

names(short_control_small) <- c(
  "gene",
  "Short_vs_Control_logFC",
  "Short_vs_Control_FDR"
)


# Long Washout vs Control

long_control_small <- long_control[
  ,
  c(
    "gene",
    "logFC",
    "adj.P.Val"
  )
]

names(long_control_small) <- c(
  "gene",
  "Long_vs_Control_logFC",
  "Long_vs_Control_FDR"
)


# Long Washout vs DTP

long_dtp_small <- long_dtp[
  ,
  c(
    "gene",
    "logFC",
    "adj.P.Val"
  )
]

names(long_dtp_small) <- c(
  "gene",
  "Long_vs_DTP_logFC",
  "Long_vs_DTP_FDR"
)


# ------------------------------------------------------------
# 6. Merge all comparisons into one gene trajectory table
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 7. Check whether DTP and washout changes have same direction
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 8. Calculate how much of the DTP change remains after washout
# ------------------------------------------------------------

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


# ------------------------------------------------------------
# 9. Measure how different Long Washout is from DTP
# ------------------------------------------------------------

trajectory$distance_from_DTP <- abs(
  trajectory$Long_vs_DTP_logFC
)


# ------------------------------------------------------------
# 10. Save complete trajectory table
# ------------------------------------------------------------

write.csv(
  
  trajectory,
  
  file.path(
    base_dir,
    "H1975_gene_trajectory_table.csv"
  ),
  
  row.names = FALSE
  
)


# ------------------------------------------------------------
# 11. Basic checks
# ------------------------------------------------------------

cat(
  "\nTrajectory table saved with",
  nrow(trajectory),
  "genes.\n"
)


cat(
  "\nGenes with same DTP and Short Washout direction:",
  sum(
    trajectory$short_same_direction,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Genes with same DTP and Long Washout direction:",
  sum(
    trajectory$long_same_direction,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nMedian Short Washout retention ratio:",
  median(
    trajectory$short_retention_ratio,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "Median Long Washout retention ratio:",
  median(
    trajectory$long_retention_ratio,
    na.rm = TRUE
  ),
  "\n"
)


cat(
  "\nH1975 trajectory analysis finished.\n"
)


# ------------------------------------------------------------
# 12. Focus on genes that truly changed in DTP
# ------------------------------------------------------------

dtp_candidates_05 <- trajectory[
  trajectory$DTP_vs_Control_FDR < 0.05 &
    abs(trajectory$DTP_vs_Control_logFC) >= 0.5,
]

dtp_candidates_1 <- trajectory[
  trajectory$DTP_vs_Control_FDR < 0.05 &
    abs(trajectory$DTP_vs_Control_logFC) >= 1,
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


# ------------------------------------------------------------
# 13. Retention summary for DTP-responsive genes
# ------------------------------------------------------------

cat("\n--- Threshold |logFC| >= 0.5 ---\n")

cat(
  "Median short retention:",
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

cat("\nShort retention quartiles:\n")

print(
  quantile(
    dtp_candidates_05$short_retention_ratio,
    probs = c(0.25, 0.5, 0.75),
    na.rm = TRUE
  )
)

cat("\nLong retention quartiles:\n")

print(
  quantile(
    dtp_candidates_05$long_retention_ratio,
    probs = c(0.25, 0.5, 0.75),
    na.rm = TRUE
  )
)


cat("\n--- Threshold |logFC| >= 1 ---\n")

cat(
  "Median short retention:",
  median(
    dtp_candidates_1$short_retention_ratio,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Median long retention:",
  median(
    dtp_candidates_1$long_retention_ratio,
    na.rm = TRUE
  ),
  "\n"
)

# ------------------------------------------------------------
# 14. Visualize DTP-to-washout persistence
# ------------------------------------------------------------

dir.create(
  "figures/differential_expression",
  recursive = TRUE,
  showWarnings = FALSE
)

# Use DTP-responsive genes
plot_data <- dtp_candidates_05


# Short washout vs DTP
png(
  "figures/differential_expression/H1975_DTP_vs_ShortWashout.png",
  width = 1000,
  height = 900
)

plot(
  plot_data$DTP_vs_Control_logFC,
  plot_data$Short_vs_Control_logFC,
  pch = 16,
  cex = 0.5,
  xlab = "DTP vs Control logFC",
  ylab = "Short Washout vs Control logFC",
  main = "H1975: DTP vs Short Washout"
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


# Long washout vs DTP
png(
  "figures/differential_expression/H1975_DTP_vs_LongWashout.png",
  width = 1000,
  height = 900
)

plot(
  plot_data$DTP_vs_Control_logFC,
  plot_data$Long_vs_Control_logFC,
  pch = 16,
  cex = 0.5,
  xlab = "DTP vs Control logFC",
  ylab = "Long Washout vs Control logFC",
  main = "H1975: DTP vs Long Washout"
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


# ------------------------------------------------------------
# 15. Correlation of DTP and washout effects
# ------------------------------------------------------------

short_cor <- cor(
  plot_data$DTP_vs_Control_logFC,
  plot_data$Short_vs_Control_logFC,
  method = "spearman",
  use = "complete.obs"
)

long_cor <- cor(
  plot_data$DTP_vs_Control_logFC,
  plot_data$Long_vs_Control_logFC,
  method = "spearman",
  use = "complete.obs"
)

cat(
  "\nSpearman correlation - DTP vs Short Washout:",
  round(short_cor, 3),
  "\n"
)

cat(
  "Spearman correlation - DTP vs Long Washout:",
  round(long_cor, 3),
  "\n"
)