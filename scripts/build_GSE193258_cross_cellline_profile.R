# ============================================================
# GSE193258
# Cross-cell-line persistence profile
# ============================================================

cell_lines <- c(
  "H1975",
  "PC9",
  "HCC827",
  "HCC2935"
)

dir.create(
  "results/persistence",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Function to read and prepare one cell line
# ============================================================

prepare_cell_line <- function(cell_line) {
  
  file_path <- file.path(
    "results/differential_expression",
    cell_line,
    paste0(
      cell_line,
      "_gene_trajectory_table.csv"
    )
  )
  
  df <- read.csv(file_path)
  
  # ----------------------------------------
  # Is the gene genuinely DTP-responsive?
  # ----------------------------------------
  
  df$DTP_responsive <- (
    df$DTP_vs_Control_FDR < 0.05 &
      abs(df$DTP_vs_Control_logFC) >= 0.5
  )
  
  
  # ----------------------------------------
  # Is the change still significant
  # after long washout?
  # ----------------------------------------
  
  df$Long_significant <- (
    df$Long_vs_Control_FDR < 0.05
  )
  
  
  # ----------------------------------------
  # Exploratory persistence evidence
  #
  # Requirements:
  # 1. Changed in DTP
  # 2. Same direction after long washout
  # 3. Still significant vs control
  # 4. At least 50% of DTP magnitude retained
  #
  # This is NOT the final Persistence Score.
  # ----------------------------------------
  
  df$Persistence_support <- (
    df$DTP_responsive &
      df$long_same_direction &
      df$Long_significant &
      df$long_retention_ratio >= 0.5
  )
  
  
  # Keep only useful columns
  
  out <- df[
    ,
    c(
      "gene",
      "DTP_vs_Control_logFC",
      "DTP_vs_Control_FDR",
      "Long_vs_Control_logFC",
      "Long_vs_Control_FDR",
      "long_retention_ratio",
      "long_same_direction",
      "DTP_responsive",
      "Persistence_support"
    )
  ]
  
  
  # Add cell-line name to column names
  
  names(out)[-1] <- paste0(
    names(out)[-1],
    "_",
    cell_line
  )
  
  return(out)
}


# ============================================================
# Read all four cell lines
# ============================================================

cellline_tables <- lapply(
  cell_lines,
  prepare_cell_line
)


# ============================================================
# Merge all genes across all cell lines
# ============================================================

combined <- Reduce(
  function(x, y) {
    merge(
      x,
      y,
      by = "gene",
      all = TRUE
    )
  },
  cellline_tables
)


# ============================================================
# Count DTP responsiveness across cell lines
# ============================================================

dtp_columns <- paste0(
  "DTP_responsive_",
  cell_lines
)

combined$n_DTP_responsive <- rowSums(
  combined[, dtp_columns],
  na.rm = TRUE
)


# ============================================================
# Count persistence support across cell lines
# ============================================================

persistence_columns <- paste0(
  "Persistence_support_",
  cell_lines
)

combined$n_persistence_support <- rowSums(
  combined[, persistence_columns],
  na.rm = TRUE
)


# ============================================================
# Median long-washout retention across cell lines
# ============================================================

retention_columns <- paste0(
  "long_retention_ratio_",
  cell_lines
)

combined$median_long_retention <- apply(
  combined[, retention_columns],
  1,
  function(x) {
    
    if (all(is.na(x))) {
      return(NA)
    }
    
    median(
      x,
      na.rm = TRUE
    )
  }
)


# ============================================================
# Simple recurrence category
# ============================================================

combined$persistence_recurrence <- ifelse(
  combined$n_persistence_support == 4,
  "4_of_4",
  ifelse(
    combined$n_persistence_support == 3,
    "3_of_4",
    ifelse(
      combined$n_persistence_support == 2,
      "2_of_4",
      ifelse(
        combined$n_persistence_support == 1,
        "1_of_4",
        "0_of_4"
      )
    )
  )
)


# ============================================================
# Sort genes:
# strongest cross-cell-line evidence first
# ============================================================

combined <- combined[
  order(
    -combined$n_persistence_support,
    -combined$n_DTP_responsive,
    -combined$median_long_retention
  ),
]


# ============================================================
# Save complete table
# ============================================================

write.csv(
  combined,
  "results/persistence/GSE193258_cross_cellline_persistence_profile.csv",
  row.names = FALSE
)


# ============================================================
# Save compact summary table
# ============================================================

compact <- combined[
  ,
  c(
    "gene",
    "n_DTP_responsive",
    "n_persistence_support",
    "median_long_retention",
    "persistence_recurrence"
  )
]

write.csv(
  compact,
  "results/persistence/GSE193258_cross_cellline_persistence_compact.csv",
  row.names = FALSE
)


# ============================================================
# Print summary
# ============================================================

cat("\nCross-cell-line profile created.\n\n")

cat(
  "Total genes:",
  nrow(combined),
  "\n"
)

cat(
  "Persistent-like in 4/4 cell lines:",
  sum(
    combined$n_persistence_support == 4,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Persistent-like in 3/4 cell lines:",
  sum(
    combined$n_persistence_support == 3,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Persistent-like in 2/4 cell lines:",
  sum(
    combined$n_persistence_support == 2,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "Persistent-like in 1/4 cell lines:",
  sum(
    combined$n_persistence_support == 1,
    na.rm = TRUE
  ),
  "\n"
)

cat(
  "\nFiles saved in results/persistence/\n"
)