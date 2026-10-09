# ============================================================
# GSE193258
# Cell-line persistence summary
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

summary_list <- list()


for (cell_line in cell_lines) {
  
  # ---------------------------------------
  # Read trajectory table
  # ---------------------------------------
  
  file_path <- file.path(
    "results/differential_expression",
    cell_line,
    paste0(
      cell_line,
      "_gene_trajectory_table.csv"
    )
  )
  
  trajectory <- read.csv(
    file_path
  )
  
  
  # ---------------------------------------
  # Keep DTP-responsive genes
  # ---------------------------------------
  
  dtp_candidates <- trajectory[
    trajectory$DTP_vs_Control_FDR < 0.05 &
      abs(
        trajectory$DTP_vs_Control_logFC
      ) >= 0.5,
  ]
  
  
  # ---------------------------------------
  # Retention
  # ---------------------------------------
  
  median_short <- median(
    dtp_candidates$short_retention_ratio,
    na.rm = TRUE
  )
  
  median_long <- median(
    dtp_candidates$long_retention_ratio,
    na.rm = TRUE
  )
  
  
  # ---------------------------------------
  # Correlations
  # ---------------------------------------
  
  short_cor <- cor(
    dtp_candidates$DTP_vs_Control_logFC,
    dtp_candidates$Short_vs_Control_logFC,
    method = "spearman",
    use = "complete.obs"
  )
  
  long_cor <- cor(
    dtp_candidates$DTP_vs_Control_logFC,
    dtp_candidates$Long_vs_Control_logFC,
    method = "spearman",
    use = "complete.obs"
  )
  
  
  # ---------------------------------------
  # Save row
  # ---------------------------------------
  
  summary_list[[cell_line]] <- data.frame(
    
    cell_line = cell_line,
    
    DTP_candidates = nrow(
      dtp_candidates
    ),
    
    median_short_retention = median_short,
    
    median_long_retention = median_long,
    
    short_spearman = short_cor,
    
    long_spearman = long_cor
    
  )
}


# ============================================================
# Combine all cell lines
# ============================================================

summary_table <- do.call(
  rbind,
  summary_list
)

rownames(summary_table) <- NULL


# ============================================================
# Save CSV
# ============================================================

write.csv(
  summary_table,
  "results/persistence/GSE193258_cellline_persistence_summary.csv",
  row.names = FALSE
)


# ============================================================
# Print result
# ============================================================

print(summary_table)

cat(
  "\nPersistence summary saved successfully.\n"
)