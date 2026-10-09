library(edgeR)
library(limma)


# ============================================================
# Paths
# ============================================================

counts_file <- "data/processed/GSE255958_core_expected_counts.tsv.gz"

metadata_file <- "data/metadata/GSE255958_core_metadata.csv"

output_root <- "results/differential_expression/GSE255958"

figure_root <- "figures/differential_expression/GSE255958"

dir.create(
  output_root,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  figure_root,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Read data
# ============================================================

counts <- read.delim(
  counts_file,
  check.names = FALSE,
  row.names = 1
)

metadata <- read.csv(
  metadata_file
)

cat(
  "Count matrix:",
  nrow(counts),
  "genes x",
  ncol(counts),
  "samples\n"
)

cat(
  "Metadata samples:",
  nrow(metadata),
  "\n"
)


# ============================================================
# Check sample matching
# ============================================================

stopifnot(
  all(
    metadata$sample_id %in%
      colnames(counts)
  )
)

counts <- counts[
  ,
  metadata$sample_id
]

stopifnot(
  all(
    colnames(counts) ==
      metadata$sample_id
  )
)


# ============================================================
# Cell lines
# ============================================================

cell_lines <- c(
  "PC9",
  "H3122",
  "H358"
)


# ============================================================
# Loop over cell lines
# ============================================================

for (cell_line in cell_lines) {
  
  cat(
    "\n====================================\n"
  )
  
  cat(
    "Cell line:",
    cell_line,
    "\n"
  )
  
  cat(
    "====================================\n"
  )
  
  
  # ----------------------------------------------------------
  # Select metadata
  # ----------------------------------------------------------
  
  meta_sub <- metadata[
    metadata$cell_line == cell_line,
  ]
  
  sample_names <- meta_sub$sample_id
  
  
  # ----------------------------------------------------------
  # Define groups
  # ----------------------------------------------------------
  
  group <- factor(
    meta_sub$treatment_state,
    levels = c(
      "control",
      "acute",
      "DTP",
      "resistant"
    )
  )
  
  print(
    table(group)
  )
  
  
  # ----------------------------------------------------------
  # Select counts
  # ----------------------------------------------------------
  
  counts_sub <- counts[
    ,
    sample_names
  ]
  
  
  # ----------------------------------------------------------
  # Output folders
  # ----------------------------------------------------------
  
  output_dir <- file.path(
    output_root,
    cell_line
  )
  
  dir.create(
    output_dir,
    recursive = TRUE,
    showWarnings = FALSE
  )
  
  
  # ----------------------------------------------------------
  # edgeR object
  # ----------------------------------------------------------
  
  y <- DGEList(
    counts = counts_sub,
    group = group
  )
  
  
  # ----------------------------------------------------------
  # Filter low-expression genes
  # ----------------------------------------------------------
  
  keep <- filterByExpr(
    y,
    group = group
  )
  
  y <- y[
    keep,
    ,
    keep.lib.sizes = FALSE
  ]
  
  cat(
    "Genes retained after filtering:",
    nrow(y),
    "\n"
  )
  
  
  # ----------------------------------------------------------
  # Normalize
  # ----------------------------------------------------------
  
  y <- normLibSizes(y)
  
  
  # ----------------------------------------------------------
  # Design matrix
  # ----------------------------------------------------------
  
  design <- model.matrix(
    ~ 0 + group
  )
  
  colnames(design) <- levels(group)
  
  write.csv(
    design,
    file.path(
      output_dir,
      "design_matrix.csv"
    )
  )
  
  
  # ----------------------------------------------------------
  # voom
  # ----------------------------------------------------------
  
  png(
    file.path(
      figure_root,
      paste0(
        cell_line,
        "_voom.png"
      )
    ),
    width = 1000,
    height = 800
  )
  
  v <- voom(
    y,
    design,
    plot = TRUE
  )
  
  dev.off()
  
  
  # ----------------------------------------------------------
  # Linear model
  # ----------------------------------------------------------
  
  fit <- lmFit(
    v,
    design
  )
  
  
  # ----------------------------------------------------------
  # Contrasts
  # ----------------------------------------------------------
  
  contrasts <- makeContrasts(
    
    Acute_vs_Control =
      acute - control,
    
    DTP_vs_Control =
      DTP - control,
    
    DTP_vs_Acute =
      DTP - acute,
    
    Resistant_vs_Control =
      resistant - control,
    
    Resistant_vs_DTP =
      resistant - DTP,
    
    levels = design
  )
  
  
  # ----------------------------------------------------------
  # Apply contrasts
  # ----------------------------------------------------------
  
  fit2 <- contrasts.fit(
    fit,
    contrasts
  )
  
  fit2 <- eBayes(
    fit2,
    robust = TRUE
  )
  
  
  # ----------------------------------------------------------
  # Save results
  # ----------------------------------------------------------
  
  for (comparison in colnames(contrasts)) {
    
    result <- topTable(
      fit2,
      coef = comparison,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )
    
    result$gene <- rownames(result)
    
    result <- result[
      ,
      c(
        "gene",
        setdiff(
          colnames(result),
          "gene"
        )
      )
    ]
    
    write.csv(
      result,
      file.path(
        output_dir,
        paste0(
          comparison,
          ".csv"
        )
      ),
      row.names = FALSE
    )
    
    
    # --------------------------------------------------------
    # Summary counts
    # --------------------------------------------------------
    
    fdr_only <- sum(
      result$adj.P.Val < 0.05,
      na.rm = TRUE
    )
    
    fdr_fc05 <- sum(
      result$adj.P.Val < 0.05 &
        abs(result$logFC) >= 0.5,
      na.rm = TRUE
    )
    
    fdr_fc1 <- sum(
      result$adj.P.Val < 0.05 &
        abs(result$logFC) >= 1,
      na.rm = TRUE
    )
    
    
    cat(
      "\n",
      comparison,
      "\n"
    )
    
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
  
  
  cat(
    "\n",
    cell_line,
    "finished.\n"
  )
}


cat(
  "\nAll GSE255958 core differential expression analyses completed.\n"
)