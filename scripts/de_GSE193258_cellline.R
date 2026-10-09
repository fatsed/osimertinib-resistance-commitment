library(edgeR)
library(limma)

# ============================================================
# Choose cell line
# ============================================================

#cell_line <- "PC9"
#cell_line <- "HCC827"
cell_line <- "HCC2935"
# ============================================================
# Output folders
# ============================================================

output_dir <- file.path(
  "results/differential_expression",
  cell_line
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/differential_expression",
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# Read estimated counts
# ============================================================

counts <- read.delim(
  "data/raw/GSE193258/GSE193258_RNAseq_estimated_counts.tsv.gz",
  check.names = FALSE
)

rownames(counts) <- counts$gene
counts$gene <- NULL


# ============================================================
# Read metadata
# ============================================================

metadata <- read.csv(
  "data/metadata/GSE193258_metadata.csv"
)


# ============================================================
# Select samples for chosen cell line
# ============================================================

meta_sub <- metadata[
  metadata$cell_line == cell_line,
]

group <- factor(
  meta_sub$treatment_state,
  levels = c(
    "control",
    "acute",
    "DTP",
    "short_washout",
    "long_washout"
  )
)

sample_names <- meta_sub$expression_sample_name

counts_sub <- counts[, sample_names]


# ============================================================
# Check sample alignment
# ============================================================

stopifnot(
  all(colnames(counts_sub) == sample_names)
)

cat(
  "\nCell line:",
  cell_line,
  "\n"
)

cat(
  "Number of samples:",
  length(sample_names),
  "\n"
)

print(table(group))


# ============================================================
# edgeR object
# ============================================================

y <- DGEList(
  counts = counts_sub,
  group = group
)


# ============================================================
# Filter low-expression genes
# ============================================================

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


# ============================================================
# Normalize
# ============================================================

y <- normLibSizes(y)


# ============================================================
# Design matrix
# ============================================================

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


# ============================================================
# voom
# ============================================================

png(
  file.path(
    "figures/differential_expression",
    paste0(cell_line, "_voom.png")
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


# ============================================================
# Linear model
# ============================================================

fit <- lmFit(
  v,
  design
)


# ============================================================
# Contrasts
# ============================================================

contrasts <- makeContrasts(
  
  Acute_vs_Control =
    acute - control,
  
  DTP_vs_Control =
    DTP - control,
  
  DTP_vs_Acute =
    DTP - acute,
  
  ShortWashout_vs_Control =
    short_washout - control,
  
  LongWashout_vs_Control =
    long_washout - control,
  
  LongWashout_vs_DTP =
    long_washout - DTP,
  
  levels = design
)


# ============================================================
# Apply contrasts
# ============================================================

fit2 <- contrasts.fit(
  fit,
  contrasts
)

fit2 <- eBayes(
  fit2,
  robust = TRUE
)


# ============================================================
# Save results
# ============================================================

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
  "differential expression analysis finished.\n"
)