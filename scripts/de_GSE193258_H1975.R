library(edgeR)
library(limma)

# -----------------------------------
# 1. Create output folders
# -----------------------------------

dir.create(
  "results/differential_expression/H1975",
  recursive = TRUE,
  showWarnings = FALSE
)

dir.create(
  "figures/differential_expression",
  recursive = TRUE,
  showWarnings = FALSE
)


# -----------------------------------
# 2. Read estimated counts
# -----------------------------------

counts <- read.delim(
  "data/raw/GSE193258/GSE193258_RNAseq_estimated_counts.tsv.gz",
  check.names = FALSE
)

rownames(counts) <- counts$gene
counts$gene <- NULL


# -----------------------------------
# 3. Read metadata
# -----------------------------------

metadata <- read.csv(
  "data/metadata/GSE193258_metadata.csv"
)


# -----------------------------------
# 4. Select H1975 samples
# -----------------------------------

meta_h1975 <- metadata[
  metadata$cell_line == "H1975",
]

group <- factor(
  meta_h1975$treatment_state,
  levels = c(
    "control",
    "acute",
    "DTP",
    "short_washout",
    "long_washout"
  )
)

sample_names <- meta_h1975$expression_sample_name

counts_h1975 <- counts[, sample_names]


# -----------------------------------
# 5. Check sample alignment
# -----------------------------------

stopifnot(
  all(colnames(counts_h1975) == sample_names)
)

print("Samples matched successfully.")
print(table(group))


# -----------------------------------
# 6. Create edgeR object
# -----------------------------------

y <- DGEList(
  counts = counts_h1975,
  group = group
)


# -----------------------------------
# 7. Remove very low-expression genes
# -----------------------------------

keep <- filterByExpr(
  y,
  group = group
)

y <- y[keep, , keep.lib.sizes = FALSE]

cat(
  "Genes retained after filtering:",
  nrow(y),
  "\n"
)


# -----------------------------------
# 8. Normalize libraries
# -----------------------------------

y <- calcNormFactors(y)


# -----------------------------------
# 9. Build design matrix
# -----------------------------------

design <- model.matrix(
  ~ 0 + group
)

colnames(design) <- levels(group)

print("Design matrix:")
print(design)

write.csv(
  design,
  "results/differential_expression/H1975/design_matrix.csv"
)


# -----------------------------------
# 10. Apply voom transformation
# -----------------------------------

png(
  "figures/differential_expression/H1975_voom.png",
  width = 1000,
  height = 800
)

v <- voom(
  y,
  design,
  plot = TRUE
)

dev.off()


# -----------------------------------
# 11. Fit linear model
# -----------------------------------

fit <- lmFit(
  v,
  design
)


# -----------------------------------
# 12. Define biological comparisons
# -----------------------------------

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


# -----------------------------------
# 13. Apply contrasts
# -----------------------------------

fit2 <- contrasts.fit(
  fit,
  contrasts
)

fit2 <- eBayes(
  fit2,
  robust = TRUE
)


# -----------------------------------
# 14. Save differential expression results
# -----------------------------------

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
    paste0(
      "results/differential_expression/H1975/",
      comparison,
      ".csv"
    ),
    row.names = FALSE
  )
  
  significant <- sum(
    result$adj.P.Val < 0.05,
    na.rm = TRUE
  )
  
  cat(
    comparison,
    ":",
    significant,
    "genes with FDR < 0.05\n"
  )
}

print("H1975 differential expression analysis finished.")