# Disentangling Generic Drug Persistence from Osimertinib-Associated Resistance Commitment in EGFR-Mutant NSCLC Using Cross-Dataset Interpretable Machine Learning

## Overview

This project investigates the transition from drug-tolerant persistence to acquired osimertinib resistance in EGFR-mutant non-small cell lung cancer (NSCLC).

The central goal is to distinguish:

- broadly shared and potentially reversible drug-persistence programs
- persistent molecular programs that remain after drug withdrawal and may later be associated with acquired osimertinib resistance

The project integrates multiple public transcriptomic datasets representing different treatment states and will ultimately use cross-dataset and interpretable machine-learning approaches to identify reproducible molecular programs.

## Research Question

Can cross-dataset transcriptomic analysis distinguish broadly shared drug-persistence programs from molecular programs that persist after drug withdrawal and recur or strengthen during acquired osimertinib resistance?

## Working Hypothesis

DTP-associated transcriptional changes contain at least two components:

1. a broadly shared and predominantly reversible drug-persistence program
2. a more persistent osimertinib-associated program that remains after drug withdrawal and may be associated with progression toward acquired resistance

## Analysis Strategy

The planned analysis follows the biological sequence:

Control → Acute Treatment → DTP → Drug Withdrawal → Acquired Resistance → Independent Validation → Clinical Validation

The main computational stages include:

- sample metadata construction and validation
- expression-data quality control
- PCA and exploratory transcriptomic analysis
- differential expression analysis
- DTP-to-washout trajectory analysis
- cross-cell-line persistence analysis
- generic persistence analysis across datasets
- resistance-progression analysis
- independent resistant-model validation
- patient-level clinical validation
- interpretable machine-learning analysis

## Datasets

| Dataset | Role in Project |
|---|---|
| GSE193258 | Core DTP and drug-withdrawal analysis |
| GSE255958 | Generic drug-persistence analysis |
| GSE254485 | Resistance-progression time-course analysis |
| GSE202859 | Independent acquired-resistance validation |
| GSE262582 | Independent HCC827 resistance validation |
| GSE272162 | Independent H1975 resistance validation |
| GSE289619 | Clinical validation in osimertinib-treated patients |

Raw and processed expression files are obtained from NCBI GEO and are not stored in this repository.

## Current Analysis: GSE193258

GSE193258 contains four EGFR-mutant NSCLC cell-line models:

- PC9
- H1975
- HCC827
- HCC2935

Each model contains samples representing:

Control, acute osimertinib treatment, DTP, short drug washout, and long drug washout.

### Initial Quality Control

The expression matrix contained:

- 19,712 genes
- 60 samples
- no missing values
- no duplicated gene identifiers

All 60 metadata samples were successfully matched to the corresponding expression-matrix columns.

### PCA

Initial PCA showed that cell-line identity is a major source of transcriptomic variation.

PCA was therefore repeated separately within each cell line to examine treatment-associated changes independently of baseline cell-line differences.

### Differential Expression

Differential expression analysis is performed using the limma-voom framework.

The main comparisons include:

- Acute vs Control
- DTP vs Control
- DTP vs Acute
- Short Washout vs Control
- Long Washout vs Control
- Long Washout vs DTP

### DTP Persistence Analysis

For DTP-responsive genes, transcriptional persistence after drug withdrawal is evaluated using:

- direction consistency
- short-washout retention
- long-washout retention
- DTP-to-washout effect-size correlation

The four cell lines show different degrees of transcriptomic persistence.

| Cell Line | Median Short Retention | Median Long Retention | Short Spearman | Long Spearman |
|---|---:|---:|---:|---:|
| H1975 | 0.819 | 1.021 | 0.961 | 0.978 |
| PC9 | 0.781 | 0.541 | 0.964 | 0.941 |
| HCC827 | 0.950 | 0.937 | 0.981 | 0.975 |
| HCC2935 | 0.843 | 0.745 | 0.947 | 0.927 |

These results suggest that persistence is not binary and varies substantially across cellular backgrounds.

### Cross-Cell-Line Persistence

An exploratory persistence-support definition was applied to identify recurrent DTP-associated changes across the four cell lines.

Current results:

- 366 genes show persistence support in 4/4 cell lines
- 1,538 genes show persistence support in 3/4 cell lines
- 3,286 genes show persistence support in 2/4 cell lines
- 5,058 genes show persistence support in 1/4 cell lines

The current persistence-support threshold is exploratory and will be evaluated through sensitivity analyses before final use.

Importantly, these genes are currently considered persistent DTP-associated candidates and are not yet interpreted as acquired-resistance genes.

## Next Step

The next stage will analyze GSE255958 to determine which persistent programs identified in GSE193258 are broadly shared across different drug-tolerant states.

This will help separate generic drug-persistence programs from molecular signals that may be more specifically associated with the osimertinib resistance trajectory.

## Repository Structure

```text
data/
    metadata/

figures/
    differential_expression/
    pca_by_cell_line/

results/
    differential_expression/
    pca/
    pca_distances/
    persistence/

scripts/
```

Large raw datasets and local processed-data files are excluded from version control.

## Tools

Python is used for data handling, metadata processing, quality control, PCA, cross-dataset analysis, visualization, and later machine-learning analyses.

R is currently used where methodologically appropriate for differential expression analysis with limma-voom.

## Status

Work in progress.

The current repository contains the completed first-stage analysis of GSE193258. Additional datasets and validation stages will be added progressively.
