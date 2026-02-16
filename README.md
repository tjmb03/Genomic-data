# Genomic Data Analysis in R

Foundational genomic data analysis workflows in R, focused on RNA-seq differential expression, dimensionality reduction, and Bioconductor/statistical fundamentals.

---

## Contents

### 📄 DESeq2 examples.R  
RNA-seq differential expression workflow using **DESeq2** (design setup, normalization, DE testing, results extraction).

### 📄 PCA.R  
PCA workflow for expression data QC and sample-level structure (variance exploration, clustering, visualization).

### 📄 bioconductor_*.R  
Exercises exploring core **Bioconductor** data structures and interoperability (e.g., common container objects and annotation handling).

### 📄 stat_*.R  
Statistical foundations used in genomics (distributions, hypothesis testing, multiple testing concepts, interpretation).

### 📄 package_related.R  
Notes/utilities related to R package usage for genomic workflows.

---

## Requirements

- R ≥ 4.2  
- Bioconductor (DESeq2)

Install DESeq2:

```r
if (!require("BiocManager")) install.packages("BiocManager")
BiocManager::install("DESeq2")
