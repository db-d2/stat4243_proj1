---
layout: default
title: Home
---

# Comparative Analysis of scRNA-seq Clustering Methods

**Systematic evaluation of clustering algorithms, denoising techniques, and batch integration for single-cell RNA sequencing**

[View Code](https://github.com/db-d2/stat4243_proj1/blob/main/code/ROC_xenopus_colab_final.ipynb) | [Open in Colab](https://colab.research.google.com/github/db-d2/stat4243_proj1/blob/main/code/ROC_xenopus_colab_final.ipynb) | [Write-up Markdown version](./Writeup.md) | [Write-up PDF version](./Writeup.pdf)

## Abstract

This project systematically evaluates computational methods for identifying Regeneration-Organizing Cells (ROCs) in *Xenopus* tail tissue using single-cell RNA sequencing data from 13,199 cells. I implemented two clustering algorithms (Walktrap and Leiden), two marker selection methods (logistic regression and Wilcoxon rank-sum test), two denoising techniques (PCA reconstruction and k-nearest neighbor smoothing), and two batch integration approaches (Harmony and BBKNN). The analysis demonstrates that method selection substantially impacts ROC identification accuracy, with k-NN smoothing denoising providing 54% improvement in marker validation and batch integration achieving 8-10x clustering quality improvements.

## Key Results

### Methods Comparison
- **Clustering Performance**: Leiden outperforms Walktrap (0.132 vs 0.046 silhouette, ARI=0.330)
- **Denoising Impact**: k-NN smoothing achieves 54% marker improvement (+14 genes), PCA reconstruction shows minimal benefit (+2 genes)
- **Batch Integration**: Harmony 10.3x and BBKNN 8.0x clustering improvements, but no marker validation benefit
- **Optimal Pipeline**: Leiden + k-NN smoothing + Harmony for best ROC detection

### Biological Validation
- **ROC populations** identified through enrichment analysis in key clusters (23, 33, 6, 41, 38)
- **19-26 genes** baseline overlap with published markers depending on selection threshold
- **40 markers** validated after optimal preprocessing (k-NN smoothing)

## Main Figures

### ROC Identification
![ROC Identification UMAP](./figures/figure1_roc_identification_improved.png)
*UMAP visualization showing 214 computationally identified ROC cells (red) among 13,199 total cells*

### Marker Validation
![Venn Diagram](./figures/figure_venn_validation.png)
*Overlap between computational markers and published ROC signatures (p=7.99×10⁻⁴⁰)*

### Method Optimization
![Denoising Comparison](./figures/denoising_comparison.png)
*Comparison of clustering quality across denoising methods*

### Gene Expression
![Gene Expression Heatmap](./figures/gene_expression_heatmap.png)
*Differential expression of top ROC markers*

## Documentation

**Technical Write-up**: [Markdown version](./Writeup.md) | [PDF version](./Writeup.pdf)

Complete methods, results, and analysis including all equations, tables, and detailed methodology.

**Supplementary Materials**: [View all figures and visualizations](./Supplementary_Materials.pdf)

Comprehensive collection of clustering visualizations, method comparisons, ROC identification plots, and validation figures generated during the analysis.

## Data Availability

- **Processed dataset**: `cleaned_processed_frogtail.h5ad` (1.2GB)
- **Published markers**: `aav9996_tables3.xlsx` (Aztekin et al. 2019)
- **Analysis notebook**: `ROC_xenopus_colab_final.ipynb`
- **Technical report**: [Writeup_Updated.pdf](./Writeup.pdf)

## Citation

Original data from: Aztekin et al. (2019) Science 364:653

---
