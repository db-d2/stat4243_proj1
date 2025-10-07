# Comparative Analysis of scRNA-seq Clu## Quick Start

### View Analysis
- **Website**: Visit the [GitHub Pages site](https://db-d2.github.io/stat4243_proj1/)
- **Technical Write-up**: [Markdown version](./Writeup.md) | [PDF version](./Writeup.pdf)ing Methods

[![GitHub Pages](https://img.shields.io/badge/docs-GitHub%20Pages-blue)](https://db-d2.github.io/stat4243_proj1/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Systematic comparison of clustering algorithms, denoising techniques, and batch integration methods for single-cell RNA sequencing analysis, applied to Regeneration-Organizing Cell (ROC) identification in *Xenopus laevis* tail regeneration.

## Key Results

### Methods Comparison
- **Walktrap vs Leiden**: Substantial agreement (ARI=0.637, Rand=0.944)
- **Denoising impact**: PCA reconstruction and k-NN smoothing evaluated
- **Batch integration**: Harmony and BBKNN methods compared
- **11-fold improvement** in clustering quality through optimization

### Biological Application
- **214 ROC cells** identified from 13,199 total cells
- **21 genes** validated against published markers (p=7.99×10⁻⁴⁰)

## Repository Contents

```
roc-analysis/
├── README.md                               # This file
├── index.md                                # GitHub Pages homepage
├── _config.yml                             # Jekyll configuration
├── Writeup.md                              # Complete analysis (Markdown)
├── Writeup.pdf                             # Complete analysis (PDF)
├── code/
│   └── ROC_xenopus_colab_final.ipynb     # Analysis notebook
├── data/
│   ├── aav9996_tables3.xlsx              # Published markers (Aztekin et al. 2019)
│   └── cleaned_processed_frogtail.h5ad   # Processed dataset
└── figures/
    ├── denoising_comparison.png
    ├── figure1_roc_identification_improved.png
    ├── figure_venn_validation.png
    └── gene_expression_heatmap.png
```

## Quick Start

### View Analysis
Visit the [GitHub Pages site](https://db-d2.github.io/stat4243_proj1/) or read the [technical write-up](./Writeup.md).

### Run Analysis
Open the notebook in Google Colab:
[ROC_xenopus_colab_final.ipynb](https://colab.research.google.com/github/db-d2/stat4243_proj1/blob/main/code/ROC_xenopus_colab_final.ipynb)

## Methods Summary

### Core Analysis Pipeline
1. **Normalization**: CP10K scaling with log2 transformation
2. **Feature Selection**: Fano factor-based HVG identification (2,308 genes)
3. **Dimensionality Reduction**: PCA (50 components) + UMAP visualization

### Comparative Methods Evaluated
4. **Clustering Algorithms**: 
   - Walktrap (graph-based, fuzzy simplicial sets)
   - Leiden (modularity optimization on kNN graph)
5. **Denoising Approaches**:
   - PCA reconstruction (20 components)
   - k-NN smoothing (k=10 in PCA space)
6. **Batch Integration**:
   - Harmony (PCA-based correction)
   - BBKNN (batch-balanced k-nearest neighbors)
7. **Quality Metrics**: Silhouette scores, ARI, Rand index, modularity

## Citation

Original data: Aztekin et al. (2019) Science 364:653