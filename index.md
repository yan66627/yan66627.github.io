---
layout: single
title: "Yan Jin"
author_profile: true
---

Welcome!

I am Yan Jin, a recent M.S. graduate in Biostatistics from Washington University in St. Louis. I received my B.S. in Biological Sciences from the University of California, Irvine in 2023.

My research interests include:

- Aging and longevity genomics
- Spatial transcriptomics
- Single-cell multi-omics
- Epigenomics & regulatory genomics
- Statistical & computational genomics

I apply statistical modeling, probabilistic inference, and computational genomics tools to study biological mechanisms and promote meaningful biological discoveries.

🔬 **Research Interests**

- Statistical Genetics & Genomics
- Spatial transcriptomics / scRNA-seq
- Multi-omics integration
- EM algorithms & mixture models
- High-dimensional data analysis

📂 **Selected Projects**

**EM Algorithm for Genotype Frequency Estimation**
Custom R implementation of an EM algorithm that estimates minor allele frequency (MAF) under Hardy–Weinberg equilibrium. Implements full E-step/M-step iterations with log-likelihood convergence tracking.
`R` · Statistical Genetics · [`code`](https://github.com/yan66627/yan66627.github.io/blob/master/_portfolio/em_algorithm.R)

---

**scRNA-seq Integration & Comparative Analysis**
End-to-end single-cell RNA-seq pipeline comparing PCKO and PE2F1DKO lung samples. Includes doublet detection (scDblFinder), SCTransform normalization, UMAP/clustering, cell-type annotation via Azimuth, and cell proportion testing.
`R` · Seurat · scDblFinder · ComplexHeatmap · scProportionTest

---

**Spatial Transcriptomics — Visium HD**
Analysis of high-definition spatial transcriptomics data (8 µm / 16 µm bin resolution) from mouse lung CKO samples. Implements sketch-based downsampling for large datasets and cell-type deconvolution with UCell marker scoring.
`R` · Seurat · UCell · Visium HD

---

**Transcription Factor Motif Enrichment Pipeline**
Batch GO enrichment analysis pipeline for TF motif target genes. Reduces redundant GO terms via semantic similarity (rrvgo), generates heatmap/scatter/treemap visualizations, and builds Cytoscape-ready network tables.
`R` · rrvgo · org.Mm.eg.db · [`code`](https://github.com/yan66627/yan66627.github.io/blob/master/_portfolio/go_redundancy_utils.R)

---

**Genomic Motif Coordinate Mapping**
Utility functions to map HOMER/FIMO motif hits back to precise genomic coordinates using peak center + offset, with strand-aware endpoint adjustment and IRanges-based overlap deduplication.
`R` · IRanges · Bioconductor · [`code`](https://github.com/yan66627/yan66627.github.io/blob/master/_portfolio/motif_location_utils.R)

📄 **CV**  
My latest CV is available in the CV tab on the navigation bar.

📫 **Contact**  
Email: yjin66627@gmail.com  
GitHub: yan66627  
LinkedIn: http://linkedin.com/in/yan-jin-080246236 

Feel free to reach out for research discussions or collaboration!
