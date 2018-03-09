Asthma Transcriptomic Analysis
======

Authors: Mengyuan Kan, Maya Shumyatcher, Blanca Himes

## Introduction
We integrated asthma-related publicly available datasets to investigate global and cell-specific gene expression signatures.

25 asthma-related datasets (23 microarray and 2 RNA-Seq datasets) were downloaded from the Gene Expression Omnibus (GEO) or the Sequence Read Archive (SRA). A total of 17 asthma vs. non-asthma and 13 glucocorticoid vs. control comparisons from the transcriptomic studies were selected for integration analyses and analyzed by [RAVED](https://github.com/HimesGroup/raved) (https://github.com/HimesGroup/raved). We integrated the differential expression results from individual studies using three summary statistics-based methods. We integrated studies across all tissue and cell types to identify genes that globally differentially expressed across asthma or glucocorticoid exposure conditions. We also performed separate analyses for blood and structural cells to identify cell-specific expression patterns.

## Cell-specific integration
![](<./figs/GC_bloodcell.png>)

## Description of files in repository
`integration_results folder` contains integration results for significant genes in asthma and glucocorticoid response.
