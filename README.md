Asthma Transcriptomic Analysis
======

Authors: Mengyuan Kan, Maya Shumyatcher, Blanca Himes

## Introduction
We integrated asthma-related publicly available datasets to investigate global and cell-specific gene expression signatures.

25 asthma-related datasets (23 microarray and 2 RNA-Seq datasets) were downloaded from the Gene Expression Omnibus (GEO) or the Sequence Read Archive (SRA). A total of 17 asthma vs. non-asthma and 13 glucocorticoid vs. control comparisons from the transcriptomic studies were selected for integration analyses and analyzed by [RAVED](https://github.com/HimesGroup/raved) (https://github.com/HimesGroup/raved). We integrated the differential expression results from individual studies using three summary statistics-based methods (i.e effect size-, p-value-, and rank-based methods). We integrated studies across all tissue and cell types to identify genes that globally differentially expressed across asthma or glucocorticoid exposure conditions. We also performed separate analyses for blood and structural cells to identify cell-specific expression patterns.

## Description of files in repository
`integration_results folder` contains integration results for significant genes in asthma and glucocorticoid response. For each phenotype and cell/tissue group of interest, genes with q-value <0.05 obtained by at least two integration methods were considered significant. The files share the same format. Here are the descriptions of the header:

Field | Description
--- | ---
Gene | gene symbol
N | number of studies
es_method.pval | p-value from effect size-based method
es_method.qval | q-value from effect size-based method
es_method.eff.size | combined effect size from effect size-based method
es_method.std.err | standard error from effect size-based method
p_method.pval-p | value from p-value-based method
p_method.qval | q-value from p-value-based method
rank_method.pval | p-value from rank-based method
rank_method.qval | q-value from rank-based method
