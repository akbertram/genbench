protein
==========
Benchmark for protein data (RPPA), and unsupervised clustering as found in the base stats package.


Data
-----------
<a href="http://www.cell.com/cell/abstract/S0092-8674(14)00876-9">Hoadley et al</a> were cool enough to not only publish their input data, they also deposit the output classifications for comparison.

- <a href="https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/">Data at TCGA</a>

Code
-----------
- data loading
- compute distance matrix as described in paper
- several methods of unsupervised clustering
- comparison of results to published and/or expected results

Aims
-----------
- "Chocolate" R code, using as little additional data containers as possible

- "Native" R code, written for clarity and function (a la Analyst) rather than optimization and performance

- clock timings for data loading, distance matrix calculation and various clustering methods

- output comparisons to Hoadley et al
