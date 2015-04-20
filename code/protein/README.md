# Protein benchmarks
#### boilerplate disclaimer :)
__As for other benchmarks, we are deliberately avoiding *"(pre-)processing"* steps, instead focussing on statistic analyses typical for this datatype. We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance. Suggestions for datasets or methods are welcome.__

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
