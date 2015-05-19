# Protein benchmarks
#### boilerplate disclaimer :)
__As for other benchmarks, we are deliberately avoiding *"(pre-)processing"* steps, instead focussing on statistic analyses typical for this datatype. We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance. Suggestions for datasets or methods are welcome.__

Benchmark for protein data, including RPPA and Mass Spectrometry.

Data
-----------
- RPPA
 - <a href="http://www.cell.com/cell/abstract/S0092-8674(14)00876-9">Hoadley et al</a> were cool enough to not only publish their input data, they also deposit the output classifications for comparison.

 - <a href="https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/">Data at TCGA</a>

- Mass Spectrometry
 - Quantitative (iTRAQ) Data from [Bantscheff et al](http://www.ncbi.nlm.nih.gov/pubmed/17721511), sourced from the [PRIDE archvive](http://www.ebi.ac.uk/pride/archive/) (dataset [PRD000032](http://www.ebi.ac.uk/pride/archive/projects/PRD000032))

Code
-----------
- RPPA
 - data loading
 - compute distance matrix as described in paper
 - several methods of unsupervised clustering as found in the base stats package
 
- MS (Thanks to [Yann Abraham](https://github.com/yannabraham) for suggestions and contributions)
 - data loading
 - analysis roughly based on the [Bioconductor proteomics workflow](http://www.bioconductor.org/help/workflows/proteomics/)
 
- TO DO
 - comparison of results to published and/or expected results

Aims
-----------
- "Chocolate" R code, using as little additional data containers as possible
- "Native" R code, written for clarity and function (a la Analyst) rather than optimization and performance
- clock timings for data loading, distance matrix calculation and various clustering methods
- summarised results for checking stability of computations across mulitple versions of R, Renjin, etc or across multiple computation runs.
