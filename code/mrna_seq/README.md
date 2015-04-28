# mRNAseq benchmarks
#### boilerplate disclaimer :)
__As for other benchmarks, we are deliberately avoiding *"(pre-)processing"* steps, instead focussing on statistic analyses typical for this datatype. We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance. Suggestions for datasets or methods are welcome.__

Benchmark for expression data (mRNAseq), differential expression.

Data
-----------
? https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/

Code
-----------
- data loading (edgeR)[http://www.bioconductor.org/packages/release/bioc/html/edgeR.html]
- processing as described in (limma user guide)[http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf]

Aims
-----------
- "Chocolate" R code, using as little additional data containers as possible

- "Native" R code, written for clarity and function (a la Analyst) rather than optimization and performance

- clock timings for data loading, distance matrix calculation and various clustering methods
