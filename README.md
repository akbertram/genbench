genbench
=======

GenBench: Realworld Genomics Benchmarks for R, forked from <a href= "https://github.com/hannesmuehleisen/genbase">GenBase</a> and inspired by [the original GenBase](https://github.com/mitdbg/genbase).


Why do we need a new benchmark:
===============================

NOTE: This benchmark does not cover all possible operations performed in genomics. In particular, we chose to not focus on processing of raw sequence data and instead focus on higher level processing.

We also chose to use real-world data, to allow not just benchmarking of clock speeds upon completion of a given analysis, but also to allow for testing of "correctness" and to ensure ecological validity of tests.

The following boilerplate disclaimer applies to everything here :)

__As for other benchmarks, we are deliberately avoiding *"(pre-)processing"* steps, instead focussing on statistic analyses typical for this datatype. We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance.__ 

__Suggestions for datasets or methods are welcome.__


Benchmark:
==========

This benchmark was developed in collaboration with <a href= "https://www.bedatadriven.com">BeDataDriven</a>.

Data:
-----

This benchmark focuses on genomics data, and specifically on statistical analysis, as opposed to processing (e.g. alignment):

All data is publicly available, sourced from:
- <a href="https://tcga-data.nci.nih.gov/docs/publications/">TCGA</a>
- <a href="http://www.ncbi.nlm.nih.gov/gds">GEO</a>
- "simulated" (this is retained from the original <a href="https://github.com/mitdbg/genbase">GenBase</a>)
- [public datasets detailed here](https://github.com/caesar0301/awesome-public-datasets)

(a) Gene Expression Data
- 3' Microarray
- mRNAseq

(b) Protein Expression Data
 - [RPPA platform data](http://www.cell.com/cell/abstract/S0092-8674(14)00876-9)

(c) Genetics Data
- [Population studies]()
- [Pedigree studies]()


(d) Simulated matrices (to allow for testing scale, as in orginal <a href="https://github.com/mitdbg/genbase">GenBase</a> Project 

Ideally, studies are chosen for scale, and for having validated and replicable results.

Benchmarks:
--------

(a) Gene Expression Data
- 3' Microarray

RMA, MAS150

- mRNAseq

Voom/LIMMA, DEseq

(b) Protein Expression Data

Focus on classification:

- RPPA data from <a href="http://www.cell.com/cell/abstract/S0092-8674(14)00876-9">Hoadley et al</a>, <a href="https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/">here</a>

unpsupervised clustering: hierarchical, kmeans, random forrest, bayesian


(c) Genetics Data

TBD

(d) Simulated data

Focus on the linear algebra and stats operations below: 

- Linear Regression: build regression model to predict drug response from expression data

- Covariance: determine which pairs of genes have expression values that are correlated

- SVD: reduce the dimensionality of the problem to the top 50 components

- Biclustering: simultaneously cluster rows and columns in the expression matrix to find related genes

- Statistics: determine if certain sets of genes are highly expressed compared to the entire set of genes


Systems:
--------

As part of this work, we have tested the benchmark queries on a variety of systems, specifically:
- GNU R (version R 2.14.2)
- Renjin (forked from GNU R 2.14.2)


<a href=""></a>
