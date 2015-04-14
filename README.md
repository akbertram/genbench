genbench
=======

GenBench: Realworld Genomics Benchmarks for R, forked from and inspired by <a href= "https://github.com/hannesmuehleisen/genbase">GenBase</a>.

Why genomics:
=============

Genomics is quickly becoming a major source of big data due to advances in sequencing technology. It is now possible to sequence over 2000 people per day at a sequencing facility. At 3 GB/genome, such a facility produces 6TB of data per day. Furthermore, it has gotten exponentially cheaper to sequence full genomes. In spite of the large amount of data available, we are unable to analyze it at scale.

Why do we need a new benchmark:
===============================

NOTE: This benchmark does not cover all possible operations performed in genomics. In particular, we chose to not focus on processing of raw sequence data and instead focus on higher level processing.

We also chose to use real-world data, to allow not just benchmarking of clock speeds upon completion of a given analysis, but also to allow for testing of "correctness" and to ensure ecological validity of tests.


Benchmark:
==========

This benchmark was developed in collaboration with <a href= "https://www.bedatadriven.com">BeDataDriven</a>.

Data:
-----

This benchmark focuses on genomics data, and specifically on statistical analysis, as opposed to processing (e.g. alignment):

(a) Gene Expression Data
- 3' Microarray
- mRNAseq

(b) Protein Expression Data
 - <a href=""></a>http://www.cell.com/cell/abstract/S0092-8674(14)00876-9

(c) Genetics Data

(d) Simulated matrices (to allow for testing scale, as in orginal <a href="https://github.com/mitdbg/genbase">GenBase</a> Project 

The above will be sourced from:
- <a href="https://tcga-data.nci.nih.gov/docs/publications/">TCGA</a>
- <a href="http://www.ncbi.nlm.nih.gov/gds">GEO</a>
- "simulated" (this is retained from the original <a href="https://github.com/mitdbg/genbase">GenBase</a>)

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