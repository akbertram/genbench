# genbench

GenBench: Realworld Genomics Benchmarks for R, forked from <a href= "https://github.com/hannesmuehleisen/genbase">GenBase</a> and inspired by [the original GenBase](https://github.com/mitdbg/genbase).


### Why do we need new benchmarks?

NOTE: This benchmark does not cover all possible operations performed in genomics. In particular, we chose to not focus on processing of raw sequence data and instead focus on higher level processing.

We also chose to use real-world data, to allow not just benchmarking of clock speeds upon completion of a given analysis, but also to allow for testing of "correctness" and to ensure ecological validity of tests.

The following boilerplate disclaimer applies to everything here :)

__As for other benchmarks, we are deliberately avoiding *"(pre-)processing"* steps, instead focussing on statistic analyses typical for this datatype. We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance.__ 

__Suggestions for datasets or methods are welcome.__


## Benchmarks

This benchmark was developed in collaboration with <a href= "https://www.bedatadriven.com">BeDataDriven</a>.

### Data:

This benchmark focuses on genomics data, and specifically on statistical analysis, as opposed to processing (e.g. alignment).

All data is publicly available, sourced from various places (for further information please see the individual README accompanying each benchmark, and please let us know if something is incorrect/missing/...):

(a) Gene Expression Data
  - 3' Microarray
    - [Ramsey and Fontes, 2013](http://www.ncbi.nlm.nih.gov/pubmed/23954399)
      - [dataset](http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5070)
      - [data](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45417)
      - [full article](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3783546/)

  - mRNAseq
    - Data sourced from [ReCount](http://bowtie-bio.sourceforge.net/recount/)
      - Specifically, the "gilad" study, [PMID: 20009012](http://www.ncbi.nlm.nih.gov/pubmed?term=20009012)

(b) Protein Expression Data
  - RPPA platform
    - Data was sourced from [Hoadley et al](http://www.cell.com/cell/abstract/S0092-8674(14)00876-9), using the [TCGA portal](https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/)
  - Mass Spectrometry is not yet implemented

(c) Genetics Data
  - Population studies were simulated using [data](http://tcga-data.nci.nih.gov/docs/publications/laml_2012/) from the [AML paper](http://www.nejm.org/doi/full/10.1056/NEJMoa1301689) authored by the TCGA consortium
  - Pedigree studies were simulated using data obtained from the nice people at [Genomes Unzipped](http://genomesunzipped.org/members), who make [their own genomic data](http://genomesunzipped.org/data) publicly available

(d) Simulated matrices (to allow for testing scale, as in orginal [GenBase](https://github.com/mitdbg/genbase) project 

(e) Clinical data and Data Integration
  - Data sourced from package [ncvreg](http://cran.r-project.org/web/packages/ncvreg/ncvreg.pdf), further references below:

    - heart dataset
      - Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning. Springer. 
      - Rousseauw, J., et al. (1983). Coronary risk factor screening in three rural communities. South African Medical Journal, 64, 430-436.
    - prostate dataset
      - Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning. Springer. 
      - Stamey, T., et al. (1989). Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate. II. Radical prostatectomy treated patients. Journal of Urology, 16: 1076-1083.
    - lung dataset
      - package [survival](http://CRAN.R-project.org/package=survival)
      - Kalbfleisch D and Prentice RL (1980), The Statistical Analysis of Failure Time Data. Wiley, New York.

  - [gene RIFs](http://www.ncbi.nlm.nih.gov/gene/about-generif) provided interaction data to be used for graphical modelling with [igraph](http://igraph.sourceforge.net/)
  - We also reproduce aspects of integrative analysis carried out in the Human Liver Cohort project:
    - [synapse entry](https://www.synapse.org/#!Synapse:syn299418)
    - [Schadt et al, 2008](http://journals.plos.org/plosbiology/article?id=10.1371/journal.pbio.0060107)

### The Code

Scripts can be run individually, or all together using `run_benchmarks.R` (see below for more details). All generated output, including timings for individual microbenchmark blocks and results calculations, will end up in `/generated/*`.

See individual `/code/*/README.md` files for details on the individual benchmarks, but here is a general overview of topics covered so far:

(a) Gene Expression Data
  - 3' Microarray
Methods typical for microarrays including RMA and MAS150 normalisation, differential expression with limma, and some gene-set tests.

  - mRNAseq
Methods typical for mRNAseq expression analyses, including the Voom/edgeR/limma approach (The "DEseq" approach is to come).

(b) Protein Expression Data
Focussing on classification, RPPA [data](https://tcga-data.nci.nih.gov/docs/publications/TCGApancan_2014/) from [Hoadley et al](http://www.cell.com/cell/abstract/S0092-8674(14)00876-9) was used a range of unpsupervised clustering methods including hierarchical, kmeans, random forrest, bayesian

(c) Genetics Data
This section is still quite under-developed, focussing on summarisation of allele frequencies, and some sliding window methods. Implementation of "proper" genetics analysis methods is a burning priority!

(d) Simulated data
Focus on the linear algebra and stats operations below: 

  - Linear Regression: build regression model to predict drug response from expression data

  - Covariance: determine which pairs of genes have expression values that are correlated

  - SVD: reduce the dimensionality of the problem to the top 50 components

  - Biclustering: simultaneously cluster rows and columns in the expression matrix to find related genes

  - Statistics: determine if certain sets of genes are highly expressed compared to the entire set of genes

(e) Clinical data and Data Integration
Alongside some graphical methods using iGraph, we focus on machine learning approaches including clustering, and prediction using naive Bayes and robust linear model approaches

### `run_benchmarks.R`
This wrapper script will install all required packages, and run all benchmark scripts in `/code/*/`. Several commandline arguments can be passed:
* Controlling the number of times each benchmark is run:
  * `run_benchmarks.R --args [an integer]`
* Resetting the locally cached timings:
  * `run_benchmarks.R --args --reset`
  
### `examine_benchmarks.R`
This script produces some overview plots of current data. Either from a database or from locally cached files in `/generated/timings/`. Several commandline arguments can be passed, to activate the use of a database (local files are used by default):
* Controlling the number of times each benchmark is run:
  * `examine_benchmarks.R --args --use-db --usr=foo --pwd=baz --conn=bar`

Uploading locally cached files to a database is done using the code in `/db/upload_benchmarks.R`.

### Systems

As part of this work, so far, we have tested the benchmarks on a variety of systems, specifically:
- GNU R (version R-2.14.2, R-2.15.3, R-3.0.3, R-3.1.3, R-2.0)
- Renjin (forked from GNU R 2.14.2)


