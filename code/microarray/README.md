# microarray benchmarks
#### boilerplate disclaimer :)
__We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may be very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance.__

__Suggestions for datasets or methods are welcome.__

Benchmark for expression data (microarray), differential expression. Unlike other benchmarks, we are including some *"(pre-)processing"* steps for this datatype, as they are so standardised and relatively stable.

Data
-----------
Data was sourced from the [GEO repository](http://www.ncbi.nlm.nih.gov/geo/), specifically the data for [Ramsey and Fontes, 2013](http://www.ncbi.nlm.nih.gov/pubmed/23954399), as this was a relatively small and simple experimental design (2x2 design matrix, 3 samples per condition):
- [dataset](http://www.ncbi.nlm.nih.gov/sites/GDSbrowser?acc=GDS5070)
- [data](http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE45417)
- [full article](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3783546/)


Code
-----------
- data loading
- processing as according to [limma manual](http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf)
    - scale and normalise
    - filter
    - fitting linear model
    - extraction of top differentially expressed genes
- Gene set enirchment testing
 - random examples, as found in [limma manual](http://www.bioconductor.org/packages/release/bioc/manuals/limma/man/limma.pdf)
 - [CAMERA](http://nar.oxfordjournals.org/content/early/2012/05/24/nar.gks461.full)


Aims
-----------
- "Chocolate" R code, using as little additional data containers as possible

- "Native" R code, written for clarity and function (a la Analyst) rather than optimization and performance

- clock timings for data loading, distance matrix calculation and various clustering methods
