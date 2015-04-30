# mRNAseq benchmarks
#### boilerplate disclaimer :)
__As for other benchmarks, we are deliberately avoiding *"(pre-)processing"* steps, instead focussing on statistic analyses typical for this datatype. We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance. Suggestions for datasets or methods are welcome.__

Benchmark for expression data (mRNAseq), differential expression.

Data
-----------
Data sourced from (ReCount)[http://bowtie-bio.sourceforge.net/recount/], which provide:
~"...an online resource consisting of RNA-seq gene count datasets built using the raw data from 18 different studies. The raw sequencing data (.fastq files) were processed with (Myrna)[http://bowtie-bio.sourceforge.net/myrna/index.shtml] to obtain tables of counts for each gene. For ease of statistical analysis, we combined each count table with sample phenotype data to form an R object of class ExpressionSet. The count tables, ExpressionSets, and phenotype tables are ready to use and freely available here. By taking care of several preprocessing steps and combining many datasets into one easily-accessible website, we make finding and analyzing RNA-seq data considerably more straightforward."~


We used thesmallest dataset at ReCount, which has at least 3 replicates per experimental group:

- the __"gilad"__ study 
- (PMID: 20009012)[http://www.ncbi.nlm.nih.gov/pubmed?term=20009012]
- human data
- 6 samples
- 41,356,738 total aligned reads

Code
-----------
__edgeR_voom.R__
- data loading to DGEList instance from (edgeR)[http://www.bioconductor.org/packages/release/bioc/html/edgeR.html]
- processing as described in the (limma user guide)[http://www.bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf], (limma case study)[http://bioinf.wehi.edu.au/RNAseqCaseStudy/] and the (edgeR user guide)[http://www.bioconductor.org/packages/release/bioc/vignettes/edgeR/inst/doc/edgeRUsersGuide.pdf]

Aims
-----------
- "Chocolate" R code, using as little additional data containers as possible

- "Native" R code, written for clarity and function (a la Analyst) rather than optimization and performance

- clock timings for data loading, distance matrix calculation and various clustering methods
