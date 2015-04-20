# Genetics benchmarks
#### boilerplate disclaimer :)
__As for other benchmarks, we are deliberately avoiding *"(pre-)processing"* steps, instead focussing on statistic analyses typical for this datatype. We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance. Suggestions for datasets or methods are welcome.__

Benchmarks for typical genetics data, focussing on the two main types of studies, familial/pedigree studies and population studies.

## Population study

Using data from the [AML paper](http://www.nejm.org/doi/full/10.1056/NEJMoa1301689), authored by the TCGA consortium, we aim to reproduce part of the analyses, focussing on the genetics data (i.e. figure 1 in the paper).

### data source
- [publication archive](http://tcga-data.nci.nih.gov/docs/publications/laml_2012/)

- [maf (mutations and annotations)](http://tcga-data.nci.nih.gov/docs/publications/laml_2012/genome.wustl.edu_LAML.IlluminaGA_DNASeq.Level_2.2.12.0.tar.gz)

- [patient meta data](http://tcga-data.nci.nih.gov/docs/publications/laml_2012/clinical_patient_laml.tsv)

- [paper]

### analysis
1. import data
2. mutation rate per gene, per disease stage


## familial monogenic data

### data source

Data was obtained from the nice people at [Genomes Unzipped](http://genomesunzipped.org/members), who make [their own genomic data](http://genomesunzipped.org/data) publicly available. 

The dataset we are using comes from the [23andme v2](https://www.23andme.com) sequencing service.

though not related, this data can still be used to perform some typical tests carried out on pedigree studies, such as determining "relatedness" between individuals.

member | dataset id | link
--- | --- | ---
Daniel MacArthur | DGM001 | http://s3.amazonaws.com/gnz.genotypes/DGM001_genotypes.zip
Luke Jostins | LXJ001 | http://s3.amazonaws.com/gnz.genotypes/LXJ001_genotypes.zip
Dan Vorhaus | DBV001 | http://s3.amazonaws.com/gnz.genotypes/DBV001_genotypes.zip
Caroline Wright | CFW001 | http://s3.amazonaws.com/gnz.genotypes/CFW001_genotypes.zip
Kate Morley  | KIM001 | http://s3.amazonaws.com/gnz.genotypes/KIM001_genotypes.zip
Vincent Plagnol | VXP001 | http://s3.amazonaws.com/gnz.genotypes/VXP001_genotypes.zip
Jeff Barrett | JCB001 | http//s3.amazonaws.com/gnz.genotypes/JCB001_genotypes.zip
Jan Aerts  | JXA001  | http://s3.amazonaws.com/gnz.genotypes/JXA001_genotypes.zip
Joe Pickrell | JKP001 | http://s3.amazonaws.com/gnz.genotypes/JKP001_genotypes.zip
Don Conrad | DFC001 | http://s3.amazonaws.com/gnz.genotypes/JKP001_genotypes.zip
Carl Anderson | CAA001 | http://s3.amazonaws.com/gnz.genotypes/CAA001_genotypes.zip
Ilana Fisher | IPF001 | http://s3.amazonaws.com/gnz.genotypes/IPF001_genotypes.zip

### Analysis
- import data, remove autosomes, remove non-common SNPs
- summarize allele frequencies
- identity by descent:

- break genome into blocks
- 

