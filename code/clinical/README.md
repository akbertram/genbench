# Genetics benchmarks
#### boilerplate disclaimer :)
__As for other benchmarks, we are deliberately avoiding *"(pre-)processing"* steps, instead focussing on statistic analyses typical for this datatype. We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance. Suggestions for datasets or methods are welcome.__

Benchmarks for supervised and unsupervised approaches on "clinical" datasets. Replicating examples in [Elements of Statistical Learning](http://www.e-booksdirectory.com/details.php?ebook=3267)

### Data sources
Data sourced from package [ncvreg](http://cran.r-project.org/web/packages/ncvreg/ncvreg.pdf), further references below:

- heart dataset
 - Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning. Springer. 
 - Rousseauw, J., et al. (1983). Coronary risk factor screening in three rural communities. South African Medical Journal, 64, 430-436.
 ```{r}
 data(heart)
 ```
- prostate dataset
 - Hastie, T., Tibshirani, R., and Friedman, J. (2001). The Elements of Statistical Learning. Springer. 
 - Stamey, T., et al. (1989). Prostate specific antigen in the diagnosis and treatment of adenocarcinoma of the prostate. II. Radical prostatectomy treated patients. Journal of Urology, 16: 1076-1083.
 ```{r}
 data(prostate)
 ```
- lung dataset
 - pacakge [survival](http://CRAN.R-project.org/package=survival)
 - Kalbfleisch D and Prentice RL (1980), The Statistical Analysis of Failure Time Data. Wiley, New York.
 ```{r}
 data(Lung)
 ```
 
 ### Analysis
 
 todo
 
 ###