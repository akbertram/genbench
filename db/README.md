# genbench/db/

GenBench: Realworld Genomics Benchmarks for R, forked from <a href= "https://github.com/hannesmuehleisen/genbase">GenBase</a> and inspired by [the original GenBase](https://github.com/mitdbg/genbase).

The following boilerplate disclaimer applies to everything here :)

__As for other benchmarks, we are deliberately avoiding *"(pre-)processing"* steps, instead focussing on statistic analyses typical for this datatype. We do not endorse any of the methods used as *"standards"* or *"recommended"*, in fact, because we aim to start simple and avoid as far as possible non-essential packages, methods may very much not recommended. Future updates will implement more advanced methods, i.e. code and datasets are simply intended to represent good, ecologically viable tests of performance.__ 

__Suggestions for datasets or methods are welcome.__

Developed in collaboration with <a href= "https://www.bedatadriven.com">BeDataDriven</a>.

## `upload_benchmarks.R`
This script uploads locally cached files to a database from locally cached files in `/generated/timings/`. Several commandline arguments can be passed:
* Connection details **must** be passed:
  * `upload_benchmarks.R --args --usr=foo --pwd=bar --conn=baz`
* Drop any existing tables in the database and recreate schema:
  * `upload_benchmarks.R --args --create_new`
