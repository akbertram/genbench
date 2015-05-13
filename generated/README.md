# generated
Generated results and timings for each set of benchmarks, organised by assay type.

## contents
- generated/results
 - each benchmarks script should write out results to here, <assay type>/<script name>.results.tsv
 - output should match files in "expected", i.e. <assay type>/<script name>.expected.tsv
- generated/timings
 - timestamped clock timings for each script, organised by assay type: <assay type>/<script name>.<timestamp>.tsv
 - expects a dataframe (each row one timing as produced by system.time(), code block names as rownames)