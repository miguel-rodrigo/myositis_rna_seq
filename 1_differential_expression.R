library(data.table)
library(DESeq2)

source("utils_get_differential_expression.R")

analyze.by.autoantibody <- create.rnaseq.analysis(aggregation.level = "autoantibody")
analyze.by.disease <- create.rnaseq.analysis(aggregation.level = "disease")

results.by.autoantibody <- analyze.by.autoantibody()
results.by.disease <- analyze.by.disease()


