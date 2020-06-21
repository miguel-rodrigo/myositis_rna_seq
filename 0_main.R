library(data.table)
library(DESeq2)

source("1_differential_expression.R")
source("2_build_results_table.R")
source("3_plots.R")

#### 1. Analisis por enfermedad y filtrar los que me salen + los que salen en literatura ####
analyze.by.disease <- create.rnaseq.analysis(aggregation.level = "disease")
results.by.disease <- analyze.by.disease()

table.by.disease <- createResultsTable(results.by.disease)

mirna.genes <- fread("Datos/all_mir_genes.csv")$mirgenes
onlymirs.disease <- table.by.disease[gene.names %in% mirna.genes]

set.re.adjusted.pvalues(onlymirs.disease)
getGenesQualityColumns(onlymirs.disease)

# We take those which are relevant for at least 1 disease
relevant.mirs <- onlymirs.disease[n_relevant >= 1, gene.names]
literature.mirs <- fread("Datos/improved_literature_list.csv")$gene_name

mirs.in.common.disease <- intersect(relevant.mirs, literature.mirs)


# which are the relevant ones
# --- union ---
# which are the literature ones
filtered.disease <- onlymirs.disease[gene.names %in% union(relevant.mirs, literature.mirs)]

set.re.adjusted.pvalues(filtered.disease)
getGenesQualityColumns(filtered.disease)

plot.disease <- create.results.heatmap(filtered.disease)


#### 2. Como el (1) pero con autoanticuerpos ####
analyze.by.autoantibody <- create.rnaseq.analysis(aggregation.level = "autoantibody")
results.by.autoantibody <- analyze.by.autoantibody()

table.by.aas <- createResultsTable(results.by.autoantibody)

mirna.genes <- fread("Datos/all_mir_genes.csv")$mirgenes
onlymirs.aas <- table.by.aas[gene.names %in% mirna.genes]  # Check where mirrgenes is coming from

# which are the relevant ones
# --- union ---
# which are the literature ones
set.re.adjusted.pvalues(onlymirs.aas)
getGenesQualityColumns(onlymirs.aas)

# We take those which are relevant for at least 1 autoantibody
relevant.mirs <- onlymirs.aas[n_relevant >= 1, gene.names]
literature.mirs <- fread("Datos/literature_mir_genes.csv")$gene_name

mirs.in.common.aas <- intersect(relevant.mirs, literature.mirs)

filtered.aas <- onlymirs.aas[gene.names %in% union(relevant.mirs, literature.mirs)]
set.re.adjusted.pvalues(filtered.aas)
getGenesQualityColumns(filtered.aas)

plot.aas <- create.results.heatmap(filtered.aas)

# TODO: Reshape into markdown to show things more easily
