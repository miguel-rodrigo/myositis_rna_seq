#### 1. Analisis por enfermedad y filtrar los que me salen + los que salen en literatura ####
source("2_build_results_table.R")
table.by.disease <- createResultsTable(results.by.disease)
onlymirs.disease <- table.by.disease[gene.names %in% mirgenes]

# which are the relevant ones
# --- union ---
# which are the literature ones
set.re.adjusted.pvalues(onlymirs.disease)
getGenesQualityColumns(onlymirs.disease)

# We take those which are relevant for at least 1 disease
relevant.mirs <- onlymirs.disease[n_relevant >= 1, gene.names]
literature.mirs <- fread("Datos/literature_mir_genes.csv")$gene_name

mirs.in.common.aas <- intersect(relevant.mirs, literature.mirs)

filtered.disease <- onlymirs.disease[gene.names %in% union(relevant.mirs, literature.mirs)]
set.re.adjusted.pvalues(filtered.disease)
getGenesQualityColumns(onlymirs.disease)


#### 2. Como el (1) pero con autoanticuerpos ####
table.by.aas <- createResultsTable(results.by.autoantibody)
onlymirs.aas <- table.by.aas[gene.names %in% mirgenes]  # Check where mirrgenes is coming from

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
getGenesQualityColumns(onlymirs.aas)

# TODO: Create the plots (2nd easiest but almost 1st)
# TODO: Reshape into markdown to show things more easily (a bit harder but also easy)
# TODO: Reshape 1_... to contain the functions and call them here (easiest)
