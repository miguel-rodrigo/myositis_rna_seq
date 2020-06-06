library(data.table)
library(DESeq2)

# -- Load unnormalized counts
count_matrix = read.csv("Datos/anonymized_gene_counts.csv",
                        row.names = 1)
count_matrix = as.matrix(count_matrix)

# -- Remove low count rows --
COUNT.THRESHOLD <- 10
low_counts <- rowSums2(count_matrix) <= COUNT.THRESHOLD

count_matrix <- count_matrix[!low_counts, ]


# -- Create colData with sample description
samples_summary <- data.table(
  sample = colnames(count_matrix)
)

samples_summary[startsWith(sample, 'NT'), condition := "Normal_muscle"]
samples_summary[startsWith(sample, 'IBM'), condition := "Inclusion_body"]
samples_summary[startsWith(sample, 'HMGCR'), condition := "Necrotizing_myopathy_Anti-HMGCR"]
samples_summary[startsWith(sample, 'SRP'), condition := "Necrotizing_myopathy_Anti-SRP"]
samples_summary[startsWith(sample, 'Mi2'), condition := "Dermatomyositis_Anti-Mi2"]
samples_summary[startsWith(sample, 'NXP2'), condition := "Dermatomyositis_Anti-NXP2"]
samples_summary[startsWith(sample, 'MDA5'), condition := "Dermatomyositis_Anti-MDA5"]
samples_summary[startsWith(sample, 'TIF1'), condition := "Dermatomyositis_Anti-TIF1g"]
samples_summary[startsWith(sample, 'Jo1'), condition := "Antisynthetase_Syndrom_Anti-Jo1"]

samples_summary[, group := fifelse(startsWith(sample, "NT"), "control", "treated")]


# ---- RESULTS BY EACH DIFFERENT TYPE OF AUTOANTIBODY ----
# -- Create DESeqData object
dds = DESeqDataSetFromMatrix(countData = count_matrix,
                             colData = samples_summary,
                             design = ~ condition)

# -- Perform differential expression analysis
analysis = DESeq(dds, parallel=TRUE)

sample_groups = samples_summary[, unique(condition)]
sample_groups = setdiff(sample_groups, "Normal_muscle")

results.by.autoantibody <- lapply(sample_groups, function(group){
  return(results(analysis, contrast=c("condition", "Normal_muscle", group)))
})
names(results.by.autoantibody) <- sample_groups


# ---- RESULTS BY EACH DIFFERENT TYPE OF DISEASE
samples_summary[condition == "Normal_muscle", grouped.condition := "Normal_muscle"]
samples_summary[startsWith(condition, "Inclusion_body"), grouped.condition := "Inclusion_body"]
samples_summary[startsWith(condition, "Necrotizing"), grouped.condition := "Necrotizing"]
samples_summary[startsWith(condition, "Dermatomyositis"), grouped.condition := "Dermatomyositis"]
samples_summary[startsWith(condition, "Antisynthetase"), grouped.condition := "Antisynthetase"]


dds = DESeqDataSetFromMatrix(countData = count_matrix,
                             colData = samples_summary,
                             design = ~ grouped.condition)
analysis = DESeq(dds, parallel=TRUE)

sample_groups <- samples_summary[, unique(grouped.condition)]
sample_groups <- setdiff(sample_groups, "Normal_muscle")

results.by.disease <- lapply(sample_groups, function(group){
  return(results(analysis, contrast=c("grouped.condition", "Normal_muscle", group)))
})
names(results.by.disease) <- sample_groups
