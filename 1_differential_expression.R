#' Create RNASeq analysis function
#' 
#' This function acts as a closure to create functions to perform differential expression
#' analysis from the count matrix. This way, the count matrix is read and processed only
#' once, and only once. It is intentionally hard to do it wrong. The path to the count matrix
#' is hardcoded in the function itself.
#' 
#' @param aggregation.level atomic character vector indicating the desired aggregation level.
#' At the moment, "autoantibody" and "disease" are supported
#' 
#' @return function to be used to get the DESeq2 analysis for each disease/autoantibody/...
create.rnaseq.analysis <- function(aggregation.level){
  require(DESeq2)
  require(data.table)
  
  # -- Assert that aggregation.level has a valid value
  if(!aggregation.level %in% c("autoantibody", "disease") | length(aggregation.level) != 1){
    stop("aggregation.level must be a character atomic vector with value either autoantibody or disease")
  }
  
  # -- Load unnormalized counts
  # If either count_matrix or samples_summary is not in the attributes, we create them
  #  and add the to the attributes. Else, we just read them to prevent unneeded reloading
  if( (!"count_matrix" %in% names(attributes(create.rnaseq.analysis))) |
      (!"samples_summary" %in% names(attributes(create.rnaseq.analysis)))
  ){
    count_matrix = read.csv("Datos/anonymized_gene_counts.csv",
                            row.names = 1)
    count_matrix = as.matrix(count_matrix)
    
    # -- Remove low count rows --
    COUNT.THRESHOLD <- 10
    low_counts <- rowSums2(count_matrix) <= COUNT.THRESHOLD
    
    count_matrix <- count_matrix[!low_counts, ]
    
    # -- Create colData with sample description to satisfy DESeqDataSet needs
    samples_summary <- data.table(
      sample = colnames(count_matrix)
    )
    
    # TODO: Move this samples summary creation outside of this closure and remove the closure all together
    
    # Experiment names for autoantibody aggregation level
    samples_summary[startsWith(sample, 'NT'), autoantibody := "Normal_muscle"]
    samples_summary[startsWith(sample, 'IBM'), autoantibody := "Inclusion_body"]
    samples_summary[startsWith(sample, 'HMGCR'), autoantibody := "Necrotizing_myopathy_Anti.HMGCR"]
    samples_summary[startsWith(sample, 'SRP'), autoantibody := "Necrotizing_myopathy_Anti.SRP"]
    samples_summary[startsWith(sample, 'Mi2'), autoantibody := "Dermatomyositis_Anti.Mi2"]
    samples_summary[startsWith(sample, 'NXP2'), autoantibody := "Dermatomyositis_Anti.NXP2"]
    samples_summary[startsWith(sample, 'MDA5'), autoantibody := "Dermatomyositis_Anti.MDA5"]
    samples_summary[startsWith(sample, 'TIF1'), autoantibody := "Dermatomyositis_Anti.TIF1g"]
    samples_summary[startsWith(sample, 'Jo1'), autoantibody := "Antisynthetase_Syndrom_Anti.Jo1"]
    
    # Experiment names for disease aggregation level
    samples_summary[autoantibody == "Normal_muscle", disease := "Normal_muscle"]
    samples_summary[startsWith(autoantibody, "Inclusion_body"), disease := "Inclusion_body"]
    samples_summary[startsWith(autoantibody, "Necrotizing"), disease := "Necrotizing"]
    samples_summary[startsWith(autoantibody, "Dermatomyositis"), disease := "Dermatomyositis"]
    samples_summary[startsWith(autoantibody, "Antisynthetase"), disease := "Antisynthetase"]
    
    # TODO: Temporary fix for plots
    samples_summary[startsWith(sample, 'NT'), autoantibody := "Normal_muscle"]
    samples_summary[startsWith(sample, 'IBM'), autoantibody := "Inclusion_body"]
    samples_summary[startsWith(sample, 'HMGCR'), autoantibody := "Anti.HMGCR"]
    samples_summary[startsWith(sample, 'SRP'), autoantibody := "Anti.SRP"]
    samples_summary[startsWith(sample, 'Mi2'), autoantibody := "Anti.Mi2"]
    samples_summary[startsWith(sample, 'NXP2'), autoantibody := "Anti.NXP2"]
    samples_summary[startsWith(sample, 'MDA5'), autoantibody := "Anti.MDA5"]
    samples_summary[startsWith(sample, 'TIF1'), autoantibody := "Anti.TIF1g"]
    samples_summary[startsWith(sample, 'Jo1'), autoantibody := "Syndrom_Anti.Jo1"]
    
    # Which samples are control group and which ones are not?
    samples_summary[, group := fifelse(startsWith(sample, "NT"), "control", "treated")]
    
    # Add them to the attributes for re-useability
    attr(create.rnaseq.analysis, "count_matrix", count_matrix)
    attr(create.rnaseq.analysis, "samples_summary", samples_summary)
    
  } else {
    count_matrix <- attr(create.rnaseq.analysis, "count_matrix", exact=TRUE)
    samples_summary <- attr(create.rnaseq.analysis, "samples_summary", exact=TRUE)
  }
  
  
  # Customize actions based on selected aggregation level
  design.formula <- formula(paste('~', aggregation.level))
  sample_groups <- unname(unlist((unique(samples_summary[, ..aggregation.level]))))
  sample_groups <- setdiff(sample_groups, "Normal_muscle")
  
  # Create the function to perform the anaysis
  analyze <- function(){
    dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                                  colData = samples_summary,
                                  design = design.formula)
    
    analysis = DESeq(dds, parallel=TRUE)
    
    # Extract results information from the analysis object
    results <- lapply(sample_groups, function(group){
      return(results(analysis, contrast=c(aggregation.level, "Normal_muscle", group)))
    })
    
    # Give meaningful names to the results
    names(results) <- sample_groups
    
    return(results)
  }
  
  return(analyze)
}

