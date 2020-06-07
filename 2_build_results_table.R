#' Create results table
#' 
#' It receives the list with the results for each group (each disease, each autoantibody...)
#' and it create a column per group and value needed. One for the mean log2 fold change,
#' one for the p-value and one for the adjusted p-value. This, for each disease, autoantibody...
#' 
#' @param results.list list containing the results from deseq2 as created on step 1
#' @param gene.names character vector containing the names of the genes, used for filtering
#' the desired genes in next steps on the analysis
#'
#' @return results.table
createResultsTable <- function(results.list){
  fold_change_table <- data.table()
  invisible(lapply(seq_along(results.list), function(i){
    this.name <- names(results.list)[i]
    this.group <- results.list[[i]]
    
    fold_change_table[, as.character(this.name) := this.group$log2FoldChange]
    fold_change_table[, paste0(this.name, "_pvalue") := this.group$pvalue]
    fold_change_table[, paste0(this.name, "_padj") := this.group$padj]
  }))
  
  gene.names <- results.list[[1]]@rownames
  fold_change_table[, gene.names := gene.names]
  setcolorder(fold_change_table, "gene.names")
  
  return(fold_change_table)
}


#' Get genes quality columns
#' 
#' Adds to the results table columns to account for some basic quality checks, like how
#' many times padj or pvalues are NA
#' 
#' @param results.table: results table computed using the above function. It assumes only
#' padj, pvalues and fold change columns are present.
#' 
#' @return data.table with the quality columns added
getGenesQualityColumns <- function(results.table){
  # First, we need to clean-up old quality columns (match everything starting by n_)
  quality.cols <- grep("^n_.+$", colnames(results.table), value=TRUE)
  results.table[, (quality.cols) := NULL]
  
  # If it matches nothing, it must be a fold change column:
  fold.change.cols <- grep("padj|pvalue|gene.names", colnames(results.table),
                           value=T,
                           invert=T)
  
  pvalue.cols <- grep("pvalue", colnames(results.table),
                      value=T)
  padj.cols <- grep("padj", colnames(results.table),
                    value=T)
  
  results.table[, n_nas_value := Reduce("+", lapply(.SD, is.na)), .SDcols=fold.change.cols]
  results.table[, n_nas_pvalue := Reduce("+", lapply(.SD, is.na)), .SDcols=pvalue.cols]
  results.table[, n_nas_padj := Reduce("+", lapply(.SD, is.na)), .SDcols=padj.cols]
  results.table[, n_relevant := Reduce("+", lapply(.SD, `<=`, 0.05)), .SDcols=padj.cols]
  
  invisible(results.table)
}


#' Filters table to keep only the desired genes
#' 
#' Additionally, it is necesary to re-compute p-values adjustment using BH (the method used)
#' by DeSEQ2. This is because these adjusted values depend on the whole set of p-values that
#' we have, and we now have a smaller subset.
#' 
#' @param results.table
#' @param gene.names.to.keep character vector containing the genes that we want to keep
filterResultsTable <- function(results.table, gene.names.to.keep){
  filtered <- results.table[gene.names %in% gene.names.to.keep]
  
  filtered[, padj := p.adjust(pvalue, "BH")]
  
  return(filtered)
}


#' Re-adjust p-values for a subset of RNAs
#' 
#' DeSEQ2 adjust p-values of the experiments using B-H based on the vector of p-values of
#' the whole list of genes. Since we are interested in only a fraction of this list, we
#' need to perform this same adjustment using only the vector of p-values associated with
#' the pertinent genes.
#' 
#' @param results.table: Data.table containing at least one column of p-values for each of
#' the experiments, with one row for each gene. These columns must include "pvalue" in
#' their names.
#' 
#' @return Silently returns the modified element. This function uses data.table's set*
#' syntax, and therefore the element is modified in-place.
#'
set.re.adjusted.pvalues <- function(results.table){
  pvalues.cols <- grep("pvalue", colnames(results.table), value=T)
  padj.cols <- grep("padj", colnames(results.table), value=T)
  
  for(i in seq_along(pvalues.cols)){
    this.pvalue.col <- pvalues.cols[[i]]
    this.padj.col <- padj.cols[[i]]
    
    set(results.table,
        j = this.padj.col,
        value = p.adjust(`[[`(results.table, this.pvalue.col), method="BH"))
  }
  
  invisible(results.table)
}


# # -- Get miRNA gene names
# mir_genes <- fold_change_table[startsWith(gene.names, "MIR"), unique(gene.names)]
# 
# # -- Check if there is any miRNA not starting by "MIR"
# possible_other_mir_genes <- fold_change_table[grep("MIR", gene.names), unique(gene.names)]
# setequal(mir_genes, possible_other_mir_genes)
# # --> Conclusion: All microRNAS start by "MIR"
# 
# # -- Filter-out non-miRNA
# only.mirna <- fold_change_table[gene.names %in% mir_genes]
# 
# 
# 
# # Filter if either all padj or all change values are NAs
# only.mirna <- only.mirna[n_nas_value != 8 & n_nas_padj != 8]
# 
# only.mirna[, `:=`(
#   n_nas_value = NULL,
#   n_nas_padj = NULL
# )]
# 
# # -- Filter-out rows with no significance for any disease
# only.mirna[, n_significant := rowSums(only.mirna[,lapply(.SD, `<`, 0.05), .SDcols=padj.cols], na.rm=TRUE)]


# 1. Filter by all microRNAs and plot using the function from utils.R
# 2. Same but using the list of microRNAs from literature only
# 3. Using literature + microRNAs found to have a lot of confidence to be good (i.e (1) + (2))
# (REMEMBER WE HAVE SAVED_RESULTS TO SPEED-UP)



