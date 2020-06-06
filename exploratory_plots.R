# ---- Histogram of adjusted p-values ----
hist(results[[1]]$pvalue,
     col = "grey", border = "white", xlab = "", ylab="",
     main = "frequencies of p-values")
hist(results[[1]]$padj,
     col = "grey", border = "white", xlab = "", ylab="",
     main = "frequencies of adjusted p-values")

# ---- MA plots (log fold change vs. mean expression) ----
plotMA(results[[1]], alpha = 0.05, main = "asdf", ylim = c(-7, 7))

# ---- Heatmaps ----
library(NMF)

results1.sorted <- results[[1]][order(results[[1]]$padj), ]
genes <- rownames(subset(results1.sorted, padj < 0.05))

# TODO: adapt these declarations to my code:
counts.sf_normalized <- counts(dds , normalized = FALSE)
log.norm.counts <- log2(counts.sf_normalized + 1)
hm.mat_genes <- log.norm.counts
aheatmap ( hm.mat_genes ,
           Rowv = TRUE , Colv = TRUE , # add dendrograms to rows and columns
           distfun = "euclidean" , hclustfun = "average" )
