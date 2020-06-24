#' Create results heatmap
#' 
#' Draws a rectangular matrix with columns each of the diseases (or autoantibodies) and
#' rows for each microRNA. Each cell uses a colorscale for the log2foldchange and it's 
#' shadowed if its correspondent padj is greater than the threshold
#' 
#' @param results.table: Data.table containing at least one column of p-values for each of
#' the experiments, as well as one log2foldchange column for each, with one row for each
#' gene. These columns must include "pvalue" in and "log2foldchange" their names.
#'
#' @return ggplot2 object
create.results.heatmap <- function(results.table){
  require(data.table)
  require(ggplot2)
  
  # Clean-up quality cols to avoid errors
  quality.cols <- grep("^n_.+$", colnames(results.table), value=T)
  results.table <- data.table::copy(results.table)
  results.table[, (quality.cols) := NULL]
  
  # Get which columns contain each type of information
  info.cols <- c("gene.names", quality.cols)
  pvalue.cols <- grep("pvalue", colnames(results.table), value=T)
  padj.cols <- grep("padj", colnames(results.table), value=T)
  
  foldchange.cols <- setdiff(colnames(results.table), c(pvalue.cols, padj.cols, info.cols))
  
  # Each rectangle in the drawing will be a row in the drawing.data
  num.rows.drawing <- nrow(results.table)
  num.cols.drawing <- length(foldchange.cols)
  num.rectangles <- num.cols.drawing * num.rows.drawing

  # Create a data.table with the x and y positions of each rectangle and the color value
  rect.width = 1
  rect.height = 1
  
  xmin = rep(0:(num.cols.drawing-1), each=num.rows.drawing) * rect.width
  xmax = (xmin/rect.width + 1) * rect.width
  ymin = rep(0:(num.rows.drawing-1), times=num.cols.drawing) * rect.height
  ymax = (ymin/rect.height + 1) * rect.height

  
  # Stack vertically each foldchange column
  stacked.values <- unname(unlist(results.table[, ..foldchange.cols]))
  stacked.padjs <- unname(unlist(results.table[, ..padj.cols]))
  
  
  # Create table for drawing the rectangles
  drawing.data <- data.table(
    xmin = xmin,
    xmax = xmax,
    ymin = ymin,
    ymax = ymax,
    value=stacked.values
  )
  
  bad.data <- drawing.data[stacked.padjs >= 0.05]
  
  
  # Draw!
  plot <- ggplot() +
    # Axes
    scale_x_continuous("Disease", expand=c(0,0)) +
    scale_y_continuous("Gene",
                       trans="reverse",
                       breaks=seq_along(results.table[, unique(gene.names)]),
                       labels=results.table[, unique(gene.names)]) +
    # Plot values
    geom_rect(data = drawing.data,
              mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax, fill=value)) +
    # Set values with high padj as mustard-ish color
    geom_rect(data = bad.data,
              mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), alpha = 1,
              fill="#E69F00") +
    # Remove ugly elements from the graph
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank())
  
  # Add name of the diseases
  for(i in seq_along(foldchange.cols)){
    disease.name <- foldchange.cols[[i]]
    plot <- plot + annotate("text",
                            x=(i-0.5)*rect.width,
                            y=-rect.height,
                            hjust=0.5,
                            vjust=0,
                            label=disease.name)
  }

  return(plot)
}


#' Plot fold change of every gene in the result table, classified by group (IBM, DM,...)
#' @param results.table
#' 
#' @return ggplot2 plot object
plot.fold.change <- function(results.table){
  require(data.table)
  require(ggplot2)
  require(ggrepel)
  
  # TODO: Make the column names extraction an external function somewhere (like utils.R)
  # Clean-up quality cols to avoid errors
  quality.cols <- grep("^n_.+$", colnames(results.table), value=T)
  results.table <- data.table::copy(results.table)
  results.table[, (quality.cols) := NULL]
  
  # Get which columns contain each type of information
  info.cols <- c("gene.names", quality.cols)
  pvalue.cols <- grep("pvalue", colnames(results.table), value=T)
  padj.cols <- grep("padj", colnames(results.table), value=T)
  foldchange.cols <- setdiff(colnames(results.table), c(pvalue.cols, padj.cols, info.cols))
  
  # TODO: Make this casting process into an external funcion somewhere (like utils.R)
  by.group <- lapply(seq_along(foldchange.cols), function(i){
    data.table(gene = results.table[['gene.names']],
               disease = foldchange.cols[[i]],
               fold.change = results.table[[foldchange.cols[[i]]]],
               pvalue = results.table[[pvalue.cols[[i]]]])
  })
  casted.results <- rbindlist(by.group)
  
  # TODO: Make this ranking and re-adjusting into an external function, maybe the same function as the casting above
  top.genes <- casted.results[, rank := frankv(pvalue), by=disease]
  top.genes[, padj := p.adjust(pvalue, "BH"), by=disease]
  
  plot <- ggplot(data=top.genes[rank<=15], mapping=aes(x=disease, y=fold.change)) +
    geom_jitter(data=top.genes[rank>15], mapping=aes(color=rank)) +
    geom_jitter(position = position_jitter(seed=123, width=0.2), color="#132C44") +
    geom_text_repel( aes(label=gene), position = position_jitter(seed=123, width=0.2), min.segment.length = 0.1)
  
  return(plot)
}


#' Plot counts for every sample of every group
