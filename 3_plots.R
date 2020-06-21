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


