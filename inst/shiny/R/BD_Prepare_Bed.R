Prepare_Bed <- function(bed, rank_, gene, peaks_width = 0.3, peak_col = "blue", margin = 100) {
  # Assign each group_name to a vertical row index.
  # rev(rank_) makes the first item in rank_ appear at the TOP of the plot.
  bed$rank    <- match(bed$group_name, rev(rank_))

  # Define the vertical band (row) for each group
  bed$y_start <- peaks_width*bed$rank - peaks_width
  bed$y_end   <- peaks_width*bed$rank

  # Single fill color for all peaks on this bed row
  bed$col     <- peak_col

  # X position for the left-side labels: margin bp to the left of the gene start.
  bed$xpos    <- min(gene$start) - margin


  # Return bed with added columns: rank, y_start, y_end, col, xpos
  bed
}
