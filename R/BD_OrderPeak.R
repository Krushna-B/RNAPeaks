#' Order Protein/Target Peaks for Display
#'
#' Determines the display order of protein tracks in the peak visualization
#' based on various criteria.
#'
#' @param bed A BED data frame with a `group_name` column.
#' @param order_by Ordering method:
#'   \describe{
#'     \item{"Count"}{Order by number of peaks (most peaks first)}
#'     \item{"Target"}{Order alphabetically by target name}
#'     \item{"Region"}{Order by genomic clustering similarity}
#'   }
#'
#' @return A character vector of target names in the desired display order.
#'
#'
#'
#' @examples
#' \dontrun{
#'   # Order by peak count
#'   order <- OrderPeak(bed_data, order_by = "Count")
#'
#'   # Order alphabetically
#'   order <- OrderPeak(bed_data, order_by = "Target")
#' }
OrderPeak <- function(bed, order_by) {
  if (order_by == "Region") {
    rank <- By_region(bed)
  } else if (order_by == "Target") {
    rank <- By_name(bed)
  } else if (order_by == "Count") {
    rank <- By_count(bed)
  } else {
    stop("Please provide acceptable feature to rank by: Target, Count")
  }
  return(rank)
}


# Order by name
By_name <- function(bed) {
  return(sort(unique(bed$group_name)))
}


# Order by number of peaks
By_count <- function(bed) {
  df <- names(sort(table(bed$group_name), decreasing = TRUE))
  return(df)
}


# Cluster by most similar
By_region <- function(bed) {
  bins <- IRanges::IRanges(
    start = min(bed$start):max(bed$end),
    width = 1,
    names = min(bed$start):max(bed$end)
  )
  regions_by_target <- list()
  for (i in 1:nrow(bed)) {
    ran <- IRanges::IRanges(start = bed$start[i]:bed$end[i], width = 1)
    if (!(bed$group_name[i] %in% names(regions_by_target))) {
      regions_by_target[[bed$group_name[i]]] <- ran
    } else if (bed$group_name[i] %in% names(regions_by_target)) {
      regions_by_target[[bed$group_name[i]]] <- c(
        regions_by_target[[bed$group_name[i]]],
        ran
      )
    }
    bins_df <- data.frame("Region" = names(bins))
    for (i in 1:length(names(regions_by_target))) {
      m <- rep(0, length(bins))
      m[S4Vectors::queryHits(IRanges::findOverlaps(
        bins,
        regions_by_target[[names(regions_by_target)[i]]]
      ))] <- 1
      bins_df <- cbind(bins_df, m)
    }
    rownames(bins_df) <- bins_df$Region
    bins_df <- bins_df[, -1]
    colnames(bins_df) <- names(regions_by_target)
  }
}
