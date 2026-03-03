
#'
#' Filters a BED data frame to retain only peaks within specified genomic
#' coordinates and optionally merges nearby peaks.
#'
#' @param bed A validated BED data frame (from `checkBed()`).
#' @param chr Chromosome to filter on.
#' @param start Start position of the region.
#' @param end End position of the region.
#' @param strand Strand to filter on ("+" or "-").
#' @param omit Character vector of tag/target names to exclude.
#' @param collapse Minimum gap width (bp) for merging nearby peaks.
#'   Set to 0 for no merging.
#'
#' @return A filtered data frame containing only peaks within the specified
#'   region, with nearby peaks optionally merged.
#'
#' @noRd
#' @examples
#' \dontrun{
#'
#'   filtered_bed <- FilterBed(
#'     bed = bed_data,
#'     chr = "17",
#'     start = 7565097,
#'     end = 7590856,
#'     strand = "-",
#'     omit = c("IgG"),
#'     collapse = 10
#'   )
#' }
FilterBed <- function(bed, chr, start, end, strand, omit = c(), collapse) {
  bed <- bed[which(as.character(bed$chr) == chr &
                   bed$start >= start &
                   bed$end <= end &
                   bed$strand == strand), ]

  # Removing Targets if given
  bed <- bed[which(!(bed$tag %in% omit)), ]

  # Converting into Granges list to merge
  bed_gr <- GenomicRanges::makeGRangesListFromDataFrame(
    bed,
    keep.extra.columns = TRUE,
    split.field = "tag"
  )
  bed_gr <- GenomicRanges::reduce(bed_gr, min.gapwidth = collapse)

  # Converting back to df
  red_bed <- data.frame(bed_gr)
  if (nrow(red_bed) == 0) {
    stop("No Peaks Found!")
  }
  return(red_bed)
}
