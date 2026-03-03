# Internal function to build region structure
Build_Region_Structure <- function(
    gtf = NULL,
    Chr,
    Start,
    End,
    Strand = "+",
    geneID = NULL,
    TxID = NULL,
    levels = NULL,
    species = "Human"
) {

  # Build GRanges Window from inputted Start Stop and End
  Range <- GenomicRanges::GRanges(
    seqnames = Chr,
    ranges = IRanges::IRanges(start = Start, end = End),
    strand = Strand
  )

  # Gets unclipped region which contains all longest transcripts for each gene in the desired region
  unclipped_region <- getRegion(
    gtf = gtf,
    Chr = Chr,
    Start = Start,
    End = End,
    geneID = geneID,
    TxID = TxID,
    species = species
  )

  # Clip any overlapping features exactly to the window
  clipped_region <- clip_overlapping_parts(unclipped_region, Range, Strand = Strand)

  # If empty return empty data frame
  if (!nrow(clipped_region)) {
    stop(sprintf("No rows in GTF found"))
  }
  return(clipped_region)
}


# ---------HELPER FUNCTIONS----------

# Clips the parts of exons, introns, UTRS not included in exact given region
clip_overlapping_parts <- function(gtf_rows_df, range_gr, Strand = NULL) {

  # If empty return it
  if (!nrow(gtf_rows_df)) {
    return(gtf_rows_df)
  }

  # Convert incoming rows to GRanges and keep all metadata columns
  gr <- GenomicRanges::makeGRangesFromDataFrame(gtf_rows_df, keep.extra.columns = TRUE)

  # Find overlaps of features with the region; if none, return empty of same shape
  overlaps <- GenomicRanges::findOverlaps(gr, range_gr, ignore.strand = FALSE)
  if (!length(overlaps)) {
    return(gtf_rows_df[0, , drop = FALSE])
  }

  # Intersect overlapping features with the region
  clip <- GenomicRanges::pintersect(
    gr[S4Vectors::queryHits(overlaps)],
    range_gr[S4Vectors::subjectHits(overlaps)],
    drop.nohit = TRUE
  )

  # Preserve original metadata
  S4Vectors::mcols(clip) <- S4Vectors::mcols(gr[S4Vectors::queryHits(overlaps)])

  # Return back data frame
  as.data.frame(clip)
}
