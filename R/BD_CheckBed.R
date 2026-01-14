#' Validate and Normalize BED Format Data
#'
#' Validates a data frame as BED format and normalizes column names and
#' chromosome identifiers. The function expects columns in standard BED order:
#' chr, start, end, tag/name, score, strand.
#'
#' @param df A data frame with at least 3 columns representing BED data.
#'   Columns are mapped by position: 1=chr, 2=start, 3=end, 4=tag, 5=score, 6=strand.
#'
#' @return A validated data frame with normalized column names (chr, start, end,
#'   tag, score, strand) and chromosome names without "chr" prefix.
#'
#' @details
#' The function performs the following validations:
#' \itemize{
#'   \item Requires at least 3 columns
#'   \item Checks that chromosome is character type
#'   \item Verifies end >= start for all rows
#'   \item Validates strand values are "+" or "-"
#'   \item Removes "chr" prefix from chromosome names
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   bed <- read.table("peaks.bed", header = FALSE)
#'   bed <- checkBed(bed)
#' }
checkBed <- function(df) {
  # Make data frame and make marker values lowercase
  df <- as.data.frame(df, stringsAsFactors = FALSE)
  colnames(df) <- tolower(colnames(df))

  # Require at least the 3 key columns to exist (by POSITION below)
  stopifnot(ncol(df) >= 3)

  # Map the FIRST up-to-6 columns BY POSITION to canonical BED-like names:
  # 1->chr, 2->start, 3->end, 4->tag, 5->score, 6->strand
  # NOTE: This assumes your input columns are already in this order

  col_names <- c("chr", "start", "end", "tag", "score", "strand")
  k <- min(ncol(df), length(col_names))

  for (i in seq_len(k)) {
    colnames(df)[i] <- col_names[i]
  }

  # Validate basic types/values:
  # - chr must be character
  # - end >= start for all rows
  # - strand must be character and only "+" or "-"
  if (is.character(df$chr) &
      all(as.numeric(df$end) >= as.numeric(df$start)) &
      is.character(df$strand) &
      all(df$strand %in% c("+", "-"))) {
    # Normalize chromosome names by removing a leading "chr"

    colchr_tmp <- tolower(df$chr)
    df$chr <- sub("chr", "", colchr_tmp)
    return(df)

  } else {
    stop("Check Bed File!")
  }
}
