
#' Validates a data frame as BED format and normalizes column names and
#' chromosome identifiers. The function expects columns in standard BED order:
#' chr, start, end, tag/name, score, strand.
#'
#' @param df A data frame with at least 3 columns representing BED data.
#'   Columns are mapped by position: 1=chr, 2=start, 3=end, 4=tag, 5=score, 6=strand.
#'   If fewer than 6 columns are provided, missing columns are filled with defaults:
#'   tag="peak", score=0, strand="+".
#'
#' @return A validated data frame with normalized column names (chr, start, end,
#'   tag, score, strand) and chromosome names without "chr" prefix.
#'
#' @details
#' The function performs the following:
#' \itemize{
#'   \item Requires at least 3 columns (chr, start, end)
#'   \item Maps columns by position to canonical BED names
#'   \item Fills missing columns with defaults (tag="peak", score=0, strand="+")
#'   \item Validates chromosome is character type
#'   \item Verifies end >= start for all rows
#'   \item Validates strand values are "+" or "-"
#'   \item Removes "chr" prefix from chromosome names
#'   \item Converts chromosome names to uppercase
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Standard 6-column BED
#'   bed <- read.table("peaks.bed", header = FALSE)
#'   bed <- checkBed(bed)
#'
#'   # 3-column BED (will add defaults for tag, score, strand)
#'   bed3 <- data.frame(chr = "chr1", start = 100, end = 200)
#'   bed3 <- checkBed(bed3)
#' }
checkBed <- function(df) {
  # Convert to data frame

  df <- as.data.frame(df, stringsAsFactors = FALSE)

  # Require at least 3 columns (chr, start, end)

if (ncol(df) < 3) {
    stop("BED file must have at least 3 columns (chr, start, end)")
  }

  # Map columns by position to canonical BED names
  col_names <- c("chr", "start", "end", "tag", "score", "strand")
  k <- min(ncol(df), length(col_names))
  colnames(df)[1:k] <- col_names[1:k]

  # Add missing columns with defaults
  if (!"tag" %in% colnames(df)) {
    df$tag <- "peak"
  }
  if (!"score" %in% colnames(df)) {
    df$score <- 0
  }
  if (!"strand" %in% colnames(df)) {
    df$strand <- "+"
    message("No strand column found, defaulting to '+'")
  }

  # Validate chromosome
  if (!is.character(df$chr)) {
    df$chr <- as.character(df$chr)
  }

  # Validate coordinates
  df$start <- as.numeric(df$start)
  df$end <- as.numeric(df$end)

  if (any(is.na(df$start)) || any(is.na(df$end))) {
    stop("BED file contains non-numeric start/end coordinates")
  }

  if (any(df$end < df$start)) {
    bad_rows <- which(df$end < df$start)
    stop(sprintf("BED file has end < start at %d rows (first: row %d)",
                 length(bad_rows), bad_rows[1]))
  }

  # Validate strand
  df$strand <- as.character(df$strand)
  invalid_strand <- !df$strand %in% c("+", "-")
  if (any(invalid_strand)) {
    bad_vals <- unique(df$strand[invalid_strand])
    stop(sprintf("Invalid strand values found: %s (must be '+' or '-')",
                 paste(bad_vals, collapse = ", ")))
  }

  # Normalize chromosome names: remove "chr" prefix, uppercase
df$chr <- toupper(df$chr)
  df$chr <- sub("^chr", "", df$chr, ignore.case = TRUE)

  message("BED file validated successfully")
  return(df)
}
