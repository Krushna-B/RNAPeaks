#' Restrict a validated BED to a chromosome / strand / coordinate window,
#' drop omitted targets, and merge nearby peaks per target.
#'
#' @param bed A BED data frame returned by [checkBed()].
#' @param chr Chromosome to retain; `"chr"` prefix and case are normalized.
#' @param start,end Region bounds (numeric, `start <= end`, both `>= 0`).
#' @param strand Strand to retain: `"+"` or `"-"`.
#' @param omit Optional character vector of `target` values to drop.
#' @param collapse Min gap width (bp) for merging adjacent peaks
#'
#' @return Data frame of reduced peaks with a `group_name` column
#'   (= source `target`).
#'
#' @noRd
#' @family bed
filter_bed <- function(bed, chr, start, end, strand, omit = NULL, collapse) {

    #Validate parameters
    if (length(chr) != 1L || is.na(chr)) {
      abort_invalid_arg("{.arg chr} must be a single non-NA value.")
    }

    if (!is.numeric(start) || length(start) != 1L || is.na(start) || start < 0) {
      abort_invalid_arg("{.arg start} must be a single non-negative number.")
    }

    if (!is.numeric(end) || length(end) != 1L || is.na(end) || end < start) {
     abort_invalid_arg(c(
      "{.arg end} must be a single number {.code >= start}.",
      "x" = "Got {.val {start}} and {.val {end}}."
      ))
    }

    if (!identical(strand, "+") && !identical(strand, "-")) {
         abort_invalid_arg(c(
          "{.arg strand} must be {.val +} or {.val -}.",
          "x" = "You supplied {.val {strand}}."
        ))
    }

    if (!is.null(omit) && (!is.character(omit) || anyNA(omit))) {
        abort_invalid_arg("{.arg omit} must be a character vector or {.code NULL}.")
    }

    if (!is.numeric(collapse) || length(collapse) != 1L ||
        is.na(collapse) || collapse < 0) {
        abort_invalid_arg("{.arg collapse} must be a single non-negative number.")
    }

    #Normalize chr to ensure match with BED df
    chr <- sub("^CHR", "", toupper(as.character(chr)))

    #Keep rows that overlap with specified region
    keep <- as.character(bed$chr) == chr &
      bed$start  >= start &
      bed$end    <= end   &
      bed$strand == strand
    bed <- bed[which(keep), , drop = FALSE]

    #Remove rows that should be omitted
    if (length(omit) > 0L) {
        bed <- bed[!bed$target %in% omit, , drop = FALSE]
    }

    if (nrow(bed) == 0L) {
      abort_not_found(c(
        "No peaks found in {.field {chr}:{start}-{end}} ({strand}).",
        "i" = if (length(omit) > 0L)
        "After omitting target{?s}: {.field {omit}}." else NULL
        ))
    }

    #Merge peaks per target and reduce
    bed_gr <- tryCatch({
      gr <- GenomicRanges::makeGRangesListFromDataFrame(
        bed,
        keep.extra.columns = TRUE,
        split.field        = "target"
      )
      GenomicRanges::reduce(gr, min.gapwidth = collapse)
    }, error = function(e) {
      cli::cli_abort(c(
        "Failed to merge peaks via {.pkg GenomicRanges}.",
        "x" = "{conditionMessage(e)}"
      ))
    })

  cli::cli_alert_success("BED filtered successfully")
  return(as.data.frame(bed_gr))
}
