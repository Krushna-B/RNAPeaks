#' Prepare a BED for plotting: filter peaks, choose display order, lay out tracks.

#' @param bed Validated BED data frame
#' @param filter Named list: `chr`, `start`, `end`, `strand`, `omit`, `collapse`.
#' @param order  Named list: `by` (`"Count"` or `"Alphabetical"`), `in_`,
#'   `max_proteins`.
#' @param track_height Vertical height (y units) of each protein track row.
#'
#' @return Data frame with columns for plotting peaks: `rank`, `y_start`, `y_end`.
#' @noRd
#' @family bed
prepare_bed <- function(bed,
                        filter = list(chr = NULL, start = NULL, end = NULL,
                                      strand = NULL, omit = NULL, collapse = 0),
                        order  = list(by = "Count", in_ = NULL, max_proteins = 40),
                        track_height = 0.3) {

  #Validate Params
  if (!is.data.frame(bed) || nrow(bed) == 0L) {
    abort_invalid_bed("{.arg bed} must be a non-empty data frame.")
  }
  if (!is.list(filter)) abort_invalid_arg("{.arg filter} must be a named list.")
  if (!is.list(order))  abort_invalid_arg("{.arg order} must be a named list.")

  if (!is.numeric(order$max_proteins) || length(order$max_proteins) != 1L ||
      is.na(order$max_proteins) || order$max_proteins < 1) {
    abort_invalid_arg("{.arg order$max_proteins} must be a positive integer.")
  }
  if (!is.numeric(track_height) || length(track_height) != 1L ||
      track_height <= 0) {
    abort_invalid_arg("{.arg track_height} must be a positive number.")
  }

  #Filter Bed
  bed <- tryCatch(
    filter_bed(bed, filter$chr, filter$start, filter$end,
               filter$strand, filter$omit, filter$collapse),
    error = function(e) cli::cli_abort(
      "Filtering stage failed.",
      parent = e,
      class  = c("rnapeaks_prepare_filter_bed", "rnapeaks_error"),
      call   = NULL
    )
  )

  # No peaks left after filtering return out
  if (is.null(bed)) return(NULL)

  #Create the order of proteins
  if (!is.null(order$in_)) {
    if (!is.character(order$in_) || anyNA(order$in_)) {
      abort_invalid_arg("{.arg order$in_} must be a character vector with no NAs.")
    }
    missing_targets <- setdiff(order$in_, unique(bed$group_name))
    if (length(missing_targets)) {
      cli::cli_warn(c(
        "{length(missing_targets)} target{?s} in {.arg order$in_} not present in BED:",
        "x" = "{.field {missing_targets}}"
      ))
    }
    rank_ <- order$in_
  } else {
    rank_ <- switch(order$by,
      Count        = names(sort(table(bed$group_name), decreasing = TRUE)),
      Alphabetical = sort(unique(bed$group_name)),
      abort_invalid_arg(
        "{.arg order$by} must be {.val Count} or {.val Alphabetical}; got {.val {order$by}}."
      )
    )
  }

  #Remove proteins above the max protein number after the ordering
  if (length(rank_) > order$max_proteins) {
    cli::cli_warn(c(
      "Showing top {order$max_proteins} of {length(rank_)} proteins (by {order$by}).",
      "i" = "Raise {.arg order$max_proteins} to show more."
    ))
    rank_ <- rank_[seq_len(order$max_proteins)]
  }
  bed <- bed[bed$group_name %in% rank_, , drop = FALSE]
  if (nrow(bed) == 0L) {
    abort_not_found("No peaks remain after restricting to the ordered target set.")
  }

  # Building peaks df with order and y axis values.
  # rev(rank_) makes the first item in rank_ appear at the TOP of the plot.
  bed$rank    <- match(bed$group_name, rev(rank_))
  bed$y_start <- track_height * bed$rank - track_height
  bed$y_end   <- track_height * bed$rank
  return(bed)
}

#' Restrict a validated BED to a chromosome / strand / coordinate window,
#' drop omitted targets, and merge nearby peaks per target.
#'
#' @param bed A BED data frame returned by [checkBed()].
#' @param chr Chromosome to retain; `"chr"` prefix and case are normalized.
#' @param start,end Region bounds (numeric, `start <= end`, both `>= 0`).
#' @param strand Strand to retain: `"+"` or `"-"`.
#' @param omit Optional character vector of `target` values to drop.
#' @param collapse Min gap width (bp) for merging adjacent peaks.
#'
#' @return Data frame of reduced peaks with a `group_name` column
#'   (= source `target`).
#'
#' @noRd
#' @family bed
filter_bed <- function(bed, chr, start, end, strand, omit = NULL, collapse) {

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

  # Normalize chr (strip "chr" prefix, uppercase)
  chr <- normalize_chr(chr)

  keep <- as.character(bed$chr) == chr &
    bed$start  >= start &
    bed$end    <= end   &
    bed$strand == strand
  bed <- bed[which(keep), , drop = FALSE]

  if (length(omit) > 0L) {
    bed <- bed[!bed$target %in% omit, , drop = FALSE]
  }

  if (nrow(bed) == 0L) {
    cli::cli_alert_info(c(
      "No peaks found in {.field {chr}:{start}-{end}} ({strand}).",
      "i" = if (length(omit) > 0L)
        "After omitting target{?s}: {.field {omit}}." else NULL
    ))
    return(NULL)
  }

  bed_gr <- tryCatch({
    gr <- GenomicRanges::makeGRangesListFromDataFrame(
      bed,
      keep.extra.columns = TRUE,
      split.field        = "target"
    )
    GenomicRanges::reduce(gr, min.gapwidth = collapse)
  }, error = function(e) {
    cli::cli_abort(
      "Failed to merge peaks via {.pkg GenomicRanges}.",
      parent = e,
      class  = c("rnapeaks_prepare_merge_peaks", "rnapeaks_error"),
      call   = NULL
    )
  })

  cli::cli_alert_success("BED filtered successfully")
  return (as.data.frame(bed_gr))
}
