#' Internal: shared helpers for splicing / sequence map entry points
#'
#' @keywords internal
#' @name sm_entry_helpers


#Read, validate, and turn a BED into a reduced GRanges of peaks.
.peaks_to_granges <- function(bed_file) {
  bed <- if (is.character(bed_file)) utils::read.table(bed_file) else bed_file
  bed <- check_bed(bed)
  gr  <- GenomicRanges::makeGRangesFromDataFrame(
    bed,
    seqnames.field     = "chr",
    start.field        = "start",
    end.field          = "end",
    strand.field       = "strand",
    keep.extra.columns = TRUE
  )
  GenomicRanges::reduce(gr)
}


#Resolve `genome` to a BSgenome object. Accepts NULL (hg38 default), one of
#"hg38" / "mm10" / "mm39"..
.resolve_genome <- function(genome) {
  if (is.null(genome))                 genome <- "hg38"
  if (methods::is(genome, "BSgenome")) return(genome)

  pkg <- switch(genome,
    hg38 = "BSgenome.Hsapiens.UCSC.hg38",
    mm10 = "BSgenome.Mmusculus.UCSC.mm10",
    mm39 = "BSgenome.Mmusculus.UCSC.mm39",
    NULL
  )
  if (is.null(pkg)) {
    abort_invalid_arg(c(
      "{.arg genome} must be {.val NULL}, one of {.or {.val hg38} {.val mm10} {.val mm39}}, or a BSgenome object.",
      "x" = "Got {.cls {class(genome)[1]}}."
    ))
  }
  if (!requireNamespace(pkg, quietly = TRUE)) {
    abort_not_found(c(
      "Required BSgenome package {.pkg {pkg}} is not installed.",
      "i" = "Install with: {.code BiocManager::install({.val {pkg}})}."
    ))
  }
  BSgenome::getBSgenome(pkg)
}


#Entry-level shape checks for splicing/sequence map functions.
#
#Only validates arguments whose downstream check would surface with a
#different name (`bed_file` -> "bed" in check_bed) or that aren't validated
#downstream at all (`title`). Everything else (events, sequence, genome,
#opts, style) already aborts with the right arg name from its own validator,
#so re-checking here would just be duplication.
.assert_sm_entry_args <- function(title, bed_file = NULL, has_bed = FALSE) {
  check_string(title, "title")

  if (has_bed) {
    ok_bed <- is.character(bed_file) || is.data.frame(bed_file) ||
              (is.list(bed_file) && length(bed_file) > 0L &&
               all(vapply(bed_file, is.data.frame, logical(1))))
    if (!ok_bed) {
      abort_invalid_arg(c(
        "{.arg bed_file} must be a file path, BED data frame, or list of BED data frames.",
        "x" = "Got {.cls {class(bed_file)[1]}}.",
        "i" = "Did the argument order get mixed up?"
      ))
    }
  }
}


#Error wrapper for the splicing/sequence map entry points.
.run_sm_entry <- function(map_name, body) {
  tryCatch(
    body,
    error = function(cnd) {
      msg <- paste0("Failed to generate ", map_name, ".")
      if (inherits(cnd, "rnapeaks_error")) {
        cli::cli_abort(msg, parent = cnd)
      } else {
        cli::cli_abort(
          c(msg, "x" = "An unexpected error occurred."),
          parent = cnd
        )
      }
    }
  )
}


#Normalize motifs: trim, uppercase, U to T
.normalize_motifs <- function(sequence) {
  if (!is.character(sequence) || length(sequence) == 0L ||
      anyNA(sequence) || any(nchar(trimws(sequence)) == 0L)) {
    abort_invalid_arg(c(
      "{.arg sequence} must be a non-empty character vector of motifs.",
      "x" = "Got {.cls {class(sequence)[1]}}."
    ))
  }


  has_u    <- grepl("U", sequence, fixed = TRUE)
  if (any(has_u)) {
    cli::cli_inform("Converting U to T in motif{?s}: {.val {sequence[has_u]}}.")
    sequence <- gsub("U", "T", sequence, fixed = TRUE)
  }
  sequence
}
