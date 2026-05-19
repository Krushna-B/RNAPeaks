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
