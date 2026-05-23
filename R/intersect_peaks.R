#' Intersect peaks with transcripts
#'
#' R replacement for `bedtools intersect -a peaks -b transcripts
#' -f <fraction> -wa -wb [-s]`.
#' For each peak, finds every transcript that
#' overlaps it under the given fraction-and-strand constraints.
#'
#' @param peaks Peak input. One of:
#'     file path to a BED file (>= 6 columns),
#'     a data.frame in BED layout (cols 1=chr, 2=start, 3=end,
#'     4=name, 5=score, 6=strand).
#'
#' @param transcripts Transcript input. One of:
#'     a GTF-shaped data.frame
#'     a BED data.frame (6 cols);
#'     a file path to a BED file;
#'     a `GRanges` with a `transcript_name` or `name` `mcols` column.
#'
#' @param output Optional file path. If supplied, writes a tab-separated BED
#'   and returns the path invisibly.
#'
#' @param fraction Minimum fraction of each peak that must overlap a
#'   transcript for the hit to be kept. Defaults to `1.0` (peak fully
#'   contained).
#'
#' @param same_strand If `TRUE` (default), only same-strand overlaps are
#'   kept.
#'
#' @return Either a data.frame of peak+transcript rows, or the output file
#'   path invisibly when `output` is supplied.
#'
#' @export
#' @family analysis
intersect_peaks <- function(peaks, transcripts, output = NULL,
                            fraction = 1.0, same_strand = TRUE) {
  #Validate Params
  check_scalar_number(fraction, "fraction", min = 0, max = 1)
  check_flag(same_strand, "same_strand")
  if (!is.null(output)) check_string(output, "output")

  #Normalize Values
  cli::cli_progress_step("Reading peaks and transcripts")
  peaks_df <- .normalize_peaks(peaks)
  tx_df    <- .normalize_transcripts(transcripts)

  #Create GRanges
  cli::cli_progress_step(
    "Building GRanges ({nrow(peaks_df)} peaks, {nrow(tx_df)} transcripts)"
  )
  peaks_gr <- .bed_to_gr(peaks_df, "peaks")
  tx_gr    <- .bed_to_gr(tx_df,    "transcripts")

  # Restrict both objects to chromosomes present in both.
  common <- intersect(GenomeInfoDb::seqlevels(peaks_gr),
                      GenomeInfoDb::seqlevels(tx_gr))
  if (length(common) == 0L) {
    abort_not_found(c(
      "No shared chromosomes between peaks and transcripts.",
      "i" = "Check naming convention (UCSC {.val chr1} vs Ensembl {.val 1})."
    ))
  }
  peaks_keep <- as.character(GenomicRanges::seqnames(peaks_gr)) %in% common
  tx_keep    <- as.character(GenomicRanges::seqnames(tx_gr))    %in% common
  dropped_peaks <- sum(!peaks_keep)
  dropped_tx    <- sum(!tx_keep)
  if (dropped_peaks > 0L || dropped_tx > 0L) {
    cli::cli_alert_warning(c(
      "Dropping ranges on non-shared chromosomes: ",
      "{dropped_peaks} peak{?s}, {dropped_tx} transcript{?s}."
    ))
  }
  peaks_gr <- GenomeInfoDb::keepSeqlevels(peaks_gr, common,
                                          pruning.mode = "coarse")
  tx_gr    <- GenomeInfoDb::keepSeqlevels(tx_gr, common,
                                          pruning.mode = "coarse")

  #Find Hits
  cli::cli_progress_step(
    "Finding overlaps ({.field fraction = {fraction}, same_strand = {same_strand}})"
  )
  hits <- if (fraction >= 1.0) {
    GenomicRanges::findOverlaps(
      peaks_gr, tx_gr, type = "within",
      ignore.strand = !isTRUE(same_strand)
    )
  } else {
    raw <- GenomicRanges::findOverlaps(
      peaks_gr, tx_gr, type = "any",
      ignore.strand = !isTRUE(same_strand)
    )
    pq <- peaks_gr[S4Vectors::queryHits(raw)]
    ts <- tx_gr[S4Vectors::subjectHits(raw)]
    ov <- IRanges::pintersect(pq, ts, ignore.strand = !isTRUE(same_strand))
    frac <- BiocGenerics::width(ov) / BiocGenerics::width(pq)
    raw[frac >= fraction]
  }

  #Select Hits
  q <- S4Vectors::queryHits(hits)
  s <- S4Vectors::subjectHits(hits)
  if (length(q) == 0L) {
    cli::cli_progress_done()
    abort_not_found(c(
      "No peak/transcript overlaps found.",
      "i" = "Tried fraction = {fraction}, same_strand = {same_strand}."
    ))
  }

  #Create output
  cli::cli_progress_step(
    "Assembling output ({length(q)} peak-transcript pair{?s})"
  )
  # Column-wise list assembly
  combined <- c(lapply(peaks_df, `[`, q),
                lapply(tx_df,    `[`, s))
  names(combined) <- paste0("V", seq_along(combined))
  out <- as.data.frame(combined, stringsAsFactors = FALSE)
  rownames(out) <- NULL

  if (!is.null(output)) {
    cli::cli_progress_step("Writing {.path {output}}")
    utils::write.table(out, file = output, sep = "\t",
                       quote = FALSE, row.names = FALSE, col.names = FALSE)
    cli::cli_progress_done()
    return(invisible(output))
  }
  cli::cli_progress_done()
  out
}



# Internal normalizers
# Always returns a character-typed data.frame in BED layout (0-based start).
.normalize_peaks <- function(x) {
  if (is.character(x) || is.data.frame(x)) {
    df <- .read_bed_df(x, "peaks")
    if (ncol(df) < 6L) {
      abort_invalid_bed(c(
        "{.arg peaks} must have at least 6 columns (BED6).",
        "x" = "Got {ncol(df)}."
      ))
    }
    return(df)
  }
  if (inherits(x, "GRanges")) {
    nm <- if (!is.null(S4Vectors::mcols(x)$name)) {
      as.character(S4Vectors::mcols(x)$name)
    } else {
      rep(".", length(x))
    }
    sc <- if (!is.null(S4Vectors::mcols(x)$score)) {
      as.character(S4Vectors::mcols(x)$score)
    } else {
      rep("0", length(x))
    }
    df <- data.frame(
      V1 = as.character(GenomicRanges::seqnames(x)),
      V2 = as.character(as.integer(BiocGenerics::start(x)) - 1L),
      V3 = as.character(as.integer(BiocGenerics::end(x))),
      V4 = nm,
      V5 = sc,
      V6 = as.character(GenomicRanges::strand(x)),
      stringsAsFactors = FALSE
    )
    return(df)
  }
  abort_invalid_arg(c(
    "{.arg peaks} must be a file path, data frame, or GRanges.",
    "x" = "Got {.cls {class(x)[1]}}."
  ))
}

# Always returns a 6-column character-typed BED data.frame (0-based start).
.normalize_transcripts <- function(x) {
  if (is.character(x)) {
    df <- .read_bed_df(x, "transcripts")
    if (ncol(df) < 6L) {
      abort_invalid_bed(c(
        "{.arg transcripts} BED must have at least 6 columns.",
        "x" = "Got {ncol(df)}."
      ))
    }
    return(df[, 1:6])
  }
  if (is.data.frame(x)) {
    # GTF data.frame? Detected by presence of a `type` column.
    if ("type" %in% colnames(x)) {
      tx <- x[as.character(x$type) == "transcript", , drop = FALSE]
      if (nrow(tx) == 0L) {
        abort_not_found("{.arg transcripts}: no rows with {.field type} == {.val transcript}.")
      }
      name_col <- if ("transcript_name" %in% colnames(tx)) {
        "transcript_name"
      } else if ("transcript_id" %in% colnames(tx)) {
        "transcript_id"
      } else {
        abort_invalid_arg("{.arg transcripts} GTF needs {.field transcript_name} or {.field transcript_id}.")
      }
      df <- data.frame(
        V1 = as.character(tx$seqnames),
        V2 = as.character(as.integer(tx$start) - 1L),  # GTF 1-based -> BED 0-based
        V3 = as.character(as.integer(tx$end)),
        V4 = as.character(tx[[name_col]]),
        V5 = ".",
        V6 = as.character(tx$strand),
        stringsAsFactors = FALSE
      )
      return(df)
    }
    df <- .read_bed_df(x, "transcripts")
    if (ncol(df) < 6L) {
      abort_invalid_bed(c(
        "{.arg transcripts} BED must have at least 6 columns.",
        "x" = "Got {ncol(df)}."
      ))
    }
    return(df[, 1:6])
  }
  if (inherits(x, "GRanges")) {
    nm <- S4Vectors::mcols(x)$transcript_name %||% S4Vectors::mcols(x)$name
    if (is.null(nm)) {
      abort_invalid_arg("{.arg transcripts} GRanges needs a {.field transcript_name} or {.field name} mcols column.")
    }
    df <- data.frame(
      V1 = as.character(GenomicRanges::seqnames(x)),
      V2 = as.character(as.integer(BiocGenerics::start(x)) - 1L),
      V3 = as.character(as.integer(BiocGenerics::end(x))),
      V4 = as.character(nm),
      V5 = ".",
      V6 = as.character(GenomicRanges::strand(x)),
      stringsAsFactors = FALSE
    )
    return(df)
  }
  abort_invalid_arg(c(
    "{.arg transcripts} must be a file path, BED data frame, GTF data frame, or GRanges.",
    "x" = "Got {.cls {class(x)[1]}}."
  ))
}

# Read a tab-separated BED file, or turn a data.frame, into a character-typed
# BED data.frame.
.read_bed_df <- function(x, kind) {
  if (is.character(x)) {
    if (length(x) != 1L || is.na(x) || !nzchar(x)) {
      abort_invalid_arg(c(
        "{.arg {kind}} must be a single non-empty file path.",
        "x" = "Got {.cls character} of length {length(x)}."
      ))
    }
    if (!file.exists(x)) {
      abort_not_found(c(
        "{.arg {kind}} file does not exist.",
        "x" = "{.path {x}}"
      ))
    }
    info <- file.info(x)
    if (isTRUE(info$isdir)) {
      abort_invalid_arg(c(
        "{.arg {kind}} path is a directory, not a file.",
        "x" = "{.path {x}}"
      ))
    }
    if (isTRUE(info$size == 0)) {
      abort_invalid_bed(c(
        "{.arg {kind}} file is empty.",
        "x" = "{.path {x}}"
      ))
    }
    df <- tryCatch(
      data.table::fread(
        x,
        header       = FALSE,
        sep          = "\t",
        quote        = "",
        fill         = TRUE,
        colClasses   = "character",
        data.table   = FALSE,
        showProgress = FALSE
      ),
      error = function(e) {
        abort_invalid_bed(c(
          "Could not read {.arg {kind}} as a tab-separated BED file.",
          "x" = "Reader error: {conditionMessage(e)}",
          "i" = "Expected a TSV with chr, start, end, name, score, strand (>=6 cols)."
        ))
      }
    )
    if (!is.data.frame(df) || nrow(df) == 0L) {
      abort_invalid_bed(c(
        "{.arg {kind}} file parsed to zero rows.",
        "x" = "{.path {x}}",
        "i" = "Is the file tab-separated with no header line?"
      ))
    }
    return(df)
  }

  if (is.data.frame(x)) {
    if (nrow(x) == 0L) {
      abort_invalid_bed("{.arg {kind}} data frame has no rows.")
    }
    df <- as.data.frame(
      lapply(x, function(col) as.character(col)),
      stringsAsFactors = FALSE,
      optional         = TRUE
    )
    names(df) <- names(x)
    return(df)
  }

  abort_invalid_arg(c(
    "{.arg {kind}} must be a file path (character) or a data frame.",
    "x" = "Got {.cls {class(x)[1]}}."
  ))
}

# Parse a BED coordinate column to integer with a clear error pointing at the
.coerce_coord <- function(x, kind, what) {
  vals <- suppressWarnings(as.integer(x))
  bad  <- is.na(vals)
  if (any(bad)) {
    i <- which(bad)[1L]
    abort_invalid_bed(c(
      "{.arg {kind}} {.field {what}} column has non-integer or missing values.",
      "x" = "First bad row: {i}, value {.val {as.character(x)[i]}}.",
      "i" = "BED coordinates must be non-negative integers."
    ))
  }
  if (any(vals < 0L)) {
    i <- which(vals < 0L)[1L]
    abort_invalid_bed(c(
      "{.arg {kind}} {.field {what}} column has negative values.",
      "x" = "First bad row: {i}, value {.val {vals[i]}}."
    ))
  }
  vals
}

# Normalize a BED strand column: accept "+", "-", ".", "*"; map "." to "*".
.coerce_strand <- function(x, kind) {
  s <- as.character(x)
  s[is.na(s) | s == "."] <- "*"
  valid <- s %in% c("+", "-", "*")
  if (!all(valid)) {
    i <- which(!valid)[1L]
    abort_invalid_bed(c(
      "{.arg {kind}} {.field strand} must be {.val +}, {.val -}, {.val .}, or {.val *}.",
      "x" = "First bad row: {i}, value {.val {as.character(x)[i]}}."
    ))
  }
  s
}

# Build a GRanges from a character-typed BED data.frame.
.bed_to_gr <- function(df, kind) {
  starts <- .coerce_coord(df[[2L]], kind, "start")
  ends   <- .coerce_coord(df[[3L]], kind, "end")
  if (any(ends < starts)) {
    i <- which(ends < starts)[1L]
    abort_invalid_bed(c(
      "{.arg {kind}} has rows where {.field end} < {.field start}.",
      "x" = "First bad row: {i} ({starts[i]}..{ends[i]})."
    ))
  }
  strands <- .coerce_strand(df[[6L]], kind)
  GenomicRanges::GRanges(
    seqnames = as.character(df[[1L]]),
    ranges   = IRanges::IRanges(start = starts + 1L, end = ends),
    strand   = strands
  )
}
