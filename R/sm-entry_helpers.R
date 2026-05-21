#' Internal: shared helpers for splicing / sequence map entry points
#'
#' @keywords internal
#' @name sm_entry_helpers


#Read, validate, and turn a BED into a reduced GRanges of peaks.
.peaks_to_granges <- function(bed_file) {
  if (is.character(bed_file)) {
    if (!file.exists(bed_file)) {
      abort_not_found(c(
        "BED file does not exist.",
        "x" = "Path: {.path {bed_file}}."
      ))
    }
    bed <- utils::read.table(bed_file)
  } else {
    bed <- bed_file
  }
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


#Resolve `genome` to a BSgenome object. Shape (NULL, string, or BSgenome) is
#already validated by validate_sm_inputs() upstream; we only check that a
#string value is one of the supported genome keys.
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
      "{.arg genome} must be one of {.or {.val hg38} {.val mm10} {.val mm39}}, or a BSgenome object.",
      "x" = "Got {.val {genome}}."
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


#Shape-check every argument supplied to a splicing / sequence map entry
#point. Optional args (`bed_file`, `sequence`, `genome`) are validated only
#when supplied at the call site. Every failure ends with a "Did the
#argument order get mixed up?" hint so misordered calls (e.g. forgetting
#`genome =` and having `opts` interpreted as the genome) surface clearly.
validate_sm_inputs <- function(events, opts, style, title,
                                bed_file, sequence, genome, motif_mode) {
  check_string(title, "title")

  if (!is.data.frame(events)) {
    abort_invalid_arg(c(
      "{.arg events} must be a data frame of rMATS events.",
      "x" = "Got {.cls {class(events)[1]}}.",
      "i" = "Did the argument order get mixed up?"
    ))
  }

  if (!is.list(opts) || is.null(opts$width_exon) ||
      is.null(opts$moving_average)) {
    abort_invalid_arg(c(
      "{.arg opts} must be a {.fn splicing_options} result.",
      "x" = "Got {.cls {class(opts)[1]}}.",
      "i" = "Did the argument order get mixed up?"
    ))
  }

  if (!is.list(style) || is.null(style$show_significance)) {
    abort_invalid_arg(c(
      "{.arg style} must be a {.fn splicing_style} result.",
      "x" = "Got {.cls {class(style)[1]}}.",
      "i" = "Did the argument order get mixed up?"
    ))
  }

  if (!missing(bed_file)) {
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

  if (!missing(sequence)) {
    if (!is.character(sequence) || length(sequence) == 0L ||
        anyNA(sequence) || any(nchar(trimws(sequence)) == 0L)) {
      abort_invalid_arg(c(
        "{.arg sequence} must be a non-empty character vector of motifs.",
        "x" = "Got {.cls {class(sequence)[1]}} of length {length(sequence)}.",
        "i" = "Did the argument order get mixed up?"
      ))
    }
  }

  if (!missing(genome)) {
    if (!is.null(genome) && !is.character(genome) &&
        !methods::is(genome, "BSgenome")) {
      abort_invalid_arg(c(
        "{.arg genome} must be {.val NULL}, a string ({.val hg38} / {.val mm10} / {.val mm39}), or a BSgenome object.",
        "x" = "Got {.cls {class(genome)[1]}}.",
        "i" = "Did the argument order get mixed up?"
      ))
    }
  }

  if (!missing(motif_mode)) {
    valid_motif_modes <- c("combined", "individual")
    if (!is.character(motif_mode) || length(motif_mode) != 1L ||
        is.na(motif_mode) || !motif_mode %in% valid_motif_modes) {
      abort_invalid_arg(c(
        "{.arg motif_mode} must be one of {.or {.val {valid_motif_modes}}}.",
        "x" = "Got {.val {motif_mode}}.",
        "i" = "Did the argument order get mixed up?"
      ))
    }
  }
}


#Wrap an entry-point body with unified error context, mirroring the pattern
#in plot_gene() / plot_region(): classed rnapeaks_error conditions are
#re-raised with just the context line, everything else gets an extra
#"unexpected error" note. `call` defaults to the entry function's frame so
#the outer "Error in `<fn>()`" line names the public API instead of this
#helper.
wrap_sm_errors <- function(map_name, body, call = parent.frame()) {
  tryCatch(
    body,
    error = function(cnd) {
      msg <- paste0("Failed to generate ", map_name, ".")
      if (inherits(cnd, "rnapeaks_error")) {
        cli::cli_abort(msg, parent = cnd, call = call)
      } else {
        cli::cli_abort(c(msg, "x" = "An unexpected error occurred."),
                       parent = cnd, call = call)
      }
    }
  )
}


#Align regions to the genome's seqlevels and extract sequences once.
#`extension` adds bp at the transcript 3' end of each region so motifs can
#start in the last (motif_len - 1) positions. Returned `regions` keep their
#original widths.
.extract_region_seqs <- function(regions_gr, genome, extension = 0L) {
  empty <- list(regions = GenomicRanges::GRanges(),
                seqs    = Biostrings::DNAStringSet())
  if (length(regions_gr) == 0L) return(empty)

  aligned <- .align_seqnames_to_genome(regions_gr, genome)
  if (length(aligned) == 0L) return(empty)

  to_extract <- aligned
  if (extension > 0L) {
    is_minus <- as.character(GenomicRanges::strand(aligned)) == "-"
    seqlens  <- GenomeInfoDb::seqlengths(genome)
    chr_lens <- seqlens[as.character(GenomeInfoDb::seqnames(aligned))]

    new_start <- pmax(GenomicRanges::start(aligned) -
                        ifelse(is_minus, extension, 0L), 1L)
    new_end   <- pmin(GenomicRanges::end(aligned) +
                        ifelse(is_minus, 0L, extension), chr_lens)
    GenomicRanges::ranges(to_extract) <- IRanges::IRanges(
      start = new_start, end = new_end
    )
  }

  seqs <- tryCatch(
    BSgenome::getSeq(genome, to_extract),
    error = function(e) {
      cli::cli_warn(c(
        "Sequence extraction failed.",
        "x" = "{e$message}",
        "i" = "Check that {.arg genome} matches the chromosome naming of your events."
      ))
      Biostrings::DNAStringSet()
    }
  )

  list(regions = aligned, seqs = seqs)
}


.prepare_sequence_map_prep <- function(events, schema, opts, genome, motifs) {
  extension <- max(nchar(motifs)) - 1L
  prep <- prepare_event_map(events, schema, opts)
  prep$genome <- genome

  seqs_by_group <- vector("list", length(prep$groups_idx))
  names(seqs_by_group) <- names(prep$groups_idx)
  for (g in names(prep$groups_idx)) {
    rs <- .extract_region_seqs(prep$regions_by_group[[g]], genome,
                                extension = extension)
    prep$regions_by_group[[g]] <- rs$regions
    seqs_by_group[[g]]         <- rs$seqs
  }
  prep$seqs_by_group <- seqs_by_group
  prep
}


#Branches out to combined / individual motif modes. Shared by SE and RI
#sequence-map entry points.
.run_sequence_map <- function(motifs, motif_mode, prep,
                              schema, opts, style, title) {
  one_run <- function(motif_set, run_title) {
    scorer <- function(regions_gr, n_regions, region_width, group_name) {
      motif_scorer(regions_gr, prep$seqs_by_group[[group_name]],
                   motif_set, n_regions, region_width)
    }
    event_map_pipeline(
      schema  = schema,
      scorer  = scorer,
      opts    = opts,
      style   = style,
      plot_fn = plot_event_map,
      title   = run_title,
      prep    = prep
    )
  }

  if (motif_mode == "combined") {
    return(one_run(motifs, title))
  }

  results <- lapply(motifs, function(m) {
    sub_title <- if (nzchar(title)) paste0(title, " - ", m) else m
    one_run(m, sub_title)
  })
  names(results) <- motifs
  results
}


#Uppercase, trim, convert U->T, then check every motif consists only of
#IUPAC nucleotide codes. Shape (non-empty character, no NAs, no blanks) is
#already validated by validate_sm_inputs() upstream; this guards against
#characters that would otherwise blow up inside Biostrings with messages
#like "key 70 (char 'F') not in lookup table".
.normalize_motifs <- function(sequence) {
  sequence <- toupper(trimws(sequence))

  has_u <- grepl("U", sequence, fixed = TRUE)
  if (any(has_u)) {
    cli::cli_inform("Converting U to T in motif{?s}: {.val {sequence[has_u]}}.")
    sequence <- gsub("U", "T", sequence, fixed = TRUE)
  }

  invalid <- !grepl("^[ACGTRYSWKMBDHVN]+$", sequence)
  if (any(invalid)) {
    bad_motifs <- sequence[invalid]
    abort_invalid_arg(c(
      "{.arg sequence} contains invalid IUPAC motif{?s}: {.val {bad_motifs}}.",
      "i" = "Allowed: {.val A C G T U} (canonical) and {.val R Y S W K M B D H V N} (ambiguity codes)."
    ))
  }

  sequence
}
