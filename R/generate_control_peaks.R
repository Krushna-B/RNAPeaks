#' Generate strand- and region-matched control peaks
#'
#' For each peak, finds the transcripts that fully contain it (via
#' [intersect_peaks()] with bedtools `-f 1 -wa -wb -s` semantics) and samples
#' a control region matched in length, strand, and region type
#' (UTR3 / UTR5 / CDS / exon / intron).
#'
#' @param raw_peaks Peak BED. File path or data.frame with >= 6 BED columns.
#' @param anno Per-region annotation BED. File path or data.frame. Either
#'   has named columns `chr, start, end, name, strand, gene` (with optional
#'   `transcript`, `region`), or is positional with col 1=chr, 2=start,
#'   3=end, 4=name (`{transcript}_{region}_*`), 6=strand, 7=gene.
#' @param gene Gene BED. Named columns `chr, start, end, gene`, or positional
#'   cols 1=chr, 2=start, 3=end, 4=gene.
#' @param transcripts Transcripts BED passed to [intersect_peaks()].
#' @param threads Number of worker threads (>= 1).
#' @param seed RNG seed.
#'
#' @return A data.frame with columns `chr`, `peak_start`, `peak_end`, `name`,
#'   `strand`, `control_start`, `control_end`.
#'
#' @export
#' @family control_peaks
generate_control_peaks <- function(raw_peaks, anno, gene, transcripts,
                                   threads = 1L, seed = 1234) {
  #Validate Params
  check_scalar_int(threads, "threads", min = 1L)
  check_scalar_int(seed, "seed")
  # Required inputs have no default; report a forgotten one as a clear
  # validation error instead of a base-R "missing argument" error.
  if (missing(raw_peaks)) {
    abort_invalid_arg(c("{.arg raw_peaks} is required.",
                        "i" = "Supply a file path or a data.frame."))
  }
  if (missing(anno)) {
    abort_invalid_arg(c("{.arg anno} is required.",
                        "i" = "Supply a file path or a data.frame."))
  }
  if (missing(gene)) {
    abort_invalid_arg(c("{.arg gene} is required.",
                        "i" = "Supply a file path or a data.frame."))
  }
  if (missing(transcripts)) {
    abort_invalid_arg(c("{.arg transcripts} is required.",
                        "i" = "Supply a file path or a data.frame."))
  }
  for (nm in c("raw_peaks", "anno", "gene", "transcripts")) {
    val <- get(nm)
    if (is.null(val) || (!is.character(val) && !is.data.frame(val))) {
      abort_invalid_arg(c(
        "{.arg {nm}} must be a file path or data.frame.",
        "x" = "Got {.cls {class(val)[1]}}."
      ))
    }
  }

  cli::cli_progress_step("Intersecting peaks with transcripts")
  peaks_df <- intersect_peaks(raw_peaks, transcripts,
                              fraction = 1.0, same_strand = TRUE)
  if (ncol(peaks_df) < 12L) {
    abort_invalid_bed(c(
      "intersect_peaks() returned fewer than 12 columns.",
      "x" = "Got {ncol(peaks_df)}.",
      "i" = "Expected 8 peak columns + 6 transcript columns = 14."
    ))
  }

  #Create input tables
  cli::cli_progress_step("Reading annotation")
  anno_dt <- read_annotation(anno)

  cli::cli_progress_step("Reading gene table")
  genes_dt <- read_genes(gene)

  #Group peaks for the engine
  peaks_grouped <- .group_peaks_for_engine(peaks_df)

  #Identify common chromosomes
  chromosomes <- intersect(unique(peaks_grouped$chr), unique(anno_dt$chr))
  if (length(chromosomes) == 0L) {
    cli::cli_progress_done()
    cli::cli_alert_warning("No chromosomes shared between peaks and annotation.")
    return(.empty_control_df())
  }

  cli::cli_progress_step(
    "Generating controls across {length(chromosomes)} chromosome{?s}"
  )

  # Build per-chromosome payload list (main thread only - touches Rcpp).
  # Each chromosome's RNG is seeded with (seed + k) so output is deterministic
  # regardless of thread scheduling.
  per_chrom <- lapply(seq_along(chromosomes), function(k) {
    chrom <- chromosomes[k]
    pk <- peaks_grouped[chr == chrom]
    an <- anno_dt[chr == chrom]
    ge <- genes_dt[chr == chrom]
    if (nrow(pk) == 0L || nrow(an) == 0L) return(NULL)
    list(
      anno_transcript  = an$transcript,
      anno_region      = an$region,
      anno_start       = an$start,
      anno_end         = an$end,
      anno_strand      = an$strand,
      anno_gene        = an$gene,
      gene_name        = ge$gene,
      gene_start       = ge$start,
      gene_end         = ge$end,
      peak_start       = pk$start,
      peak_end         = pk$end,
      peak_FC          = pk$FC,
      peak_transcripts = pk$transcripts,
      seed             = as.integer(seed) + as.integer(k),
      chrom            = chrom
    )
  })

  keep <- !vapply(per_chrom, is.null, logical(1L))
  per_chrom <- per_chrom[keep]
  chrom_names <- vapply(per_chrom, function(x) x$chrom, character(1L))

  # One call into C++. Engine errors are re-raised as rnapeaks_error.
  raw <- tryCatch(
    process_chromosomes_threaded_cpp(per_chrom, as.integer(threads)),
    rnapeaks_error = function(e) stop(e),
    error = function(e) {
      abort_invalid_arg(c(
        "Control peak engine failed.",
        "x" = conditionMessage(e)
      ))
    }
  )

  result_list <- lapply(seq_along(raw), function(i) {
    r <- raw[[i]]
    if (length(r$peak_start) == 0L) return(NULL)
    data.frame(chr = chrom_names[i], r, stringsAsFactors = FALSE)
  })

  #Combine results
  total <- data.table::rbindlist(result_list, fill = TRUE)
  cli::cli_progress_done()
  cli::cli_alert_info("Generated {nrow(total)} control peak{?s}.")

  if (nrow(total) == 0L) return(.empty_control_df())

  data.frame(
    chr           = total$chr,
    peak_start    = as.integer(total$peak_start),
    peak_end      = as.integer(total$peak_end),
    name          = total$name,
    strand        = total$strand,
    control_start = as.integer(total$control_start),
    control_end   = as.integer(total$control_end),
    stringsAsFactors = FALSE
  )
}

# Read the per-region annotation table. Accepts a file path or a data.frame.
# Returns a data.table with columns:
read_annotation <- function(anno) {
  if (!is.character(anno) && !is.data.frame(anno)) {
    abort_invalid_arg(c(
      "{.arg anno} must be a file path or data.frame.",
      "x" = "Got {.cls {class(anno)[1]}}."
    ))
  }
  df <- .read_bed_df(anno, "anno")
  has_named <- all(c("chr", "start", "end", "name", "strand", "gene") %in% colnames(df))
  if (has_named) {
    dt <- data.table::as.data.table(df)
  } else {
    if (ncol(df) < 7L) {
      abort_invalid_bed(c(
        "{.arg anno} positional BED needs at least 7 columns.",
        "x" = "Got {ncol(df)}.",
        "i" = "Expected 1=chr, 2=start, 3=end, 4=name, 6=strand, 7=gene."
      ))
    }
    dt <- data.table::data.table(
      chr    = as.character(df[[1L]]),
      start  = as.integer(df[[2L]]),
      end    = as.integer(df[[3L]]),
      name   = as.character(df[[4L]]),
      strand = as.character(df[[6L]]),
      gene   = as.character(df[[7L]])
    )
  }
  dt[, start := as.integer(start)]
  dt[, end   := as.integer(end)]

  if (!"transcript" %in% colnames(dt) || !"region" %in% colnames(dt)) {
    parts <- data.table::tstrsplit(dt$name, "_", fixed = TRUE)
    if (length(parts) < 2L) {
      abort_invalid_bed(c(
        "{.arg anno}: cannot parse transcript/region from {.field name}.",
        "i" = "Expected {.val transcript_region_*} (e.g. ENST..._CDS_1)."
      ))
    }
    if (!"transcript" %in% colnames(dt)) dt[, transcript := parts[[1L]]]
    if (!"region"     %in% colnames(dt)) dt[, region     := parts[[2L]]]
  }
  dt[]
}

# Read the gene table. Accepts a file path or a data.frame. Returns a
# data.table with columns: chr, start, end, gene.
read_genes <- function(gene) {
  if (!is.character(gene) && !is.data.frame(gene)) {
    abort_invalid_arg(c(
      "{.arg gene} must be a file path or data.frame.",
      "x" = "Got {.cls {class(gene)[1]}}."
    ))
  }
  df <- .read_bed_df(gene, "gene")
  has_named <- all(c("chr", "start", "end", "gene") %in% colnames(df))
  if (has_named) {
    dt <- data.table::as.data.table(df)
    dt <- dt[, .(chr  = as.character(chr),
                 start = as.integer(start),
                 end   = as.integer(end),
                 gene  = as.character(gene))]
  } else {
    if (ncol(df) < 4L) {
      abort_invalid_bed(c(
        "{.arg gene} positional BED needs at least 4 columns.",
        "x" = "Got {ncol(df)}.",
        "i" = "Expected 1=chr, 2=start, 3=end, 4=gene."
      ))
    }
    dt <- data.table::data.table(
      chr   = as.character(df[[1L]]),
      start = as.integer(df[[2L]]),
      end   = as.integer(df[[3L]]),
      gene  = as.character(df[[4L]])
    )
  }
  dt[]
}

# Empty result with the final column shape.
.empty_control_df <- function() {
  data.frame(
    chr           = character(0L),
    peak_start    = integer(0L),
    peak_end      = integer(0L),
    name          = character(0L),
    strand        = character(0L),
    control_start = integer(0L),
    control_end   = integer(0L),
    stringsAsFactors = FALSE
  )
}

# Regroup intersect_peaks()'s 14-col output into the (chr, peak_range) shape
# Preserves original peak file order via first_line.
.group_peaks_for_engine <- function(df) {
  dt <- data.table::as.data.table(df)
  dt[, line_id := .I]
  dt[, chr        := as.character(V1)]
  dt[, start      := as.integer(V2)]
  dt[, end        := as.integer(V3)]
  dt[, FC         := as.character(V7)]
  dt[, transcript := as.character(V12)]
  dt[, peak_range := paste0(start, "-", end)]

  grouped <- dt[
    ,
    .(first_line  = min(line_id),
      start       = start[1L],
      end         = end[1L],
      FC          = FC[1L],
      transcripts = list(transcript)),
    by = .(chr, peak_range)
  ]
  data.table::setorder(grouped, first_line)
  grouped[]
}


