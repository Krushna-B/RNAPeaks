# Generate strand and region-matched control peaks.

#' Annotates each peak with the transcripts that fully contain it via
#' [intersect_peaks()] (bedtools `-f 1 -wa -wb -s` semantics), then chooses
#' a control region per peak matched in length, strand, and transcript-region
#' type (UTR3 / UTR5 / CDS / exon / intron).
#'
#' @param raw_peaks Peak BED. File path or data.frame in BED layout (>= 6 cols).
#
#' @param anno      Per-region annotation BED. File path or data.frame.
#'   Either has named columns `chr, start, end, name, strand, gene` (with
#'   optional `transcript`, `region`) or is positional with col 1=chr,
#'   2=start, 3=end, 4=name (`{transcript}_{region}_*`), 6=strand, 7=gene.
#' @param gene      Gene BED. File path or data.frame with named columns
#'   `chr, start, end, gene` or positional cols 1=chr, 2=start, 3=end, 4=gene.
#' @param transcripts Transcripts BED passed to [intersect_peaks()] (file
#'   path or data.frame, same format used in the bedtools intersect example).
#' @param pool      Parallel cores (>= 1).
#' @param seed      RNG seed. Default `1234L`.
#'
#' @return A data.frame with columns `chr`, `peak_start`, `peak_end`, `name`,
#'   `strand`, `control_start`, `control_end`.
#'
#' @export
#' @family control_peaks
generate_control_peaks <- function(raw_peaks, anno, gene, transcripts,
                                   pool = 1L, seed = 1234L) {
  #Validate Params
  check_scalar_int(pool, "pool", min = 1L)
  check_scalar_int(seed, "seed")

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

  # Dispatch each chromosome to the C++ engine. process_chromosome_cpp()
  # reseeds std::mt19937 with (seed + chrom_index) on entry.
  result_list <- lapply(seq_along(chromosomes), function(k) {
    .run_chrom_cpp(chromosomes[k], k, peaks_grouped, anno_dt, genes_dt, seed)
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

# Run the C++ engine for one chromosome. Slices anno / genes / peaks for the
# chromosome, hands them to process_chromosome_cpp, and attaches `chr` to the
# returned rows. Returns NULL when there is nothing to process.
.run_chrom_cpp <- function(chrom, k, peaks_dt, anno_dt, genes_dt, seed) {
  pk <- peaks_dt[chr == chrom]
  an <- anno_dt[chr == chrom]
  ge <- genes_dt[chr == chrom]

  if (nrow(pk) == 0L || nrow(an) == 0L) return(NULL)

  res <- process_chromosome_cpp(
    anno_line_id    = an$line_id,
    anno_transcript = an$transcript,
    anno_region     = an$region,
    anno_start      = an$start,
    anno_end        = an$end,
    anno_strand     = an$strand,
    anno_gene       = an$gene,
    gene_name       = ge$gene,
    gene_start      = ge$start,
    gene_end        = ge$end,
    peak_start       = pk$start,
    peak_end         = pk$end,
    peak_FC          = pk$FC,
    peak_range       = pk$peak_range,
    peak_transcripts = pk$transcripts,
    seed             = as.integer(seed) + as.integer(k)
  )

  if (length(res$peak_start) == 0L) return(NULL)
  data.frame(chr = chrom, res, stringsAsFactors = FALSE)
}

