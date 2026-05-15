#' Plot RBP peaks over a multi-gene genomic region
#'
#' Selects every gene whose body overlaps the window `chr:start-end` on
#' `strand` (one transcript per gene, longest by default), stacks their
#' structures vertically, and overlays peaks from `bed`.
#'
#' @param bed Validated BED data frame or named list of BED data frames.
#' @param chr Chromosome (with or without `"chr"` prefix; normalized).
#' @param start,end Region bounds (bp, `start <= end`).
#' @param strand `"+"` or `"-"`.
#' @param gtf Optional GTF data frame, or path to a local GTF file.
#'   If `NULL`, the bundled annotation for `species` is used.
#' @param species One of `"hg38"`, `"mm10"`, or `"mm39"`. Ignored when
#'   `gtf` is supplied.
#' @param peaks_opts Output of [peaks_options()]: BED filtering / ordering options.
#' @param style Output of [peaks_plot_style()]: visual settings.
#'
#' @return A ggplot object.
#' @export
#' @family plot
plot_region <- function(bed,
                        chr,
                        start,
                        end,
                        strand,
                        gtf           = NULL,
                        species       = "hg38",
                        peaks_opts    = peaks_options(),
                        style         = peaks_plot_style()) {
  tryCatch(
    {
      # 0. Normalize inputs
      chr     <- normalize_chr(chr)
      strand  <- normalize_str(strand)
      species <- normalize_str(species)
      start   <- normalize_coord(start, "start")
      end     <- normalize_coord(end,   "end")

      # 1. Resolve annotation source
      gtf <- get_GTF(species = species, file = gtf)

      # 2. Select all transcripts overlapping the window one per gene
      txs <- select_region(gtf, chr = chr, start = start, end = end, strand = strand)

      # 3. Start Peaks Plotting Pipeline
      plot_peaks_pipeline(
        transcripts = txs,
        bed         = bed,
        is_region   = TRUE,
        peaks_opts  = peaks_opts,
        style       = style
      )
    },
    error = function(cnd) {
      if (inherits(cnd, "rnapeaks_error")) {
        cli::cli_abort("Failed to generate plot.", parent = cnd)
      } else {
        cli::cli_abort(
          c("Failed to generate plot.",
            "x" = "An unexpected error occurred."),
          parent = cnd
        )
      }
    }
  )
}
