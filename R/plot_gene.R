#' Plot RBP peaks over a single transcript structure
#'
#' Picks one transcript from `gtf` and overlays peaks from `bed`.
#'
#' @param bed Validated BED data frame or named list of BED data frames.
#' @param gene Gene symbol or Ensembl gene id.
#' @param transcript Optional Ensembl transcript id or transcript name.
#'   Defaults to the longest transcript of `gene`.
#' @param gtf Optional GTF data frame, or path to a local GTF file.
#'   If `NULL`, the bundled annotation for `species` is used.
#' @param species One of `"hg38"`, `"mm10"`, or `"mm39"`.
#' @param peaks_opts Output of [peaks_options()]: BED filtering / ordering options.
#' @param style Output of [peaks_plot_style()]: visual settings.
#'
#' @return A ggplot object.
#' @export
#' @family plot
plot_gene <- function(bed,
                      gene,
                      transcript    = NULL,
                      gtf           = NULL,
                      species       = "hg38",
                      peaks_opts    = peaks_options(),
                      style         = peaks_plot_style()) {
  tryCatch(
    {
      # 1. Determine annotation source
      gtf <- get_GTF(species = species, file = gtf)

      # 2. Select a single transcript
      tx <- select_transcript(gtf, geneID = gene, TxID = transcript)

      # 3. Start Peaks Plotting Pipeline
      plot_peaks_pipeline(
        transcripts = tx,
        bed         = bed,
        is_region   = FALSE,
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
