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
#' @param bam_files Optional named character vector of BAM file paths drawn as
#'   coverage tracks above the gene structure. Names become track labels; if
#'   unnamed, the filename (without extension) is used. Each BAM must be
#'   sorted and have a `.bai` index next to it.
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
                      bam_files     = NULL,
                      peaks_opts    = peaks_options(),
                      style         = peaks_plot_style()) {
  tryCatch(
    {
      # 0. Normalize string inputs
      gene       <- normalize_str(gene)
      transcript <- normalize_str(transcript)
      species    <- normalize_str(species)

      # 1. Determine annotation source
      gtf <- get_GTF(species = species, file = gtf)

      # 2. Select a single transcript
      tx <- select_transcript(gtf, geneID = gene, TxID = transcript)

      # 3. Start Peaks Plotting Pipeline
      plot_peaks_pipeline(
        transcripts = tx,
        bed         = bed,
        is_region   = FALSE,
        bam_files   = bam_files,
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
