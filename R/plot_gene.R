#' Plot RBP peaks over a single transcript structure
#'
#' Picks one transcript from `gtf` and overlays peaks from `bed`.
#'
#' @param bed Validated BED data frame or named list of BED data frames.
#' @param gene Gene symbol or Ensembl gene id.
#' @param transcript Optional Ensembl transcript id or transcript name.
#'   Defaults to the longest transcript of `gene`.
#' @param gtf Optional path to a local GTF file. If `NULL`, the bundled
#'   annotation for `species` is used.
#' @param species One of `"hg38"`, `"mm10"`, or `"mm39"`. Ignored when
#'   `gtf` is supplied.
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
      # Required args have no default; catch them here so a forgotten
      # argument reports as a clear validation error instead of falling
      # through to the generic "unexpected error" branch below.
      if (missing(bed)) {
        abort_invalid_arg(c(
          "{.arg bed} is required.",
          "i" = "Supply a BED data frame, a file path, or a named list of either."
        ))
      }
      if (missing(gene)) {
        abort_invalid_arg(c(
          "{.arg gene} is required.",
          "i" = "Supply a gene symbol (e.g. {.val APOE}) or an Ensembl gene id."
        ))
      }

      # 0. Normalize string inputs
      gene       <- normalize_str(gene)
      transcript <- normalize_str(transcript)
      species    <- normalize_str(species)

      # 1. Determine annotation source
      cli::cli_progress_step("Loading GTF")
      gtf <- get_GTF(species = species, file = gtf)

      # 2. Select a single transcript
      cli::cli_progress_step("Selecting transcript {.val {gene}}")
      tx <- select_transcript(gtf, geneID = gene, TxID = transcript)
      cli::cli_progress_done()

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
