#' Get GTF gene annotation for a supported species or from a local file.
#'
#' Returns a GTF annotation as a data frame, either from the bundled
#' \code{gtf_hg38}, \code{gtf_mm10}, or \code{gtf_mm39} datasets, or by
#' importing a user-supplied GTF file.
#'
#' @param species One of \code{"hg38"}, \code{"mm10"}, or \code{"mm39"}.
#'   Ignored when \code{file} is supplied.
#'
#' @param file Optional path to a local GTF file.
#'
#' @return A data frame containing GTF annotation with at least the columns
#'   \code{seqnames, start, end, strand, type, gene_id, gene_name,
#'   transcript_id}.
#' @export
#' @family gene
#'
getGTF <- function(species = "hg38", file = NULL) {
  if (!is.null(file)) {
    # File arg must be a single, non-NA character path
    if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file)) {
      abort_invalid_arg(c(
        "{.arg file} must be a single non-empty file path or {.code NULL}.",
        "x" = "You supplied {.cls {class(file)[1]}} of length {length(file)}."
      ))
    }
    if (!file.exists(file)) {
      abort_not_found(c(
        "{.arg file} does not exist.",
        "x" = "Path: {.path {file}}."
      ))
    }

    gtf <- tryCatch(
      data.frame(rtracklayer::import(
        file,
        colnames = c("type", "gene_id", "gene_name", "transcript_id", "gene_biotype")
      )),
      error = function(e) {
        abort_invalid_gtf(c(
          "Failed to import {.arg file} as GTF.",
          "x" = "{conditionMessage(e)}",
          "i" = "Path: {.path {file}}."
        ))
      }
    )
    verify_gtf(gtf)
    cli::cli_alert_success("GTF imported and validated successfully")
    return(gtf)
  }

  if (!is.character(species) || length(species) != 1L || is.na(species)) {
    abort_invalid_arg(c(
      "{.arg species} must be a single string.",
      "x" = "You supplied {.cls {class(species)[1]}} of length {length(species)}."
    ))
  }
  valid <- c("hg38", "mm10", "mm39")
  if (!species %in% valid) {
    abort_invalid_arg(c(
      "{.arg species} must be one of {.val {valid}}.",
      "x" = "You supplied {.val {species}}."
    ))
  }

  dataset <- switch(species,
    hg38 = "gtf_hg38",
    mm10 = "gtf_mm10",
    mm39 = "gtf_mm39"
  )
  gtf <- tryCatch(
    get(dataset, envir = asNamespace("RNAPeaks"), inherits = FALSE),
    error = function(e) {
      abort_not_found(c(
        "Bundled annotation {.val {dataset}} could not be loaded.",
        "x" = "{conditionMessage(e)}",
        "i" = "Pass a local GTF via {.arg file} instead."
      ))
    }
  )
  verify_gtf(gtf)
  return(gtf)
}

# Verify a GTF data frame has the columns required by downstream functions
verify_gtf <- function(gtf) {
  if (!is.data.frame(gtf)) {
    abort_invalid_gtf(c(
      "GTF must be a data frame.",
      "x" = "Got {.cls {class(gtf)[1]}}."
    ))
  }
  required <- c("seqnames", "start", "end", "strand", "type",
                "gene_id", "gene_name", "transcript_id")
  missing <- setdiff(required, colnames(gtf))
  if (length(missing)) {
    abort_invalid_gtf(c(
      "GTF is missing required column{?s}: {.field {missing}}.",
      "i" = "Expected columns: {.field {required}}."
    ))
  }
  if (nrow(gtf) == 0L) {
    abort_invalid_gtf("GTF has zero rows.")
  }
  invisible(TRUE)
}
