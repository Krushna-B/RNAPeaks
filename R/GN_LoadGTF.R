
# Loads GTF gene annotation data from AnnotationHub for Human or Mouse.
#' Can be called once and the result passed to \code{PlotGene()} or
#' \code{PlotRegion()} to avoid repeated downloads.
#'
#' @param species Species to load annotation for: "Human" or "Mouse".
#' @param file Optional file path to a local GTF file. If provided, imports
#'   directly without connecting to AnnotationHub.
#'
#' @return A data frame containing GTF annotation with columns including
#'   seqnames, start, end, strand, type, gene_id, gene_name, transcript_id, etc.
#'
#' @details
#' Loading annotations from AnnotationHub can take time on first use.
#' By calling this function separately, you can:
#' \itemize{
#'   \item Load the annotation once and reuse it across multiple plots
#'   \item Save the annotation to disk for faster future sessions
#' }
#'
#' Human annotations use Ensembl GTF (AH110867).
#' Mouse annotations use Ensembl GTF (AH47076).
#'
#' @export
#'
#' @examples
#' \dontrun{
#'   # Load human GTF once
#'   gtf <- LoadGTF(species = "Human")
#'
#'   # Optionally save for future sessions
#'   saveRDS(gtf, "human_gtf.rds")
#'
#'   # Use in multiple plots without reloading
#'   PlotGene(bed = bed, geneID = "TP53", gtf = gtf)
#'   PlotGene(bed = bed, geneID = "BRCA1", gtf = gtf)
#'
#'   # Load from saved file in future sessions
#'   gtf <- readRDS("human_gtf.rds")
#' }
LoadGTF <- function(species = "Human", file = NULL) {
  if (!is.null(file)) {
    gtf <- data.frame(rtracklayer::import(file))
    return(gtf)
  }

  if (!species %in% c("Human", "Mouse"))
    stop("Species must be 'Human' or 'Mouse'")

  # Skip network calls if cache is already populated in deployment setting
  cache_exists <- tryCatch({
    cache_dir <- AnnotationHub::getAnnotationHubOption("CACHE")
    file.exists(cache_dir) && length(list.files(cache_dir)) > 0
  }, error = function(e) FALSE)

  ah <- AnnotationHub::AnnotationHub(ask = FALSE, localHub = cache_exists)

  if (species == "Human") {
    # Ensembl Human GTF
    gtf <- data.frame(ah[["AH110867"]])
  } else if (species == "Mouse") {
    # Ensembl Mouse GTF
    gtf <- data.frame(ah[["AH47076"]])
  }

  return(gtf)
}
