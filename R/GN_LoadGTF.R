#' Load GTF Annotation from AnnotationHub
#'
#' Loads GTF gene annotation data from AnnotationHub for Human or Mouse.
#' This function can be called once to load the annotation, which can then
#' be passed to other functions like `PlotGene()` or `PlotRegion()` to avoid
#' repeated downloads.
#'
#' @param species Species to load annotation for: "Human" or "Mouse".
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
LoadGTF <- function(species = "Human") {
  if (!species %in% c("Human", "Mouse"))

    stop("Species must be 'Human' or 'Mouse'")

  ah <- AnnotationHub::AnnotationHub(ask = FALSE)

  if (species == "Human") {
    # Ensembl Human GTF
    gtf <- data.frame(ah[["AH110867"]])
  } else if (species == "Mouse") {
    # Ensembl Mouse GTF
    gtf <- data.frame(ah[["AH47076"]])
  }

  return(gtf)
}
