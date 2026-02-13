

#' Retrieves gene or transcript annotation from a GTF data frame or
#' automatically loads annotations from AnnotationHub for Human or Mouse.
#'
#' @param geneID Gene identifier (gene symbol like "TP53" or Ensembl ID
#'   starting with "ENSG").
#' @param species Species for automatic annotation loading: "Human" or "Mouse".
#' @param TxID Optional specific transcript ID. If NA, the longest transcript
#'   is automatically selected.
#' @param gtf Optional pre-loaded GTF annotation data frame. If NULL,
#'   annotations are loaded from AnnotationHub.
#'
#' @return A data frame containing GTF annotation rows for the specified
#'   gene/transcript, including exons, UTRs, and CDS features.
#'
#' @details
#' If both `geneID` and `TxID` are provided, the function verifies that the

#' transcript belongs to the specified gene. If only `geneID` is provided,
#' the longest transcript is selected.
#'
#'
#'
#' @noRd
#' @examples
#' \dontrun{
#'   # Get TP53 annotation (longest transcript)
#'   gene_data <- GetGene(geneID = "TP53", species = "Human")
#'
#'   # Get specific transcript
#'   gene_data <- GetGene(
#'     geneID = "TP53",
#'     TxID = "ENST00000269305",
#'     species = "Human"
#'   )
#' }
GetGene <- function(geneID, species, TxID, gtf) {
  if (is.null(geneID) & is.null(TxID)) {
    stop("Enter Gene or Transcript ID!")
  }

  if (is.null(gtf)) {
    gtf <- LoadGTF(species = species)
  }

  df <- gtf

  if (!is.na(TxID) & is.null(geneID)) {
    # Gets transcript if only transcript id is provided
    gtf <- gtf[which(gtf$transcript_id == TxID), ]
    if (is.null(gtf) | nrow(gtf) == 0) {
      stop("Check Transcript ID")
    }
  } else if (is.na(TxID) & !is.null(geneID)) {
    if (grepl("ENSG", geneID)) {
      gtf <- gtf[which(gtf$gene_id == geneID), ]
    } else {
      # Gets geneID if only geneID is provided
      gtf <- gtf[which(gtf$gene_name == geneID), ]
    }
    if (is.null(gtf) | nrow(gtf) == 0) {
      stop("Check Gene ID")
    }
    transcripts_all <- unique(gtf[which(gtf$type == "transcript"), c("transcript_id", "width")])
    TxID <- transcripts_all$transcript_id[which.max(transcripts_all$width)]
    gtf <- gtf[which(gtf$transcript_id == TxID), ]
  } else {
    # Check if both geneID and transcript ID are provided
    all_possible_transcripts <- unique(gtf[which(gtf$type == "transcript"), c("transcript_id", "width")])
    if (!TxID %in% all_possible_transcripts$transcript_id) {
      stop("Transcript ID not in Gene")
    }
    gtf <- gtf[which(gtf$transcript_id == TxID), ]
  }
  return(gtf)
}
