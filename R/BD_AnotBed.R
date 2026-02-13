
#' Annotates each peak in a BED file with its genomic region type (CDS, UTR,
#' intron, splice site, etc.) using Ensembl database annotations.
#'
#' @param peak A BED data frame or GRanges object containing peaks to annotate.
#' @param enDb Optional pre-loaded EnsDb object. If NULL, loaded from AnnotationHub.
#' @param species Species for annotation: "Human" or "Mouse".
#'
#' @return A data frame with original peak coordinates plus a `Region` column
#'   indicating the genomic feature type:
#'   \describe{
#'     \item{miRNA}{microRNA regions}
#'     \item{CDS}{Coding sequence (exons)}
#'     \item{UTR5}{5' untranslated region}
#'     \item{UTR3}{3' untranslated region}
#'     \item{NCExon}{Non-coding exons (lncRNA)}
#'     \item{SpliceSite5}{5' splice site (50bp flanking)}
#'     \item{SpliceSite3}{3' splice site (50bp flanking)}
#'     \item{Intron}{Intronic regions}
#'     \item{Others}{Intergenic or unclassified}
#'   }
#'
#' @details
#' Peaks are assigned to regions hierarchically - if a peak overlaps multiple
#' region types, it is assigned to the first matching category in the order
#' listed above. A peak must have >50% overlap with a region to be assigned.
#'
#' @noRd
#' @examples
#' \dontrun{
#'   # Annotate peaks from a BED file
#'   bed <- read.table("peaks.bed", header = FALSE)
#'   annotated <- AnotBed(peak = bed, species = "Human")
#'
#'   # Check region distribution
#'   table(annotated$Region)
#' }
AnotBed <- function(peak = NULL, enDb = NULL, species = "Human") {

  if (is.null(peak)) {
    stop("Provide peak bed file.")
  }
  if (is.null(enDb)) {
    ah <- AnnotationHub::AnnotationHub(ask = FALSE)
    if (species == "Human") {
      # query(ah, pattern = c("Homo sapiens", "EnsDb"))
      enDb <- ah[["AH109606"]]
    } else if (species == "Mouse") {
      # query(ah, pattern = c("Mus musculus", "EnsDb"))
      enDb <- ah[["AH89211"]]
    } else {
      stop("The species is not supported.")
    }
  }

  if (is.data.frame(peak)) {
    peak <- GenomicRanges::makeGRangesFromDataFrame(
      peak,
      seqnames.field = colnames(peak)[1],
      start.field = colnames(peak)[2],
      end.field = colnames(peak)[3],
      strand.field = colnames(peak)[6],
      starts.in.df.are.0based = TRUE
    )
  }
  GenomeInfoDb::seqlevelsStyle(peak) <- "NCBI"

  # Getting regions
  miRNA <- unlist(ensembldb::exonsBy(enDb, by = "gene")[
    ensembldb::genes(enDb)$gene_id[which(ensembldb::genes(enDb)$gene_biotype == "miRNA")]
  ])
  Coding_Exons <- unlist(ensembldb::cdsBy(enDb, by = "gene"))
  Introns <- unlist(ensembldb::intronsByTranscript(enDb))
  Utr5 <- unlist(ensembldb::fiveUTRsByTranscript(enDb))
  Utr3 <- unlist(ensembldb::threeUTRsByTranscript(enDb))
  NCexons <- unlist(ensembldb::exonsBy(enDb, by = "gene")[
    ensembldb::genes(enDb)$gene_id[which(ensembldb::genes(enDb)$gene_biotype == "lncRNA")]
  ])

  # Getting Splice Site annotation
  all_cds <- ensembldb::cdsBy(enDb, by = "gene")
  all_ncexons <- ensembldb::exonsBy(enDb, by = "gene")[
    ensembldb::genes(enDb)$gene_id[which(ensembldb::genes(enDb)$gene_biotype == "lncRNA")]
  ]

  all_cds <- all_cds[-which(lapply(all_cds, length) == 1)]
  all_ncexons <- all_ncexons[-which(lapply(all_ncexons, length) == 1)]

  # Getting 5' Splice site
  all_cds_fivess <- unlist(lapply(all_cds, Get_Junctions_five))
  all_ncexons_fivess <- unlist(lapply(all_ncexons, Get_Junctions_five))

  # Getting 3' Splice site
  all_cds_threes <- lapply(all_cds, Get_Junctions_three)
  all_ncexons_threes <- lapply(all_ncexons, Get_Junctions_three)

  # Converting them to GRanges
  all_cds_fivess <- unlist(GenomicRanges::GRangesList(all_cds_fivess))
  all_ncexons_fivess <- unlist(GenomicRanges::GRangesList(all_ncexons_fivess))
  all_cds_threes <- unlist(GenomicRanges::GRangesList(all_cds_threes))
  all_ncexons_threes <- unlist(GenomicRanges::GRangesList(all_ncexons_threes))

  splicesite5 <- c(all_cds_fivess, all_ncexons_fivess)
  splicesite3 <- c(all_cds_threes, all_ncexons_threes)

  # Getting the peaks
  # miRNA
  miRNA_peaks <- peak[Get_overlaps(peak = peak, anot = miRNA)]
  miRNA_peaks$Region <- "miRNA"
  peak <- peak[-Get_overlaps(peak = peak, anot = miRNA)]

  Coding_Exons_peaks <- peak[Get_overlaps(peak = peak, anot = Coding_Exons)]
  Coding_Exons_peaks$Region <- "CDS"
  peak <- peak[-Get_overlaps(peak = peak, anot = Coding_Exons)]

  Utr5_peaks <- peak[Get_overlaps(peak = peak, anot = Utr5)]
  Utr5_peaks$Region <- "UTR5"
  peak <- peak[-Get_overlaps(peak = peak, anot = Utr5_peaks)]

  Utr3_peaks <- peak[Get_overlaps(peak = peak, anot = Utr3)]
  Utr3_peaks$Region <- "UTR3"
  peak <- peak[-Get_overlaps(peak = peak, anot = Utr3)]

  NCExons_peaks <- peak[Get_overlaps(peak = peak, anot = NCexons)]
  NCExons_peaks$Region <- "NCExon"
  peak <- peak[-Get_overlaps(peak = peak, anot = NCexons)]

  ss3_peaks <- peak[Get_overlaps(peak = peak, anot = splicesite3)]
  ss3_peaks$Region <- "SpliceSite3"
  peak <- peak[-Get_overlaps(peak = peak, anot = splicesite3)]

  ss5_peaks <- peak[Get_overlaps(peak = peak, anot = splicesite5)]
  ss5_peaks$Region <- "SpliceSite5"
  peak <- peak[-Get_overlaps(peak = peak, anot = splicesite5)]

  Introns_peaks <- peak[Get_overlaps(peak = peak, anot = Introns)]
  Introns_peaks$Region <- "Intron"
  peak <- peak[-Get_overlaps(peak = peak, anot = Introns)]

  peak$Region <- "Others"

  comb <- c(
    miRNA_peaks, Coding_Exons_peaks, Utr5_peaks, Utr3_peaks,
    NCExons_peaks, ss3_peaks, ss5_peaks, Introns_peaks, peak
  )
  return(data.frame(comb))
}

# Internal helper: Get overlapping peak indices
Get_overlaps <- function(peak, anot) {
  ov <- GenomicRanges::findOverlaps(peak, anot)
  overlaps <- GenomicRanges::pintersect(
    peak[S4Vectors::queryHits(ov)],
    anot[S4Vectors::subjectHits(ov)]
  )
  percentOverlap <- BiocGenerics::width(overlaps) / BiocGenerics::width(peak[S4Vectors::queryHits(ov)])
  m <- S4Vectors::queryHits(ov)[which(percentOverlap > 0.5)]
  return(unique(m))
}


# Internal helper: Get 5' splice junctions
Get_Junctions_five <- function(x) {
  if (length(x) == 1) {
    return(NULL)
  }
  strand <- BiocGenerics::strand(x)@values[1]
  if (strand == "+") {
    x_resize <- GenomicRanges::flank(x[1:(length(x) - 1)], width = 50, start = FALSE, both = TRUE)
  } else if (strand == "-") {
    x_resize <- GenomicRanges::flank(x[-1], width = 50, start = FALSE, both = TRUE)
  }
  return(x_resize)
}


# Internal helper: Get 3' splice junctions
Get_Junctions_three <- function(x) {
  if (length(x) == 1) {
    return(NULL)
  }
  strand <- BiocGenerics::strand(x)@values[1]
  if (strand == "+") {
    x_resize <- GenomicRanges::flank(x[-1], width = 50, start = TRUE, both = TRUE)
  } else if (strand == "-") {
    x_resize <- GenomicRanges::flank(x[1:(length(x) - 1)], width = 50, start = TRUE, both = TRUE)
  }
  return(x_resize)
}
