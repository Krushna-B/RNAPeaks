#' Sample K562 RBP Binding Peaks
#'
#' A dataset containing RNA-binding protein (RBP) peaks from K562 cell line
#' This data can be used to test and demonstrate the
#' plotting functions in RNAPeaks.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{chr}{Chromosome (without "chr" prefix)}
#'   \item{start}{Start position (0-based)}
#'   \item{end}{End position}
#'   \item{tag}{Protein/peak identifier}
#'   \item{score}{Peak score}
#'   \item{strand}{Strand (+ or -)}
#' }
#'
#'
#'
#' @examples
#' # Load the data
#' data(sample_bed)
#'
#' # View first few rows
#' head(sample_bed)
#'
#' # Use with PlotGene
#' \dontrun{
#'   gtf <- LoadGTF("Human")
#'   PlotGene(bed = sample_bed,geneID ="GAPDH", gtf=gtf)
#' }
"sample_bed"

#' Sample SE.MATS Output
#'
#' A dataset containing skipped exon (SE) alternative splicing events.
#' This data can be used to test and demonstrate the splicing map functions in RNAPeaks.
#'
#' @format A data frame with 87,736 observations and 23 variables:
#' \describe{
#'   \item{ID}{Event ID}
#'   \item{GeneID}{Ensembl gene ID}
#'   \item{geneSymbol}{Gene symbol}
#'   \item{chr}{Chromosome}
#'   \item{strand}{Strand (+ or -)}
#'   \item{exonStart_0base}{Skipped exon start position (0-based)}
#'   \item{exonEnd}{Skipped exon end position}
#'   \item{upstreamES}{Upstream exon start position}
#'   \item{upstreamEE}{Upstream exon end position}
#'   \item{downstreamES}{Downstream exon start position}
#'   \item{downstreamEE}{Downstream exon end position}
#'   \item{ID.1}{Duplicate ID column}
#'   \item{IJC_SAMPLE_1}{Inclusion junction counts for sample 1}
#'   \item{SJC_SAMPLE_1}{Skipping junction counts for sample 1}
#'   \item{IJC_SAMPLE_2}{Inclusion junction counts for sample 2}
#'   \item{SJC_SAMPLE_2}{Skipping junction counts for sample 2}
#'   \item{IncFormLen}{Inclusion form length}
#'   \item{SkipFormLen}{Skipping form length}
#'   \item{PValue}{P-value for differential splicing}
#'   \item{FDR}{False discovery rate adjusted p-value}
#'   \item{IncLevel1}{Inclusion levels for sample 1}
#'   \item{IncLevel2}{Inclusion levels for sample 2}
#'   \item{IncLevelDifference}{Difference in inclusion levels between samples}
#' }
#'
#'
#' @examples
#' # Load the data
#' data(sample_se.mats)
#'
#' # View first few rows
#' head(sample_se.mats)
#'
#' # Use with createSplicingMap
#' \dontrun{
#'   createSplicingMap(bed_file = sample_bed, SEMATS = sample_se.mats)
#'   createSequenceMap(SEMATS = sample_se.mats, sequence = "CCCC")
#' }
"sample_se.mats"
