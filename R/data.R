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
#'   gtf <- LoadGTF("human")
#'   PlotGene(bed = sample_bed,geneID ="GAPDH", gtf=gtf)
#' }
"sample_bed"
