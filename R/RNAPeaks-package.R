#' RNAPeaks: Visualize RNA-Binding Protein Peaks on Gene Structures
#'
#' @description
#' RNAPeaks creates publication-quality visualizations of RNA-binding protein
#' (RBP) peaks overlaid on gene structures. It supports single-gene plots,
#' multi-gene genomic region plots, and integrates with GTF annotations
#' from Ensembl via AnnotationHub.
#'
#' @section Main Functions:
#' The two primary functions for generating plots are:
#'
#' \describe{
#'   \item{\code{\link{PlotGene}}}{Plot RBP peaks on a single gene structure.
#'     Use this when you want to visualize peaks on one specific gene.}
#'   \item{\code{\link{PlotRegion}}}{Plot RBP peaks across a genomic region
#'     containing multiple genes. Use this for broader regional views.}
#' }
#'
#' @section Workflow:
#' A typical workflow involves:
#' \enumerate{
#'   \item Load your BED file with peak data
#'   \item Validate it with \code{\link{checkBed}}
#'   \item Load GTF annotation once with \code{\link{LoadGTF}} and store locally
#'     (e.g., \code{gtf <- LoadGTF("Human")}; save with \code{saveRDS(gtf, "gtf.rds")}
#'     for future sessions)
#'   \item Call \code{\link{PlotGene}} or \code{\link{PlotRegion}}, passing the
#'     stored \code{gtf} object
#' }
#'
#' @section Helper Functions:
#' \describe{
#'   \item{\code{\link{checkBed}}}{Validate and standardize BED file format}
#'   \item{\code{\link{LoadGTF}}}{Load GTF annotation from AnnotationHub. Call once
#'     and store the result in a local variable or save to disk with
#'     \code{saveRDS()} to avoid repeated downloads. Supports "Human" and "Mouse".}
#' }
#'
#' @examples
#' \dontrun{
#' # Basic workflow
#' library(RNAPeaks)
#'
#' # Load and validate BED data
#' bed <- read.table("peaks.bed", header = FALSE)
#' bed <- checkBed(bed)
#'
#' # Load GTF once and store locally (do this once per session)
#' gtf <- LoadGTF(species = "Human")
#'
#' # Optionally save for future sessions to avoid re-downloading
#' saveRDS(gtf, "human_gtf.rds")
#' # In future sessions: gtf <- readRDS("human_gtf.rds")
#'
#' # Plot peaks on a single gene
#' PlotGene(
#'   bed = bed,
#'   geneID = "TP53",
#'   gtf = gtf
#' )
#'
#' # Plot peaks across a genomic region
#' PlotRegion(
#'   bed = bed,
#'   Chr = "17",
#'   Start = 7565097,
#'   End = 7590856,
#'   Strand = "-",
#'   gtf = gtf
#' )
#' }
#'
"_PACKAGE"

## usethis namespace: start
#' @import GenomicRanges
#' @import IRanges
#' @importFrom ggplot2 aes element_blank element_rect element_text geom_rect
#'   geom_text ggplot ggtitle scale_x_continuous theme theme_classic
#'   geom_segment geom_line geom_vline geom_hline annotate coord_cartesian
#'   scale_fill_identity scale_linetype_identity scale_y_continuous labs
#'   theme_minimal margin element_line arrow
#' @importFrom dplyr filter mutate group_by ungroup bind_rows case_when arrange
#' @importFrom S4Vectors mcols queryHits subjectHits DataFrame
#' @importFrom BiocGenerics strand width
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<-
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom scales label_comma scientific
#' @importFrom grid unit textGrob gpar convertWidth grobWidth
#' @importFrom magrittr %>%
#' @importFrom stats aggregate na.omit
#' @importFrom utils write.csv
## usethis namespace: end
NULL

# Global variables for non-standard evaluation (NSE) in ggplot2 and dplyr
utils::globalVariables(c(
  # ggplot2 aes variables
  "dir_start", "dir_end", "mid_y", "y_start", "y_end", "x", "xend", "y", "yend",
  "x_start", "x_end", "group_name", "xpos", "label_y", "Label", "X", "Y",
  "label_x", "label", "xintercept", "xmin", "xmax", "ymin", "ymax", "fill",
  "linetype", "type", "schematic_position", "moving_avg", "group",
  # dplyr variables
  "global_position", "bin", "frequency", "event_id", "bin_index", "position",
  "PValue", "FDR", "IncLevelDifference", "BEDFILE",
  # Other
  "bed_df", "MASTER_FILE"
))
