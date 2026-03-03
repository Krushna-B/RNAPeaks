#' RNAPeaks: Visualize RNA-Binding Protein Peaks on Gene Structures
#'
#' @description
#' RNAPeaks creates publication-quality visualizations of RNA-binding protein
#' (RBP) peaks overlaid on gene structures. It supports single-gene plots,
#' multi-gene genomic region plots, splicing maps, and sequence motif analysis
#' around splice junctions. Integrates with GTF annotations from Ensembl via
#' AnnotationHub.
#'
#' @section Main Functions:
#' The primary functions for generating plots are:
#'
#' \describe{
#'   \item{\code{\link{PlotGene}}}{Plot RBP peaks on a single gene structure.
#'     Use this when you want to visualize peaks on one specific gene.}
#'   \item{\code{\link{PlotRegion}}}{Plot RBP peaks across a genomic region
#'     containing multiple genes. Use this for broader regional views.}
#'   \item{\code{\link{createSplicingMap}}}{Analyze protein binding frequency
#'     around splice junctions. Shows where RBPs bind relative to exon/intron
#'     boundaries in retained vs excluded splicing events.}
#'   \item{\code{\link{createSequenceMap}}}{Analyze sequence motif frequency
#'     around splice junctions. Identifies position-specific enrichment of
#'     motifs in different splicing event categories.}
#' }
#'
#' @section Included Data:
#' The package includes sample data for testing:
#' \describe{
#'   \item{\code{\link{sample_bed}}}{K562 cell line RBP binding peaks, ready to use}
#'   \item{\code{\link{sample_se.mats}}}{Sample SE.MATS output for splicing analysis}
#' }
#'
#' @section Workflow:
#' A typical workflow involves:
#' \enumerate{
#'   \item Use the included \code{sample_bed} data, OR load your own BED file
#'     and validate it with \code{\link{checkBed}}
#'   \item Load GTF annotation once with \code{\link{LoadGTF}} and store locally
#'     (e.g., \code{gtf <- LoadGTF("Human")}; save with \code{saveRDS(gtf, "gtf.rds")}
#'     for future sessions)
#'   \item Call \code{\link{PlotGene}} or \code{\link{PlotRegion}}, passing the
#'     stored \code{gtf} object
#' }
#'
#' For splicing analysis:
#' \enumerate{
#'   \item Prepare SE.MATS output from rMATS or use \code{sample_se.mats}
#'   \item Call \code{\link{createSplicingMap}} with BED and SE.MATS data
#'   \item Call \code{\link{createSequenceMap}} to analyze sequence motifs
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
#' library(RNAPeaks)
#'
#' # Load GTF once and store locally (do this once per session)
#' gtf <- LoadGTF(species = "Human")
#'
#' # Optionally save for future sessions to avoid re-downloading
#' saveRDS(gtf, "human_gtf.rds")
#' # In future sessions: gtf <- readRDS("human_gtf.rds")
#'
#' # ----- Gene-level visualization -----
#' result <- PlotGene(
#'   bed = sample_bed,
#'   geneID = "GAPDH",
#'   gtf = gtf
#' )
#' result$plot
#'
#' # ----- Region-level visualization -----
#' PlotRegion(
#'   bed = sample_bed,
#'   Chr = "12",
#'   Start = 56000000,
#'   End = 56050000,
#'   Strand = "+",
#'   gtf = gtf
#' )
#'
#' # ----- Splicing map analysis -----
#' createSplicingMap(
#'   bed_file = sample_bed,
#'   SEMATS = sample_se.mats
#' )
#'
#' # ----- Sequence motif analysis -----
#' library(BSgenome.Hsapiens.UCSC.hg38)
#' createSequenceMap(
#'   SEMATS = sample_se.mats,
#'   sequence = "YCAY"
#' )
#' }
#'
"_PACKAGE"

## usethis namespace: start
#' @import GenomicRanges
#' @import IRanges
#' @importFrom AnnotationHub AnnotationHub
#' @importFrom BiocGenerics strand width start
#' @importFrom GenomeInfoDb seqlevelsStyle seqlevelsStyle<- seqlengths
#' @importFrom S4Vectors DataFrame mcols queryHits subjectHits
#' @importFrom dplyr filter mutate group_by ungroup bind_rows case_when arrange
#'   select rename left_join if_else summarise group_modify
#' @importFrom ggplot2 aes element_blank element_rect element_text geom_rect
#'   geom_text ggplot ggtitle scale_x_continuous theme theme_classic
#'   geom_segment geom_line geom_vline geom_hline annotate coord_cartesian
#'   scale_fill_identity scale_linetype_identity scale_y_continuous labs
#'   theme_minimal margin element_line arrow geom_ribbon scale_color_manual
#' @importFrom grid unit textGrob gpar convertWidth grobWidth
#' @importFrom magrittr %>%
#' @importFrom scales label_comma scientific
#' @importFrom stats aggregate na.omit setNames
#' @importFrom utils write.csv read.table
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
  # Splicing/Sequence map variables
  "moving_avg_sd", "control_sd", "grp_freq", "control_mean", "z_score",
  "Inc_1", "ribbon_fill", "start_pos", "end_pos", "max_y",
  "schematic_start", "schematic_end", "bar_y", "overlap_count", "match_count",
  "position_in_bin",
  # Other
  "bed_df", "MASTER_FILE"
))
