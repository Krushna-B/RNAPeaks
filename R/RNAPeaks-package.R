#' RNAPeaks: Visualize RNA-Binding Protein Peaks on Gene Structures
#'
#' @description
#' RNAPeaks creates publication-quality visualizations of RNA-binding protein
#' (RBP) peaks overlaid on gene structures. It supports single-gene and
#' multi-gene region plots, optional BAM coverage tracks, splicing maps, and
#' sequence motif analysis around splice junctions.
#' @section Plotting:
#' \describe{
#'   \item{[plot_gene()]}{RBP peaks over a single transcript.}
#'   \item{[plot_region()]}{RBP peaks across every gene overlapping a
#'     genomic window.}
#'   \item{[peaks_options()]}{BED filtering / ordering options.}
#'   \item{[peaks_plot_style()]}{Visual settings (colors, sizes, BAM track
#'     layout, highlight band, junction lines).}
#' }
#'
#' @section Splicing analysis:
#' \describe{
#'   \item{[createSplicingMap()]}{RBP binding frequency around skipped-exon
#'     splice junctions with bootstrap significance testing.}
#'   \item{[createSequenceMap()]}{Position-specific motif enrichment around
#'     skipped-exon splice junctions.}
#'   \item{[createRetainedIntronSplicingMap()]}{Binding frequency around
#'     retained-intron junctions.}
#'   \item{[createRetainedIntronSequenceMap()]}{Motif enrichment around
#'     retained-intron junctions.}
#' }
#'
#' @section Helpers:
#' \describe{
#'   \item{[check_bed()]}{Validate and normalize a BED data frame.}
#' }
#'
#' @section Bundled data:
#' \describe{
#'   \item{[gtf_hg38], [gtf_mm10], [gtf_mm39]}{GENCODE basic annotations.}
#'   \item{[K562_bed], [HepG2_bed]}{ENCODE eCLIP peak calls for testing.}
#'   \item{[sample_se.mats]}{rMATS skipped-exon events for splicing analysis.}
#' }
#'
#' @examples
#' \dontrun{
#' library(RNAPeaks)
#'
#' # Single-gene plot using bundled annotation + bundled peaks
#' plot_gene(bed = K562_bed, gene = "GAPDH", species = "hg38")
#'
#' # Multi-gene region plot
#' plot_region(
#'   bed    = K562_bed,
#'   chr    = "12", start = 56000000, end = 56050000, strand = "+",
#'   gtf    = gtf_hg38
#' )
#'
#' # With BAM coverage tracks above the gene
#' plot_gene(
#'   bed       = K562_bed,
#'   gene      = "GAPDH",
#'   gtf       = gtf_hg38,
#'   bam_files = c("Control" = "ctrl.bam", "Treated" = "treat.bam")
#' )
#'
#' # Splicing + sequence analysis
#' createSplicingMap(bed_file = K562_bed, SEMATS = sample_se.mats)
#' createSequenceMap(SEMATS = sample_se.mats, sequence = "CCCC")
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
  "pos", "y_bot", "y_top", "track", "text",
  # dplyr variables
  "global_position", "bin", "frequency", "event_id", "bin_index", "position",
  "PValue", "FDR", "IncLevelDifference", "BEDFILE",
  # Splicing/Sequence map variables
  "moving_avg_sd", "control_sd", "grp_freq", "control_mean", "z_score", "p_adjusted",
  "Inc_1", "ribbon_fill", "start_pos", "end_pos", "max_y",
  "schematic_start", "schematic_end", "bar_y", "overlap_count", "match_count",
  "position_in_bin", "position_in_region", "region_idx", "frequency_sd",
  "n_events", "pvalue_adj", "significant", "plot_group",
  # Other
  "bed_df", "MASTER_FILE"
))
