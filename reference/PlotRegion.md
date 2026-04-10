# Plot RNA-Binding Protein Peaks on a Genomic Region

Creates a quality visualization of RNA-binding protein peaks overlaid on
multiple gene structures within a specified genomic region. Unlike
[`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md),
this function can display multiple genes that fall within the specified
coordinates.

## Usage

``` r
PlotRegion(
  Chr = NULL,
  Start = NULL,
  End = NULL,
  Strand = NULL,
  geneID = NULL,
  gtf = NULL,
  bed = NULL,
  species = "Human",
  TxID = NA,
  Target_col = NULL,
  omit = c(),
  order_by = "Count",
  order_in = NULL,
  merge = 0,
  peaks_width = 0.3,
  utr_col = "darkgray",
  peak_col = "Blue",
  exon_width = 0.5,
  utr_width = 0.3,
  exon_col = "navy",
  total_arrows = 12,
  max_per_intron = 5,
  five_to_three = FALSE,
  bam_files = NULL,
  bam_fill_col = "navy",
  bam_fill_alpha = 0.75,
  bam_label_size = 9,
  bam_axis_text_size = 8,
  bam_ylim = NULL,
  bam_track_height = 1,
  RNA_Peaks_File_Path = "~/Desktop/RNAPeaks.pdf",
  Bed_File_Path = "~/Desktop/BEDFILE_PEAKS.csv",
  ...
)
```

## Arguments

- Chr:

  Chromosome name (e.g., "1", "X", "chr1").

- Start:

  Start position of the genomic region (bp).

- End:

  End position of the genomic region (bp).

- Strand:

  Strand to display ("+" or "-").

- geneID:

  Optional gene identifier to focus on a specific gene.

- gtf:

  Optional pre-loaded GTF annotation data frame. If NULL, annotations
  are loaded from AnnotationHub based on species.

- bed:

  A data frame containing BED-format peak data.

- species:

  Species for annotation lookup. Either "Human" or "Mouse".

- TxID:

  Optional transcript ID to plot a specific transcript isoform.

- Target_col:

  Column name in bed containing the protein/target identifiers.

- omit:

  Character vector of target names to exclude from the plot.

- order_by:

  Method for ordering protein tracks: "Count" (default), "Target", or
  "Region".

- order_in:

  Optional character vector specifying exact order of targets.

- merge:

  Minimum gap width (bp) for merging nearby peaks.

- peaks_width:

  Vertical height of each peak track row.

- utr_col:

  Color for UTR regions.

- peak_col:

  Color for peak rectangles.

- exon_width:

  Vertical height of exon rectangles.

- utr_width:

  Vertical height of UTR rectangles.

- exon_col:

  Color for exon/CDS regions.

- total_arrows:

  Total number of directional arrows drawn across all introns to
  indicate transcription direction. Default is 12.

- max_per_intron:

  Maximum number of directional arrows drawn per intron. Default is 5.

- five_to_three:

  Logical. If TRUE and Strand is "-", flips the x-axis so 5' is on the
  left. Default FALSE.

- bam_files:

  Optional. A named character vector of BAM file paths to display as
  coverage tracks above the gene structure. Names are used as track
  labels on the left-hand side of each panel. If unnamed, the filename
  (without extension) is used as the label. BAM files must be sorted and
  indexed (a `.bai` file must exist alongside each BAM). Example:
  `c("Sample A" = "/path/to/a.bam", "Sample B" = "/path/to/b.bam")`

- bam_fill_col:

  Fill color for BAM coverage tracks. A single color applied to all
  tracks, or a character vector the same length as `bam_files` for
  per-track colours. Default `"navy"`.

- bam_fill_alpha:

  Opacity of BAM track fill. Default `0.75`.

- bam_label_size:

  Font size of the BAM track name label on the left. Default `9`.

- bam_axis_text_size:

  Font size of the 0 and max coverage values. Default `8`.

- bam_ylim:

  Optional global y-axis limits `c(min, max)` applied to all BAM tracks.
  If `NULL`, all tracks share a common scale derived from the maximum
  coverage across all BAMs.

- bam_track_height:

  Relative height of each BAM panel compared to the gene plot panel
  (which is always 4 units). Default `1`.

- RNA_Peaks_File_Path:

  File path to save the output PDF plot.

- Bed_File_Path:

  File path to save the filtered BED data as CSV.

- ...:

  Additional styling arguments passed to internal plotting functions.
  See Styling Parameters section below.

## Value

A named list containing:

- plot:

  A ggplot2 object of the peak visualization

- csv:

  The filtered BED data frame used for plotting

Access with `result$plot` and `result$csv`.

## Styling Parameters

The following parameters can be passed via `...` to customize the plot
appearance:

**Gene Structure Colors:**

- exon_fill:

  Fill color for exon/CDS regions. Default: "navy"

- utr_fill:

  Fill color for UTR regions. Default: "lightgray"

- intron_color:

  Color for intron lines. Default: "gray60"

- intron_linewidth:

  Line width for introns. Default: 0.9

- intron_arrow_len_in:

  Length of intron direction arrows in inches. Default: 0.15

**Peak Styling:**

- peak_alpha:

  Opacity of peak rectangles. Default: 0.95

- peak_border_color:

  Border color for peaks. Default: NA (no border)

- peak_border_linewidth:

  Border line width for peaks. Default: 0.4

**Background Bands:**

- band_even_fill:

  Fill color for even-numbered protein track bands. Default: "#F7F8FA"

- band_odd_fill:

  Fill color for odd-numbered protein track bands. Default: "#FFFFFF"

- band_sep_color:

  Color for band separator lines. Default: "#E5E7EB"

- band_sep_linewidth:

  Line width for band separators. Default: 0.4

**Labels:**

- label_size:

  Font size for protein labels. Default: 5

- label_color:

  Color for protein labels. Default: "black"

- strand_label_size:

  Font size for 5'/3' strand labels. Default: 5

- strand_label_color:

  Color for strand labels. Default: "black"

- protein_label_x_offset:

  Horizontal offset for protein labels in bp. Default: 100

**Title and Axes:**

- title_size:

  Font size for plot title. Default: 25

- title_color:

  Color for plot title. Default: "black"

- subtitle_size:

  Font size for plot subtitle. Default: 12

- subtitle_color:

  Color for plot subtitle. Default: "black"

- subtitle_sep:

  Separator between gene name and coordinates in subtitle. Default: ": "

- axis_title_size:

  Font size for axis titles. Default: 11

- axis_text_size:

  Font size for axis tick labels. Default: 9

**Axis and Layout:**

- x_lims:

  Custom x-axis limits as c(min, max). Default: NULL (auto)

- axis_pad_bp:

  Padding in base pairs added to each side of the plot. Default: 500

- axis_breaks_n:

  Number of axis tick breaks. Default: 5

- max_proteins:

  Maximum number of protein tracks to display. Default: 40

**Plot Margins:**

- plot_right_margin:

  Right margin in points. Default: 50

- plot_top_margin:

  Top margin in points. Default: 30

- plot_bottom_margin:

  Bottom margin in points. Default: 30

- plot_left_margin:

  Left margin in points. Default: NULL (auto)

**Region Plot Specific:**

- gene_label_x_offset:

  Horizontal offset for gene labels as fraction of axis_pad_bp. Default:
  0.25

- gene_label_size:

  Font size for gene name labels. Default: 5

- gene_label_color:

  Color for gene name labels. Default: "black"

**Highlighted Region:**

- highlighted_region_start:

  Start position of region to highlight. Default: NULL

- highlighted_region_stop:

  End position of region to highlight. Default: NULL

- highlighted_region_color:

  Color for highlighted region. Default: "pink"

- highlighted_region_opacity:

  Opacity for highlighted region. Default: 0.30

**Junction Lines:**

- show_junctions:

  Logical. If TRUE, draws vertical dashed lines at exon/intron
  boundaries. Default: FALSE

- junction_color:

  Color for junction lines. Default: "gray40"

- junction_linetype:

  Line type for junction lines. Default: "dashed"

- junction_linewidth:

  Line width for junction lines. Default: 0.4

- junction_alpha:

  Opacity for junction lines. Default: 0.7

## Examples

``` r
if (FALSE) { # \dontrun{
  # Load GTF annotation (do this once, takes time on first call)
  gtf <- LoadGTF(species = "Human")

  # ----- Using included sample data -----
  # sample_bed is included with the package and ready to use
  result <- PlotRegion(
    bed = sample_bed,
    gtf = gtf,
    Chr = "12",
    Start = 56000000,
    End = 56050000,
    Strand = "+"
  )

  # Access results
  result$plot
  result$csv

  # ----- Using your own BED file -----
  # 1. Read your BED file
  my_bed <- read.table("my_peaks.bed", header = FALSE, sep = "\t")

  # 2. Check Bed file
  my_bed <- checkBed(my_bed)

  # 3. Plot peaks on region
  result <- PlotRegion(
    bed = my_bed,
    gtf = gtf,
    Chr = "12",
    Start = 56000000,
    End = 56050000,
    Strand = "+"
  )

  # Access results
  result$plot
  result$csv
} # }
```
