# Create Retained Intron Sequence Map

Analyzes the frequency of a target sequence motif across retained intron
junction regions. Compares motif frequency between Retained, Excluded,
and Control events to identify position-specific enrichment patterns
around the upstream exon/intron and intron/downstream exon boundaries.

## Usage

``` r
createRetainedIntronSequenceMap(
  RIMATS,
  sequence,
  genome = NULL,
  moving_average = 40,
  WidthIntoExon = 50,
  WidthIntoIntron = 250,
  p_valueRetainedAndExclusion = 0.05,
  p_valueControls = 0.95,
  retained_IncLevelDifference = 0.1,
  exclusion_IncLevelDifference = -0.1,
  Min_Count = 50,
  groups = c("Retained", "Excluded", "Control"),
  control_multiplier = 2,
  control_iterations = 20,
  z_threshold = 1.96,
  min_consecutive = 10,
  one_sided = TRUE,
  use_fdr = TRUE,
  fdr_threshold = 0.05,
  show_significance = TRUE,
  return_data = FALSE,
  return_diagnostics = FALSE,
  verbose = TRUE,
  progress_callback = NULL,
  title = "",
  retained_col = "blue",
  excluded_col = "red",
  control_col = "black",
  line_width = 0.8,
  line_alpha = 0.7,
  ribbon_alpha = 0.3,
  title_size = 20,
  title_color = "black",
  axis_text_size = 11,
  boundary_col = "gray70",
  exon_col = "navy",
  legend_position = "bottom",
  ylab = "Frequency"
)
```

## Arguments

- RIMATS:

  A data frame containing rMATS output with columns: chr, strand,
  upstreamES, upstreamEE, downstreamES, downstreamEE, GeneID, PValue,
  FDR, IncLevelDifference, IJC_SAMPLE_1, SJC_SAMPLE_1, IJC_SAMPLE_2,
  SJC_SAMPLE_2, IncLevel1, IncLevel2

- sequence:

  Character string of the target sequence motif to search for (e.g.,
  "CCCC", "YGCY"). Supports IUPAC ambiguity codes.

- genome:

  A BSgenome object. Default uses BSgenome.Hsapiens.UCSC.hg38.

- moving_average:

  Integer specifying the window size for moving average smoothing. Set
  to NULL or 0 to disable smoothing. Default is 40.

- WidthIntoExon:

  Integer specifying how many bp to extend into exons. Default is 50.

- WidthIntoIntron:

  Integer specifying how many bp to extend into introns. Default is 250.

- p_valueRetainedAndExclusion:

  P-value threshold for retained/excluded events. Default is 0.05.

- p_valueControls:

  P-value threshold for control events. Default is 0.95.

- retained_IncLevelDifference:

  Inclusion level difference threshold for retained events. Default is
  0.1.

- exclusion_IncLevelDifference:

  Inclusion level difference threshold for excluded events. Default is
  -0.1.

- Min_Count:

  Minimum read count threshold. Default is 50.

- groups:

  Character vector specifying which event groups to process. Options are
  "Retained", "Excluded", and/or "Control". Default is c("Retained",
  "Excluded", "Control") to process all groups.

- control_multiplier:

  Numeric multiplier for control sample size. Default is 2.0.

- control_iterations:

  Integer number for sampling iterations for control sampling. Default
  is 20.

- z_threshold:

  Z-score threshold for significance testing. Default is 1.96. Only used
  when use_fdr = FALSE.

- min_consecutive:

  Minimum number of consecutive significant positions required to form a
  significant region. Default is 10.

- one_sided:

  Logical. If TRUE (default), only test for enrichment.

- use_fdr:

  Logical. If TRUE, use FDR-corrected p-values. Default is TRUE.

- fdr_threshold:

  FDR threshold for significance when use_fdr = TRUE. Default is 0.05.

- show_significance:

  Logical. If TRUE (default), displays colored bars above the plot
  indicating significant regions.

- return_data:

  Logical. If TRUE, returns the frequency data frame instead of a plot.
  Default is FALSE.

- return_diagnostics:

  Logical. If TRUE, returns a list containing the frequency data and
  bootstrap diagnostics. Default is FALSE.

- verbose:

  Logical. If TRUE, prints progress messages. Default is TRUE.

- progress_callback:

  Optional function to report progress. Default is NULL.

- title:

  Character string for the plot title. Default is "".

- retained_col:

  Color for the Retained group line. Default is "blue".

- excluded_col:

  Color for the Excluded group line. Default is "red".

- control_col:

  Color for the Control group line. Default is "black".

- line_width:

  Numeric line width for the frequency lines. Default is 0.8.

- line_alpha:

  Numeric alpha for the frequency lines. Default is 0.7.

- ribbon_alpha:

  Numeric alpha for the SD ribbon around Control. Default is 0.3.

- title_size:

  Numeric font size for the plot title. Default is 20.

- title_color:

  Color for the plot title text. Default is "black".

- axis_text_size:

  Numeric font size for y-axis tick labels. Default is 11.

- boundary_col:

  Color for the dashed vertical boundary lines. Default is "gray70".

- exon_col:

  Unused parameter kept for API consistency. Default is "navy".

- legend_position:

  Position of the legend. Default is "bottom".

- ylab:

  Label for the y-axis. Default is "Frequency".

## Value

A ggplot object showing sequence motif frequency across the 2 regions
for Retained, Excluded, and Control groups. The bottom schematic shows
two exon boxes connected by a single intron line. Returns a data frame
if return_data = TRUE.

## Details

The function divides each retained intron event into 2 regions of
(WidthIntoExon + WidthIntoIntron) bp each:

- Region 1 (UE-RI5): Upstream exon end to retained intron

- Region 2 (RI3-DE): Retained intron end to downstream exon start

At each position, the function checks if the target sequence starts
there. The frequency is calculated as: (events with motif at position) /
(total events)

## Examples

``` r
if (FALSE) { # \dontrun{
library(BSgenome.Hsapiens.UCSC.hg38)
rimats <- read.table("RI.MATS.JC.txt", header = TRUE)

# Basic usage
createRetainedIntronSequenceMap(RIMATS = rimats, sequence = "CCCC")

# Search for YCAY motif (Y = C or T)
createRetainedIntronSequenceMap(RIMATS = rimats, sequence = "YCAY")

# Return data instead of plot
freq_data <- createRetainedIntronSequenceMap(RIMATS = rimats,
                                              sequence = "GGGG",
                                              return_data = TRUE)
} # }
```
