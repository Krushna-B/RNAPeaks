# Analyzes the frequency of a target sequence motif across splicing junction regions. Compares motif frequency between Retained, Excluded, and Control splicing events to identify position-specific enrichment patterns.

Analyzes the frequency of a target sequence motif across splicing
junction regions. Compares motif frequency between Retained, Excluded,
and Control splicing events to identify position-specific enrichment
patterns.

## Usage

``` r
createSequenceMap(
  SEMATS,
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
  cores = 1,
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

- SEMATS:

  A data frame containing SE.MATS output with columns: chr, strand,
  upstreamES, upstreamEE, exonStart_0base, exonEnd, downstreamES,
  downstreamEE, GeneID, PValue, FDR, IncLevelDifference, IJC_SAMPLE_1,
  SJC_SAMPLE_1, IJC_SAMPLE_2, SJC_SAMPLE_2, IncLevel1, IncLevel2

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
  "Excluded", "Control") to process all groups. Use c("Retained",
  "Excluded") to skip the Control group (which can be large).

- control_multiplier:

  Numeric multiplier for control sample size. The number of control
  events sampled per iteration is (n_retained + n_excluded) \*
  control_multiplier. Default is 1.0.

- control_iterations:

  Integer number for sampling iterations for control sampling. The final
  control frequency is the mean across iterations, with standard
  deviation shown as a shaded band. Default is 20.

- cores:

  Integer number of cores for parallel processing. Default is 1
  (sequential). Set higher for faster processing on multi-core systems.

- z_threshold:

  Z-score threshold for significance testing. Default is 1.96
  (corresponds to p \< 0.05 two-tailed). Only used when use_fdr = FALSE.

- min_consecutive:

  Minimum number of consecutive significant positions required to form a
  significant region. Default is 10. Helps reduce false positives from
  noise.

- one_sided:

  Logical. If TRUE (default), only test for enrichment (frequency \>
  control). If FALSE, test for both enrichment and depletion.

- use_fdr:

  Logical. If TRUE, use FDR-corrected p-values (Benjamini-Hochberg) for
  significance testing. If FALSE (default), use z_threshold directly.

- fdr_threshold:

  FDR threshold for significance when use_fdr = TRUE. Default is 0.05.

- show_significance:

  Logical. If TRUE (default), displays colored bars above the plot
  indicating regions where Retained/Excluded differ significantly from
  Control based on z-test.

- return_data:

  Logical. If TRUE, returns the frequency data frame instead of a plot.
  Default is FALSE.

- return_diagnostics:

  Logical. If TRUE, returns a list containing the frequency data, raw
  bootstrap iteration results (for normality testing), and significance
  results. Useful for validating bootstrap assumptions. Default is
  FALSE.

- verbose:

  Logical. If TRUE, prints progress messages. Default is TRUE.

- progress_callback:

  Optional function to report progress. Called with two arguments:
  current iteration number and total iterations. Used by Shiny app for
  progress display. Default is NULL (no callback).

- title:

  Character string for the plot title. Default is "" (no title).

- retained_col:

  Color for the Retained group line. Default is "blue".

- excluded_col:

  Color for the Excluded group line. Default is "red".

- control_col:

  Color for the Control group line. Default is "black".

- line_width:

  Numeric line width for the frequency lines. Default is 0.8.

- line_alpha:

  Numeric alpha (opacity) for the frequency lines. Default is 0.7.

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

  Fill color for the skipped (middle) exon in the schematic. Default is
  "navy".

- legend_position:

  Position of the legend. Default is "bottom".

- ylab:

  Label for the y-axis. Default is "Frequency".

## Value

A ggplot object showing sequence frequency across the 4 regions for
Retained, Excluded, and Control groups. Significant regions (z-test vs
Control) are shown as colored bars above the plot. Returns a data frame
if return_data = TRUE.

## Details

The function divides each splicing event into 4 regions of
(WidthIntoExon + WidthIntoIntron) bp each:

- Region 1: Upstream exon end to first intron

- Region 2: First intron end to middle (skipped) exon start

- Region 3: Middle exon end to second intron

- Region 4: Second intron end to downstream exon start

Events are filtered into three groups:

- Retained: Significant events (PValue \< threshold) with positive
  inclusion

- Excluded: Significant events (PValue \< threshold) with negative
  inclusion

- Control: Non-significant events

At each position, the function checks if the target sequence starts
there. The frequency is calculated as: (events with motif at position) /
(total events)

## Examples

``` r
if (FALSE) { # \dontrun{
library(BSgenome.Hsapiens.UCSC.hg38)

# Basic usage
createSequenceMap(SEMATS = sample_se.mats, sequence = "CCCC")

# Search for YCAY motif (Y = C or T)
createSequenceMap(SEMATS = sample_se.mats, sequence = "YCAY")

# Return data instead of plot
freq_data <- createSequenceMap(SEMATS = sample_se.mats,
                                sequence = "GGGG",
                                return_data = TRUE)
} # }
```
