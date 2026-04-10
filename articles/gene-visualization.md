# PlotGene

## Overview

[`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)
draws all RBP peaks that overlap a single gene of interest, overlaid on
the annotated transcript structure — exons, UTRs, introns, and CDS — at
accurate genomic scale.

``` r
library(RNAPeaks)
data(gtf_human)
```

------------------------------------------------------------------------

## Basic usage

``` r
result <- PlotGene(
    bed    = sample_bed,
    geneID = "GAPDH",
    gtf    = gtf_human
)
result$plot   # ggplot2 object
result$csv    # data frame of peaks used in the figure
```

------------------------------------------------------------------------

## Selecting a transcript isoform

By default the longest transcript is used. To target a specific isoform,
supply its Ensembl transcript ID:

``` r
result <- PlotGene(
    bed    = sample_bed,
    geneID = "TP53",
    gtf    = gtf_human,
    TxID   = "ENST00000269305"
)
```

Alternatively, pass an Ensembl gene ID directly:

``` r
PlotGene(bed = sample_bed, geneID = "ENSG00000141510", gtf = gtf_human)
```

------------------------------------------------------------------------

## Ordering protein tracks

The `order_by` argument controls the vertical stacking order of RBP
tracks:

``` r
# Order by peak count — most peaks at bottom (default)
PlotGene(bed = sample_bed, geneID = "GAPDH", gtf = gtf_human, order_by = "Count")

# Alphabetical by protein name
PlotGene(bed = sample_bed, geneID = "GAPDH", gtf = gtf_human, order_by = "Target")

# Explicit custom order
PlotGene(
    bed      = sample_bed,
    geneID   = "GAPDH",
    gtf      = gtf_human,
    order_in = c("PROTEIN_A", "PROTEIN_B", "PROTEIN_C")
)
```

------------------------------------------------------------------------

## Merging nearby peaks

Use `merge` to collapse peaks within a given distance (bp) into a single
rectangle, reducing visual noise from fragmented binding sites:

``` r
PlotGene(
    bed    = sample_bed,
    geneID = "GAPDH",
    gtf    = gtf_human,
    merge  = 50     # merge peaks within 50 bp
)
```

------------------------------------------------------------------------

## 5′ → 3′ orientation

For negative-strand genes the x-axis runs right-to-left by default
(genomic coordinate order). Set `five_to_three = TRUE` to flip so 5′ is
always on the left:

``` r
PlotGene(
    bed           = sample_bed,
    geneID        = "ACTB",
    gtf           = gtf_human,
    five_to_three = TRUE
)
```

------------------------------------------------------------------------

## Highlighting a region

Draw a shaded rectangle over any genomic interval to call attention to a
region of interest:

``` r
PlotGene(
    bed                        = sample_bed,
    geneID                     = "GAPDH",
    gtf                        = gtf_human,
    highlighted_region_start   = 6643000,
    highlighted_region_stop    = 6643500,
    highlighted_region_color   = "gold",
    highlighted_region_opacity = 0.4
)
```

------------------------------------------------------------------------

## Exon-intron junction lines

Add vertical dashed lines at every exon-intron boundary:

``` r
PlotGene(
    bed            = sample_bed,
    geneID         = "GAPDH",
    gtf            = gtf_human,
    show_junctions = TRUE,
    junction_color = "steelblue"
)
```

------------------------------------------------------------------------

## Saving output

Both a PDF and a CSV are written by default. Set paths explicitly or
pass `NULL` to suppress file output:

``` r
result <- PlotGene(
    bed                 = sample_bed,
    geneID              = "GAPDH",
    gtf                 = gtf_human,
    RNA_Peaks_File_Path = "figures/GAPDH_peaks.pdf",
    Bed_File_Path       = "data/GAPDH_peaks.csv"
)

# Suppress file output
result <- PlotGene(
    bed                 = sample_bed,
    geneID              = "GAPDH",
    gtf                 = gtf_human,
    RNA_Peaks_File_Path = NULL,
    Bed_File_Path       = NULL
)

# Save manually with custom dimensions
ggplot2::ggsave("GAPDH.pdf", result$plot, width = 18, height = 10)
```
