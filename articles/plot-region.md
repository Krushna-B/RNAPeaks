# PlotRegion

## Overview

[`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)
renders all genes whose features overlap a specified genomic window,
stacking each gene’s transcript structure and overlapping RBP peaks
vertically. It is useful for visualizing binding patterns across a locus
or comparing multiple genes in a region simultaneously.

``` r
library(RNAPeaks)
data(gtf_human)
```

------------------------------------------------------------------------

## Basic usage

``` r
result <- PlotRegion(
    bed    = sample_bed,
    gtf    = gtf_human,
    Chr    = "12",
    Start  = 56000000,
    End    = 56050000,
    Strand = "+"
)
result$plot
```

------------------------------------------------------------------------

## Focusing on a single gene within a region

Supply `geneID` to include only one gene while still using the
coordinate window for peak filtering:

``` r
PlotRegion(
    bed    = sample_bed,
    gtf    = gtf_human,
    Chr    = "12",
    Start  = 56000000,
    End    = 56050000,
    Strand = "+",
    geneID = "GAPDH"
)
```

------------------------------------------------------------------------

## Saving output

``` r
result <- PlotRegion(
    bed                 = sample_bed,
    gtf                 = gtf_human,
    Chr                 = "12",
    Start               = 56000000,
    End                 = 56050000,
    Strand              = "+",
    RNA_Peaks_File_Path = "figures/region_peaks.pdf",
    Bed_File_Path       = "data/region_peaks.csv"
)

# Suppress file output
PlotRegion(
    bed                 = sample_bed,
    gtf                 = gtf_human,
    Chr                 = "12",
    Start               = 56000000,
    End                 = 56050000,
    Strand              = "+",
    RNA_Peaks_File_Path = NULL,
    Bed_File_Path       = NULL
)
```
