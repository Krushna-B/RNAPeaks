# Splicing Maps and Sequence Motif Analysis

## Overview

[`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md)
and
[`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)
quantify signal — protein binding or sequence motif frequency — at every
position across a four-region window that spans the exon-intron
junctions of skipped-exon (SE) alternative splicing events.

Both functions accept SE.MATS output from
[rMATS](https://rnaseq-mats.sourceforge.net/) and classify events into
three groups for comparison:

- **Retained** — exon inclusion significantly increases in condition 1
- **Excluded** — exon skipping significantly increases in condition 1
- **Control** — events with stable inclusion levels across conditions

``` r
library(RNAPeaks)
```

------------------------------------------------------------------------

## The four-region window

Each splicing event is decomposed into four equal-length windows of
`(WidthIntoExon + WidthIntoIntron)` positions:

     Upstream exon        Intron 1       Skipped exon      Intron 2     Downstream exon
    |──────|──────────────|──────|────────────────|──────|──────────────|──────|
      UE    ←Region 1→    UI3    ←──Region 2──→   EX5    ←─Region 3─→   DI3    ←Region 4→
           (UE-UI5)              (UI3-EX3)               (EX5-DI5)             (DI3-DE)

The x-axis of the output plot shows positions within this linearized
window, with vertical dashed lines at region boundaries and a schematic
of the skipped exon drawn beneath the x-axis.

------------------------------------------------------------------------

## Splicing map

### Basic usage

``` r
createSplicingMap(
    bed_file = sample_bed,
    SEMATS   = sample_se.mats
)
```

### Event classification thresholds

Control which events fall into each group:

``` r
createSplicingMap(
    bed_file                        = sample_bed,
    SEMATS                          = sample_se.mats,
    p_valueRetainedAndExclusion     = 0.05,   # significance threshold for R/E
    p_valueControls                 = 0.95,   # min p-value for controls
    retained_IncLevelDifference     = 0.1,    # min ΔPSI for Retained
    exclusion_IncLevelDifference    = -0.1,   # max ΔPSI for Excluded
    Min_Count                       = 50,     # min junction read count
    use_fdr                         = TRUE    # use FDR instead of raw p-value
)
```

### Window size

Adjust how far into exons and introns the analysis extends:

``` r
createSplicingMap(
    bed_file        = sample_bed,
    SEMATS          = sample_se.mats,
    WidthIntoExon   = 75,    # 75 bp into each flanking exon
    WidthIntoIntron = 400    # 400 bp into each flanking intron
)
```

### Smoothing

A moving average is applied to reduce positional noise. Adjust the
window size or disable smoothing entirely:

``` r
# Wider smoothing window
createSplicingMap(bed_file = sample_bed, SEMATS = sample_se.mats, moving_average = 100)

# No smoothing
createSplicingMap(bed_file = sample_bed, SEMATS = sample_se.mats, moving_average = NULL)
```

### Parallel processing

For large datasets, use multiple cores to speed up the binding frequency
calculation:

``` r
createSplicingMap(
    bed_file = sample_bed,
    SEMATS   = sample_se.mats,
    cores    = 4
)
```

### Bootstrap control and significance testing

The control distribution is built by bootstrapping: at each of
`control_iterations` replicates, a random sample of
`(n_retained + n_excluded) × control_multiplier` control events is drawn
and their binding frequency computed. The mean ± SD across replicates
forms the control ribbon.

Significance is then assessed at each position using a z-score:

    z = (observed_frequency − control_mean) / control_SD

Contiguous runs of ≥ `min_consecutive` significant positions are
reported as enriched regions and drawn as colored bars above the plot.

``` r
createSplicingMap(
    bed_file           = sample_bed,
    SEMATS             = sample_se.mats,
    control_multiplier = 2.0,   # sample_size = (R + E) * 2
    control_iterations = 30,    # more iterations = tighter control band
    z_threshold        = 1.96,  # p < 0.05 (used when use_fdr = FALSE)
    min_consecutive    = 10,    # minimum run length to call a region
    one_sided          = TRUE,  # test enrichment only (not depletion)
    show_significance  = TRUE
)
```

### Validating bootstrap normality

The z-score test assumes that bootstrap replicates are approximately
normally distributed. Verify this with `return_diagnostics = TRUE`:

``` r
diag <- createSplicingMap(
    bed_file           = sample_bed,
    SEMATS             = sample_se.mats,
    return_diagnostics = TRUE
)

# Fraction of positions where bootstrap distribution passes Shapiro-Wilk
mean(
  apply(diag$bootstrap_matrix, 1, function(x) shapiro.test(x)$p.value) > 0.05
)
```

### Returning data

``` r
freq_df <- createSplicingMap(
    bed_file    = sample_bed,
    SEMATS      = sample_se.mats,
    return_data = TRUE
)
head(freq_df)
```

------------------------------------------------------------------------

## Sequence motif analysis

[`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)
accepts the same SE.MATS input and parameters as
[`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md),
but instead of protein binding sites it measures the per-position
frequency of a target nucleotide motif in the genomic sequence.

### Basic usage

``` r
library(BSgenome.Hsapiens.UCSC.hg38)

createSequenceMap(
    SEMATS   = sample_se.mats,
    sequence = "YCAY"    # NOVA binding motif; Y = C or T
)
```

### IUPAC ambiguity codes

Any standard IUPAC ambiguity code is supported:

| Code | Matches    |
|------|------------|
| `R`  | A, G       |
| `Y`  | C, T       |
| `S`  | G, C       |
| `W`  | A, T       |
| `K`  | G, T       |
| `M`  | A, C       |
| `B`  | C, G, T    |
| `D`  | A, G, T    |
| `H`  | A, C, T    |
| `V`  | A, C, G    |
| `N`  | A, C, G, T |

``` r
# YGCY motif (MBNL binding site)
createSequenceMap(SEMATS = sample_se.mats, sequence = "YGCY")

# Poly-C run
createSequenceMap(SEMATS = sample_se.mats, sequence = "CCCC")
```

### Supplying a custom genome

By default
[`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)
uses `BSgenome.Hsapiens.UCSC.hg38`. Supply any BSgenome object to use a
different assembly:

``` r
library(BSgenome.Mmusculus.UCSC.mm10)

createSequenceMap(
    SEMATS   = sample_se.mats,
    sequence = "YCAY",
    genome   = BSgenome.Mmusculus.UCSC.mm10
)
```

------------------------------------------------------------------------

## Appearance options

Both functions share a common set of visual parameters:

| Parameter         | Default       | Description                        |
|-------------------|---------------|------------------------------------|
| `retained_col`    | `"blue"`      | Retained group line color          |
| `excluded_col`    | `"red"`       | Excluded group line color          |
| `control_col`     | `"black"`     | Control group line color           |
| `line_width`      | `0.8`         | Line width                         |
| `ribbon_alpha`    | `0.3`         | Opacity of control SD ribbon       |
| `boundary_col`    | `"gray70"`    | Region boundary line color         |
| `exon_col`        | `"navy"`      | Skipped exon fill in the schematic |
| `legend_position` | `"bottom"`    | Legend position                    |
| `ylab`            | `"Frequency"` | Y-axis label                       |
| `title`           | `""`          | Plot title                         |

------------------------------------------------------------------------

## Session info

``` r
sessionInfo()
#> R version 4.5.3 (2026-03-11)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.4 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] RNAPeaks_1.2.0
#> 
#> loaded via a namespace (and not attached):
#>  [1] KEGGREST_1.50.0      gtable_0.3.6         xfun_0.57           
#>  [4] bslib_0.10.0         ggplot2_4.0.2        httr2_1.2.2         
#>  [7] Biobase_2.70.0       vctrs_0.7.2          tools_4.5.3         
#> [10] generics_0.1.4       stats4_4.5.3         curl_7.0.0          
#> [13] tibble_3.3.1         AnnotationDbi_1.72.0 RSQLite_2.4.6       
#> [16] blob_1.3.0           pkgconfig_2.0.3      RColorBrewer_1.1-3  
#> [19] dbplyr_2.5.2         S7_0.2.1             desc_1.4.3          
#> [22] S4Vectors_0.48.0     lifecycle_1.0.5      farver_2.1.2        
#> [25] compiler_4.5.3       textshaping_1.0.5    Biostrings_2.78.0   
#> [28] Seqinfo_1.0.0        GenomeInfoDb_1.46.2  htmltools_0.5.9     
#> [31] sass_0.4.10          yaml_2.3.12          pillar_1.11.1       
#> [34] pkgdown_2.2.0        crayon_1.5.3         jquerylib_0.1.4     
#> [37] cachem_1.1.0         AnnotationHub_4.0.0  tidyselect_1.2.1    
#> [40] digest_0.6.39        dplyr_1.2.1          BiocVersion_3.22.0  
#> [43] fastmap_1.2.0        grid_4.5.3           cli_3.6.5           
#> [46] magrittr_2.0.5       filelock_1.0.3       UCSC.utils_1.6.1    
#> [49] scales_1.4.0         rappdirs_0.3.4       bit64_4.6.0-1       
#> [52] rmarkdown_2.31       XVector_0.50.0       httr_1.4.8          
#> [55] bit_4.6.0            ragg_1.5.2           png_0.1-9           
#> [58] memoise_2.0.1        evaluate_1.0.5       knitr_1.51          
#> [61] GenomicRanges_1.62.1 IRanges_2.44.0       BiocFileCache_3.0.0 
#> [64] rlang_1.1.7          glue_1.8.0           DBI_1.3.0           
#> [67] BiocManager_1.30.27  BiocGenerics_0.56.0  jsonlite_2.0.0      
#> [70] R6_2.6.1             systemfonts_1.3.2    fs_2.0.1
```
