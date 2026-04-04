# Gene-Level Visualization

## Overview

This article covers
[`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)
and
[`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)
in depth — the two functions for producing peak-on-gene-structure
diagrams.

``` r
library(RNAPeaks)
data(gtf_human)
```

------------------------------------------------------------------------

## PlotGene()

### Basic usage

``` r
result <- PlotGene(
    bed    = sample_bed,
    geneID = "GAPDH",
    gtf    = gtf_human
)
result$plot
```

### Selecting a transcript isoform

By default, RNAPeaks selects the longest transcript. To target a
specific isoform, supply its Ensembl transcript ID:

``` r
result <- PlotGene(
    bed    = sample_bed,
    geneID = "TP53",
    gtf    = gtf_human,
    TxID   = "ENST00000269305"
)
result$plot
```

### Ordering protein tracks

The `order_by` argument controls the vertical stacking order of RBP
tracks:

``` r
# Order by peak count (most peaks at the bottom — default)
PlotGene(bed = sample_bed, geneID = "GAPDH", gtf = gtf_human, order_by = "Count")

# Alphabetical order
PlotGene(bed = sample_bed, geneID = "GAPDH", gtf = gtf_human, order_by = "Target")

# Explicit custom order
PlotGene(
    bed      = sample_bed,
    geneID   = "GAPDH",
    gtf      = gtf_human,
    order_in = c("PROTEIN_A", "PROTEIN_B", "PROTEIN_C")
)
```

### Merging nearby peaks

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

### 5′ → 3′ orientation

For negative-strand genes, the x-axis runs right-to-left by default
(genomic coordinate order). Set `five_to_three = TRUE` to flip the axis
so 5′ is always on the left:

``` r
PlotGene(
    bed           = sample_bed,
    geneID        = "ACTB",
    gtf           = gtf_human,
    five_to_three = TRUE
)
```

### Highlighting a region

Draw a shaded rectangle over any genomic interval to call attention to a
region of interest:

``` r
PlotGene(
    bed                       = sample_bed,
    geneID                    = "GAPDH",
    gtf                       = gtf_human,
    highlighted_region_start  = 6643000,
    highlighted_region_stop   = 6643500,
    highlighted_region_color  = "gold",
    highlighted_region_opacity = 0.4
)
```

### Exon-intron junction lines

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

## PlotRegion()

[`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)
renders all genes (using their longest transcript) whose features
overlap the specified genomic window, stacked vertically.

### Basic usage

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

### Focusing on a single gene within a region

Supply `geneID` to include only one gene while still using the specified
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

## Styling reference

Both
[`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)
and
[`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)
accept a shared set of styling arguments via `...`. The most commonly
used are listed below.

### Gene structure

| Parameter          | Default       | Description           |
|--------------------|---------------|-----------------------|
| `exon_fill`        | `"navy"`      | Exon / CDS fill color |
| `utr_fill`         | `"lightgray"` | UTR fill color        |
| `intron_color`     | `"gray60"`    | Intron line color     |
| `intron_linewidth` | `0.9`         | Intron line width     |

### Peaks

| Parameter           | Default    | Description                      |
|---------------------|------------|----------------------------------|
| `peak_col`          | `"purple"` | Peak rectangle fill color        |
| `peak_alpha`        | `0.95`     | Peak opacity                     |
| `peak_border_color` | `NA`       | Peak border color (`NA` = none)  |
| `peaks_width`       | `0.3`      | Vertical height of each peak row |

### Labels

| Parameter     | Default   | Description             |
|---------------|-----------|-------------------------|
| `label_size`  | `5`       | Protein label font size |
| `label_color` | `"black"` | Protein label color     |
| `title_size`  | `25`      | Plot title font size    |

### Layout

| Parameter       | Default | Description                                 |
|-----------------|---------|---------------------------------------------|
| `axis_pad_bp`   | `500`   | Padding (bp) added to each side of the plot |
| `axis_breaks_n` | `5`     | Number of x-axis tick marks                 |
| `max_proteins`  | `40`    | Maximum number of protein tracks to show    |
| `x_lims`        | `NULL`  | Custom x-axis limits as `c(min, max)`       |

### Background bands

| Parameter        | Default     | Description               |
|------------------|-------------|---------------------------|
| `band_even_fill` | `"#F7F8FA"` | Even-row band fill        |
| `band_odd_fill`  | `"#FFFFFF"` | Odd-row band fill         |
| `band_sep_color` | `"#E5E7EB"` | Band separator line color |

------------------------------------------------------------------------

## Saving output

Both functions write a PDF and a CSV by default. Set paths explicitly or
pass `NULL` to suppress file output:

``` r
result <- PlotGene(
    bed                = sample_bed,
    geneID             = "GAPDH",
    gtf                = gtf_human,
    RNA_Peaks_File_Path = "figures/GAPDH_peaks.pdf",
    Bed_File_Path       = "data/GAPDH_peaks.csv"
)

# Suppress file output — work with the ggplot object directly
result <- PlotGene(
    bed                = sample_bed,
    geneID             = "GAPDH",
    gtf                = gtf_human,
    RNA_Peaks_File_Path = NULL,
    Bed_File_Path       = NULL
)

# Save manually with custom dimensions
ggplot2::ggsave("GAPDH.pdf", result$plot, width = 18, height = 10)
```

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
#> [46] magrittr_2.0.4       filelock_1.0.3       UCSC.utils_1.6.1    
#> [49] scales_1.4.0         rappdirs_0.3.4       bit64_4.6.0-1       
#> [52] rmarkdown_2.31       XVector_0.50.0       httr_1.4.8          
#> [55] bit_4.6.0            ragg_1.5.2           png_0.1-9           
#> [58] memoise_2.0.1        evaluate_1.0.5       knitr_1.51          
#> [61] GenomicRanges_1.62.1 IRanges_2.44.0       BiocFileCache_3.0.0 
#> [64] rlang_1.1.7          glue_1.8.0           DBI_1.3.0           
#> [67] BiocManager_1.30.27  BiocGenerics_0.56.0  jsonlite_2.0.0      
#> [70] R6_2.6.1             systemfonts_1.3.2    fs_2.0.1
```
