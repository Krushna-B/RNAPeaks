# Getting Started with RNAPeaks

## Overview

RNAPeaks is an R package for producing publication-quality
visualizations of RNA-binding protein (RBP) peaks overlaid on annotated
gene structures. Starting from a BED file of protein binding sites,
RNAPeaks draws each peak as a rectangle above the corresponding gene
model — showing exons, UTRs, and introns at accurate genomic scale.

The package provides four analysis modes:

| Mode                     | Function                                                                                     | Input             |
|--------------------------|----------------------------------------------------------------------------------------------|-------------------|
| Single-gene peak plot    | [`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)                   | BED + gene ID     |
| Genomic region peak plot | [`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)               | BED + coordinates |
| Splicing map             | [`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md) | BED + SE.MATS     |
| Sequence motif map       | [`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md) | SE.MATS + motif   |

------------------------------------------------------------------------

## Installation

``` r
# Bioconductor dependencies
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "GenomicRanges", "IRanges", "S4Vectors", "BiocGenerics",
    "GenomeInfoDb", "AnnotationHub", "Biostrings", "BSgenome"
))

# Install RNAPeaks from GitHub
devtools::install_github("Krushna-B/RNAPeaks")
```

------------------------------------------------------------------------

## Loading the package

``` r
library(RNAPeaks)
```

------------------------------------------------------------------------

## Loading GTF annotation

Gene structure diagrams are built from an Ensembl GTF annotation. Load
it once per session with
[`LoadGTF()`](https://krushna-b.github.io/RNAPeaks/reference/LoadGTF.md)
and optionally cache it to disk to avoid re-downloading on subsequent
sessions.

``` r
# First call downloads from AnnotationHub (~1–2 min)
gtf <- LoadGTF(species = "Human")

# Cache locally for instant loading next time
saveRDS(gtf, "human_gtf.rds")
```

``` r
# The package ships with a pre-loaded GTF — no download required
data(gtf_human)
gtf <- gtf_human
```

Supported species: `"Human"` (Ensembl GRCh38, AH110867) and `"Mouse"`
(Ensembl GRCm38, AH47076).

------------------------------------------------------------------------

## Preparing a BED file

Use
[`checkBed()`](https://krushna-b.github.io/RNAPeaks/reference/checkBed.md)
to validate and normalize any BED-format data frame before passing it to
a plotting function. Columns are mapped **by position**:

| Position | Name     | Default  |
|----------|----------|----------|
| 1        | `chr`    | required |
| 2        | `start`  | required |
| 3        | `end`    | required |
| 4        | `tag`    | `"peak"` |
| 5        | `score`  | `0`      |
| 6        | `strand` | `"+"`    |

``` r
my_bed <- read.table("my_peaks.bed", header = FALSE, sep = "\t")
my_bed <- checkBed(my_bed)
```

The package includes a ready-to-use sample dataset:

``` r
head(sample_bed)
#>   chr    start      end           tag score strand
#> 1  21  8401778  8401840 AARS_K562_IDR  1000      +
#> 2  19  4035770  4035869 AARS_K562_IDR  1000      +
#> 3  11 93721690 93721719 AARS_K562_IDR  1000      +
#> 4  19  2270224  2270300 AARS_K562_IDR  1000      +
#> 5   m    12167    12208 AARS_K562_IDR  1000      +
#> 6  16  2760062  2760143 AARS_K562_IDR  1000      +
```

------------------------------------------------------------------------

## Gene-level visualization

[`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)
plots all RBP peaks that overlap a single gene of interest.

``` r
result <- PlotGene(
    bed    = sample_bed,
    geneID = "GAPDH",
    gtf    = gtf
)

result$plot   # display the plot
result$csv    # filtered peaks used in the figure
```

**Key parameters:**

| Parameter       | Default   | Description                                               |
|-----------------|-----------|-----------------------------------------------------------|
| `geneID`        | —         | Gene symbol or Ensembl ID (`"ENSG..."`)                   |
| `TxID`          | `NA`      | Specific transcript isoform; defaults to longest          |
| `species`       | `"Human"` | `"Human"` or `"Mouse"`                                    |
| `order_by`      | `"Count"` | Peak track ordering: `"Count"`, `"Target"`, or `"Region"` |
| `merge`         | `0`       | Merge peaks within this many bp of each other             |
| `five_to_three` | `FALSE`   | Orient plot 5′ → 3′ regardless of strand                  |

------------------------------------------------------------------------

## Region-level visualization

[`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)
displays all genes and their overlapping peaks within a specified
genomic window.

``` r
result <- PlotRegion(
    bed    = sample_bed,
    gtf    = gtf,
    Chr    = "12",
    Start  = 56000000,
    End    = 56050000,
    Strand = "+"
)
result$plot
```

------------------------------------------------------------------------

## Splicing map analysis

[`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md)
quantifies RBP binding frequency at every position across a four-region
window spanning exon-intron junctions of skipped-exon (SE) events from
rMATS output.

``` r
createSplicingMap(
    bed_file  = sample_bed,
    SEMATS    = sample_se.mats
)
```

The function classifies events into three groups based on
`IncLevelDifference` and `PValue` / `FDR`:

- **Retained** — significantly more inclusion in condition 1
- **Excluded** — significantly more skipping in condition 1
- **Control** — non-significant, stable inclusion

A bootstrap procedure samples a matched number of control events across
`control_iterations` replicates to build a null distribution. Positions
where Retained or Excluded frequencies exceed the control z-score
threshold are highlighted as significant regions.

``` r
createSplicingMap(
    bed_file           = sample_bed,
    SEMATS             = sample_se.mats,
    cores              = 4,           # parallel processing
    control_iterations = 20,          # bootstrap replicates
    show_significance  = TRUE,        # overlay significance bars
    WidthIntoExon      = 50,          # bp into flanking exons
    WidthIntoIntron    = 300          # bp into flanking introns
)
```

------------------------------------------------------------------------

## Sequence motif analysis

[`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)
works identically to
[`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md)
but instead of protein binding, it measures the per-position frequency
of a target sequence motif. Full IUPAC ambiguity codes are supported.

``` r
library(BSgenome.Hsapiens.UCSC.hg38)

# YCAY motif (Y = C or T; classic NOVA binding site)
createSequenceMap(
    SEMATS   = sample_se.mats,
    sequence = "YCAY"
)
```

Common IUPAC codes:

| Code | Bases    |
|------|----------|
| `R`  | A or G   |
| `Y`  | C or T   |
| `S`  | G or C   |
| `W`  | A or T   |
| `N`  | any base |

------------------------------------------------------------------------

## Citation

If you use RNAPeaks in your research, please cite:

> Bhanushali K, Giri G, et al. (2025). *RNAPeaks: publication-quality
> visualization of RNA-binding protein peaks on gene structures.*
> GitHub: <https://github.com/Krushna-B/RNAPeaks>

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
