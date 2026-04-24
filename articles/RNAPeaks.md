# Getting Started with RNAPeaks

## Overview

RNAPeaks is an R package for producing publication-quality
visualizations of RNA-binding protein (RBP) peaks overlaid on annotated
gene structures. Starting from a BED file of protein binding sites,
RNAPeaks draws each peak as a rectangle above the corresponding gene
model — showing exons, UTRs, and introns at accurate genomic scale.

The package provides four analysis modes:

| Mode                     | Function                                                                                                                 | Input             |
|--------------------------|--------------------------------------------------------------------------------------------------------------------------|-------------------|
| Single-gene peak plot    | [`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)                                               | BED + gene ID     |
| Genomic region peak plot | [`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)                                           | BED + coordinates |
| Splicing map (SE)        | [`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md)                             | BED + SE.MATS     |
| Splicing map (RI)        | [`createRetainedIntronSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSplicingMap.md) | BED + RI.MATS     |
| Sequence motif map (SE)  | [`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)                             | SE.MATS + motif   |
| Sequence motif map (RI)  | [`createRetainedIntronSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSequenceMap.md) | RI.MATS + motif   |

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

The package includes ready-to-use sample datasets from the ENCODE eCLIP
project:

``` r
data(K562_bed)
head(K562_bed)
#>      V1       V2       V3   V4        V5 V6       V7          V8
#> 1 chr16 88894182 88894282 AATF   CBFA2T3  - 9.394340 5.68023e-10
#> 2 chr16 88894282 88894382 AATF   CBFA2T3  - 9.074450 5.68023e-10
#> 3 chr16 88894382 88894482 AATF   CBFA2T3  - 8.939215 5.68023e-10
#> 4 chr20 26209844 26209944 AATF MIR663AHG  - 8.637735 5.68023e-10
#> 5  chr3 13377662 13377762 AATF    NUP210  - 8.551355 5.68023e-10
#> 6 chr16   571060   571160 AATF      PIGQ  + 8.434275 5.68023e-10
```

------------------------------------------------------------------------

## Gene-level visualization

[`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)
plots all RBP peaks that overlap a single gene of interest.

``` r
result <- PlotGene(
    bed    = K562_bed,
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
    bed    = K562_bed,
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
    bed_file  = K562_bed,
    SEMATS    = sample_se.mats
)
```

The function classifies events into three groups based on
`IncLevelDifference` and `PValue` / `FDR`:

- **Retained** — `IncLevelDifference < -0.1`: exon inclusion
  significantly higher in condition 1
- **Excluded** — `IncLevelDifference > 0.1`: exon skipping significantly
  higher in condition 1
- **Control** — non-significant, stable inclusion (PValue \> 0.95)

A bootstrap procedure samples a matched number of control events across
`control_iterations` replicates to build a null distribution. Positions
where Retained or Excluded frequencies exceed the control z-score
threshold are highlighted as significant regions.

``` r
createSplicingMap(
    bed_file           = K562_bed,
    SEMATS             = sample_se.mats,
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

# YCAY motif (Y = C or T)
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

## Retained intron splicing map

[`createRetainedIntronSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSplicingMap.md)
quantifies RBP binding frequency around retained-intron junctions: the
upstream exon/intron boundary (UE-RI5) and the downstream intron/exon
boundary (RI3-DE).

``` r
createRetainedIntronSplicingMap(
    bed_file = K562_bed,
    RIMATS   = sample_se.mats
)
```

------------------------------------------------------------------------

## Retained intron sequence map

[`createRetainedIntronSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSequenceMap.md)
measures per-position motif frequency at retained intron junctions using
the same approach as
[`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md).

``` r
library(BSgenome.Hsapiens.UCSC.hg38)

createRetainedIntronSequenceMap(
    RIMATS   = sample_se.mats,
    sequence = "YCAY"
)
```
