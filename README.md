# RNAPeaks <img src="man/figures/logo.png" align="right" height="120" alt="" />

**Publication-quality visualization of RNA-binding protein peaks on gene structures**

[![R-CMD-check](https://github.com/Krushna-B/RNAPeaks/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Krushna-B/RNAPeaks/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)

RNAPeaks is an R package for visualizing RNA-binding protein (RBP) peaks
overlaid on annotated gene structures. It accepts BED files of protein
binding sites and Ensembl GTF annotations to produce publication-ready
figures showing where RBPs bind relative to exons, introns, and UTRs.
Splicing map and sequence motif analyses quantify position-specific binding
frequency around splice junctions with bootstrap-based significance testing.

---

## Documentation

Full function reference and articles are available at the package website:

**<https://krushna-b.github.io/RNAPeaks/>**

---

## Features

| Feature | Function |
|---|---|
| Single-gene RBP peak visualization | `PlotGene()` |
| Multi-gene genomic region visualization | `PlotRegion()` |
| RBP binding frequency around splice junctions | `createSplicingMap()` |
| Sequence motif frequency around splice junctions | `createSequenceMap()` |
| Interactive web interface (no R required) | `launchApp()` |
| GTF annotation loading (Human / Mouse) | `LoadGTF()` |
| BED file validation and normalization | `checkBed()` |

Additional capabilities:

- Bootstrap-based z-score significance testing for enriched binding regions
- IUPAC ambiguity code support in sequence motif search
- Multi-core parallel processing (`cores` parameter) for large datasets
- Extensive ggplot2-based styling options
- Curated sample datasets for immediate exploration

---

## Installation

### Prerequisites

Install Bioconductor dependencies first:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "GenomicRanges", "IRanges", "S4Vectors", "BiocGenerics",
    "GenomeInfoDb", "AnnotationHub", "Biostrings", "BSgenome"
))
```

### Install RNAPeaks

```r
# install.packages("devtools")
devtools::install_github("Krushna-B/RNAPeaks")
```

For sequence motif analysis, also install the reference genome:

```r
BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
```

---

## Quick Start

### Gene-level visualization

```r
library(RNAPeaks)

# Load GTF annotation once per session (or load from a saved .rds file)
gtf <- LoadGTF(species = "Human")
saveRDS(gtf, "human_gtf.rds")   # cache for future sessions

# Plot RBP peaks on a single gene using the included sample data
result <- PlotGene(
    bed    = sample_bed,
    geneID = "GAPDH",
    gtf    = gtf
)

result$plot   # ggplot2 object
result$csv    # filtered BED data used in the plot
```

### Region-level visualization

```r
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

### Splicing map

Quantify RBP binding frequency around exon-intron boundaries across
skipped-exon splicing events:

```r
createSplicingMap(
    bed_file  = sample_bed,
    SEMATS    = sample_se.mats
)
```

With significance testing and parallel processing:

```r
createSplicingMap(
    bed_file           = sample_bed,
    SEMATS             = sample_se.mats,
    cores              = 4,
    control_iterations = 20,
    show_significance  = TRUE
)
```

### Sequence motif analysis

Identify enrichment of a sequence motif (IUPAC codes supported) around
splice junctions:

```r
library(BSgenome.Hsapiens.UCSC.hg38)

createSequenceMap(
    SEMATS   = sample_se.mats,
    sequence = "YCAY"           # Y = C or T
)
```

---

## Input data format

### BED file

Columns are mapped **by position** to the canonical BED fields. Use
`checkBed()` to validate and normalize your file before plotting.

| Column | Name | Description | Required |
|--------|------|-------------|----------|
| 1 | `chr` | Chromosome (e.g., `"chr1"` or `"1"`) | Yes |
| 2 | `start` | Start coordinate (0-based) | Yes |
| 3 | `end` | End coordinate | Yes |
| 4 | `tag` | Protein / peak identifier | No (default: `"peak"`) |
| 5 | `score` | Peak score | No (default: `0`) |
| 6 | `strand` | Strand (`"+"` or `"-"`) | No (default: `"+"`) |

```r
my_bed <- read.table("my_peaks.bed", header = FALSE, sep = "\t")
my_bed <- checkBed(my_bed)   # validates, normalizes, strips "chr" prefix
```

### SE.MATS file

For `createSplicingMap()` and `createSequenceMap()`, provide the
`SE.MATS.JC.txt` output from
[rMATS](https://rnaseq-mats.sourceforge.net/). Required columns:

`chr`, `strand`, `exonStart_0base`, `exonEnd`, `upstreamES`,
`upstreamEE`, `downstreamES`, `downstreamEE`, `PValue`, `FDR`,
`IncLevelDifference`, `IJC_SAMPLE_1`, `SJC_SAMPLE_1`,
`IJC_SAMPLE_2`, `SJC_SAMPLE_2`, `IncLevel1`, `IncLevel2`

---

## Included sample data

| Dataset | Description |
|---------|-------------|
| `sample_bed` | K562 cell line RBP binding peaks (eCLIP) |
| `sample_se.mats` | 87,736 skipped-exon splicing events from rMATS |
| `gtf_human` | Pre-loaded Ensembl GRCh38 GTF (no download required) |

```r
head(sample_bed)
head(sample_se.mats)
```

---

## Key parameters

### Splicing / Sequence maps

| Parameter | Default | Description |
|-----------|---------|-------------|
| `WidthIntoExon` | `50` | bp to extend into flanking exons |
| `WidthIntoIntron` | `300` | bp to extend into flanking introns |
| `moving_average` | `50` | smoothing window size (positions) |
| `control_multiplier` | `2.0` | control sample size = (retained + excluded) × multiplier |
| `control_iterations` | `20` | bootstrap iterations for control distribution |
| `z_threshold` | `1.96` | z-score threshold (p < 0.05) when `use_fdr = FALSE` |
| `min_consecutive` | `10` | minimum run of significant positions to report a region |
| `show_significance` | `TRUE` | overlay significance bars on the plot |
| `cores` | `1` | cores for parallel processing |
| `return_data` | `FALSE` | return the frequency data frame instead of a plot |

---

## Web interface

A browser-based interface is available for users who prefer not to write R code:

```r
launchApp()
```

The app supports all four analysis modes, file upload, and export to PDF or CSV.

---

## Citation

If you use RNAPeaks in your research, please cite:

> Bhanushali K, Giri G, et al. (2025). *RNAPeaks: publication-quality
> visualization of RNA-binding protein peaks on gene structures.*
> GitHub: <https://github.com/Krushna-B/RNAPeaks>

---

## Contributing

Bug reports and feature requests are welcome on the
[issue tracker](https://github.com/Krushna-B/RNAPeaks/issues).

---

## License

MIT License — see [LICENSE](LICENSE) for details.
