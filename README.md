# RNAPeaks

**Visualize RNA-Binding Protein Peaks on Gene Structures**

[![R-CMD-check](https://github.com/Krushna-B/RNAPeaks/workflows/R-CMD-check/badge.svg)](https://github.com/Krushna-B/RNAPeaks/actions)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

RNAPeaks is an R/Bioconductor package for creating publication-quality visualizations of RNA-binding protein (RBP) peaks overlaid on gene structure diagrams. It takes BED files containing protein binding sites and overlays them on gene annotations from Ensembl GTF files, enabling researchers to visualize where RBPs bind relative to exons, introns, and UTRs.

## Features

- **Single-gene visualization**: Plot RBP peaks on individual genes with `PlotGene()`
- **Region-based visualization**: View multiple genes within a genomic window with `PlotRegion()`
- **Splicing maps**: Analyze RBP binding frequency around splice junctions with `createSplicingMap()`
- **Sequence motif analysis**: Identify motif enrichment patterns around splice sites with `createSequenceMap()`
- **Automatic GTF integration**: Fetches Ensembl annotations via AnnotationHub (Human/Mouse)
- **Customizable styling**: Extensive options for colors, sizes, labels, and layout

## Installation

### From Bioconductor (when available)

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("RNAPeaks")
```

### Development Version from GitHub

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "GenomicRanges", "IRanges", "S4Vectors", "BiocGenerics",
    "GenomeInfoDb", "AnnotationHub", "Biostrings", "BSgenome"
))

# Install from GitHub
devtools::install_github("Krushna-B/RNAPeaks")
```

## Quick Start

### Load the Package

```r
library(RNAPeaks)
```

### Gene-Level Visualization

```r
# Load GTF annotations (do once per session)
gtf <- LoadGTF(species = "Human")
saveRDS(gtf, "human_gtf.rds")  # Cache for future sessions

# Plot RBP peaks on a gene using included sample data
result <- PlotGene(
    bed = sample_bed,
    geneID = "GAPDH",
    gtf = gtf
)
result$plot
```

### Region-Level Visualization

```r
PlotRegion(
    bed = sample_bed,
    gtf = gtf,
    Chr = "12",
    Start = 56000000,
    End = 56050000,
    Strand = "+"
)
```

### Splicing Map Analysis

Analyze where RBPs bind relative to exon/intron boundaries in splicing events:

```r
# Using included sample data
createSplicingMap(
    bed_file = sample_bed,
    SEMATS = sample_se.mats
)
```

### Sequence Motif Analysis

Identify sequence motif enrichment patterns around splice sites:
```r
library(BSgenome.Hsapiens.UCSC.hg38)

# Search for YCAY motif (Y = C or T)
createSequenceMap(
    SEMATS = sample_se.mats,
    sequence = "YCAY"
)

# Return data for custom analysis
freq_data <- createSequenceMap(
    SEMATS = sample_se.mats,
    sequence = "CCCC",
    return_data = TRUE
)
```

## Input Data Format

### BED File Requirements

Your BED file should contain peak/binding site data with the following columns:

| Column | Description |
|--------|-------------|
| 1 | Chromosome (e.g., "chr1" or "1") |
| 2 | Start position (0-based) |
| 3 | End position |
| 4 | Name/Target (protein identifier) |
| 5 | Score (optional) |
| 6 | Strand ("+" or "-") |

The `checkBed()` function validates and standardizes your BED file format.

### SE.MATS Format (for Splicing Analysis)

For `createSplicingMap()` and `createSequenceMap()`, provide SE.MATS output from rMATS with columns including: chr, strand, exonStart_0base, exonEnd, upstreamES, upstreamEE, downstreamES, downstreamEE, PValue, FDR, IncLevelDifference.

## Included Sample Data

- `sample_bed`: K562 cell line RBP binding peaks
- `sample_se.mats`: Skipped exon alternative splicing events

```r
# View sample data
head(sample_bed)
head(sample_se.mats)
```

## Documentation

See the package vignette for detailed usage:

```r
browseVignettes("RNAPeaks")
```

Full function documentation:

```r
?PlotGene
?PlotRegion
?createSplicingMap
?createSequenceMap
```

## Contributing

Report issues and contribute at [GitHub](https://github.com/Krushna-B/RNAPeaks/issues).

## License

MIT License. See [LICENSE](LICENSE) for details.
