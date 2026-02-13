# RNAPeaks

**Visualize RNA-Binding Protein Peaks on Gene Structures**

[![Status: Beta](https://img.shields.io/badge/Status-Beta-yellow.svg)]()
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

RNAPeaks is an R package for creating quality visualizations of RNA-binding protein (RBP) peaks overlaid on gene structure diagrams, sequence maps, and splicing maps. It takes BED files containing protein binding sites and overlays them on gene annotations from Ensembl GTF files, enabling researchers to visualize where RBPs bind relative to exons, introns, and UTRs.

## Features

- **Single-gene visualization**: Plot RBP peaks on individual genes with `PlotGene()`
- **Region-based visualization**: View multiple genes within a genomic window with `PlotRegion()`
- **Splicing maps**: Analyze RBP binding frequency around splice junctions with `createSplicingMap()`
- **Sequence motif analysis**: Identify motif enrichment patterns around splice sites with `createSequenceMap()`
- **Statistical significance testing**: Bootstrap-based z-score testing to identify significant enrichment regions
- **Automatic GTF integration**: Fetches Ensembl annotations via AnnotationHub (Human)
- **Parallel processing**: Multi-core support for faster analysis of large datasets
- **Customizable styling**: Extensive options for colors, sizes, labels, and layout

## Installation

### Development Version from GitHub

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
    "GenomicRanges", "IRanges", "S4Vectors", "BiocGenerics",
    "GenomeInfoDb", "AnnotationHub", "Biostrings", "BSgenome"
))

# Install from GitHub
# Install CRAN dependencies
install.packages(c("ggplot2", "dplyr", "scales", "magrittr"))
```

### Install RNAPeaks

```r
# Install from GitHub
# install.packages("devtools")
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
result <- PlotRegion(
    bed = sample_bed,
    gtf = gtf,
    Chr = "12",
    Start = 56000000,
    End = 56050000,
    Strand = "+"
)
result$plot
```

### Splicing Map Analysis
---

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

The `checkBed()` function automatically validates your BED file format.

Analyze where RBPs bind relative to exon/intron boundaries in splicing events:

```r
# Using included sample data
createSplicingMap(
    bed_file = sample_bed,
    SEMATS = sample_se.mats
)

# With parallel processing and custom parameters
createSplicingMap(
    bed_file = sample_bed,
    SEMATS = sample_se.mats,
    cores = 4,
    control_iterations = 20,      # Bootstrap iterations for control sampling
    control_multiplier = 2.0,     # Sample size = (retained + excluded) * multiplier
    show_significance = TRUE      # Show significant enrichment regions
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

# With significance testing
createSequenceMap(
    SEMATS = sample_se.mats,
    sequence = "YGCY",
    control_iterations = 20,
    z_threshold = 1.96,          
    min_consecutive = 10,         # Minimum positions for significant region
    show_significance = TRUE
)

# Return data for custom analysis
freq_data <- createSequenceMap(
    SEMATS = sample_se.mats,
    sequence = "CCCC",
    return_data = TRUE
)
```

### Validating Bootstrap Normality

The bootstrap sampling assumes normality. You can verify this:

```r
# Get diagnostics including raw bootstrap samples
diag <- createSequenceMap(
    SEMATS = sample_se.mats,
    sequence = "YGCY",
    return_diagnostics = TRUE
)

# Test normality at all positions
mean(apply(diag$bootstrap_matrix, 1, function(x) shapiro.test(x)$p.value) > 0.05)
```

## API Reference

### Main Functions

| Function | Description |
|----------|-------------|
| `PlotGene()` | Plot RBP peaks on a single gene structure |
| `PlotRegion()` | Plot RBP peaks across a genomic region with multiple genes |
| `createSplicingMap()` | Analyze protein binding frequency around splice junctions |
| `createSequenceMap()` | Analyze sequence motif frequency around splice junctions |

### Helper Functions

| Function | Description |
|----------|-------------|
| `LoadGTF()` | Load GTF annotation from AnnotationHub (Human) |
| `checkBed()` | Validate and normalize BED file format |

### Key Parameters for Splicing/Sequence Maps

| Parameter | Default | Description |
|-----------|---------|-------------|
| `control_multiplier` | 2.0 | Sample size = (retained + excluded) * multiplier |
| `control_iterations` | 20 | Number of bootstrap iterations for control group |
| `z_threshold` | 1.96 | Z-score threshold for significance (1.96 = p < 0.05) |
| `min_consecutive` | 10 | Minimum consecutive positions for significant region |
| `show_significance` | TRUE | Display significance bars on plot |
| `cores` | 1 | Number of cores for parallel processing |
| `return_data` | FALSE | Return data frame instead of plot |
| `return_diagnostics` | FALSE | Return bootstrap samples for normality testing |

## Input Data Format

### BED File Requirements

Your BED file should contain peak/binding site data. Columns are mapped by position:

| Column | Name | Description | Required |
|--------|------|-------------|----------|
| 1 | chr | Chromosome (e.g., "chr1" or "1") | Yes |
| 2 | start | Start position (0-based) | Yes |
| 3 | end | End position | Yes |
| 4 | tag | Name/Target (protein identifier) | No (default: "peak") |
| 5 | score | Score | No (default: 0) |
| 6 | strand | Strand ("+" or "-") | No (default: "+") |

The `checkBed()` function validates and standardizes your BED file, filling in defaults for missing columns.

### SE.MATS Format (for Splicing Analysis)

For `createSplicingMap()` and `createSequenceMap()`, provide SE.MATS output from rMATS with columns including:

- `chr`, `strand` - Genomic location
- `exonStart_0base`, `exonEnd` - Skipped exon coordinates
- `upstreamES`, `upstreamEE` - Upstream exon coordinates
- `downstreamES`, `downstreamEE` - Downstream exon coordinates
- `PValue`, `FDR`, `IncLevelDifference` - Statistical results
- `IJC_SAMPLE_1`, `SJC_SAMPLE_1`, `IJC_SAMPLE_2`, `SJC_SAMPLE_2` - Junction counts
- `IncLevel1`, `IncLevel2` - Inclusion levels

## Included Sample Data

- `sample_bed`: K562 cell line RBP binding peaks
- `sample_se.mats`: Skipped exon alternative splicing events (87,736 events)

```r
# View sample data
head(sample_bed)
head(sample_se.mats)
```

## Documentation

Full function documentation:

```r
?PlotGene
?PlotRegion
?createSplicingMap
?createSequenceMap
?checkBed
?LoadGTF
```

## Contributing

Report issues and contribute at [GitHub](https://github.com/Krushna-B/RNAPeaks/issues).

## License

MIT License. See [LICENSE](LICENSE) for details.
