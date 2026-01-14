# RNAPeaks

**Visualize RNA-Binding Protein Peaks on Gene Structures**

[![Beta](https://img.shields.io/badge/status-beta-yellow.svg)](https://github.com/Krushna-B/RNAPeaks)
[![R Version](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

RNAPeaks is an R package for creating quality visualizations of RNA-binding protein (RBP) peaks overlaid on gene structure diagrams. It takes BED files containing protein binding sites and overlays them on gene annotations from Ensembl GTF files, enabling researchers to visualize where RBPs bind relative to exons, introns, and UTRs.

---

## Features

- **Single-gene visualization**: Plot RBP peaks on individual genes with `PlotGene()`
- **Region-based visualization**: View multiple genes within a genomic window with `PlotRegion()`
- **Automatic GTF integration**: Fetches Ensembl annotations via AnnotationHub (Human)
- **Customizable styling**: Extensive options for colors, sizes, labels, and layout
- **Region highlighting**: Highlight specific genomic regions of interest

---

## Installation

### Prerequisites

RNAPeaks requires R 4.0 or higher and several Bioconductor packages.

```r
# Install Bioconductor if not already installed
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

# Install Bioconductor dependencies
BiocManager::install(c(
    "GenomicRanges",
    "IRanges",
    "S4Vectors",
    "BiocGenerics",
    "GenomeInfoDb",
    "AnnotationHub"
))

# Install CRAN dependencies
install.packages(c("ggplot2", "dplyr", "scales", "magrittr"))
```

### Install RNAPeaks

```r
# Install from GitHub (recommended for beta)
# install.packages("devtools")
devtools::install_github("Krushna-B/RNAPeaks")
```

Or install from a local clone:

```r
devtools::install("/path/to/RNAPeaks")
```

---

## Quick Start

### 1. Load the Package

```r
library(RNAPeaks)
```

### 2. Load Your BED File

```r
# Load BED file with RBP peaks
bed <- read.table("path/to/your/peaks.bed", header = FALSE, sep = "\t")

# Validate the BED format
bed <- checkBed(bed)
Success will be marked by printing bedfile in console
```

An example BED file is provided in the repository under `Testing Bed/`:

```r
# If working from the cloned repository
bed <- read.table(
    "Testing Bed/Sample_Bed_(All_K562_peaks).bed",
    header = FALSE
)
```

This sample contains K562 cell line peaks and can be used to test the package functionality.

### 3. Load GTF Annotations (One-Time Setup)

Loading GTF annotations from AnnotationHub can take several minutes on first run. We recommend saving the result for future sessions:

```r
# Load GTF (do this once per species)
gtf <- LoadGTF(species = "Human")

# Save for future sessions
saveRDS(gtf, "human_gtf.rds")

# In future R sessions, load directly:
# gtf <- readRDS("human_gtf.rds")
```

### 4. Generate Plots

#### Plot peaks on a single gene:

```r
result <- PlotGene(
    bed = bed,
    geneID = "BRCA1",
    gtf = gtf
)

# Access the plot
result[[1]]

# Access the filtered BED data
result[[2]]
```

#### Plot peaks across a genomic region:

```r
result <- PlotRegion(
    bed = bed,
    gtf = gtf,
    Chr = "17",
    Start = 7565097,
    End = 7590856,
    Strand = "-"
)
```

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

The `checkBed()` function automatically validates and standardizes your BED file format.

### GTF Annotations

RNAPeaks fetches GTF annotations from Ensembl via AnnotationHub. Currently supported:

- **Human**: Ensembl GRCh38 annotations
- **Mouse**: Ensembl GRCm39 annotations

---

## Customization Options

RNAPeaks provides extensive styling options. Pass these as additional arguments:

### Gene Structure Appearance

```r
PlotGene(
    ...,
    exon_fill = "navy",           # Exon/CDS color
    utr_fill = "lightgray",       # UTR color
    intron_color = "gray60",      # Intron line color
    intron_linewidth = 0.9        # Intron line width
)
```

### Peak Styling

```r
PlotGene(
    ...,
    peak_col = "purple",          # Peak fill color
    peak_alpha = 0.95,            # Peak opacity
    peak_border_color = "black"   # Peak border color
)
```

### Labels and Text

```r
PlotGene(
    ...,
    title_size = 25,
    subtitle_size = 12,
    label_size = 5,               # Protein labels
    axis_text_size = 9
)
```

### Region Highlighting

```r
PlotGene(
    ...,
    highlighted_region_start = 1000000,
    highlighted_region_stop = 1005000,
    highlighted_region_color = "pink",
    highlighted_region_opacity = 0.30
)
```

See `?PlotGene` or `?PlotRegion` for complete documentation of all styling parameters.

---

## Example Workflow

```r
library(RNAPeaks)

# Step 1: Load and validate BED data
bed <- read.table("peaks.bed", header = FALSE, sep = "\t")
bed <- checkBed(bed)

# Step 2: Load GTF (or use cached version)
gtf <- LoadGTF(species = "Human")
saveRDS(gtf, "human_gtf.rds")  # Cache for next time

# Step 3: Generate single-gene plot
result <- PlotGene(
    bed = bed,
    geneID = "TP53",
    gtf = gtf,
    peak_col = "steelblue",
    title_size = 20,
    RNA_Peaks_File_Path = "TP53_peaks.pdf"
)

# Step 4: Generate region plot for broader context
result <- PlotRegion(
    bed = bed,
    gtf = gtf,
    Chr = "17",
    Start = 7560000,
    End = 7600000,
    Strand = "-",
    RNA_Peaks_File_Path = "chr17_region.pdf"
)
```

---

## Troubleshooting

### Common Issues

**"No rows in GTF found"**
- Ensure your gene ID matches Ensembl gene symbols (e.g., "TP53", "BRCA1")
- Check that the chromosome naming matches (with or without "chr" prefix)
- Verify the strand orientation

**GTF loading is slow**
- This is expected on first run. Save with `saveRDS()` and reload with `readRDS()` for subsequent sessions.

**Empty plot or no peaks shown**
- Verify your BED file coordinates overlap with the gene region
- Check that the strand in your BED file matches the gene strand
- Use `checkBed()` to validate your BED format

---
## Contributing

RNAPeaks is currently in beta testing. We welcome feedback and bug reports:

- Report issues on [GitHub Issues](https://github.com/Krushna-B/RNAPeaks/issues)
---

## License

MIT License. See [LICENSE](LICENSE) for details.
