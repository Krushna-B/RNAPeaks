# RNAPeaks: Visualize RNA-Binding Protein Peaks on Gene Structures

RNAPeaks creates publication-quality visualizations of RNA-binding
protein (RBP) peaks overlaid on gene structures. It supports single-gene
plots, multi-gene genomic region plots, splicing maps, and sequence
motif analysis around splice junctions. Integrates with GTF annotations
from Ensembl via AnnotationHub.

## Main Functions

The primary functions for generating plots are:

- [`PlotGene`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md):

  Plot RBP peaks on a single gene structure. Use this when you want to
  visualize peaks on one specific gene.

- [`PlotRegion`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md):

  Plot RBP peaks across a genomic region containing multiple genes. Use
  this for broader regional views.

- [`createSplicingMap`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md):

  Analyze protein binding frequency around splice junctions. Shows where
  RBPs bind relative to exon/intron boundaries in retained vs excluded
  splicing events.

- [`createSequenceMap`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md):

  Analyze sequence motif frequency around splice junctions. Identifies
  position-specific enrichment of motifs in different splicing event
  categories.

- [`createRetainedIntronSequenceMap`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSequenceMap.md):

  Analyze sequence motif frequency around retained intron splice
  junctions. Identifies position-specific enrichment of motifs in
  different splicing event categories.

- [`createRetainedIntronSplicingMap`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSplicingMap.md):

  Analyze protein binding frequency around retained intron splice
  junctions. Shows where RBPs bind relative to exon/intron boundaries in
  retained vs excluded splicing events.

## Included Data

The package includes sample data for testing:

- [`sample_se.mats`](https://krushna-b.github.io/RNAPeaks/reference/sample_se.mats.md):

  Sample SE.MATS output for splicing analysis

## Workflow

A typical workflow involves:

1.  Load your own BED file and validate it with
    [`checkBed`](https://krushna-b.github.io/RNAPeaks/reference/checkBed.md)

2.  Load GTF annotation once with
    [`LoadGTF`](https://krushna-b.github.io/RNAPeaks/reference/LoadGTF.md)
    and store locally (e.g., `gtf <- LoadGTF("Human")`; save with
    `saveRDS(gtf, "gtf.rds")` for future sessions), or use the gtf_human
    that is pre-loaded

3.  Call
    [`PlotGene`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)
    or
    [`PlotRegion`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md),
    passing the stored `gtf` object

For splicing analysis:

1.  Prepare SE.MATS file or use `sample_se.mats`

2.  Call
    [`createSplicingMap`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md)
    with BED and SE.MATS data to analyze protein binding

3.  Call
    [`createSequenceMap`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)
    to analyze sequence motifs

4.  Call
    [`createRetainedIntronSplicingMap`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSplicingMap.md)
    with BED and SE.MATS data to analyze protein binding for retained
    introns

5.  Call
    [`createRetainedIntronSequenceMap`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSequenceMap.md)
    to analyze sequence motifs for retained introns

## Helper Functions

- [`checkBed`](https://krushna-b.github.io/RNAPeaks/reference/checkBed.md):

  Validate and standardize BED file format

- [`LoadGTF`](https://krushna-b.github.io/RNAPeaks/reference/LoadGTF.md):

  Load GTF annotation from AnnotationHub. Call once and store the result
  in a local variable or save to disk with
  [`saveRDS()`](https://rdrr.io/r/base/readRDS.html) to avoid repeated
  downloads.

## See also

Useful links:

- <https://krushna-b.github.io/RNAPeaks/>

- <https://github.com/Krushna-B/RNAPeaks>

- Report bugs at <https://github.com/Krushna-B/RNAPeaks/issues>

## Author

**Maintainer**: Krushna Bhanushali <Krushna.Bhanushali@unc.edu>

Authors:

- Gilbert Giri <Ggiri@unc.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
library(RNAPeaks)

# Load GTF once and store locally (do this once per session)
gtf <- LoadGTF(species = "Human")

# Optionally save for future sessions to avoid re-downloading
saveRDS(gtf, "human_gtf.rds")
# In future sessions: gtf <- readRDS("human_gtf.rds")

# ----- Gene-level visualization -----
result <- PlotGene(
  bed = your_bed,
  geneID = "GAPDH",
  gtf = gtf
)
result$plot

# ----- Region-level visualization -----
PlotRegion(
  bed = your_bed,
  Chr = "12",
  Start = 56000000,
  End = 56050000,
  Strand = "+",
  gtf = gtf
)

# ----- Splicing map analysis -----
createSplicingMap(
  bed_file = your_bed,
  SEMATS = sample_se.mats
)

# ----- Sequence motif analysis -----
library(BSgenome.Hsapiens.UCSC.hg38)
createSequenceMap(
  SEMATS = sample_se.mats,
  sequence = "CCCC"
)
} # }
```
