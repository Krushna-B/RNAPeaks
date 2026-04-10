# Can be called once and the result passed to `PlotGene()` or `PlotRegion()` to avoid repeated downloads.

Can be called once and the result passed to
[`PlotGene()`](https://krushna-b.github.io/RNAPeaks/reference/PlotGene.md)
or
[`PlotRegion()`](https://krushna-b.github.io/RNAPeaks/reference/PlotRegion.md)
to avoid repeated downloads.

## Usage

``` r
LoadGTF(species = "Human", file = NULL)
```

## Arguments

- species:

  Species to load annotation for: "Human" or "Mouse".

- file:

  Optional file path to a local GTF file. If provided, imports directly
  without connecting to AnnotationHub.

## Value

A data frame containing GTF annotation with columns including seqnames,
start, end, strand, type, gene_id, gene_name, transcript_id, etc.

## Details

Loading annotations from AnnotationHub can take time on first use. By
calling this function separately, you can:

- Load the annotation once and reuse it across multiple plots

- Save the annotation to disk for faster future sessions

Human annotations use Ensembl GTF (AH110867). Mouse annotations use
Ensembl GTF (AH47076).

## Examples

``` r
if (FALSE) { # \dontrun{
  # Load human GTF once
  gtf <- LoadGTF(species = "Human")

  # Optionally save for future sessions
  saveRDS(gtf, "human_gtf.rds")

  # Use in multiple plots without reloading
  PlotGene(bed = bed, geneID = "TP53", gtf = gtf)
  PlotGene(bed = bed, geneID = "BRCA1", gtf = gtf)

  # Load from saved file in future sessions
  gtf <- readRDS("human_gtf.rds")
} # }
```
