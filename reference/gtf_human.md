# Human GTF Gene Annotation

Pre-loaded Ensembl GTF annotation for human genes (GRCh38). This bundled
dataset eliminates the need to download from AnnotationHub, enabling
offline use and faster loading.

## Usage

``` r
gtf_human
```

## Format

A data frame containing GTF annotation with the following columns:

- seqnames:

  Chromosome

- start:

  Feature start position

- end:

  Feature end position

- width:

  Feature width in bp

- strand:

  Strand (+ or -)

- source:

  Annotation source

- type:

  Feature type (gene, transcript, exon, CDS, UTR, etc.)

- score:

  Annotation score

- phase:

  CDS phase (0, 1, or 2)

- gene_id:

  Ensembl gene ID

- gene_version:

  Ensembl gene version

- gene_name:

  Gene symbol

- gene_source:

  Gene annotation source

- gene_biotype:

  Gene biotype (protein_coding, lncRNA, etc.)

- transcript_id:

  Ensembl transcript ID

- transcript_version:

  Ensembl transcript version

- transcript_name:

  Transcript name

- transcript_source:

  Transcript annotation source

- transcript_biotype:

  Transcript biotype

- transcript_support_level:

  Transcript support level (TSL)

- exon_number:

  Exon number within the transcript

- exon_id:

  Ensembl exon ID

- exon_version:

  Ensembl exon version

- protein_id:

  Ensembl protein ID

- protein_version:

  Ensembl protein version

- ccds_id:

  CCDS identifier

- tag:

  Feature tag (e.g., basic, Ensembl_canonical)

## Source

Ensembl via AnnotationHub (AH110867)

## Examples

``` r
if (FALSE) { # \dontrun{
  # Load the bundled GTF
  data(gtf_human)

  # Use with PlotGene (no download required)
  data(sample_bed)
  PlotGene(bed = sample_bed, geneID = "GAPDH", gtf = gtf_human)
} # }
```
