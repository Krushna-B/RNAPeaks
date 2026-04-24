# K562 RBP Binding Peaks

ENCODE eCLIP peak calls from the K562 (chronic myelogenous leukemia)
cell line, provided for testing and demonstration of RNAPeaks
visualization and analysis functions.

## Usage

``` r
K562_bed
```

## Format

A data frame with the following columns:

- chr:

  Chromosome identifier (without "chr" prefix, e.g., "1", "X")

- start:

  Peak start coordinate

- end:

  Peak end coordinate

- tag:

  RBP name or peak identifier

- score:

  Peak score or confidence value

- strand:

  Genomic strand ("+" or "-")

## Source

ENCODE Project (<https://www.encodeproject.org/>)

## Examples

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

if (FALSE) { # \dontrun{
  data(K562_bed)
  data(gtf_human)
  PlotGene(bed = K562_bed, geneID = "GAPDH", gtf = gtf_human)
} # }
```
