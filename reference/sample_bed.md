# Sample K562 RBP Binding Peaks

A curated dataset of RNA-binding protein (RBP) peak calls from the K562
chronic myelogenous leukemia cell line, provided for testing and
demonstration of the RNAPeaks visualization functions.

## Usage

``` r
sample_bed
```

## Format

A data frame with the following columns:

- chr:

  Chromosome identifier (without "chr" prefix, e.g., "1", "X")

- start:

  Peak start coordinate (0-based, BED convention)

- end:

  Peak end coordinate (1-based, BED convention)

- tag:

  Protein or peak identifier

- score:

  Peak score or confidence value

- strand:

  Genomic strand ("+" or "-")

## Source

ENCODE Project (<https://www.encodeproject.org/>)

## Examples

``` r
data(sample_bed)
head(sample_bed)
#>   chr    start      end           tag score strand
#> 1  21  8401778  8401840 AARS_K562_IDR  1000      +
#> 2  19  4035770  4035869 AARS_K562_IDR  1000      +
#> 3  11 93721690 93721719 AARS_K562_IDR  1000      +
#> 4  19  2270224  2270300 AARS_K562_IDR  1000      +
#> 5   m    12167    12208 AARS_K562_IDR  1000      +
#> 6  16  2760062  2760143 AARS_K562_IDR  1000      +

if (FALSE) { # \dontrun{
  data(gtf_human)
  PlotGene(bed = sample_bed, geneID = "GAPDH", gtf = gtf_human)
} # }
```
