# HepG2 RBP Binding Peaks

ENCODE eCLIP peak calls from the HepG2 (hepatocellular carcinoma) cell
line, provided for testing and demonstration of RNAPeaks visualization
and analysis functions.

## Usage

``` r
HepG2_bed
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
data(HepG2_bed)
head(HepG2_bed)
#>      V1        V2        V3    V4      V5 V6       V7          V8
#> 1 chr12 123510272 123510372 AGGF1  RILPL1  - 10.39385 9.13816e-09
#> 2 chr19   4036075   4036174 AGGF1   PIAS4  + 10.17775 9.13816e-09
#> 3 chr19   4035875   4035975 AGGF1   PIAS4  +  9.93753 1.35672e-08
#> 4 chr13  29558255  29558355 AGGF1  SLC7A1  -  9.91781 9.70116e-09
#> 5  chr6   2376477   2376576 AGGF1 GMDS-DT  +  9.80050 9.13816e-09
#> 6  chr2  10312403  10312503 AGGF1  HPCAL1  +  9.74314 1.45070e-08

if (FALSE) { # \dontrun{
  data(HepG2_bed)
  data(gtf_human)
  PlotGene(bed = HepG2_bed, geneID = "GAPDH", gtf = gtf_human)
} # }
```
