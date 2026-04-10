# createSequenceMap

## Overview

[`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)
measures per-position sequence motif frequency across the same
four-region window used by
[`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md),
applied to skipped-exon (SE) events. Full IUPAC ambiguity codes are
supported, enabling analysis of degenerate RBP binding motifs.

The analysis logic, event classification, bootstrap procedure, and
significance testing are identical to
[`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md).
See the [createSplicingMap
article](https://krushna-b.github.io/RNAPeaks/articles/splicing-maps.md)
for details on those parameters.

``` r
library(RNAPeaks)
```

------------------------------------------------------------------------

## Basic usage

``` r
library(BSgenome.Hsapiens.UCSC.hg38)

createSequenceMap(
    SEMATS   = sample_se.mats,
    sequence = "YCAY"    # NOVA binding motif; Y = C or T
)
```

------------------------------------------------------------------------

## IUPAC ambiguity codes

Any standard IUPAC code is supported:

| Code | Matches    |
|------|------------|
| `R`  | A, G       |
| `Y`  | C, T       |
| `S`  | G, C       |
| `W`  | A, T       |
| `K`  | G, T       |
| `M`  | A, C       |
| `B`  | C, G, T    |
| `D`  | A, G, T    |
| `H`  | A, C, T    |
| `V`  | A, C, G    |
| `N`  | A, C, G, T |

``` r
# YGCY motif (MBNL binding site)
createSequenceMap(SEMATS = sample_se.mats, sequence = "YGCY")

# Poly-C run
createSequenceMap(SEMATS = sample_se.mats, sequence = "CCCC")
```

------------------------------------------------------------------------

## Multiple motifs

Pass a character vector to `sequence` to generate one plot per motif:

``` r
plots <- createSequenceMap(
    SEMATS   = sample_se.mats,
    sequence = c("YCAY", "YGCY", "CCCC")
)
plots[["YCAY"]]
```

------------------------------------------------------------------------

## Using a custom genome

By default
[`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)
uses `BSgenome.Hsapiens.UCSC.hg38`. Supply any installed BSgenome object
to use a different assembly:

``` r
library(BSgenome.Mmusculus.UCSC.mm10)

createSequenceMap(
    SEMATS   = sample_se.mats,
    sequence = "YCAY",
    genome   = BSgenome.Mmusculus.UCSC.mm10
)
```

------------------------------------------------------------------------

## Returning data

``` r
freq_df <- createSequenceMap(
    SEMATS      = sample_se.mats,
    sequence    = "YCAY",
    return_data = TRUE
)
head(freq_df)
```
