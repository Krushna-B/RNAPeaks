# createRetainedIntronSequenceMap

## Overview

[`createRetainedIntronSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSequenceMap.md)
measures per-position sequence motif frequency at retained-intron
junctions using the same two-bin window as
[`createRetainedIntronSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSplicingMap.md).
It accepts the same parameters as
[`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md)
with `RIMATS` in place of `SEMATS`.

``` r
library(RNAPeaks)
```

------------------------------------------------------------------------

## Basic usage

``` r
library(BSgenome.Hsapiens.UCSC.hg38)

createRetainedIntronSequenceMap(
    RIMATS   = sample_se.mats,
    sequence = "YCAY"
)
```

------------------------------------------------------------------------

## Multiple motifs

``` r
plots <- createRetainedIntronSequenceMap(
    RIMATS   = sample_se.mats,
    sequence = c("YCAY", "YGCY")
)
plots[["YCAY"]]
```

------------------------------------------------------------------------

## Using a custom genome

``` r
library(BSgenome.Mmusculus.UCSC.mm10)

createRetainedIntronSequenceMap(
    RIMATS   = sample_se.mats,
    sequence = "YCAY",
    genome   = BSgenome.Mmusculus.UCSC.mm10
)
```

------------------------------------------------------------------------

## IUPAC codes and returning data

The full IUPAC code table and `return_data` behavior are identical to
[`createSequenceMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSequenceMap.md).
See the [createSequenceMap
article](https://krushna-b.github.io/RNAPeaks/articles/sequence-map.md)
for details.
