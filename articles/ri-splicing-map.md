# createRetainedIntronSplicingMap

## Overview

[`createRetainedIntronSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createRetainedIntronSplicingMap.md)
quantifies RBP binding frequency around the two junctions of
retained-intron (RI) events: the upstream exon/intron boundary and the
downstream intron/exon boundary. It accepts the same SE.MATS column
format and the same parameters as
[`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md).

``` r
library(RNAPeaks)
```

------------------------------------------------------------------------

## Basic usage

``` r
createRetainedIntronSplicingMap(
    bed_file = sample_bed,
    RIMATS   = sample_se.mats
)
```

------------------------------------------------------------------------

## Event classification

Events are classified using the same thresholds as
[`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md):

- **Retained** — `IncLevelDifference < -0.1`
- **Excluded** — `IncLevelDifference > 0.1`
- **Control** — non-significant, stable inclusion

``` r
createRetainedIntronSplicingMap(
    bed_file                     = sample_bed,
    RIMATS                       = sample_se.mats,
    retained_IncLevelDifference  = -0.1,
    exclusion_IncLevelDifference = 0.1,
    p_valueRetainedAndExclusion  = 0.05,
    p_valueControls              = 0.95,
    Min_Count                    = 50
)
```

------------------------------------------------------------------------

## Bootstrap control and significance

The same bootstrap z-score approach is used as in
[`createSplicingMap()`](https://krushna-b.github.io/RNAPeaks/reference/createSplicingMap.md).
See the [createSplicingMap
article](https://krushna-b.github.io/RNAPeaks/articles/splicing-maps.md)
for full details.

``` r
createRetainedIntronSplicingMap(
    bed_file           = sample_bed,
    RIMATS             = sample_se.mats,
    control_iterations = 30,
    use_fdr            = TRUE,
    fdr_threshold      = 0.05,
    show_significance  = TRUE
)
```

------------------------------------------------------------------------

## Returning data

``` r
freq_df <- createRetainedIntronSplicingMap(
    bed_file    = sample_bed,
    RIMATS      = sample_se.mats,
    return_data = TRUE
)
head(freq_df)
```
