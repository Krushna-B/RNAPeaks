# Validate and Normalize a BED-Format Data Frame

Validates a data frame as BED format and normalizes column names and
chromosome identifiers. The function expects columns in standard BED
order: chr, start, end, tag/name, score, strand.

## Usage

``` r
checkBed(df)
```

## Arguments

- df:

  A data frame with at least 3 columns representing BED data. Columns
  are mapped by position: 1=chr, 2=start, 3=end, 4=tag, 5=score,
  6=strand. If fewer than 6 columns are provided, missing columns are
  filled with defaults: tag="peak", score=0.

## Value

A validated data frame with normalized column names (chr, start, end,
tag, score, strand) and chromosome names without "chr" prefix.

## Details

The function performs the following:

- Requires at least 4 columns (chr, start, end, strand)

- Maps columns by position to BED names

- Fills missing columns with defaults (tag="peak", score=0)

- Validates chromosome is character type

- Verifies end \>= start for all rows

- Validates strand values are "+" or "-"

- Removes "chr" prefix from chromosome names

- Converts chromosome names to uppercase

## Examples

``` r
if (FALSE) { # \dontrun{
  # Standard 6-column BED
  bed <- read.table("peaks.bed", header = FALSE)
  bed <- checkBed(bed)

  # 4-column BED (will add defaults for tag, score)
  bed3 <- data.frame(chr = "chr1", start = 100, end = 200, strand = "+")
  bed3 <- checkBed(bed3)
} # }
```
