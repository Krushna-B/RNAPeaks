# RNAPeaks <img src="man/figures/logo.png" align="right" height="120" alt="" />

RNAPeaks is an R package for making publication quality figures of RNA binding protein (RBP) peaks on top of gene structures, together with statistical analysis of where those peaks and their motifs fall around splice junctions.

Everything here is plain R functions that return ggplot objects and tidy data frames, so you can drop the results straight into your own analysis.

## Installation

RNAPeaks is a Bioconductor style package, so the easiest path is through BiocManager. It will pull the Bioconductor dependencies for you.

```r
install.packages("BiocManager")
BiocManager::install("Krushna-B/RNAPeaks")
```

If you prefer remotes that works too:

```r
install.packages("remotes")
remotes::install_github("Krushna-B/RNAPeaks")
```

The sequence maps and the kmer enrichment pull genomic sequence, so they need a matching BSgenome. These are optional and only required if you use those functions:

```r
BiocManager::install(c(
  "BSgenome.Hsapiens.UCSC.hg38",
  "BSgenome.Mmusculus.UCSC.mm10",
  "BSgenome.Mmusculus.UCSC.mm39"
))
```

RNAPeaks ships with GENCODE annotations for hg38, mm10 and mm39, so for the visual functions you do not need to supply your own GTF unless you want to.

## Usage

### Peaks over a single gene

```r
library(RNAPeaks)

p <- plot_gene(
  bed     = K562_bed,   # a BED data frame, a file path, or a named list of them
  gene    = "SF3B1",    # gene symbol or Ensembl gene id
  species = "hg38"      # which bundled annotation to use: hg38, mm10 or mm39
)
p
```

You can pass one or more BAM files to draw coverage tracks above the gene, point `transcript` at a specific isoform, and tune how peaks are grouped and ordered through `peaks_options()`.

### Peaks over a genomic region

When you want a window rather than one gene, give coordinates instead:

```r
plot_region(
  bed    = K562_bed,
  chr    = "chr3",
  start  = 128470000,
  end    = 128500000,
  strand = "+"
)
```

### Splicing maps

A splicing map takes an rMATS event table and a set of peaks and reports how often peaks land at each position around the splice junctions of that event type. There is one function per event type.

```r
res <- skipped_exon_splicing_map(
  events   = se_mats_jc,   # rMATS SE junction count table (bundled example)
  bed_file = K562_bed,     # peaks as a BED data frame or a file path
  opts     = splicing_options(width_exon = 50, width_intron = 250),
  title    = "RBFOX2 in K562"
)

res$plot   # a ggplot you can save or theme further
res$data   # the positive and negative event tables behind the curve
```

The other three work the same way:

```r
five_prime_splicing_map(a5ss_mats_jc, K562_bed)       # alternative 5' splice site
three_prime_splicing_map(a3ss_mats_jc, K562_bed)      # alternative 3' splice site
retained_intron_splicing_map(ri_mats_jc, K562_bed)    # retained intron
```

### Sequence maps

Sequence maps are the motif version of the same idea. Instead of peaks you give one or more motifs, and RNAPeaks counts how often each motif occurs at every position.

```r
res <- skipped_exon_sequence_map(
  events     = se_mats_jc,
  sequence   = c("GCATG", "TGCATG"),  # IUPAC ambiguity codes are allowed, U is read as T
  genome     = "hg38",                # a genome key or a BSgenome object
  motif_mode = "combined"             # combined pools the motifs, individual makes one map per motif
)
```

The `five_prime_sequence_map`, `three_prime_sequence_map` and `retained_intron_sequence_map` functions cover the other event types.

### Binding across UTRs

`plot_utr_binding` averages peak density over the 5' and 3' UTRs across many genes, with one curve per BED track and an optional split by gene group.

```r
plot_utr_binding(
  bed     = list(K562 = K562_bed, HepG2 = HepG2_bed),
  species = "hg38"
)
```

### kmer enrichment

Compare kmer content between two sets. Each set can be BED peaks, a `.bed` path, or a vector of gene and transcript ids.

```r
kmer_enrichment(
  set_a = K562_bed,
  set_b = HepG2_bed,
  k     = 5
)
```

### Peak and transcript overlap, and control peaks

`intersect_peaks` is an R version of the common bedtools intersect call, and `generate_control_peaks` samples length, strand and region matched control peaks for a background set.

```r
hits <- intersect_peaks(K562_bed, gencode_v46_transcripts,
                        fraction = 1.0, same_strand = TRUE)

controls <- generate_control_peaks(
  raw_peaks   = K562_bed,
  anno        = gencode_v46_anno,
  gene        = gencode_v46_genes,
  transcripts = gencode_v46_transcripts
)
```

## Controlling the analysis and calculating significance

Most of the knobs for the splicing and sequence maps live in `splicing_options()`. The defaults are a reasonable starting point and match what we use for ENCODE eCLIP data.

Events get sorted into three groups before anything is plotted:

* Positive and Negative come from the significant events, split by the sign of `IncLevelDifference` using `event_fdr` and `psi_cutoff`.
* Control comes from the events that did not change, selected with `control_pval` and `psi_control_max`.

The Control curve is not just the raw pool. RNAPeaks draws a matched number of Control events many times over (`control_multiplier` and `control_iterations`) and uses the resampled mean and spread as the background.

Per position significance is a one sided hypergeometric (Fisher) test of the tested group against Control. You pick which Control numbers feed that test with `stat_test`:

* `"fisher-all"` (default) uses the whole raw Control pool.
* `"fisher-bootstrap"` uses the resampled Control mean, sized to the tested groups.

Set `use_fdr = TRUE` to apply Benjamini Hochberg correction across positions, and `fdr_threshold` to set the cutoff that decides which positions are drawn as significant.

## Bundled data

You get a handful of example objects so you can try things without downloading anything:

* `gtf_hg38`, `gtf_mm10`, `gtf_mm39`: GENCODE gene models used by the visual functions.
* `K562_bed`, `HepG2_bed`: example eCLIP peak tables in BED layout.
* `se_mats_jc`, `a5ss_mats_jc`, `a3ss_mats_jc`, `ri_mats_jc`: rMATS event tables, one per event type.
* `gencode_v46_anno`, `gencode_v46_genes`, `gencode_v46_transcripts`: annotation tables used by `generate_control_peaks` and `intersect_peaks`.

## Other options

A few settings you will reach for often:

* `width_exon` and `width_intron` control how many bases into the exon and the intron each region covers (defaults 50 and 250).
* `moving_average` smooths each curve with a window of this size (default 10, set to 0 to turn it off).
* `min_count` drops low coverage events before anything is scored (default 50, set to 0 to keep everything).
* `control_multiplier` sets how many Control events are drawn per iteration relative to the tested groups (default 2).

Styling for the maps lives in `splicing_style()`, and the peak plots use `peaks_options()` and `peaks_plot_style()`.

## Other notes

* The rMATS tables should be the junction count outputs (the SE.MATS, A5SS.MATS, A3SS.MATS and RI.MATS files). The bundled `*_mats_jc` objects are trimmed examples of exactly that shape.
* Sequence maps and kmer enrichment need a genome whose chromosome names line up with your events. RNAPeaks will add a `chr` prefix when it can and will warn about any chromosome it cannot place, so watch for that message if a map comes back empty.
* `generate_control_peaks` is happiest with the 8 column eCLIP peak format because it carries the fold change through, but it accepts any table with at least the 6 standard BED columns. The sampling runs through a small C++ engine and is deterministic for a fixed `seed`.
* Maps built from only a handful of events tend to look noisy. For ENCODE we like to have at least a hundred or so events in a group before reading much into the shape.

## License

MIT. See the LICENSE file.
