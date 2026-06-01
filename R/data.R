
#' K562 RBP Binding Peaks
#'
#' ENCODE eCLIP peak calls from the K562 (chronic myelogenous leukemia) cell
#' line, provided for testing and demonstration of RNAPeaks visualization and
#' analysis functions.
#'
#' @format A BED-style data frame with unnamed columns `V1`-`V8`:
#' \describe{
#'   \item{V1}{Chromosome identifier ("chr"-prefixed, e.g. "chr16")}
#'   \item{V2}{Peak start coordinate (0-based)}
#'   \item{V3}{Peak end coordinate}
#'   \item{V4}{RBP name}
#'   \item{V5}{Associated gene symbol}
#'   \item{V6}{Genomic strand ("+" or "-")}
#'   \item{V7}{Signal value (fold enrichment)}
#'   \item{V8}{Significance (p-value)}
#' }
#'
#' @source ENCODE Project (\url{https://www.encodeproject.org/})
#'
#' @examples
#' data(K562_bed)
#' head(K562_bed)
#'
#' \dontrun{
#'   data(K562_bed)
#'   plot_gene(bed = K562_bed, gene = "GAPDH", gtf = gtf_hg38)
#' }
"K562_bed"

#' HepG2 RBP Binding Peaks
#'
#' ENCODE eCLIP peak calls from the HepG2 (hepatocellular carcinoma) cell
#' line, provided for testing and demonstration of RNAPeaks visualization and
#' analysis functions.
#'
#' @format A BED-style data frame with unnamed columns `V1`-`V8`:
#' \describe{
#'   \item{V1}{Chromosome identifier ("chr"-prefixed, e.g. "chr16")}
#'   \item{V2}{Peak start coordinate (0-based)}
#'   \item{V3}{Peak end coordinate}
#'   \item{V4}{RBP name}
#'   \item{V5}{Associated gene symbol}
#'   \item{V6}{Genomic strand ("+" or "-")}
#'   \item{V7}{Signal value (fold enrichment)}
#'   \item{V8}{Significance (p-value)}
#' }
#'
#' @source ENCODE Project (\url{https://www.encodeproject.org/})
#'
#' @examples
#' data(HepG2_bed)
#' head(HepG2_bed)
#'
#' \dontrun{
#'   plot_gene(bed = HepG2_bed, gene = "GAPDH", gtf = gtf_hg38)
#' }
"HepG2_bed"

#' Human GTF Gene Annotation (GRCh38 / hg38)
#'
#' Pre-loaded GENCODE basic annotation for human genes (GRCh38, GENCODE v38).
#' This bundled dataset eliminates the need to download annotation at runtime.
#'
#' @format A data frame of GTF annotation rows with columns including
#'   `seqnames`, `start`, `end`, `width`, `strand`, `type`, `gene_id`,
#'   `gene_name`, `gene_biotype`, `transcript_id`, `transcript_name`.
#' @source GENCODE v38 basic annotation.
#' @examples
#' \dontrun{
#'   plot_gene(bed = your_bed, gene = "GAPDH", gtf = gtf_hg38)
#' }
"gtf_hg38"

#' Mouse GTF Gene Annotation (GRCm38 / mm10)
#'
#' Pre-loaded GENCODE basic annotation for mouse genes (GRCm38 assembly).
#'
#' @format See [gtf_hg38] for column layout.
#' @source GENCODE mouse basic annotation (mm10).
#' @examples
#' \dontrun{
#'   plot_gene(bed = your_bed, gene = "Gapdh", species = "mm10")
#' }
"gtf_mm10"

#' Mouse GTF Gene Annotation (GRCm39 / mm39)
#'
#' Pre-loaded GENCODE basic annotation for mouse genes (GRCm39 assembly).
#'
#' @format See [gtf_hg38] for column layout.
#' @source GENCODE mouse basic annotation (mm39).
#' @examples
#' \dontrun{
#'   plot_gene(bed = your_bed, gene = "Gapdh", species = "mm39")
#' }
"gtf_mm39"

#' Skipped-Exon Splicing Events (rMATS SE.MATS.JC)
#'
#' rMATS junction-count (JC) output for skipped-exon (SE) alternative splicing
#' events. Bundled for testing and demonstration of the splicing-map and
#' sequence-motif analysis functions.
#'
#' @format A data frame with 40,136 rows and 27 columns:
#' \describe{
#'   \item{ID}{Event ID}
#'   \item{GeneID}{Ensembl gene ID}
#'   \item{geneSymbol}{Gene symbol}
#'   \item{chr}{Chromosome}
#'   \item{strand}{Strand (`+` or `-`)}
#'   \item{exonStart_0base}{Skipped exon start (0-based)}
#'   \item{exonEnd}{Skipped exon end}
#'   \item{upstreamES, upstreamEE}{Upstream exon start / end}
#'   \item{downstreamES, downstreamEE}{Downstream exon start / end}
#'   \item{ID.1}{Duplicate ID column emitted by rMATS}
#'   \item{IJC_SAMPLE_1, SJC_SAMPLE_1}{Inclusion / skipping junction counts,
#'     sample 1}
#'   \item{IJC_SAMPLE_2, SJC_SAMPLE_2}{Inclusion / skipping junction counts,
#'     sample 2}
#'   \item{IncFormLen, SkipFormLen}{Effective inclusion / skipping form length}
#'   \item{PValue, FDR}{Differential-splicing p-value and BH-adjusted FDR}
#'   \item{IncLevel1, IncLevel2}{Per-replicate inclusion levels (PSI),
#'     comma-separated}
#'   \item{IncLevelDifference}{Mean \eqn{\Delta\Psi} (sample 1 - sample 2)}
#'   \item{upstream_to_target_count, target_to_downstream_count, target_count,
#'     upstream_to_downstream_count}{rMATS JC-specific junction counts}
#' }
#'
#' @source rMATS turbo (\url{https://github.com/Xinglab/rmats-turbo}).
#'
#' @examples
#' data(se_mats_jc)
#' head(se_mats_jc)
#'
#' \dontrun{
#'   skipped_exon_splicing_map(events = se_mats_jc, bed_file = your_bed)
#'   skipped_exon_sequence_map(events = se_mats_jc, sequence = "CCCC")
#' }
"se_mats_jc"

#' Alternative 3' Splice-Site Events (rMATS A3SS.MATS.JC)
#'
#' rMATS junction-count (JC) output for alternative 3' splice-site (A3SS)
#' events.
#'
#' @format A data frame with 6,365 rows and 27 columns. Shares the universal
#'   rMATS columns described in [se_mats_jc] (`ID`, `GeneID`, `geneSymbol`,
#'   `chr`, `strand`, `IJC_*`, `SJC_*`, `IncFormLen`, `SkipFormLen`, `PValue`,
#'   `FDR`, `IncLevel1`, `IncLevel2`, `IncLevelDifference`). A3SS-specific
#'   coordinate columns:
#' \describe{
#'   \item{longExonStart_0base, longExonEnd}{Long exon form start / end (0-based start)}
#'   \item{shortES, shortEE}{Short exon form start / end}
#'   \item{flankingES, flankingEE}{Flanking exon start / end}
#'   \item{across_short_boundary_count, long_to_flanking_count,
#'     exclusive_to_long_count, short_to_flanking_count}{rMATS JC-specific
#'     junction counts}
#' }
#'
#' @source rMATS turbo (\url{https://github.com/Xinglab/rmats-turbo}).
#'
#' @examples
#' data(a3ss_mats_jc)
#' head(a3ss_mats_jc)
"a3ss_mats_jc"

#' Alternative 5' Splice-Site Events (rMATS A5SS.MATS.JC)
#'
#' rMATS junction-count (JC) output for alternative 5' splice-site (A5SS)
#' events.
#'
#' @format A data frame with 4,282 rows and 27 columns. Same column layout as
#'   [a3ss_mats_jc] (universal rMATS columns plus `longExonStart_0base`,
#'   `longExonEnd`, `shortES`, `shortEE`, `flankingES`, `flankingEE`, and the
#'   four JC-specific junction count columns).
#'
#' @source rMATS turbo (\url{https://github.com/Xinglab/rmats-turbo}).
#'
#' @examples
#' data(a5ss_mats_jc)
#' head(a5ss_mats_jc)
"a5ss_mats_jc"

#' Retained-Intron Events (rMATS RI.MATS.JC)
#'
#' rMATS junction-count (JC) output for retained-intron (RI) events.
#'
#' @format A data frame with 5,283 rows and 27 columns. Shares the universal
#'   rMATS columns described in [se_mats_jc]. RI-specific coordinate columns:
#' \describe{
#'   \item{riExonStart_0base, riExonEnd}{Retained-intron-containing exon start
#'     / end (0-based start)}
#'   \item{upstreamES, upstreamEE}{Upstream exon start / end}
#'   \item{downstreamES, downstreamEE}{Downstream exon start / end}
#'   \item{upstream_to_intron_count, intron_to_downstream_count, intron_count,
#'     upstream_to_downstream_count}{rMATS JC-specific junction counts}
#' }
#'
#' @source rMATS turbo (\url{https://github.com/Xinglab/rmats-turbo}).
#'
#' @examples
#' data(ri_mats_jc)
#' head(ri_mats_jc)
#'
#' \dontrun{
#'   retained_intron_splicing_map(events = ri_mats_jc, bed_file = your_bed)
#'   retained_intron_sequence_map(events = ri_mats_jc, sequence = "CCCC")
#' }
"ri_mats_jc"

#' GENCODE v46 Transcript BED (human / GRCh38)
#'
#' One-row-per-transcript BED extracted from the GENCODE v46 primary-assembly
#' GTF. Used as the `transcripts` input to [intersect_peaks()] when running
#' the eCLIP control-peak pipeline against human GRCh38.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{chr}{Chromosome identifier with "chr" prefix (e.g. "chr1", "chrX").}
#'   \item{start}{Transcript start coordinate (0-based BED).}
#'   \item{end}{Transcript end coordinate (exclusive).}
#'   \item{name}{GENCODE transcript name (e.g. "DDX11L2-202").}
#'   \item{score}{BED score column; placeholder "." for this dataset.}
#'   \item{strand}{Genomic strand ("+" or "-").}
#' }
#'
#' @source GENCODE release 46, primary assembly annotation
#'   (\url{https://www.gencodegenes.org/human/release_46.html}).
#' @examples
#' \dontrun{
#'   hits <- intersect_peaks(
#'     peaks       = "K562_RBFOX2_peaks.bed",
#'     transcripts = gencode_v46_transcripts
#'   )
#' }
"gencode_v46_transcripts"

#' GENCODE v46 Gene BED (human / GRCh38)
#'
#' One-row-per-gene BED extracted from the GENCODE v46 primary-assembly GTF.
#' Used as the `genes` input to [generate_control_peaks()] to constrain
#' control regions to lie inside the same gene as the peak.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{chr}{Chromosome identifier with "chr" prefix.}
#'   \item{start}{Gene start coordinate (0-based BED).}
#'   \item{end}{Gene end coordinate (exclusive).}
#'   \item{name}{GENCODE gene symbol (e.g. "DDX11L1").}
#'   \item{score}{BED score column; placeholder "." for this dataset.}
#'   \item{strand}{Genomic strand ("+" or "-").}
#' }
#'
#' @source GENCODE release 46, primary assembly annotation.
#' @examples
#' \dontrun{
#'   generate_control_peaks(
#'     peaks      = hits,
#'     annotation = gencode_v46_anno,
#'     genes      = gencode_v46_genes,
#'     output_dir = "result"
#'   )
#' }
"gencode_v46_genes"

#' GENCODE v46 Per-Region Annotation BED (human / GRCh38)
#'
#' Per-region annotation BED derived from the GENCODE v46 GTF: one row per
#' (transcript, region) where region is one of `UTR3`, `UTR5`, `CDS`, or
#' `exon`. Introns are not stored; they are derived at runtime from the
#' exon boundaries inside [generate_control_peaks()]. Used as the
#' `annotation` input to [generate_control_peaks()].
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{chr}{Chromosome identifier with "chr" prefix.}
#'   \item{start}{Region start coordinate (0-based BED).}
#'   \item{end}{Region end coordinate (exclusive).}
#'   \item{name}{Composite `{transcript}_{region}_{key}` string, e.g.
#'     `DDX11L2-202_exon_11869`. The control-peak algorithm splits this on
#'     `_` to recover the transcript name and region type.}
#'   \item{score}{BED score column; "0" for this dataset.}
#'   \item{strand}{Genomic strand ("+" or "-").}
#'   \item{gene}{Parent gene symbol (e.g. "DDX11L2").}
#' }
#'
#' @source Derived from GENCODE release 46 (primary assembly) via the
#'   annotation generation pipeline in
#'   \url{https://github.com/jechia/control_peak} (Yue Hu, 2025).
#' @examples
#' \dontrun{
#'   generate_control_peaks(
#'     peaks      = hits,
#'     annotation = gencode_v46_anno,
#'     genes      = gencode_v46_genes,
#'     output_dir = "result"
#'   )
#' }
"gencode_v46_anno"
