
#' K562 RBP Binding Peaks
#'
#' ENCODE eCLIP peak calls from the K562 (chronic myelogenous leukemia) cell
#' line, provided for testing and demonstration of RNAPeaks visualization and
#' analysis functions.
#'
#' @format A data frame with the following columns:
#' \describe{
#'   \item{chr}{Chromosome identifier (without "chr" prefix, e.g., "1", "X")}
#'   \item{start}{Peak start coordinate}
#'   \item{end}{Peak end coordinate}
#'   \item{tag}{RBP name or peak identifier}
#'   \item{score}{Peak score or confidence value}
#'   \item{strand}{Genomic strand ("+" or "-")}
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
#' @format A data frame with the following columns:
#' \describe{
#'   \item{chr}{Chromosome identifier (without "chr" prefix, e.g., "1", "X")}
#'   \item{start}{Peak start coordinate}
#'   \item{end}{Peak end coordinate}
#'   \item{tag}{RBP name or peak identifier}
#'   \item{score}{Peak score or confidence value}
#'   \item{strand}{Genomic strand ("+" or "-")}
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

#' Sample SE.MATS Skipped-Exon Splicing Events
#'
#' A dataset of skipped exon (SE) alternative splicing event provided for testing and demonstration of the splicing map and
#' sequence motif analysis functions in RNAPeaks.
#'
#' @format A data frame with 87,736 observations and 23 variables:
#' \describe{
#'   \item{ID}{Event ID}
#'   \item{GeneID}{Ensembl gene ID}
#'   \item{geneSymbol}{Gene symbol}
#'   \item{chr}{Chromosome}
#'   \item{strand}{Strand (+ or -)}
#'   \item{exonStart_0base}{Skipped exon start position (0-based)}
#'   \item{exonEnd}{Skipped exon end position}
#'   \item{upstreamES}{Upstream exon start position}
#'   \item{upstreamEE}{Upstream exon end position}
#'   \item{downstreamES}{Downstream exon start position}
#'   \item{downstreamEE}{Downstream exon end position}
#'   \item{ID.1}{Duplicate ID column}
#'   \item{IJC_SAMPLE_1}{Inclusion junction counts for sample 1}
#'   \item{SJC_SAMPLE_1}{Skipping junction counts for sample 1}
#'   \item{IJC_SAMPLE_2}{Inclusion junction counts for sample 2}
#'   \item{SJC_SAMPLE_2}{Skipping junction counts for sample 2}
#'   \item{IncFormLen}{Inclusion form length}
#'   \item{SkipFormLen}{Skipping form length}
#'   \item{PValue}{P-value for differential splicing}
#'   \item{FDR}{False discovery rate adjusted p-value}
#'   \item{IncLevel1}{Inclusion levels for sample 1}
#'   \item{IncLevel2}{Inclusion levels for sample 2}
#'   \item{IncLevelDifference}{Difference in inclusion levels between samples
#'     (positive = more inclusion in sample 1; negative = more skipping)}
#' }
#'
#' @source Generated with
#'
#' @examples
#' data(sample_se.mats)
#' head(sample_se.mats)
#'
#' \dontrun{
#'   skipped_exon_splicing_map(events = sample_se.mats, bed_file = your_bed)
#'   skipped_exon_sequence_map(events = sample_se.mats, sequence = "CCCC")
#' }
"sample_se.mats"
