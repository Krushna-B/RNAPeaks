
#' Human GTF Gene Annotation
#'
#' Pre-loaded Ensembl GTF annotation for human genes (GRCh38).
#' This bundled dataset eliminates the need to download from AnnotationHub,
#' enabling offline use and faster loading.
#'
#' @format A data frame containing GTF annotation with the following columns:
#' \describe{
#'   \item{seqnames}{Chromosome}
#'   \item{start}{Feature start position}
#'   \item{end}{Feature end position}
#'   \item{width}{Feature width in bp}
#'   \item{strand}{Strand (+ or -)}
#'   \item{source}{Annotation source}
#'   \item{type}{Feature type (gene, transcript, exon, CDS, UTR, etc.)}
#'   \item{score}{Annotation score}
#'   \item{phase}{CDS phase (0, 1, or 2)}
#'   \item{gene_id}{Ensembl gene ID}
#'   \item{gene_version}{Ensembl gene version}
#'   \item{gene_name}{Gene symbol}
#'   \item{gene_source}{Gene annotation source}
#'   \item{gene_biotype}{Gene biotype (protein_coding, lncRNA, etc.)}
#'   \item{transcript_id}{Ensembl transcript ID}
#'   \item{transcript_version}{Ensembl transcript version}
#'   \item{transcript_name}{Transcript name}
#'   \item{transcript_source}{Transcript annotation source}
#'   \item{transcript_biotype}{Transcript biotype}
#'   \item{transcript_support_level}{Transcript support level (TSL)}
#'   \item{exon_number}{Exon number within the transcript}
#'   \item{exon_id}{Ensembl exon ID}
#'   \item{exon_version}{Ensembl exon version}
#'   \item{protein_id}{Ensembl protein ID}
#'   \item{protein_version}{Ensembl protein version}
#'   \item{ccds_id}{CCDS identifier}
#'   \item{tag}{Feature tag (e.g., basic, Ensembl_canonical)}
#' }
#'
#' @source Ensembl via AnnotationHub (AH110867)
#'
#' @examples
#' \dontrun{
#'   data(gtf_human)
#'   PlotGene(bed = your_bed, geneID = "GAPDH", gtf = gtf_human)
#' }
"gtf_human"

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
#'   createSplicingMap(bed_file = your_bed, SEMATS = sample_se.mats)
#'   createSequenceMap(SEMATS = sample_se.mats, sequence = "CCCC")
#' }
"sample_se.mats"
