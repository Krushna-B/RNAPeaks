Build_Gene_Structure <- function(gene, levels, peaks_width, exon_width, utr_width, exon_col, utr_col) {

  # Choose the appropriate gene-structure builder based on biotype
  Gene_s <- if (gene$gene_biotype[1] == "protein_coding") {
    Gene_strcuture_prot(gene, levels, peaks_width, exon_width, utr_width, exon_col, utr_col)
  } else {
    Gene_strcuture_nc(gene, levels, peaks_width, exon_width, utr_width, exon_col, utr_col)
  }

  # Separate introns from the rest of the features
  Intron_s <- Gene_s[Gene_s$type == "intron", ]
  Gene_s   <- Gene_s[Gene_s$type != "intron", ]

  # Convenience subsets for plotting layers
  if(gene$gene_biotype[1] == "protein_coding"){
    Exons    <- subset(Gene_s, type == "CDS") # coding sequence only
  } else{
    Exons    <- subset(Gene_s, type == "exon")
    }
  UTRs     <- subset(Gene_s, type %in% c("five_prime_utr", "three_prime_utr"))  # both UTR ends

  # Hand back the pieces needed by the plotting code
  list(Gene_s = Gene_s, Intron_s = Intron_s, Exons = Exons, UTRs = UTRs)
}
