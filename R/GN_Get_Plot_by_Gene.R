Get_Plot_by_Gene<-function(bed,
                           gene,
                           rank_,
                           utr_col="dark gray",
                           peak_col="purple",
                           peaks_width=0.3,
                           margin = 100,
                           exon_width=0.2,
                           utr_width=0.3,
                           exon_col="navy",
                           total_arrows = 6,
                           max_per_intron = 2,
                           five_to_three = FALSE,
                           ...){

  #Return prepared bed with labeled coordinates positions for proteins
  bed <- Prepare_Bed(bed, rank_ = rank_, gene = gene,
                        peaks_width = peaks_width, peak_col = peak_col, margin = margin)

  # Build Gene Structure and return out positions for Gene, Exons, and UTRS
  gs <- Build_Gene_Structure(gene,
                             levels     = max(bed$rank),
                             peaks_width = peaks_width,
                             exon_width  = exon_width,
                             utr_width   = utr_width,
                             exon_col    = exon_col,
                             utr_col     = utr_col)

  #Stores labels positions for gene
  labs <- Make_Strand_Labels(gs$Gene_s, offset = 100)

  #Returns individual intron gaps which is used for positioning arrows
  introns <- Compute_Intron_Positions(gs$Gene_s, gs$Intron_s)

  #Gets arrow positions
  arrow_df <- make_intron_arrows(introns$Introns_Positions, gs$Gene_s$strand[1],
                                  total_arrows = total_arrows, max_per_intron = max_per_intron)

  # Plot for Gene using ggplot2
  do.call(
    Draw_Gene_Plot,
    c(
      list(
      Gene_s = gs$Gene_s,
      Exons  = gs$Exons,
      UTRs   = gs$UTRs,
      Intron_s = introns$Intron_s,
      Arrow_df = arrow_df,
      bed = bed,
      gene = gene,
      peak_col = bed$peak_col,
      peaks_width = bed$peaks_width,
      five_to_three = five_to_three,
      gene_strand = gene$strand[1]
    ),
    list(...)  # all user-supplied extras styling
    )
  )
}


#--------Helper Functions----------
#Prepare Background
Background_Table<-function(df,Start,End){
  Background<-df[c("y_start","y_end")]
  Background<-Background[!duplicated(Background),]
  Background<-Background[order(Background$y_start),]
  Background$x_start<-Start
  Background$x_end<-End
  Background$col<-rep(c("#33333333","#11111111"),length.out = nrow(Background))
  return(Background)
}



