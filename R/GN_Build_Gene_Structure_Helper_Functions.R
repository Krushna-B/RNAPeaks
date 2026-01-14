Build_Gene_Structure <- function(gene, levels, peaks_width, exon_width, utr_width, exon_col, utr_col) {

  # Choose the appropriate gene-structure builder based on biotype
  Gene_s <- if (gene$gene_biotype[1] == "protein_coding") {
   Gene_strcuture_prot(gene, levels, peaks_width, exon_width, utr_width, exon_col, utr_col)

  } else {
    Gene_strcuture_nc(gene, levels, peaks_width, exon_width, utr_width, exon_col, utr_col)
  }

  if (is.null(Gene_s)) {
    return(NULL)
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


#-------Building Gene Structure Helper Functions-------------
#Gene Structure for non coding genes
Gene_strcuture_nc<-function(gtf,levels,peaks_width,exon_width,utr_width,exon_col,utr_col){
  gtf<-gtf[which(gtf$type %in% c("exon")),]
  exon_center<-levels*peaks_width+0.2+exon_width/2
  gtf$y_start<-exon_center-exon_width/2
  gtf$y_end<-exon_center+exon_width/2
  gtf$col<-exon_col
  gtf<-gtf[c("seqnames","type","start","end","strand","y_start","y_end","col")]
  Intron_s<-data.frame("seqnames"=gtf$seqnames[1],"type"="intron","start"=min(gtf$start),
                       "end"=max(gtf$end),"strand"=gtf$strand[1],
                       "y_start"=exon_center-utr_width/4,"y_end"=exon_center+utr_width/4,"col"="light gray")
  gtf <- rbind(gtf, Intron_s)
  return(gtf)
}


#Gene Structure for protein coding genes
Gene_strcuture_prot<-function(gtf,levels,peaks_width,exon_width,utr_width,exon_col,utr_col){

  if(!("CDS" %in% gtf$type) ){
    return(NULL)
  }

  gtf<-gtf[which(gtf$type %in% c("CDS","five_prime_utr","three_prime_utr")),]
  exon_center<-levels*peaks_width+0.2+exon_width/2

  gtf$y_start<-ifelse(gtf$type=="CDS",exon_center-exon_width/2,exon_center-utr_width/2)
  gtf$y_end<-ifelse(gtf$type=="CDS",exon_center+exon_width/2,exon_center+utr_width/2)
  gtf$col<-ifelse(gtf$type=="CDS",exon_col,utr_col)
  gtf<-gtf[c("seqnames","type","start","end","strand","y_start","y_end","col")]
  Intron_s<-data.frame("seqnames"=gtf$seqnames[1],"type"="intron","start"=min(gtf$start),
                       "end"=max(gtf$end),"strand"=gtf$strand[1],
                       "y_start"=exon_center-utr_width/4,"y_end"=exon_center+utr_width/4,"col"="light gray")
  gtf<-rbind(gtf,Intron_s)
  return(gtf)
}


#Function Returns strand labels based on positive or negative
Make_Strand_Labels <- function(Gene_s, offset = 100) {
  #Choose a vertical position for the strand tags
  ys <- min(Gene_s$y_start)

  if (Gene_s$strand[1] == "+") {
    #+ strand: transcription runs left -> right
    # Place 5′ a bit to the left of the gene, 3′ a bit to the right

    left  <- data.frame(Label="5'", X=min(Gene_s$start) - offset, Y=ys)
    right <- data.frame(Label="3'", X=max(Gene_s$end)   + offset, Y=ys)
  } else {
    # − strand: transcription runs right -> left (labels flip sides)
    # Place 3′ to the left and 5′ to the right of the gene.

    left  <- data.frame(Label="3'", X=min(Gene_s$start) - offset, Y=ys)
    right <- data.frame(Label="5'", X=max(Gene_s$end)   + offset, Y=ys)
  }

  list(left = left, right = right)
}


