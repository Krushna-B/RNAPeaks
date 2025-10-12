
#Gene Structure-protein coding genes
Gene_strcuture_prot<-function(gtf,levels,peaks_width,exon_width,utr_width,exon_col,utr_col){
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

