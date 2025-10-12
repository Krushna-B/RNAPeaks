#Gene Structure-non coding genes
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
  exon_test<<- gtf
  gtf<-rbind(gtf,Intron_s)
  return(gtf)
}
