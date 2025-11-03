
PlotRegion<-function(Chr=NULL,
                     Start=NULL,
                     End=NULL,
                     Strand=NULL,
                     geneID=NULL,
                     gtf=NULL,
                     bed=NULL,
                     species="Human",
                     TxID=NA,
                     Target_col=NULL,
                     omit=c(),
                     order_by="Count",
                     order_in=NULL,
                     merge=0,
                     peaks_width=0.3,
                     utr_col="dark gray",
                     peak_col="Blue",
                     exon_width=0.5,
                     utr_width=0.3,
                     exon_col="black",
                     ...){

  if (is.null(Chr) | is.null(Start) | is.null(End) | is.null(Strand)){
    stop("Need to provide Chr, Start, Stop and Strand parameters to visualize region.")
  }

  clipped_region <- Build_Region_Structure(
    gtf   = gtf,
    Chr   = Chr,
    Start = Start,
    End = End,
    Strand = Strand,
    geneID = geneID,
    TxID = TxID
  )
  clipped_region <- clipped_region[order(clipped_region$start, clipped_region$end), ]


  #check the bed file
  colnames(bed)[which(colnames(bed)==Target_col)]<-"target"
  bed<-checkBed(bed)

  #Filter bed
  bed<-FilterBed(bed=bed,chr=clipped_region$seqnames[1],start=min(clipped_region$start),end=max(clipped_region$end),strand=clipped_region$strand[1],omit=omit,collapse=merge)
  print(bed)

  #Option to give an order to rank proteins in
  if (!is.null(order_in)){
    Target_rank<-order_in
  } else {
    Target_rank<-OrderPeak(bed=bed,order_by=order_by)
  }

 Plot <- Get_Multi_Plot_by_Region(
    gtf   = gtf,
    bed = bed,
    Chr   = Chr,
    Start = Start,
    End = End,
    region_with_multiple_genes = clipped_region,
    Strand = Strand,
    geneID = geneID,
    rank_ = Target_rank,
    gene = clipped_region,
    TxID = TxID,
    peaks_width = peaks_width,
    exon_width = exon_width,
    utr_width = utr_width,
    exon_col = exon_col,
    utr_col = utr_col

  )

  ggsave("~/Desktop/RNAPeaks.pdf",Plot,height=12,width=16)
  return(Plot)
}

