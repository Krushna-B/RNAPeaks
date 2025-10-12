checkBed<-function(df){
  colnames(df)<-tolower(colnames(df))

  #Check if includes more than 3 most important cols
  stopifnot(ncol(df)>=3)


  col_names <- c("chr","start","end","tag","score","strand")
  k <- min(ncol(df), length(col_names))

  for(i in seq_len(k)){
    colnames(df)[i] <- col_names[i]
  }




  #if (!(all(c("chr","start","end","strand") %in% colnames(df)))) {

   # stop("Must Contain Chr,Start,End and Strand columns!")
  #}

  if(is.character(df$chr) &
             all(as.numeric(df$end)>=as.numeric(df$start)) &
             is.character(df$strand) &
             all(df$strand %in% c("+","-"))){

    #converting the version of the seqnames

    #df_gr<-GenomicRanges::makeGRangesFromDataFrame(df,keep.extra.columns = T)
   # GenomeInfoDb::seqlevelsStyle(df_gr)<-"NCBI"

   colchr_tmp = tolower(df$chr)
   df$chr = sub("chr","",colchr_tmp)
   return (df)
    #df<-data.frame(df_gr)
  } else {
    stop("Check Bed File!")
  }
}


