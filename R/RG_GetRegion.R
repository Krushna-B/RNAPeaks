
getRegion <- function(gtf,
                      Chr,
                      Start,
                      End,
                      geneID = NULL,
                      TxID = NULL,
                      species = "Human"
                      ) {

  # If there is gene or transcript, delegate to your existing GetGene()
  if (!is.null(geneID) || (!is.null(TxID) && !is.na(TxID))) {
    rows <- GetGene(
      geneID  = geneID,
      species = species,
      TxID    = if (is.null(TxID)) NA else TxID,
      gtf     = gtf
    )
    return(rows)
}

    #If GTF is not provided then load it from from AnnotationHub by species
    if (is.null(gtf) & species=="Human"){
      #AnnotationHub::query(ah, pattern = c("Homo Sapiens","gtf"))
      ah<-AnnotationHub::AnnotationHub()
      gtf<-data.frame(ah[["AH110867"]])
    } else if (is.null(gtf) & species=="Mouse"){
      #query(ah, pattern = c("Mus musculus", "gtf"))[[1]]
      ah<-AnnotationHub::AnnotationHub()
      gtf<-data.frame(ah[["AH47076"]])
    }


    #Keep any features overlapping the [Start, End] window on Chr
    rows <- subset(gtf, seqnames == Chr & end >= Start & start <= End)

    # nothing overlaps; return empty frame
    if (!nrow(rows)) return(rows)


    #Transcripts that exist in the region
    transcripts <- rows[rows$type == "transcript", c("gene_id","transcript_id","width","gene_name")]

    #If there are transcripts then keep the max transcripts for each gene
    if(nrow(transcripts)){
      transcripts_keep <- tapply(seq_len(nrow(transcripts)), transcripts$gene_id, function(idx){
        transcripts$transcript_id[idx[ which.max(transcripts$width[idx]) ]]
        })
    }

    rows <-rows[rows$transcript_id %in% transcripts_keep, ,drop = FALSE]
    transcripts <<- rows
    rows
  }


