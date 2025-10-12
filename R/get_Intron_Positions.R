get_Intron_Positions <- function(gene_df) {


    types <- c("intron","exon","CDS","five_prime_utr","three_prime_utr",
                 "five_prime_UTR","three_prime_UTR")
    gene_df <- gene_df[gene_df$type %in% types, ]

    gene_df <- gene_df[order(gene_df$start, gene_df$end),]

    print(gene_df)

    n <- nrow(gene_df)
    if (n < 2){
      return(NULL)
    }

    introns <- vector("list", n - 1)

       start <- gene_df$end[-n]
       end_prev <- gene_df$start[-1]
       gene_df$ystart
       len <- end_prev > start
       data.frame(
         seqnames = gene_df$seqnames[1],
         strand   = gene_df$strand[1],
         intron_id= seq_len(n - 1)[len],
         start    = start[len] + 1,
         end      = end_prev[len] - 1,
         length_bp= (end_prev[len] - start[len]),
         row.names = NULL,
         check.names = FALSE
       )


  }








