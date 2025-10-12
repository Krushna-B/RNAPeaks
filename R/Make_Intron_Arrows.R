make_intron_arrows <- function(
    Introns,
    strand = "+",
    total_arrows   = 12,   # target across ALL introns
    min_per_intron = 0,
    max_per_intron = 5,
    margin_buff  = 0,   # keep arrow tips away from intron ends
    seg_frac       = 0.5,  # arrow shaft ≈ this fraction of the spacing
    seg_min_bp     = 20,
    seg_max_bp     = 180
){
  #Sets Introns to dataframe and if that throws error then sets to NULL
  Introns <- tryCatch(as.data.frame(Introns), error = function(e) NULL)
  # Always return the columns ggplot needs:
  empty <- data.frame(x=numeric(0), xend=numeric(0), y=numeric(0), yend=numeric(0))

  #Checks if Introns are emtpy and exits
  if (is.null(Introns) || nrow(Introns) == 0){
    return(empty)
  }

  # Pick a y-line coordinates to draw on
  y_line <- if ("y" %in% names(Introns)){
    Introns$y
    } else if ("y_mid" %in% names(Introns)){
      Introns$y_mid
    }else if (all(c("y_start","y_end") %in% names(Introns))){
      (Introns$y_start + Introns$y_end)/2
    } else{
      0
    }

  #Raw Intron length
  len_raw <- pmax(0, Introns$end - Introns$start)


  #Find usebale length subracting margins buffers on both sides
  usable  <- pmax(0, len_raw - 2*margin_buff)

  #Total Space to draw arrows
  tot_use <- sum(usable)

  if(tot_use <400){
    total_arrows <-6
  }
  #If there is no space total introns inputted is 0 return empty
  if (tot_use == 0 || total_arrows <= 0){
    return(empty)
  }


  bp_per_arrow <- tot_use / total_arrows
  #Gets proptional number
  number_of_introns <- floor(usable / bp_per_arrow)

  #Checks that number of introsn is within bounds
  number_of_introns <- pmax(min_per_intron, pmin(max_per_intron, number_of_introns))
  if (sum(number_of_introns) == 0) {  # ensure at least 1 if something is drawable
    j <- which.max(usable)
    if (usable[j] > 0){
      number_of_introns[j] <- 1
    }
  }

  #Directionality for each arrow
  dir <- if (strand == "+"){
    1
  } else{
    -1
  }

  out <- lapply(seq_len(nrow(Introns)), function(i){
    k <- number_of_introns[i]

    #Skips arrows with 0
    if (k <= 0 || usable[i] <= 0) return(NULL)

    s2 <- Introns$start[i] + margin_buff
    e2 <- Introns$end[i]   - margin_buff
    #Divides into gaps
    spacing <- usable[i] / (k + 1)


    seg_bp  <- min(seg_max_bp, max(seg_min_bp, seg_frac * spacing))

    #X coordinates for it
    centers <- s2 + spacing * seq_len(k)
    x    <- centers - (seg_bp/2) * dir
    xend <- centers + (seg_bp/2) * dir

    # Clamp inside the intron and drop zero-length segments
    x    <- pmax(s2, pmin(e2, x))
    xend <- pmax(s2, pmin(e2, xend))

    #Error Checks for keeping
    keep <- (xend - x) * dir > 0
    if (!any(keep)) return(NULL)

    data.frame(x=x[keep], xend=xend[keep], y=y_line[i], yend=y_line[i], row.names=NULL)
  })

  out <- do.call(rbind, out)
  if (is.null(out)){
    empty
    } else{
      out
    }
}
