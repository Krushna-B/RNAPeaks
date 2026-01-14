 createSplicingMap <- function(bed_File = NULL,
              Protein = NULL,
             SEMATS = NULL,
             moving_average = 40,
             p_valueRetainedAndExclusion = 0.05,
             p_valueControls = 0.08,
            retained_IncLevelDifference = 0.1,
            exclusion_IncLevelDifference = -0.1,
            WidthIntoExon = 50,
            WidthIntoIntron = 250
             ){



#Separate out Retained, excluded, controls exon
Retained <- SEMATS %>%
  dplyr::filter(PValue < p_valueRetainedAndExclusion,
         FDR < p_valueRetainedAndExclusion,
         IncLevelDifference > retained_IncLevelDifference)

Exclusion <- SEMATS %>%
  dplyr::filter(PValue < p_valueRetainedAndExclusion,
         FDR < p_valueRetainedAndExclusion,
         IncLevelDifference < exclusion_IncLevelDifference)

Controls <- SEMATS %>% dplyr::filter(PValue > p_valueControls,abs(IncLevelDifference) <0.001)

# Retained <- Retained %>%
#   dplyr::filter((exonStart_0base - upstreamEE > 500),
#                 downstreamES - exonEnd > 500)
#
# Exclusion <- Exclusion %>%
#   dplyr::filter((exonStart_0base - upstreamEE > 500),
#                 downstreamES - exonEnd > 500)
#
# Controls <- Controls %>%
#   dplyr::filter((exonStart_0base - upstreamEE > 500),
#                 downstreamES - exonEnd > 500)


Retained_Bins_R <- make_bins_matrix(Retained,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)
Retained_Bins_E <- make_bins_matrix(Exclusion,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)
Exclusion_Bins_C <- make_bins_matrix(Controls,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)

#Make it by order and return error if incorrect file
buckets <- makeGRangesFromDataFrame(BEDFILE,
                                    seqnames.field="V1",
                                    start.field="V2",
                                    end.field="V3",
                                    strand.field="V6",
                                    keep.extra.columns=TRUE)

overlap_matrix_R <- find_overlaps(Retained_Bins_R, buckets)
overlap_matrix_E <- find_overlaps(Retained_Bins_E, buckets)
overlap_matrix_C <- find_overlaps(Exclusion_Bins_C, buckets)

total_events_Retained<- length(unique(Retained_Bins_R$event_id))
total_events_Excluded <- length(unique(Retained_Bins_E$event_id))
total_events_Controls <- length(unique(Exclusion_Bins_C$event_id))


freq_data_Retained<- calculate_overlap_frequency(overlap_matrix_R,total_events_Retained,WidthIntoExon+WidthIntoIntron+1)
freq_data_Excluded <- calculate_overlap_frequency(overlap_matrix_E,total_events_Excluded,WidthIntoExon+WidthIntoIntron+1)
freq_data_Controls <- calculate_overlap_frequency(overlap_matrix_C,total_events_Controls,WidthIntoExon+WidthIntoIntron+1)

moving_average_data_R <- calculate_moving_average(freq_data_Retained, moving_average)
moving_average_data_E <- calculate_moving_average(freq_data_Excluded, moving_average)
moving_average_data_C <- calculate_moving_average(freq_data_Controls, moving_average)

moving_average_data_R$group <- "R"
moving_average_data_E$group <- "E"
moving_average_data_C$group <- "C"

df <- rbind(moving_average_data_R,moving_average_data_E,moving_average_data_C)



ggplot(df,aes(x=global_position, y=moving_avg, color = group))+
  geom_line()+
  geom_vline(xintercept =c(51,552,653,1154),linetype="dashed")


plot_splicing_map(moving_average_data_R)

#
# Retained_P <- Retained %>%
#   dplyr::filter(strand == '+')
#
# Retained_N <- Retained %>%
#   dplyr::filter(strand == '-')
#
# Exclusion_P <- Exclusion %>%
#   dplyr::filter(strand == '+')
#
# Exclusion_N <- Exclusion %>%
#   dplyr::filter(strand == '-')
#
# Controls_P <- Controls %>%
#   dplyr::filter(strand == '+')
#
# Controls_N <- Controls %>%
#   dplyr::filter(strand == '-')
#Picking Controls
#Do later Separate Function

#Build Matrices for all 3 situations (Beware of Stands)
# Retained_Bins_P <- make_bins_matrix(Retained_P,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)
# Retained_Bins_N <- make_bins_matrix(Retained_N,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)
#
# Exclusion_Bins_P <- make_bins_matrix(Exclusion_P,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)
# Exclusion_Bins_N <- make_bins_matrix(Exclusion_N,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)
#
# Controls_Bins_P <- make_bins_matrix(Controls_P,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)
# Controls_Bins_N <- make_bins_matrix(Controls_N,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)


# #Holds regions that are overlaps
# overlap_matrix_RP <- find_overlaps(Retained_Bins_P, buckets)
# overlap_matrix_RN <- find_overlaps(Retained_Bins_N, buckets)
# overlap_matrix_EP <- find_overlaps(Exclusion_Bins_P, buckets)
# overlap_matrix_EN <- find_overlaps(Exclusion_Bins_N, buckets)
# overlap_matrix_CP <- find_overlaps(Controls_Bins_P, buckets)
# overlap_matrix_CN <- find_overlaps(Controls_Bins_N, buckets)
#
#
# Retained_Overlaps <- rbind(overlap_matrix_RP,overlap_matrix_RN)
# Exclusion_Overlaps <- rbind(overlap_matrix_EP,overlap_matrix_EN)
# Controls_Overlaps <- rbind(overlap_matrix_CP,overlap_matrix_CN)





#Total Number of events in original retained binds
# total_events_RP <- length(unique(Retained_Bins_P$event_id))
# total_events_RN <- length(unique(Retained_Bins_N$event_id))
# total_events_EP <- length(unique(Exclusion_Bins_P$event_id))
# total_events_EN <- length(unique(Exclusion_Bins_N$event_id))
# total_events_CP <- length(unique(Controls_Bins_P$event_id))
# total_events_CN <- length(unique(Controls_Bins_N$event_id))


#Frequency DataFrame
# freq_data_RP <- calculate_overlap_frequency(overlap_matrix_RP,total_events_RP,WidthIntoExon+WidthIntoIntron+1)
# freq_data_RN <- calculate_overlap_frequency(overlap_matrix_RN,total_events_RN,WidthIntoExon+WidthIntoIntron+1)
# freq_data_EP <- calculate_overlap_frequency(overlap_matrix_EP,total_events_EP,WidthIntoExon+WidthIntoIntron+1)
# freq_data_EN <- calculate_overlap_frequency(overlap_matrix_EN,total_events_EN,WidthIntoExon+WidthIntoIntron+1)
# freq_data_CP <- calculate_overlap_frequency(overlap_matrix_CP,total_events_CP,WidthIntoExon+WidthIntoIntron+1)
# freq_data_CN <- calculate_overlap_frequency(overlap_matrix_CN,total_events_CN,WidthIntoExon+WidthIntoIntron+1)

# freq_data_Retained<- calculate_overlap_frequency(Retained_Overlaps,total_events_Retained,WidthIntoExon+WidthIntoIntron+1)
# freq_data_Excluded <- calculate_overlap_frequency(Exclusion_Overlaps,total_events_Excluded,WidthIntoExon+WidthIntoIntron+1)
# freq_data_Controls <- calculate_overlap_frequency(Controls_Overlaps,total_events_Controls,WidthIntoExon+WidthIntoIntron+1)
#
# moving_average_data_R <- calculate_moving_average(freq_data_Retained, moving_average)
# moving_average_data_E <- calculate_moving_average(freq_data_Excluded, moving_average)
# moving_average_data_C <- calculate_moving_average(freq_data_Controls, moving_average)
#
#
# moving_average_data_RP <- calculate_moving_average(freq_data_RP, moving_average)
# moving_average_data_RN <- calculate_moving_average(freq_data_RN, moving_average)
# moving_average_data_EP <- calculate_moving_average(freq_data_EP, moving_average)
# moving_average_data_EN <- calculate_moving_average(freq_data_EN, moving_average)
# moving_average_data_CP <- calculate_moving_average(freq_data_CP, moving_average)
# moving_average_data_CN <- calculate_moving_average(freq_data_CN, moving_average)
#
#
# moving_average_data_RP$group <- "RP"
# moving_average_data_RN$group <- "RN"
# moving_average_data_EP$group <- "EP"
# moving_average_data_EN$group <- "EN"
# moving_average_data_CP$group <- "CP"
# moving_average_data_CN$group <- "CN"
#
# df <- rbind(moving_average_data_RP,moving_average_data_RN,moving_average_data_EP,moving_average_data_EN,moving_average_data_CP,moving_average_data_CN)
#
# ggplot(df,aes(x=global_position, y=moving_avg, color = group))+
#   geom_line()



return(plot_splicing_map(df))


}



 #-----HELPER FUNCTIONS --------

 make_bins_matrix <- function(SEMATS,WidthIntoExon,WidthIntoIntron){
    n <- nrow(SEMATS)
   strand <- SEMATS$strand

   # Normalized exon coordinates (always start < end)
   exonStart <- pmin(SEMATS$exonStart_0base, SEMATS$exonEnd)
   exonEnd   <- pmax(SEMATS$exonStart_0base, SEMATS$exonEnd)
   upS <- pmin(SEMATS$upstreamES, SEMATS$upstreamEE)
   upE <- pmax(SEMATS$upstreamES, SEMATS$upstreamEE)
   downS <- pmin(SEMATS$downstreamES, SEMATS$downstreamEE)
   downE <- pmax(SEMATS$downstreamES, SEMATS$downstreamEE)

   #Bin1
   width_into_exon = pmin(WidthIntoExon,upE-upS)

   bin1_start <- upE - width_into_exon

   width_into_intron = pmin(WidthIntoIntron,exonStart-upE)
   bin1_end <- upE +  width_into_intron

   #Bin 2
   width_into_intron = pmin(WidthIntoIntron,exonStart-upE)
   bin2_start <- exonStart -  width_into_intron


   width_into_exon = pmin(WidthIntoExon,exonEnd-exonStart)
   bin2_end <- exonStart + width_into_exon


   #Bin 3
   width_into_exon = pmin(WidthIntoExon,exonEnd-exonStart)
   bin3_start <- exonEnd - width_into_exon

   width_into_intron = pmin(WidthIntoIntron,downS-exonEnd)
   bin3_end <- exonEnd + width_into_intron

   #Bin 4
   width_into_intron = pmin(WidthIntoIntron,downS-exonEnd)
   bin4_start <- downS - width_into_intron


   width_into_exon = pmin(WidthIntoExon,downE-downS)
   bin4_end <- downS  + width_into_exon

   starts <- cbind(bin1_start, bin2_start, bin3_start, bin4_start)
   ends   <- cbind(bin1_end,   bin2_end,   bin3_end,   bin4_end)

    #Build GRanges
    GRanges(
     seqnames = rep(SEMATS$chr, each = 4),
     ranges   = IRanges(start = as.vector(t(starts)),
                        end   = as.vector(t(ends))),
     strand   = rep(SEMATS$strand, each = 4),
     geneID = rep(SEMATS$GeneID, each = 4),
     event_id = rep(seq_len(nrow(SEMATS)), each = 4)
   )
 }



 find_overlaps <- function(bins, protein){
   # Pre-filter - only check bins that actually overlap protein ranges
   overlapping_bins <- subsetByOverlaps(bins, protein)

   if(length(overlapping_bins) == 0){
     return(data.frame(
       event_id = integer(),
       bin_index = integer(),
       position = integer(),
       has_overlap = logical()
     ))
   }

   results <- lapply(seq_along(overlapping_bins), function(i) {
     bin <- overlapping_bins[i]

     chr <- as.character(seqnames(bin))

     chr_protein <- protein[seqnames(protein) == chr]

     positions <- start(bin):end(bin)

     #Every 1 bp position
     pos_gr <- GRanges(chr, IRanges(positions, width = 1), strand = strand(bin))

     overlaps <- countOverlaps(pos_gr, chr_protein, ignore.strand = TRUE) > 0
     data.frame(
       event_id = mcols(bin)$event_id,
       geneID = mcols(bin)$geneID,
       seqnames = chr,
       strand = as.character(strand(bin)),
       bin_index = ((i - 1) %% 4) + 1,  # 1-4 for each event
       position = seq_along(positions),
       genomic_position = positions,
       has_overlap = overlaps
     )
   })
     do.call(rbind, results)
   }




calculate_overlap_frequency <- function(overlap_df, total_events,bin_width) {
  # Calculate global position (1-1200) from bin_index and position
  overlap_df$global_position <- (overlap_df$bin_index - 1) * bin_width + overlap_df$position

#Flips Global Position
overlap_df <- overlap_df %>%
  group_by(event_id) %>%
  mutate(
    global_position = (bin_index - 1) * bin_width + position,
    global_position = ifelse(
      strand == "-",
      (4 * bin_width) - global_position + 1,
      global_position
    )
  ) %>%
  ungroup()

  # OVERLAPS <<- overlap_df

  # Count overlaps at each global position
  overlap_counts <- aggregate(
    has_overlap ~ global_position,
    data = overlap_df,
    FUN = sum
  )

  # Calculate frequency (overlaps / TOTAL events in dataset)
  overlap_counts$frequency <- overlap_counts$has_overlap / total_events

  # Positions 4 * bin_width are represented (fill missing with 0)
  all_positions <- data.frame(global_position = 1:(4*bin_width))
  result <- merge(all_positions, overlap_counts[, c("global_position", "frequency")],
                  by = "global_position", all.x = TRUE)
  result$frequency[is.na(result$frequency)] <- 0

  return(result)
}


calculate_moving_average <- function(freq_data, window_size = NULL) {

  # Add bin column
  freq_data <- freq_data %>%
    mutate(
      bin = case_when(
        global_position <= 301 ~ 1,
        global_position <= 602 ~ 2,
        global_position <= 903 ~ 3,
        TRUE ~ 4
      )
    ) %>%
    arrange(global_position)

  # Apply moving average only if window_size is provided
  if (!is.null(window_size) && window_size > 0) {

    freq_data <- freq_data %>%
      group_by(bin) %>%
      mutate(moving_avg = zoo::rollmean(frequency, k = window_size, fill = NA, align = "center")) %>%
      ungroup()
  } else {
    # No moving average - just copy frequency to moving_avg
    freq_data <- freq_data %>%
      mutate(moving_avg = frequency)
  }

  return(freq_data)
}



plot_splicing_map <- function(moving_average_data,
                              WidthIntoExon = 50,
                              WidthIntoIntron = 250,
                              title = "Protein Binding Overlap Frequency",
                              line_color = "red"){
  bin_width <- WidthIntoExon + WidthIntoIntron
  gap <- 80

  #Height for exons
  y_max <- max(moving_average_data$moving_avg, na.rm = TRUE)
  y_min <- min(moving_average_data$moving_avg, na.rm = TRUE)
  y_range <- y_max - y_min
  exon_height <- y_range * 0.08

  # Calculate region starts
  region_starts <- c(0, bin_width + gap, 2*bin_width + gap, 3*bin_width + 2*gap)


  boundary1 <- region_starts[1] + WidthIntoExon
  boundary2 <- region_starts[2] + WidthIntoIntron
  boundary3 <- region_starts[3] + WidthIntoExon
  boundary4 <- region_starts[4] + WidthIntoIntron

  #Dashed boundary lines locations
  boundary_lines <- data.frame(
    xintercept = c(boundary1, boundary2, boundary3, boundary4)
  )

  exon_regions <- data.frame(
    xmin = c(region_starts[1], boundary2, boundary4),
    xmax = c(boundary1, boundary3, region_starts[4] + bin_width),
    ymin = rep(y_min - exon_height, 3),
    ymax = rep(y_min, 3),
    fill = c("white", "navy", "white")
  )

  # Intron lines
  intron_y <- y_min - exon_height/2
  intron_segments <- data.frame(
    x = c(boundary1,boundary3),
    xend = c(boundary2, boundary4),

    y = rep(intron_y, 2),
    yend = rep(intron_y, 2),
    linetype = c("solid", "solid")
  )

  # // Breaks
  break_x <- c(bin_width + gap/2, 3*bin_width + gap + gap/2)

  #Aligns Data to proper region
  moving_average_data <- moving_average_data %>%
    mutate(
      position_in_bin = global_position - (bin - 1) * bin_width,
      schematic_position = case_when(
        bin == 1 ~ position_in_bin,
        bin == 2 ~ bin_width + gap + position_in_bin,

        bin == 3 ~ 2*bin_width + gap + position_in_bin,        # No extra gap!
        bin == 4 ~ 3*bin_width + 2*gap + position_in_bin
      )
    )



 plot <- ggplot(moving_average_data, aes(x = schematic_position, y = moving_avg, group = bin)) +
    geom_line(color = line_color, linewidth = 0.8) +


   # Vertical dashed lines at boundaries
   geom_vline(data = boundary_lines, aes(xintercept = xintercept),
              linetype = "dashed", color = "gray70", linewidth = 0.5) +

   # // break symbols
   annotate("text", x = break_x, y = 0, label = "//", size = 8, fontface = "bold") +

   #Zero line
   geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +

   # Exon boxes
   geom_rect(data = exon_regions,
             aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
             color = "black", linewidth = 0.5, inherit.aes = FALSE) +
   scale_fill_identity() +

   # Intron lines
   geom_segment(data = intron_segments,
                aes(x = x, xend = xend, y = y, yend = yend, linetype = linetype),
                color = "black", linewidth = 1.5, inherit.aes = FALSE) +
   scale_linetype_identity() +

   scale_y_continuous(labels = scales::scientific) +

    labs(x = NULL, y = "Overlap Frequency", title = title) +
    theme_minimal() +
    theme(
           # axis.text.x         = element_text(size = 12),
           # axis.ticks.x        = element_line(),
           # axis.title.x.top    = element_text(size = 12, face = "bold"),
           # axis.text.x.top     = element_text(size = 12),
           # axis.ticks.x.top    = element_line(),
           # axis.title.y=element_blank(),
           # axis.text.y=element_blank(),
           # axis.ticks.y=element_blank(),

           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           panel.border = element_blank(),
           axis.line = element_blank(),

           panel.background = element_rect(fill = "transparent",colour = NA),
           plot.background = element_rect(fill = "transparent",colour = NA),
           plot.title = element_text(hjust=0.5,size=20,color = "black",   face = "bold.italic")
           )
 return(plot)
}


