 createSplicingMap <- function(bed_File = NULL,
              Protein = NULL,
             SEMATS = NULL,
             moving_average = NULL,
             p_valueRetainedAndExclusion = 0.05,
             p_valueControls = 0.08,
            retained_IncLevelDifference = 0.1,
            exclusion_IncLevelDifference = -0.1,
            WidthIntoExon = 50,
            WidthIntoIntron = 250
             ){

#First row is names
names(SEMATS) <- as.character(unlist(SEMATS[1, ]))
SEMATS <- type.convert(SEMATS[-1, , drop = FALSE], as.is = TRUE)

#Remove extra ID columns
SEMATS <- SEMATS[, !duplicated(names(SEMATS)), drop = FALSE]

print(any(duplicated(names(SEMATS))))
print(names(SEMATS)[duplicated(names(SEMATS))] )
#Separate out Retained, excluded, controls exon

#TODO - Get columnames based on row1.
Retained <- SEMATS %>%
  filter(PValue < p_valueRetainedAndExclusion,
         FDR < p_valueRetainedAndExclusion,
         IncLevelDifference > retained_IncLevelDifference)
Exclusion <- SEMATS %>%
  filter(PValue < p_valueRetainedAndExclusion,
         FDR < p_valueRetainedAndExclusion,
         IncLevelDifference < exclusion_IncLevelDifference)

Controls <- SEMATS %>% filter(PValue > p_valueControls)



ENV_Retained <<- Retained
ENV_Exclusion <<- Exclusion
ENV_Controls <<- Controls


#Picking Controls
#Do later Separate Function



#Build Matrices for all 3 situations (Beware of Stands)
Retained_Bins <- make_bins_matrix(Retained,WidthIntoExon = WidthIntoExon,WidthIntoIntron = WidthIntoIntron)

Retained_Bins_Vis <<- as.data.frame(Retained_Bins)



#Buckets for specific protein
if(!is.null(Protein)){
  bed_df <- subset(bed_df,V4 == Protein)
}
buckets <<- makeGRangesFromDataFrame(bed_df,
                                    seqnames.field="V1",
                                    start.field="V2",
                                    end.field="V3",
                                    strand.field="V6",
                                    keep.extra.columns=TRUE)



# buckets_vs <<- as.data.frame(buckets)

#Find overlaps
# hits <<- findOverlaps(Retained_Bins, buckets, ignore.strand=FALSE)
# overlap_table <<- as.data.frame(hits)
# overlap_table$event_id <<- Retained_Bins$event_id[overlap_table$queryHits]
#
# overlap_counts <- table(overlap_table$queryHits)

overlap_matrix <<- make_overlap_matrix(Retained_Bins,buckets)
# overlap_index <- as.integer(names(overlap_counts))
# retained_gr$ov <- 0L
# retained_gr$ov[overlap_index] <- as.integer(overlap_counts)

#Make strand order correct
unique_events <- unique(Retained_Bins$event_id)

for (i in seq_along(unique_events)) {
  ev_id <- unique_events[i]
  # ensure same type comparison
  idx <- as.character(Retained_Bins$event_id) == as.character(ev_id)

  # skip if event not found (shouldn't happen, but safe)
  if (!any(idx)) next

  strand_i <- unique(as.character(Retained_Bins$strand[idx]))

  if (length(strand_i) > 0 && strand_i == "-") {
    overlap_matrix[i, ] <- rev(overlap_matrix[i, ])
  }
}





#Output frequency for every position value
binding_events <- colSums(overlap_matrix, na.rm = TRUE)
binding_events_proportion   <- binding_events / nrow(overlap_matrix)

plot_data <<- binding_events_proportion

#Plot Matrix
df <- data.frame(
  position   = seq_along(binding_events),
  proportion = binding_events_proportion
)

ggplot2::ggplot(df, ggplot2::aes(position, proportion)) +
  ggplot2::geom_line(linewidth = 0.7) +
  ggplot2::scale_x_continuous(expand = c(0, 0)) +
  ggplot2::labs(x = "Position", y = "Proportion (0–1)") +
  ggplot2::theme_minimal(base_size = 12)

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

   # Sign (+1 for + strand, -1 for – strand)
   sgn <- ifelse(strand == "+", 1, -1)

   #Shifting by signed distance
   shift <- function(base, dist_exon, dist_intron = NULL) {
     if (is.null(dist_intron)) dist_intron <- dist_exon
     out <- numeric(length(base))
     out[strand == "+"] <- base[strand == "+"] + dist_exon
     out[strand == "-"] <- base[strand == "-"] - dist_intron
     out
   }

   #  4 sites of interest
   site1_start <- shift(upE,  -WidthIntoExon)
   site1_end   <- shift(upE,   WidthIntoIntron)

   site2_start <- shift(exonStart, -WidthIntoIntron)
   site2_end   <- shift(exonStart,  WidthIntoExon)

   site3_start <- shift(exonEnd,  -WidthIntoExon)
   site3_end   <- shift(exonEnd,   WidthIntoIntron)

   site4_start <- shift(downS, -WidthIntoIntron)
   site4_end   <- shift(downS,  WidthIntoExon)

   # Normalize to genomic coordinates (start ≤ end)
   normalize <- function(a,b) list(start=pmin(a,b), end=pmax(a,b))
   s1 <- normalize(site1_start, site1_end)
   s2 <- normalize(site2_start, site2_end)
   s3 <- normalize(site3_start, site3_end)
   s4 <- normalize(site4_start, site4_end)

   starts <- cbind(s1$start, s2$start, s3$start, s4$start)
   ends   <- cbind(s1$end,   s2$end,   s3$end,   s4$end)


    #Build GRanges
    GRanges(
     seqnames = rep(SEMATS$chr, each = 4),
     ranges   = IRanges(start = as.vector(t(starts)),
                        end   = as.vector(t(ends))),
     strand   = rep(SEMATS$strand, each = 4),
     event_id = rep(seq_len(nrow(SEMATS)), each = 4)
   )
 }


 make_overlap_matrix <- function(event_gr, buckets_gr, window_length = 1200) {
   # For each event_id, extract its ranges and flatten into 1 vector of genomic positions
   unique_events <- unique(event_gr$event_id)
   n_events <- length(unique_events)

   overlap_matrix <- matrix(0L, nrow = n_events, ncol = window_length)

   for (i in seq_along(unique_events)) {
     e <- event_gr[event_gr$event_id == unique_events[i]]

     # Merge the 4 sites for this event into one ordered track
     merged <- reduce(e)  # merge overlapping intervals

     # Build vector of all genomic positions in these merged windows
     coords <- unlist(apply(as.data.frame(ranges(merged)), 1, function(r) seq(r[1], r[2])))
     if (length(coords) > window_length) coords <- coords[seq_len(window_length)]

     # Make GRanges of each position (1 bp)
     pos_gr <- GRanges(seqnames = seqnames(e)[1], ranges = IRanges(coords, coords), strand = strand(e)[1])

     # Find overlaps with buckets
     hits <- countOverlaps(pos_gr, buckets_gr, ignore.strand = FALSE) > 0

     overlap_matrix[i, seq_along(hits)] <- as.integer(hits)
   }

   overlap_matrix
 }
