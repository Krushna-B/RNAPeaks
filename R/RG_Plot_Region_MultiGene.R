# Internal function for multi-gene region plotting
Get_Multi_Plot_by_Region <- function(Chr,
                                     Start,
                                     End,
                                     Strand,
                                     bed,
                                     gtf,
                                     gene,
                                     region_with_multiple_genes,
                                     rank_,
                                     density = FALSE,
                                     Species = "Human",
                                     utr_col = "dark gray",
                                     peak_col = "blue",
                                     peaks_width = 0.3,
                                     exon_width = 0.5,
                                     utr_width = 0.3,
                                     exon_col = "navy",
                                     geneID = NULL,
                                     TxID = NULL,
                                     Vertical_Offset_Step = 0.7,
                                     ...) {

  # Fallback gene frame used by Prepare_Bed if none provided
  if (is.null(gene)) {
    gene <- data.frame(seqnames = Chr, start = Start, end = End, strand = Strand)
  }

  # Return prepared bed with labeled coordinates positions for proteins
  bed <- Prepare_Bed(bed, rank_ = rank_, gene = gene,
                     peaks_width = peaks_width, peak_col = peak_col)

  # If density is TRUE
  if (isTRUE(density)) {
    Dens <- Density(Chr = gene$seqnames[1], Start = min(gene$start), End = max(gene$end), df = bed)
  }

  # Error check to make sure gtf rows are not empty
  if (!nrow(region_with_multiple_genes)) {
    stop("No GTF rows overlap the requested region.")
  }

  # Gene IDs present in window one transcript for each geneID
  gene_ids <- unique(region_with_multiple_genes$gene_id)
  gene_ids <- gene_ids[!is.na(gene_ids)]
  # Make sure gene ID's exist
  stopifnot(length(gene_ids) > 0)

  # Build per-gene gene structures (CDS/exon/UTR/intron)
  gene_structs <- lapply(gene_ids, function(g_id) {
    gene_rows <- region_with_multiple_genes[region_with_multiple_genes$gene_id == g_id, , drop = FALSE]

    g_name <- gene_rows$gene_name[1]
    gs <- Build_Gene_Structure(
      gene = gene_rows,
      levels = max(bed$rank),
      peaks_width = peaks_width,
      exon_width = exon_width,
      utr_width = utr_width,
      exon_col = exon_col,
      utr_col = utr_col
    )
    if (is.null(gs)) {
      return(NULL)
    }

    add_meta <- function(df) {
      if (!is.null(df) && nrow(df) > 0) {
        df$gene_id <- g_id
        df$gene_name <- g_name
      }
      df
    }

    gs$Gene_s <- add_meta(gs$Gene_s)
    gs$Exons <- add_meta(gs$Exons)
    gs$UTRs <- add_meta(gs$UTRs)
    gs$Intron_s <- add_meta(gs$Intron_s)
    return(gs)
  })

  # Delete empty gene structs
  gene_structs <- Filter(Negate(is.null), gene_structs)

  # Combine all of them into one set of data frames
  Gene_s_all <- do.call(rbind, lapply(gene_structs, `[[`, "Gene_s"))
  Exons_all <- do.call(rbind, lapply(gene_structs, `[[`, "Exons"))
  UTRs_all <- do.call(rbind, lapply(gene_structs, `[[`, "UTRs"))
  Intron_all <- do.call(rbind, lapply(gene_structs, `[[`, "Intron_s"))

  # Apply vertical offsets
  gene_order <- unique(Gene_s_all[order(Gene_s_all$start), "gene_id"])
  Gene_s_all <- apply_vertical_offset(Gene_s_all, id_col = "gene_id", Vertical_Offset_Step, gene_order)
  Exons_all <- apply_vertical_offset(Exons_all, id_col = "gene_id", Vertical_Offset_Step, gene_order)
  UTRs_all <- apply_vertical_offset(UTRs_all, id_col = "gene_id", Vertical_Offset_Step, gene_order)
  Intron_all <- apply_vertical_offset(Intron_all, id_col = "gene_id", Vertical_Offset_Step, gene_order)

  # All Gene combined into one for reference in the plotting
  combined_gene <- data.frame(
    seqnames = Gene_s_all$seqnames[1],
    start = min(Gene_s_all$start),
    end = max(Gene_s_all$end),
    strand = "*",
    gene_name = paste(unique(Gene_s_all$gene_name), collapse = ", ")
  )

  # Arrows on all of the genes
  Intron_all$mid_y <- (Intron_all$y_start + Intron_all$y_end) / 2

  arrows_list <- lapply(seq_along(gene_structs), function(i) {
    List_of_unique_genes <- split(Gene_s_all, Gene_s_all$gene_id)
    List_of_unique_genes_introns <- split(Intron_all, Intron_all$gene_id)

    gs <- List_of_unique_genes[[i]]
    intr <- List_of_unique_genes_introns[[i]]

    introns <- Compute_Intron_Positions(gs, intr)

    if (is.null(introns) ||
        is.null(introns$Introns_Positions) ||
        nrow(introns$Introns_Positions) == 0L) {
      return(NULL)
    }

    RG_make_intron_arrows(introns$Introns_Positions, gs$strand[1], Start = Start, End = End)
  })

  Arrows_all <- dplyr::bind_rows(arrows_list)

  # Makes labels for each gene
  label_df <- if (nrow(Gene_s_all)) {
    # y per lane
    agg <- stats::aggregate(y_start ~ gene_id, Gene_s_all, mean)

    # row-wise label candidate
    lab_col <- if ("gene_name" %in% names(Gene_s_all)) {
      ifelse(!is.na(Gene_s_all$gene_name) & nzchar(Gene_s_all$gene_name),
             Gene_s_all$gene_name, Gene_s_all$gene_id)
    } else {
      Gene_s_all$gene_id
    }
    # one label per gene_id
    lab_map <- unique(data.frame(gene_id = Gene_s_all$gene_id,
                                 label = lab_col,
                                 stringsAsFactors = FALSE))
    out <- merge(agg, lab_map, by = "gene_id", all.x = TRUE)
    out$x <- Start - (End - Start) * 0.01
    out
  } else {
    data.frame(gene_id = character(0), y_start = numeric(0), label = character(0), x = numeric(0))
  }

  # Clean Gene Names to consistent values if they are empty
  clean_gene_names <- function(df) {
    df$gene_name <- as.character(df$gene_name)
    bad <- is.na(df$gene_name) |
      df$gene_name == "" |
      df$gene_name == " " |
      df$gene_name == "NA"
    df$gene_name[bad] <- NA_character_
    df
  }

  Gene_s_all <- clean_gene_names(Gene_s_all)

  # Plot Gene
  p <- Draw_Gene_Plot(
    Gene_s = Gene_s_all,
    Exons = Exons_all,
    UTRs = UTRs_all,
    Intron_s = Intron_all,
    Arrow_df = Arrows_all,
    bed = bed,
    gene = combined_gene,
    peak_col = bed$peak_col,
    peaks_width = bed$peaks_width,
    is_region_plot = TRUE,
    ...
  )
  return(p)
}


# ---------Helper Functions---------

# Calculate Density of Peaks
Density <- function(Chr, Start, End, df) {
  peaks_gr <- GenomicRanges::makeGRangesFromDataFrame(
    df,
    keep.extra.columns = TRUE,
    start.field = "Start",
    end.field = "End",
    seqnames.field = "seqnames"
  )
  bins <- GenomicRanges::GRanges(
    seqnames = Chr,
    strand = "*",
    ranges = IRanges::IRanges(start = seq(Start, End, by = 10), width = 10)
  )
  Overlaps <- GenomicRanges::countOverlaps(bins, peaks_gr)
  S4Vectors::values(bins) <- S4Vectors::DataFrame(Density = Overlaps)
  bins <- as.data.frame(bins)
  return(bins)
}


# Prepare Background
Background_Table <- function(df, Start, End) {
  Background <- df[c("y_start", "y_end")]
  Background <- Background[!duplicated(Background), ]
  Background <- Background[order(Background$y_start), ]
  Background$x_start <- Start
  Background$x_end <- End
  Background$col <- rep(c("white", "white"), length.out = nrow(Background))
  return(Background)
}

# Builds vertical offset based off pre-specified order of the genes
apply_vertical_offset <- function(df, id_col = "gene_id", offset_step = 0.8, levels_order = NULL) {
  if (is.null(levels_order)) {
    levels_order <- unique(df[[id_col]])
  }
  df[[id_col]] <- factor(df[[id_col]], levels = levels_order)
  offset_vec <- (as.numeric(df[[id_col]]) - 1) * offset_step
  df$y_start <- df$y_start + offset_vec
  df$y_end <- df$y_end + offset_vec
  df
}
