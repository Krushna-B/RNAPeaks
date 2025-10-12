Compute_Intron_Positions <- function(Gene_s, Intron_s) {
  # Compute genomic intron intervals from the ordered gene features
  # (get_Intron_Positions() finds gaps between consecutive rows in Gene_s)
  Introns_Positions <- get_Intron_Positions(Gene_s)

  Introns_Positions$y_start <- Intron_s$y_start
  Introns_Positions$y_end   <- Intron_s$y_end

  #Midline in y for drawing intron baselines/arrows
  Intron_s$mid_y <- (Intron_s$y_start + Intron_s$y_end) / 2

  list(Introns_Positions = Introns_Positions, Intron_s = Intron_s)
}
