Make_Strand_Labels <- function(Gene_s, offset = 100) {
  #Choose a vertical position for the strand tags
  ys <- min(Gene_s$y_start)

  if (Gene_s$strand[1] == "+") {
    #+ strand: transcription runs left -> right
    # Place 5′ a bit to the left of the gene, 3′ a bit to the right

    left  <- data.frame(Label="5'", X=min(Gene_s$start) - offset, Y=ys)
    right <- data.frame(Label="3'", X=max(Gene_s$end)   + offset, Y=ys)
  } else {
    # − strand: transcription runs right -> left (labels flip sides)
    # Place 3′ to the left and 5′ to the right of the gene.

    left  <- data.frame(Label="3'", X=min(Gene_s$start) - offset, Y=ys)
    right <- data.frame(Label="5'", X=max(Gene_s$end)   + offset, Y=ys)
  }

  list(left = left, right = right)
}
