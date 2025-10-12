Finding_Left_Margin <- function(bed,
                                label_size_mm = 5,              # font size (mm) used to measure labels
                                pad_pt       = 8,              # extra left padding (points)
                                pts_per_mm   = 72.27 / 25.4){

  #Converts Label Size from mm -> pt
  label_size_pt <- label_size_mm * pts_per_mm

  # pick the longest label in bed file by character count
  longest_label <- bed$group_name[ which.max(nchar(bed$group_name)) ]

  # measure width of the longest label (in points)
  w_pt <- convertWidth(
    grobWidth(textGrob(longest_label, gp = gpar(fontsize = label_size_pt))),
    "pt", valueOnly = TRUE
  )

  # add a small pad so text doesn't hug the plot area
  return (w_pt + pad_pt)
}
