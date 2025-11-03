
SM_PerformanceTesting <- function(){

#Output 1
microbenchmark(
  ouptut = createSplicingMap(bed_File = bed_df,SEMATS = MASTER_FILE, Protein = "RBFOX2_K562_IDR"),
  times = 1
)


}
