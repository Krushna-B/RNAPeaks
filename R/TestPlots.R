

testGenes <- c("APOE","BRCA1","TP53","HBA1")


print_test_Plots <- function(){
  for (name in testGenes) {
        PlotGene(bed = bed_df,geneID=name,gtf=gtf)
  }
}
