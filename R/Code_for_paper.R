



#1C
m<-plot_gene(bed=rbind(HepG2_bed,K562_bed),gene="TAF1D",
             bam_files = c("IN"="~/Desktop/RNAPeaks/Bams/HEPG2_IN.bam",
                           "BUD13"="~/Desktop/RNAPeaks/Bams/BUD13_HEPG2_IP.bam",
                           "AQR"="~/Desktop/RNAPeaks/Bams/AQR_HEPG2_IP.bam",
                           "PRPF8"="~/Desktop/RNAPeaks/Bams/PRPF8_HEPG2_IP.bam"),
             peaks_opts = peaks_options(include = c("AQR","BUD13","CDC40","EFTUD2","PPIG","PRPF4","GPKOW","PRPF8","RBM22")),
             style = peaks_plot_style(bam_label_size = 5,title_size = 16,subtitle_size = 12,bam_axis_text_size = 0,
                                      protein_label_size = 5))
ggsave("~/Desktop/RNAPeaks/Figures/Figure1C.pdf",height=8,width=10)


#1D
m<-plot_gene(bed=rbind(HepG2_bed,K562_bed),gene="ACTB",
             peaks_opts = peaks_options(max_proteins = 8),
             style = peaks_plot_style(bam_label_size = 5,protein_label_size = 5,title_size = 16,subtitle_size = 12,bam_axis_text_size = 0,
                                      highlight = c(5527600,5529660)))
m
ggsave("~/Desktop/RNAPeaks/Figures/Figure1D.pdf",height=5,width=10)


#1E
m<-plot_region(bed=rbind(HepG2_bed,K562_bed),chr = "7",start="5527600",end="5529660",strand="-",
             peaks_opts = peaks_options(max_proteins = 8,
                                        include = c("YBX3","G3BP1","RPS3","PRPF8","RBFOX2","AQR","EFTUD2")),
             style = peaks_plot_style(bam_label_size = 5,title_size = 16,subtitle_size = 12,bam_axis_text_size = 0))
m
ggsave("~/Desktop/RNAPeaks/Figures/Figure1E.pdf",height=4,width=10)


#2A
df<-read.table("~/Desktop/Dominguez_Lab/Resource/Splicing_ENCODE/All_SE_K562.txt",header=T)
df<-df[which(df$Proteins=="K562_AATF"),]
df$Group<-"ns"
df$Group[which(df$PValue<0.05 & df$FDR<0.05 & df$IncLevelDifference>0.1)]<-"up"
df$Group[which(df$PValue<0.05 & df$FDR<0.05 & df$IncLevelDifference<(-0.1))]<-"down"

cols <- c("up" = "red", "down" = "blue", "ns" = "grey")
sizes <- c("up" = 0.75, "down" = 0.75, "ns" = 0.25)
alphas <- c("up" = 0.7, "down" = 0.7, "ns" = 0.5)

volcano<-df %>%
  ggplot(aes(x = IncLevelDifference,
             y = -log10(FDR),
             fill = Group,
             color=Group,
             size = Group,
             alpha =Group)) +
  geom_point(shape = 21) +
  geom_hline(yintercept = -log10(0.05),
             linetype = "dashed") +
  geom_vline(xintercept = c(-0.1,0.1),
             linetype = "dashed")+
  scale_fill_manual(values = cols) + # Modify point colour
  scale_color_manual(values = cols) + # Modify point colour
  scale_size_manual(values = sizes) + # Modify point size
  scale_alpha_manual(values = alphas) + # Modify point transparency
  scale_x_continuous(breaks = c(-0.4,-0.2,0,0.2,0.4),
                     limits = c(-0.5, 0.5))+
  lims(y=c(0,4))+
  labs(x=bquote(Delta~PSI))+
  theme_classic()+
  theme(legend.position = "none",
        axis.text=element_text(family="Arial",size=12, color="black"),
        axis.title = element_text(family="Arial",size=12, color="black"),
        panel.background = element_blank(),
        plot.background = element_blank(),
        panel.border = element_blank(),
        legend.background = element_blank())
volcano
ggsave("~/Desktop/RNAPeaks/Figures/Figure2A.png",height=3,width=3)





