library(plyr)
library(dplyr)
library(ggplot2)
library(cowplot)

## Inputs

cblm<-"human_cblm_m6a_containing_genes_with_m6a_counts.txt"
crbm<-"human_crbm_m6a_containing_genes_with_m6a_counts.txt"

## Functions 

makeplot<-function (m6agl, fn, grp) {
  m6a_gl<-read.table(m6agl, header = FALSE, row.names = NULL)
  degl<-read.table(fn, header = TRUE, row.names = NULL, sep = "\t")
  degl<-degl[degl$P.Value<0.05,]
  m6a_1<-degl[degl$Gene.symbol %in% m6a_gl[m6a_gl$V2==1,1],]
  m6a_2_5<-degl[degl$Gene.symbol %in% m6a_gl[m6a_gl$V2>1 & m6a_gl$V2<=5,1],]
  m6a_6<-degl[degl$Gene.symbol %in% m6a_gl[m6a_gl$V2>5,1],]
  
  `%notin%` <- Negate(`%in%`)
  
  m6a_0<-degl[degl$Gene.symbol %notin% m6a_gl$V1,]
  
  newdf<-data.frame(m6a_counts = "m6a_0", lfc = m6a_0$logFC, grp = grp)
  newdf<-rbind(newdf, data.frame(m6a_counts = "m6a_1", lfc = m6a_1$logFC, grp = grp))
  newdf<-rbind(newdf, data.frame(m6a_counts = "m6a_2_5", lfc = m6a_2_5$logFC, grp = grp))
  newdf<-rbind(newdf, data.frame(m6a_counts = "m6a_6", lfc = m6a_6$logFC, grp = grp))
  newdf <- ddply(newdf, .(m6a_counts), transform, ecd = ecdf(lfc)(lfc))
  
}

pl<-function (p) {
  ggplot(p, aes(x = lfc, y = ecd, group = m6a_counts, colour = m6a_counts)) + 
    geom_line() +
    geom_vline(xintercept = 0, color="grey", alpha=0.5, linetype="dashed") +
    geom_hline(yintercept = 0.5, color="grey", alpha=0.5, linetype="dashed") +
    facet_wrap(~grp, ncol = 3) +
    xlab("Log2FoldChange") +
    ylab("Cumulative Distribution") +
    xlim(-1, 1) +
    ylim(0,1) +
    theme_classic()
}

pl1<-function (p) {
  ggplot(p, aes(x = lfc, y = ecd, group = m6a_counts, colour = m6a_counts)) + 
    geom_line() +
    geom_vline(xintercept = 0, color="grey", alpha=0.5, linetype="dashed") +
    geom_hline(yintercept = 0.5, color="grey", alpha=0.5, linetype="dashed") +
    facet_wrap(~grp, ncol = 3) +
    xlab("Log2FoldChange") +
    ylab("Cumulative Distribution") +
    xlim(-2, 2) +
    ylim(0,1) +
    theme_classic()
}

#######################################################################################################################################

## Make plots for all DEGs 
## Function takes arguements such as file with gene name and corresponding m6a sites, DEG list and group/sample/condition name as input 

p1<-makeplot(crbm, "GSE54570.txt", "DLPFC_BA9")
p2<-makeplot(crbm, "GSE54568.txt", "DLPFC_BA9_F")
p3<-makeplot(crbm, "GSE54567.txt", "DLPFC_BA9_M")
p4<-makeplot(crbm, "GSE87610_DLPFC_MDD_L3.txt", "DLPFC_L3")
p5<-makeplot(crbm, "GSE87610_DLPFC_MDD_L5.txt", "DLPFC_L5")
p6<-makeplot(crbm, "GSE54571.txt", "ACC_BA25_F")
p7<-makeplot(crbm, "GSE54572.txt", "ACC_BA25_M")
p8<-makeplot(crbm, "GSE92538_DLPFC_male", "DLPFC_M")
p9<-makeplot(crbm, "GSE92538_DLPFC_female", "DLPFC_F")
p10<-makeplot(crbm, "GSE17440_HIV_MDD.txt", "FC_HIV_M")
p11<-makeplot(crbm, "GSE54575.txt", "OVPFC_BA47")
p12<-makeplot(crbm, "GSE35977_CORT.txt", "PariCORT_Mix")
p13<-makeplot(crbm, "GSE53987_PFC_MDD.txt", "PFC_BA46")
p14<-makeplot(crbm, "GSE12654.txt", "PFC_BA10")
p15<-makeplot(crbm, "biorxiv_DLPFC_MDD.csv", "DLPFC_Jaffe_et.al.")
p16<-makeplot(crbm, "biorxiv_ACC_MDD.csv", "ACC_Jaffe_et.al.")
p17<-makeplot(crbm, "biorxiv_DLPFC_PTSD.csv", "DLPFC_Jaffe_et.al.")
p18<-makeplot(crbm, "biorxiv_ACC_PTSD.csv", "ACC_Jaffe_at.al.")
p19<-makeplot(cblm, "GSE35974_CBLM.txt", "CBLM_Mix")

## MDD CDF plot

pg<-plot_grid(pl(p2) + theme(legend.position="none"), 
              pl(p3) + theme(legend.position="none"),
              pl(p6) + theme(legend.position="none"),
              pl(p7) + theme(legend.position="none"),
              pl(p9) + theme(legend.position="none"),
              pl(p8) + theme(legend.position="none"),
              pl(p15) + theme(legend.position="none"),
              pl(p16) + theme(legend.position="none"),
              pl(p4) + theme(legend.position="none"),
              pl(p5) + theme(legend.position="none"),
              pl(p11) + theme(legend.position="none"),
              pl(p13) + theme(legend.position="none"),
              pl(p14) + theme(legend.position="none"),
              pl1(p10) + theme(legend.position="none"),
              pl(p19) + theme(legend.position="none"), ncol = 2, labels = "AUTO", label_size = 12 )

legend <- get_legend(
  pl(p2) + theme(legend.box.margin = margin(0, 0, 0, 12))
)

png("mdd_cdf_grid.png", res = 300, height = 6000, width = 3000)
plot_grid(pg, legend, rel_widths = c(3, .4))
dev.off()


## PTSD CDF plot

pgptsd<-plot_grid(pl(p17) + theme(legend.position="none"), 
                  pl(p18) + theme(legend.position="none"), ncol = 2, labels = "AUTO", label_size = 12 )

legendptsd <- get_legend(
  pl(p17) + theme(legend.box.margin = margin(0, 0, 0, 12))
)

png("ptsd_cdf_grid.png", res = 300, height = 3000, width = 6000)
plot_grid(pgptsd, legend, rel_widths = c(3, .4))
dev.off()





