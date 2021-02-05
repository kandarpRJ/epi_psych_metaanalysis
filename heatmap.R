library(ComplexHeatmap)
library(RColorBrewer)

data<-read.table("MDD_mastersheet.csv", header = TRUE, row.names = NULL, sep = ",")

data<-data[with (data, order(data$Brain.region, data$Subregion, data$Sex)),]
mat<-data[,5:ncol(data)]

brcol<-brewer.pal(n=12, name = "Paired")
brcol[13]<-brewer.pal(n=4, name = "Set3")
names(brcol)<-unique (sort (data$Brain.region))

sbrcol<-brewer.pal(n=9, name = "Set3")
sbrcol[10]<-"black"
names(sbrcol)<-unique(sort(data$Subregion))
sbrcol["L5"]<-"#E7298A"

ha = HeatmapAnnotation(df = data[,2:4], col = list(Brain.region = brcol, Subregion = sbrcol,
                       Sex=c("F"="green", "M"="yellow", "Mix"="orange", "Unknown/NA"="black")))

h<-Heatmap(as.matrix(t(mat)), cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = ha, 
           column_labels = data$Study, heatmap_legend_param = list (title="Log2FC", direction="horizontal", title_position="topcenter"),
           row_split = data.frame(c(rep("Writers", 7), rep("Readers", 9), rep("Erasers", 2))),
           column_split = data.frame(data$Brain.region), column_title = NULL)

png("figure4.png", res = 600, width = 5000, height = 5000)
draw(h, heatmap_legend_side="bottom")
dev.off()



########################################################################################################################

data<-read.table("PTSD_mastersheet.csv", header = TRUE, row.names = NULL, sep = "\t")

data<-data[with (data, order(data$Study)),]
mat<-data[,4:ncol(data)]

ha = HeatmapAnnotation(df = data[,2:3], 
                       col = list(Condition=c("PTSD_vs_NoPTSD"="gold", "PTSD_Notimproved_vs_Improved"="cyan", 
                                              "PTSD_vs_NoPTSD_ER"="darkcyan", "PTSD_vs_NoPTSD_M4"="purple")))

h<-Heatmap(as.matrix(t(mat)), cluster_rows = FALSE, cluster_columns = FALSE, top_annotation = ha, 
           column_labels = data$Study, 
           heatmap_legend_param = list (title="Log2FC", direction="horizontal", title_position="topcenter"),
           row_split = data.frame(c(rep("Writers", 11), rep("Readers", 9), rep("Erasers", 2))),
           column_split = data.frame(data$Study), column_title = NULL)

png("figure6.png", res = 600, width = 5000, height = 5000)
draw(h, heatmap_legend_side="bottom")
dev.off()


#########################################################################################################################
