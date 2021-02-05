##########################################################################################################################################
#### Acknowledgement to Malosree Maitra, PhD Student, Gustavo Turecki Lab for providing the scripts used for analysis of this dataset ####
##########################################################################################################################################

library(Matrix)
library(Seurat)
library(pheatmap)
library(ggplot2)

#Read in matrix, cell names, and gene names
mdd_matrix <- readMM("GSE144136_GeneBarcodeMatrix_Annotated.mtx") ## Available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144136
mdd_matrix <- as(mdd_matrix, "dgCMatrix")
cell_names <- read.csv("GSE144136_CellNames.csv") ## Available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144136
gene_names <- read.csv("GSE144136_GeneNames.csv") ## Available at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144136
colnames(mdd_matrix) <- cell_names$x
rownames(mdd_matrix) <- gene_names$x

#Create Seurat object and normalize data
mdd_seurat <- CreateSeuratObject(counts = mdd_matrix)
mdd_seurat <- NormalizeData(mdd_seurat, normalization.method = "LogNormalize", scale.factor = 10000)

#Create metadata and add to Seurat object
meta.data <- matrix(unlist(strsplit(colnames(mdd_seurat), split = ".", fixed = TRUE)), byrow = TRUE, ncol = 2)
meta.data <- cbind(meta.data, matrix(unlist(gsub ("_\\w*", "", meta.data[,1])))) ## extra
meta.data <- cbind(meta.data, matrix(unlist(strsplit(meta.data[,2], split ="_")), byrow = TRUE, ncol = 4))
#meta.data <- meta.data[,c(1,3,4,5)]
meta.data <- meta.data[,c(3,1,4,5,6)]
colnames(meta.data) <- c("cell_type", "cell_subtype", "subject", "condition", "batch")
rownames(meta.data) <- colnames(mdd_seurat)
mdd_seurat <- AddMetaData(mdd_seurat, as.data.frame(meta.data))

#Assign cell-type clusters
Idents(mdd_seurat) <- mdd_seurat$cell_subtype

mdd.list<-SplitObject(mdd_seurat, split.by = "condition")
mdd.list <- lapply(X = mdd.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
anchors <- FindIntegrationAnchors(object.list = mdd.list, dims = 1:20)
combined <- IntegrateData(anchorset = anchors, dims = 1:20)
DefaultAssay(combined) <- "RNA"
combined$celltype.condition <- paste(combined$condition, Idents(combined), sep = "_")
Idents(combined) <- "celltype.condition"

response <- FindMarkers(combined, ident.1 = suicide_celltype, ident.2 = control_celltype, verbose = FALSE) ## Replace suicide_celltype and control_celltype with appropriate Idents (names)

features = c("METTL3", "METTL14", "WTAP", "METTL16", "RBM15", "CBLL1", "YTHDC1", 
	   "YTHDC2", "YTHDF1", "YTHDF2", "YTHDF3", "FMR1", "IGF2BP1", "IGF2BP2", "IGF2BP3", "FTO", "ALKBH5")


AvgExp_GenesOfInterest<-AverageExpression(combined, assays = "RNA", features = features)

new_dat<-AvgExp_GenesOfInterest$RNA[,c("Control_Astros_1", "Control_Astros_2", "Suicide_Astros_2", "Control_Astros_3", "Suicide_Astros_3", "Control_Endo", "Suicide_Endo",
                                       "Control_Ex_1_L5_6", "Suicide_Ex_1_L5_6", "Control_Ex_10_L2_4", "Suicide_Ex_10_L2_4", "Control_Ex_2_L5", "Suicide_Ex_2_L5", "Control_Ex_3_L4_5",
                                       "Suicide_Ex_3_L4_5", "Control_Ex_4_L_6", "Suicide_Ex_4_L_6", "Control_Ex_5_L5", "Suicide_Ex_5_L5", "Control_Ex_6_L4_6", "Suicide_Ex_6_L4_6",
                                       "Control_Ex_7_L4_6", "Suicide_Ex_7_L4_6", "Control_Ex_8_L5_6", "Suicide_Ex_8_L5_6", "Control_Ex_9_L5_6", "Suicide_Ex_9_L5_6", "Control_Inhib_1",
                                       "Suicide_Inhib_1", "Control_Inhib_2_VIP", "Suicide_Inhib_2_VIP", "Control_Inhib_3_SST", "Suicide_Inhib_3_SST", "Control_Inhib_4_SST",
                                       "Suicide_Inhib_4_SST", "Control_Inhib_5", "Suicide_Inhib_5", "Control_Inhib_6_SST", "Suicide_Inhib_6_SST", "Control_Inhib_7_PVALB", "Suicide_Inhib_7_PVALB",
                                       "Control_Inhib_8_PVALB", "Suicide_Inhib_8_PVALB", "Control_Micro/Macro", "Suicide_Micro/Macro", "Control_Mix_1", "Suicide_Mix_1", "Control_Mix_2",
                                       "Suicide_Mix_2", "Control_Mix_3", "Suicide_Mix_3", "Control_Mix_4", "Suicide_Mix_4", "Control_Mix_5", "Suicide_Mix_5", "Control_Oligos_1", "Suicide_Oligos_1",
                                       "Control_Oligos_2", "Suicide_Oligos_2", "Control_Oligos_3", "Suicide_Oligos_3", "Control_OPCs_1", "Suicide_OPCs_1", "Control_OPCs_2", "Suicide_OPCs_2")]

cat_names<-read.table("colnames.txt", header = FALSE, row.names = NULL, sep = "\t")
colnames(cat_names)<-c("Condition", "Cell type")
ha<-HeatmapAnnotation(df=cat_names)

h<-Heatmap(as.matrix(new_dat), cluster_columns = FALSE, cluster_rows = FALSE, top_annotation = ha,
           row_split = data.frame(c(rep("Writers", 6), rep("Readers", 9), rep("Erasers", 2))),
           column_split = cat_names$`Cell type`, column_names_gp = gpar(fontsize = 9),
           heatmap_legend_param = list (title="Avg. gene expression", direction="horizontal", title_position="topcenter"))

png("figure3.png", res = 600, height = 4000, width = 7500)
draw(h, heatmap_legend_side="bottom")
dev.off()

