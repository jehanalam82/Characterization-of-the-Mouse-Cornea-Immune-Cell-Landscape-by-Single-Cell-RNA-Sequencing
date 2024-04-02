library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(ggpubr)
library(cowplot)
wd <- "H:/Manuscripts/2024_K_SCRNA-Methods/Combined_Analysis" # MODIFY
setwd(wd)
######################Data was loaded from Cellqc pipeline########################
####################k_male######################
M <- Read10X("H:/Manuscripts/2024_K_SCRNA-Methods/PF_3V3_B6_KM/filtered_feature_bc_matrix")
M <- CreateSeuratObject(counts = M, project = "K_M", min.cells = 3, min.features = 200)
M[["percent.mt"]] <- PercentageFeatureSet(M, pattern = "^MT-|^mt-")
VlnPlot(M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "nonfiltered_vlnplot_M.pdf", plot = VlnPlot(M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
plot1 <- FeatureScatter(M, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(M, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot3 <- FeatureScatter(M, feature1 = "nFeature_RNA", feature2 = "percent.mt")
ggsave(filename = "plot1.pdf", plot = plot1)
ggsave(filename = "plot2.pdf", plot = plot2)
ggsave(filename = "plot3.pdf", plot = plot3)
M <- subset(M, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
vlnplot <- VlnPlot(M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "filtered_vlnplot_M.pdf", plot = vlnplot)
scatterplot <- FeatureScatter(M, feature1 = "percent.mt", feature2 = "nFeature_RNA")
ggsave(filename = "filtered_scatterplot_M.pdf", plot = scatterplot)
scatterplot2 <- FeatureScatter(M, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(filename = "filtered_scatterplot2_M.pdf", plot = scatterplot2)
scatterplot3 <- FeatureScatter(M, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave(filename = "filtered_scatterplot3_M.pdf", plot = scatterplot3)
#Data normalization##########################
M <- NormalizeData(object = M, )
M <- FindVariableFeatures(object = M, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(M), 20)
plot4 <- VariableFeaturePlot(M)
plot4 <- LabelPoints(plot = plot4, points = top20, repel = TRUE)
ggsave(filename = "plot4.pdf", plot = plot4)
all.genes <- rownames(M)
M <- ScaleData(object = M, features = all.genes, vars.to.regress = "percent.mt")
####################k_female######################
FM <- Read10X("H:/Manuscripts/2024_K_SCRNA-Methods/PF_61353_KB6F/filtered_feature_bc_matrix")
FM <- CreateSeuratObject(counts = FM, project = "K_FM", min.cells = 3, min.features = 200)
FM[["percent.mt"]] <- PercentageFeatureSet(FM, pattern = "^MT-|^mt-")
VlnPlot(FM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(filename = "nonfiltered_vlnplot_FM.pdf", plot = VlnPlot(FM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3))
plot5 <- FeatureScatter(FM, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot6 <- FeatureScatter(FM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot7 <- FeatureScatter(FM, feature1 = "nFeature_RNA", feature2 = "percent.mt")
ggsave(filename = "plot5.pdf", plot = plot5)
ggsave(filename = "plot6.pdf", plot = plot6)
ggsave(filename = "plot7.pdf", plot = plot7)
FM <- subset(FM, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 10)
scatterplot4 <- FeatureScatter(FM, feature1 = "percent.mt", feature2 = "nFeature_RNA")
ggsave(filename = "filtered_scatterplot1_FM.pdf", plot = scatterplot4)
scatterplot5 <- FeatureScatter(FM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
ggsave(filename = "filtered_scatterplot2_FM.pdf", plot = scatterplot5)
scatterplot6 <- FeatureScatter(FM, feature1 = "nCount_RNA", feature2 = "percent.mt")
ggsave(filename = "filtered_scatterplot3_FM.pdf", plot = scatterplot6)
#Data normalization##########################
FM <- NormalizeData(object = FM, )
FM <- FindVariableFeatures(object = FM, selection.method = "vst", nfeatures = 2000)
top20 <- head(VariableFeatures(FM), 20)
plot8 <- VariableFeaturePlot(FM)
plot8 <- LabelPoints(plot = plot8, points = top20_FM, repel = TRUE)
ggsave(filename = "plot8.pdf", plot = plot4)
all.genes <- rownames(FM)
FM <- ScaleData(object = FM, features = all.genes, vars.to.regress = "percent.mt")

############################iNTEGRATION################
MFM.list <- list(M, FM)
Common.anchors <- FindIntegrationAnchors(object.list = MFM.list, dims = 1:30)
Common.combined <- IntegrateData(anchorset = Common.anchors, dims = 1:30)
DefaultAssay(Common.combined) <- "integrated"
Common.combined <-FindVariableFeatures(Common.combined, verbose = FALSE)
Common.combined <- ScaleData(Common.combined, verbose = FALSE)
Common.combined <-RunPCA(Common.combined, npcs = 30, verbose = FALSE)
Common.combined <-RunTSNE(Common.combined, dims = 1:30)
Common.combined <-RunUMAP(Common.combined, reduction = "pca", dims = 1:30)
Common.combined <-FindNeighbors(Common.combined, reduction = "pca", dims = 1:30)
resolutions <- c(1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1)
dir.create("DimPlots", showWarnings = FALSE)
for (res in resolutions) {
  Common.combined=FindClusters(Common.combined, resolution=res)
}
prefix <- "MFM_Minus_Stromal_legend"
saveRDS(Common.combined, paste0(prefix, ".rds"))
setwd("DimPlots")
####################Analysis#########################
prefix <- "rMP_SUBTYPES"
x <- readRDS(paste0(prefix, ".rds"))
DefaultAssay(x)<- ("RNA")
# Set the resolutions
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
for (res in resolutions) {
  x <- SetIdent(x, value = paste0("integrated_snn_res.", res))
  plot_filename <- file.path("Dimplots", paste0("DefaultAssay(x)_", res, ".png"))
  ggsave(filename = plot_filename, plot = DimPlot(x, split.by = "orig.ident", label = TRUE))
}
#############We will use resolution 0.4#################
x <- SetIdent(x, value = "integrated_snn_res.0.4")
DimPlot(x, label = T)+NoLegend()
TSNEPlot(x, label = T)+NoLegend()
dimplot_filename <- "Dimplots_RES0.5.pdf"
ggsave(filename = dimplot_filename, plot = DimPlot(x,  label = TRUE) + NoLegend())
###################Cells count####################
x[["my.clusters"]] <- Idents(object = x)
table <- table(x@meta.data$my.clusters, x@meta.data$orig.ident)
print(table)
write.csv(table, file.path(getwd(), "cluster_Count_legend.csv"))
############################cluster biomarkers################
x<-JoinLayers(x)
x.markers <- FindAllMarkers(x, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
x.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)
write.csv(x.markers, paste0(prefix, ".csv"))

#######################Dot-plot##################
features <- c("Timd4", "Lyve1","Folr2","Gas6","Mrc1","Cd163","Igf1","Ccr2","Cd52","H2-Eb1","H2-Ab1","Cd14","Cst3","Apoe","Mertk")

#features <- c("Krt12","Tjp1","Aldh3a1","Epcam","Fxyd3",
"Ptprc","Itgam","Mrc1","Mertk","Gas6",
"Cd3d","Cd3g","Cd3e","Cd4","Cd8a",
"Apoe","C1qa","Ccl7","Ccl2",
"Igkc","Ebf1","Cd79a","Ly6d",
"Cd209a","H2-Eb1","H2-Ab1","H2-Aa",
"Cd28","Trac","Nkg7","Cd8b1",
"S100a9","S100a8","G0s2","Hdc",
"Mmp12","Mmp13","Plau","Il1a",
"Lef1","Gimap6","Dusp10","Gramd3",
"Gzma","AW112010","Irf8","Klra8",
"Ace","Plac8","Ifitm3","Ms4a6c",
"Gzmc","Ctsw","Il2rb","Cd7",
"Il17a","Trdc","Il23r","Cxcr6",
"Il1rl1","Gata3", "Csf2","Dach2",
"Top2a","Pclaf","Birc5","Ube2c",
"Mcpt4","Cma1","Tpsb2","Cpa3")
features <- c("Ccr2","Mki67","Mrc1","Cx3cr1","Mertk","Csf1r","")
DotPlot(x, features =features, cols = c("lightGrey", "darkred"), dot.scale = 6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1))
################################################################################################################Sub-setting######
saveRDS(subset(x, ident = c("0","1","2","4","5","6","7","8","11","12", "13", "14","16", "19","21","24","25")), paste0(prefix, ".rds"))

new.cluster.ids<-c("MP","B","cDC2-1","Cd8","Neut","MMP12/13h_MP","Cd4","NK","Mono","Mono","NKT","GDT","ILC2","MP","cDC2-2","NKT","Mast")
names(new.cluster.ids) <- levels(x)
x <- RenameIdents(x, new.cluster.ids)
saveRDS(x, paste0(prefix, ".rds"))
#######################Heatmap####################
saveRDS(subset(x, ident = c("MP","cDC2-1","MMP12/13h_MP","Mono","cDC2-2")), paste0(prefix, ".rds"))

################################
#Top 500 markers between strains
################################
x<- ScaleData(x)
#Get the Differential Genes
Idents(x) <- "orig.ident"
strain.markers <- FindMarkers(x, ident.1 = "K_M", 
                              only.pos = FALSE, min.pct = 0.25, logfc.threshold = 1.0)
#Selecting Top 500 Genes 
top100 <- strain.markers %>%
  filter(p_val_adj < 0.01) %>% #significant filter
  slice_max(abs(avg_log2FC), n = 500) %>% #Using absolute for high or low genes
  arrange(avg_log2FC) #Order by value
top100genes <- rownames(top100)
#Heatmap Across All Cells
# Create the heatmap with color
DoHeatmap(x, features = top100genes, draw.lines = FALSE, slot = "data")
#Average-based Heatmap
Avg <- AverageExpression(x, return.seurat = TRUE)
DoHeatmap(Avg, features = top100genes, draw.lines = FALSE, slot = "data")

VlnPlot(x, features = c("H2-DMb1"), pt.size = 1.5, combine = T)+NoLegend()
FeaturePlot(x, features = c("Egr2"), cols = c("lightGrey", "red"), pt.size = 1.5, combine = T)

write.csv(strain.markers, paste0(prefix, ".csv"))
savehistory("R_history.txt")
# Save the session info


#####################################################################rMP-Subtypes#######
prefix <- "rMP_SUBTYPES"
x$celltype <- Idents(x)
saveRDS(subset(x, ident = c("MP")), paste0(prefix, ".rds"))
x <- readRDS(paste0(prefix, ".rds"))
DimPlot(x, label = T)+NoLegend()
DefaultAssay(x) <- "RNA"
x<- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
x<- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
x<- ScaleData(x, verbose = FALSE)
x<- RunPCA(x, npcs = 50, verbose = FALSE)
x<- RunUMAP(x, reduction = "pca", dims = 1:30)
x<- FindNeighbors(x, reduction = "pca", dims = 1:30)
x<- FindClusters(x, resolution = 0.1)
p1 <- DimPlot(x, reduction = "umap" )
p1
x.markers <- FindAllMarkers(x, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
x.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)
write.csv(x.markers, paste0(prefix, ".csv"))
FeaturePlot(x, features = c("Lyve1"), cols = c("lightGrey", "red"), pt.size = 1.5, combine = F)

###################Growth factors for epithelial cell and nerves in the cornea################
features <- c("Igf1","Vegfb","Egf","Fgf1","Fgf2", "Fgfr2","Egr2","Egr1")
DotPlot(x, features =features, cols = c("lightGrey", "darkred"), dot.scale = 6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1))
###############################Resident-MP_Subpopullation############################
wd <- "H:/Manuscripts/2024_K_SCRNA-Methods/Combined_Analysis" # MODIFY
setwd(wd)
prefix <- "rMP_SUBTYPES"
x <- readRDS(paste0(prefix, ".rds"))
DefaultAssay(x)<- ("RNA")
x<- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 10000)
x<- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
x<- ScaleData(x, verbose = FALSE)
x<- RunPCA(x, npcs = 30, verbose = FALSE)
x<- RunUMAP(x, reduction = "pca", dims = 1:30)
x<- FindNeighbors(x, reduction = "pca", dims = 1:30)
x<- FindClusters(x, resolution = 0.2)
p1 <- DimPlot(x, reduction = "umap", label = TRUE, repel = TRUE)
p1
current.cluster.ids <- c(0, 1, 2, 3, 4)
new.cluster.ids <- c("TLF-CCR2+", "TLF-CCR2-", "TLF-CCR2+", "TLF+CCR2-", "4")
x@active.ident <- plyr::mapvalues(x = x@active.ident, from = current.cluster.ids, to = new.cluster.ids)
DimPlot(x, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend() 
#######Cells-count##
x[["RNA_snn_res.0.2"]] <- Idents(object = x)
table <- table(x@meta.data$RNA_snn_res.0.2, x@meta.data$orig.ident)
print(table)
write.csv(table, file.path(getwd(), "rMP_SUBTYPES_Count.csv"))
x.markers <- FindAllMarkers(x, only.pos = T, min.pct = 0.25, logfc.threshold = 0.25)
x.markers %>%
  group_by(cluster) %>%
  slice_max(n = 4, order_by = avg_log2FC)
write.csv(x.markers, paste0(prefix, "_DEG2.csv"))
features <- c("Ccr2","Fgf13","Il7r","Mmp14","Apbb2",
              "Timd4", "Lyve1","Folr2","Fgfr1","Ednrb",
              "Ccl2","Ccl7","Hbegf","Ier3","Pf4",
              "Tyro3","Axl","Mertk","Gas6","Mrc1","Cd163","S100a4","S100a6","Ccl6",
              
              "H2-Eb1","H2-Ab1","Cxcl16","Il1b","Il10","Egr2","Igf1","Vegfb","Csf1r")
features <- c("Mki67","Top2a")
DotPlot(x, features =features, cols = c("green","red"), dot.scale = 6)+ theme(axis.text.x = element_text(angle = 90, vjust = 0.1, hjust=1))

