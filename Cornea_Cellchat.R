#####################Cell-chat##################################
####https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html
library(Seurat)
library(CellChat)
library(patchwork)
library(ggalluvial)
library(NMF)
options(stringsAsFactors = FALSE)
# reticulate::use_python("/Users/suoqinjin/anaconda3/bin/python", required=T) 
###########Create a directory to save figures#########
ptm = Sys.time()

wd <- "H:/Manuscripts/2024_K_SCRNA-Methods/Combined_Analysis" # MODIFY
setwd(wd)
B6<- readRDS("H:/Manuscripts/2024_K_SCRNA-Methods/Combined_Analysis/MFM_Minus_Stromal_legend.rds")
DimPlot(B6, label = T)
data.input <- GetAssayData(B6,assay = "RNA", layer = "counts")
B6_mod <- CreateSeuratObject(counts = data.input)
B6_mod <- NormalizeData(B6_mod)
B6$celltypes <- B6@active.ident
B6_mod$celltypes <- B6@active.ident
meta<-B6@meta.data
cellchat <- createCellChat(object = data.input, meta=meta, group.by = "celltypes")
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf(file = "K_B6_MFM_CELLCHAT", width = 8, height= 10)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
cellchat@netP$pathways
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
gg1
dev.off()
prefix <- "MFM_topactive_celltype"
saveRDS(cellchat, paste0(prefix, ".rds"))
cellchat <- readRDS(paste0(prefix, ".rds"))
saveRDS(sessionInfo(), file = "Cellchat_session_info.rds")
netVisual_bubble(cellchat, sources.use = 1, targets.use = c(1:16), remove.isolate = FALSE)
netVisual_chord_gene(cellchat, sources.use = 13, targets.use = c(5,1,3,6,9), lab.cex = 1,legend.pos.y = 20)
netVisual_chord_gene(cellchat, sources.use = 5, targets.use = c(3,1,6,9,13), signaling = c("CSF"),legend.pos.x = 8)
# show all the significant signaling pathways from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
netVisual_chord_gene(cellchat, sources.use = 1, targets.use = c(3,6,9,13), slot.name = "netP", legend.pos.x = 10)
selectK(cellchat, pattern = "outgoing")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "outgoing", k = nPatterns)
netAnalysis_river(cellchat, pattern = "outgoing")
selectK(cellchat, pattern = "incoming")
nPatterns = 3
cellchat <- identifyCommunicationPatterns(cellchat, pattern = "incoming", k = nPatterns)
netAnalysis_river(cellchat, pattern = "incoming")