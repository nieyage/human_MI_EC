library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
library(ggplot2)
library(dittoSeq)
library(dplyr)
library(tidyverse)
library(harmony)
library(cluster)
library(CellChat)
snRNA<- readRDS("./01_all_celltype/all_celltype_snRNA_add_gene_symbol.rds")
snATAC<- readRDS("./01_all_celltype/all_celltype_snATAC_add_peaks_fragment.rds")
coembed<- readRDS("./01_all_celltype/all_celltype_snRNA_snATAC_merged.rds")
# change the celltype name 
Idents(snRNA)<-snRNA$cell_type
#####further annotation########
snRNA <- RenameIdents(
  object = snRNA,
  'neuronal receptor cell' = 'Neuronal',
  'cardiac muscle myoblast' = 'vCMs',
  'smooth muscle myoblast' = 'vSMCs',
  'pericyte' = 'Pericytes',
  'lymphoid lineage restricted progenitor cell' = 'Lymphoid',
  'immature innate lymphoid cell' = 'Myeloid',
  'fibroblast of cardiac tissue' = 'Fibroblasts',
  'cardiac endothelial cell' = 'Endothelial',
  'mast cell' = 'Mast',
  'adipocyte of epicardial fat of left ventricle' = 'Adipocytes',
  'unknown' = 'Cycling cells'
  )
snRNA@meta.data$celltype<-Idents(snRNA)
table(snRNA$celltype,snRNA$major_labl)
saveRDS(snRNA,"./01_all_celltype/all_celltype_snRNA_add_gene_symbol_change_celltype.rds")

p1 <- DimPlot(snRNA,  cols=myUmapcolors,group.by = "celltype",reduction="umap" ,raster=FALSE,label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(snATAC,  cols=myUmapcolors,group.by = "cell_type",reduction="umap" , label = TRUE) + NoLegend() + ggtitle("ATAC")
pdf("./06_cellchat/01_all_celltype/RNA_ATAC_split_Umap.pdf",width=10,height=5)
p1 + p2
plot <- (p1 + p2) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))
plot
dev.off()



# Part I: Data input & processing and initialization of CellChat object
cellChat <- createCellChat(object = snRNA, group.by = "celltype", assay = "RNA")
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# set the used database in the object
cellChat@DB <- CellChatDB
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellChat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692
options(future.globals.maxSize = 2000)
# Part II: Inference of cell-cell communication network
future::plan("multisession", workers = 1) # do parallel
cellchat <- computeCommunProb(cellchat, type = "triMean")
df.net <- subsetCommunication(cellchat) # returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))
pdf("./06_cellchat/01_all_celltype/Number_of_interactions_or_interaction_strength_circleplot.pdf",width=12,height=5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

pdf("./06_cellchat/01_all_celltype/Splited_interaction_strength_circleplot.pdf",width=12,height=12)
mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

pdf("./06_cellchat/01_all_celltype/Splited_number_of_interactions_circleplot.pdf",width=12,height=12)
mat <- cellchat@net$count
par(mfrow = c(3,3), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

saveRDS(cellchat, file = "./06_cellchat/01_all_celltype/all_celltype_snRNA_cellchat.rds")




# Part III: Visualization of cell-cell communication network
# Access all the signaling pathways showing significant communications
pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)
vertex.receiver = seq(1,4)

setwd("./06_cellchat/01_all_celltype/")
#All pathways on the same graph
pdf("All_pathways.pdf",  width = 8, height = 6)
for (i in pathways.show.all) {
  circle<- netVisual_aggregate(object = cellchat, layout = "circle", signaling = as.character(i))
  chord<- netVisual_aggregate(object = cellchat, layout = "chord", signaling = as.character(i))
  h<- netVisual_heatmap(cellchat, signaling = as.character(i), color.heatmap = "Reds")
  contribution<- netAnalysis_contribution(cellchat, signaling = as.character(i));
  print(circle)
  print(chord)
  print(h)
  print(contribution)
}
dev.off()

# 
# for (i in 1:length(pathways.show.all)) {
#   # Visualize communication network associated with both signaling pathway and individual L-R pairs
#   netVisual(cellchat, signaling = pathways.show.all[i], vertex.receiver = vertex.receiver, layout = "hierarchy")
#   # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
#   gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
#   ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"), plot=gg, width = 3, height = 2, units = 'in', dpi = 300)
# }

# Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways

# (1) show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
pdf("ALL_ligand-receptors_bubbleplot.pdf",  width = 10, height = 30)
netVisual_bubble(cellchat, sources.use = c(1:7), targets.use = c(1:7), remove.isolate = FALSE)
dev.off()

pdf("EC_FBvsCM_ligand-receptors_bubbleplot.pdf",  width = 6, height = 20)
netVisual_bubble(cellchat, sources.use = c(1:5), targets.use = c(6:7), remove.isolate = FALSE)
netVisual_bubble(cellchat, sources.use = c(6:7), targets.use = c(1:5), remove.isolate = FALSE)
dev.off()

#> Comparing communications on a single object
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred inte
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
pdf("Signaling role analysis.pdf",  width = 10, height = 15)
ht1 <- netAnalysis_signalingRole_heatmap(cellchat,width = 10,height = 15, pattern = "outgoing")
ht2 <- netAnalysis_signalingRole_heatmap(cellchat,width = 10,height = 15, pattern = "incoming")
ht1 + ht2
dev.off()

# Comparison analysis of multiple datasets with different cell type compositions

library(CellChat)
library(patchwork)
library(Seurat)
library(stringr) 
library(ggplot2)
library(dplyr)
options(stringsAsFactors = FALSE)
set.seed(1234)
setwd("/data/R02/nieyg/project/ORG_CT_scRNA/01_all_celltype/")
combined<- readRDS("./02_annotation/ORG_CT_integrated_annotated.rds")
# add subtype in all object 
CM<- readRDS("./03_CM/CM_subtype_annotated.rds")

meta.data<- data.frame(barcode=c(colnames(CM),colnames(EC),colnames(FB)),
  celltype=c(as.character(CM$subtype),as.character(EC$Annotation),as.character(FB$Annotation)))
rownames(meta.data)<- meta.data[,1]
meta.data<- meta.data[colnames(combined),]
combined$celltype<- meta.data$celltype
combined$samples<- combined$orig.ident

Idents(combined)<- combined$orig.ident;
Nor<- subset(combined,idents="NOR")
DMSO<- subset(combined,idents="DMSO")
CT<- subset(combined,ident="CT")
setwd("/data/R02/nieyg/project/ORG_CT_scRNA/01_all_celltype/04_cellchat")
# Part I: Data input & processing and initialization of CellChat object
cellChat.Nor <- createCellChat(object = Nor, group.by = "celltype", assay = "RNA")
cellChat.DMSO <- createCellChat(object = DMSO, group.by = "celltype", assay = "RNA")
cellChat.CT <- createCellChat(object = CT, group.by = "celltype", assay = "RNA")

CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)
# set the used database in the object
cellChat.Nor@DB <- CellChatDB
cellChat.DMSO@DB <- CellChatDB
cellChat.CT@DB <- CellChatDB

# subset the expression data of signaling genes for saving computation cost
cellChat.Nor <- subsetData(cellChat.Nor) # This step is necessary even if using the whole database
cellChat.DMSO <- subsetData(cellChat.DMSO) # This step is necessary even if using the whole database
cellChat.CT <- subsetData(cellChat.CT) # This step is necessary even if using the whole database

future::plan("multisession", workers = 4) # do parallel
cellChat.Nor <- identifyOverExpressedGenes(cellChat.Nor)
cellChat.DMSO <- identifyOverExpressedGenes(cellChat.DMSO)
cellChat.CT <- identifyOverExpressedGenes(cellChat.CT)

cellChat.Nor <- identifyOverExpressedInteractions(cellChat.Nor)
cellChat.DMSO <- identifyOverExpressedInteractions(cellChat.DMSO)
cellChat.CT <- identifyOverExpressedInteractions(cellChat.CT)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692

# Part II: Inference of cell-cell communication network
cellChat.Nor <- computeCommunProb(cellChat.Nor, type = "triMean")
cellChat.DMSO <- computeCommunProb(cellChat.DMSO, type = "triMean")
cellChat.CT <- computeCommunProb(cellChat.CT, type = "triMean")

cellChat.Nor <- computeCommunProbPathway(cellChat.Nor)
cellChat.DMSO <- computeCommunProbPathway(cellChat.DMSO)
cellChat.CT <- computeCommunProbPathway(cellChat.CT)

cellChat.Nor <- aggregateNet(cellChat.Nor)
cellChat.DMSO <- aggregateNet(cellChat.DMSO)
cellChat.CT <- aggregateNet(cellChat.CT)

saveRDS(cellChat.Nor, file = "cellchat_NOR.rds")
saveRDS(cellChat.DMSO, file = "cellchat_DMSO.rds")
saveRDS(cellChat.CT, file = "cellchat_CT.rds")

object.list <- list(NOR = cellChat.Nor, DMSO = cellChat.DMSO,CT=cellChat.CT)
cellchat <- mergeCellChat(object.list, add.names = names(object.list), cell.prefix = TRUE)
save(object.list, file = "cellchat_object.list_ORG_CT.RData")
save(cellchat, file = "cellchat_merged__ORG_NOR_DMSO_CT.RData")

# Compare the total number of interactions and interaction strength
pdf("./cell_chat_compare/Compare the total number of interactions and interaction strength.pdf",width=10,height=5)
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2,3), measure = "weight")
gg1 + gg2 
dev.off()

# Compare the number of interactions and interaction strength among different cell populations

pdf("./cell_chat_compare/DMSOvsNOR_Compare the number of interactions and interaction strength among different cell populations.pdf",width=10,height=5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat,comparison = c(1, 2), weight.scale = T)
netVisual_diffInteraction(cellchat,comparison = c(1, 2), weight.scale = T, measure = "weight")
gg1 <- netVisual_heatmap(cellchat)
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("./cell_chat_compare/CTvsDMSO_Compare the number of interactions and interaction strength among different cell populations.pdf",width=10,height=5)
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat,comparison = c(2, 3), weight.scale = T)
netVisual_diffInteraction(cellchat,comparison = c(2, 3), weight.scale = T, measure = "weight")
gg1 <- netVisual_heatmap(cellchat,comparison = c(2, 3))
#> Do heatmap based on a merged object
gg2 <- netVisual_heatmap(cellchat,comparison = c(2, 3), measure = "weight")
#> Do heatmap based on a merged object
gg1 + gg2
dev.off()

pdf("./cell_chat_compare/Compare the number of interactions among different cell populations.pdf",width=5,height=5)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}
dev.off()

pdf("./cell_chat_compare/DMSOvsNORvsCT_Compare the interaction strength among different cell populations.pdf",width=8,height=20)
netVisual_bubble(cellchat,title.name = "CM->NCM_max_in_NOR", sources.use = 1:5, targets.use = c(6:7),  comparison = c(1, 2,3), angle.x = 45,max.dataset = 1)
netVisual_bubble(cellchat,title.name = "NCM->CM_max_in_NOR", sources.use = 6:7, targets.use = c(1:5),  comparison = c(1, 2,3), angle.x = 45,max.dataset = 1)
netVisual_bubble(cellchat,title.name = "CM->NCM_max_in_DMSO", sources.use = 1:5, targets.use = c(6:7),  comparison = c(1, 2,3), angle.x = 45,max.dataset = 2)
netVisual_bubble(cellchat,title.name = "NCM->CM_max_in_DMSO", sources.use = 6:7, targets.use = c(1:5),  comparison = c(1, 2,3), angle.x = 45,max.dataset = 2)
netVisual_bubble(cellchat,title.name = "CM->NCM_max_in_CT", sources.use = 1:5, targets.use = c(6:7),  comparison = c(1, 2,3), angle.x = 45,max.dataset = 3)
netVisual_bubble(cellchat,title.name = "NCM->CM_max_in_CT", sources.use = 6:7, targets.use = c(1:5),  comparison = c(1, 2,3), angle.x = 45,max.dataset = 3)
dev.off()
