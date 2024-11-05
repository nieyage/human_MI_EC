
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
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

EC_snRNA<- readRDS("./02_EC/snRNA_All_EC_harmony.rds")
EC_snATAC<- readRDS("./02_EC/snATAC_All_EC_harmony_predictions.rds")

# For RNA subcluster 
Idents(EC_snRNA)<- EC_snRNA$subcluster_CAV
CAV_snRNA <- subset(EC_snRNA,idents="CAV")
integrated_data <- SplitObject(CAV_snRNA, split.by = "sample")

# Calculate HVG per sample - Here we assume that batch and patient effects aren't as strong
hvg_list <- map(integrated_data, function(x) {
  DefaultAssay(x) <- "RNA"
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", 
                       nfeatures = 3000)
  x@assays[["RNA"]]@var.features 
}) %>% unlist()
hvg_list <- table(hvg_list) %>% sort(decreasing = TRUE)
gene_selection <- hvg_list[1:3000] %>% names()

# We collapse again
integrated_data <- reduce(integrated_data,
                          merge,
                          merge.data = TRUE)

# Quickly get characteristic profile of the object
integrated_data <- integrated_data %>%
  ScaleData(verbose = FALSE) %>% 
  RunPCA(features = gene_selection, 
         npcs = 30, 
         verbose = FALSE) 

# Integrate the data -----------------------
integrated_data <- RunHarmony(integrated_data, 
                              c("sample", "donor_id"), 
                              plot_convergence = TRUE,
                              assay.use = "RNA",
                              plot_convergence=TRUE,
                              max.iter.harmony = 20)

# Create the UMAP with new reduction -----------
resolutions <- seq(0.2, 2, 0.2)
integrated_data <- integrated_data %>% 
  FindNeighbors(reduction = "harmony") %>% 
  FindClusters(reduction = "harmony",
    reduction.type="cca.aligned",
    k.param=15,
    resolution = resolutions) %>% 
  RunUMAP(reduction = "harmony", dims = 1:30,min.dist=0.05,metric = 'euclidean',
          reduction.name = "umap_harmony") %>% 
  identity()

library(clustree)
pdf("./02_EC/CAV/resolutions_select_CM_subcluster_reintergrated.pdf",width=15,height=15)
clustree(integrated_data, prefix = "RNA_snn_res.")
plotlist <- lapply(resolutions, function(x){
    cols <- myUmapcolors
    p <- DimPlot(integrated_data, group.by = glue::glue("RNA_snn_res.{x}"), label = TRUE,
             reduction = "umap_harmony", shuffle = TRUE) +
    scale_color_manual(values = cols) +
    xlab("UMAP1") + ylab("UMAP2")
    p
})
p <- patchwork::wrap_plots(plotlist, ncol = 3)
p
DimPlot(object = integrated_data, group.by = "sample")
DimPlot(object = integrated_data, group.by = "patient_group")
DimPlot(object = integrated_data, group.by = "major_labl")
DimPlot(object = integrated_data, group.by = "region")
dev.off()

integrated_data$seurat_clusters<- integrated_data$RNA_snn_res.1
Idents(integrated_data)<- integrated_data$seurat_clusters

pdf("./02_EC/CAV/snRNA_CAV_Unsupervised_cluster_Umap.pdf",width=6,height=5)
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "seurat_clusters")
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "sample")
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "patient_group")
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "major_labl")
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "region")
dev.off()

table(integrated_data$seurat_clusters)
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# unsupervised clustering umap plot 
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_EC_subcluster_umap.pdf",width=12,height=10)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = c("seurat_clusters","patient_group","region","major_labl"))
dev.off();

# plot the features 
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_EC_subcluster_umap_features.pdf",width=10,height=10)
FeaturePlot(
  object = integrated_data,
  features = c("nCounts_RNA","nFeaturess_RNA","nCount_RNA","nFeature_RNA"),
  order = FALSE,
  reduction='umap_harmony' ,
  cols = c("grey", "darkred"),
  ncol = 2
) & NoLegend()
dev.off()

# make tree for cluster
object<- integrated_data
embeddings <- Embeddings(object = object, reduction = "pca")[,1:30]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_subcluster_harmony_tree.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


CAV_snRNA<- integrated_data
# The dot plot for annotation 
features <- c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","CD36","MGLL",# capillary Endo
	"FAM155A","TSHZ2","SLCO2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2","POSTN",# venous Endo
	"PCSK5","ARL15","LINC00639","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2","SLC45A4")# Arterial Endo)
features_label<-c(rep("capillary",10),rep("venous",9),rep("Arterial",10))

EC_atlas<- c("FABP4",# heart 
	"IFIT3","IFIT1","ISG15","STAT1","USP18",# interferon Endo
	"SPARCL1","APLNR","ADM","APLN","COL15A1","COL4A2",# angiogenic Endo
	"PRSS23","FXYD6","CP","PROX1","PDPN","LYVE1",# Lymphatic EC
	"FBLN5","HEY1","MECOM",# Artery 
	"CXCL12","RBP7",# C_A
	"KDR","ENDOU",# capillary
	"VCAM1","FMO1","FMO2",# C_V
	"MGP","CFH","BGN","VWF"# large vein
	)
EC_atlas_label<-c("heart_EC",rep("interferon",5),rep("angiogenic",6),
	rep("Lymphatic",6),rep("Artery",3),rep("C_A",2),
	rep("capillary",2),rep("C_V",3),rep("large_vein",4))

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_subcluster_annotation_dotplot.pdf",width=8,height=6)
p1<- DotPlot(CAV_snRNA, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
p1&scale_x_discrete(labels=features_label)
p2<- DotPlot(CAV_snRNA, features = EC_atlas,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p2
p2&scale_x_discrete(labels=EC_atlas_label)
dev.off()

saveRDS(CAV_snRNA,"./02_EC/CAV/snRNA_CAV_subcluster.rds")

CAV_snRNA<- readRDS("./02_EC/CAV/snRNA_CAV_subcluster.rds")
Idents(CAV_snRNA)<-CAV_snRNA$seurat_clusters
#####further annotation########
CAV_snRNA <- RenameIdents(
  object = CAV_snRNA,
  '0' = 'C_A',
  '1' = 'C_A',
  '2' = 'Capillary',
  '3' = 'Capillary',
  '4' = 'Capillary',
  '5' = 'Artery',
  '6' = 'Venous',
  '7' = 'Artery',
  '8' = 'Venous',
  '9' = 'C_V',
  '10' = 'Venous',
  '11' = 'Capillary',
  '12' = 'Venous',
  '13' = 'C_A'
  )
CAV_snRNA@meta.data$CAV<-Idents(CAV_snRNA)
CAV_snRNA$CAV<- factor(CAV_snRNA$CAV,levels=c("Artery","C_A","Capillary","C_V","Venous"))
table(CAV_snRNA$CAV,CAV_snRNA$major_labl)

Idents(CAV_snRNA)<-CAV_snRNA$seurat_clusters
#####further annotation########
CAV_snRNA <- RenameIdents(
  object = CAV_snRNA,
  '0' = 'C_A_1',
  '1' = 'C_A_2',
  '2' = 'Capillary_1',
  '3' = 'Capillary_2',
  '4' = 'Capillary_3',
  '5' = 'Artery_1',
  '6' = 'Venous_1',
  '7' = 'Artery_2',
  '8' = 'Venous_2',
  '9' = 'C_V',
  '10' = 'Venous_3',
  '11' = 'Capillary_4',
  '12' = 'Venous_4',
  '13' = 'C_A_3'
  )
CAV_snRNA@meta.data$CAV_detail<-Idents(CAV_snRNA)
CAV_snRNA$CAV_detail<- factor(CAV_snRNA$CAV_detail,
    levels=c("Artery_1","Artery_2","C_A_1","C_A_2","C_A_3",
        "Capillary_1","Capillary_2","Capillary_3","Capillary_4",
        "C_V","Venous_1","Venous_2","Venous_3","Venous_4"))

table(CAV_snRNA$CAV_detail,CAV_snRNA$major_labl)


pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_subcluster_annotated_umap.pdf",width=6,height=5)
DimPlot(CAV_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "CAV")
DimPlot(CAV_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "CAV_detail")
DimPlot(CAV_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "sample")
DimPlot(CAV_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "major_labl")
DimPlot(CAV_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "patient_group")
DimPlot(CAV_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "patient_group")
DimPlot(CAV_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "development_stage")
dev.off()

pdf("./02_EC/CAV/snRNA_CAV_Unsupervised_cluster_Umap.pdf",width=6,height=5)
DimPlot(object = CAV_snRNA,cols=myUmapcolors,pt.size=0.7, group.by = "seurat_clusters")
DimPlot(object = CAV_snRNA,cols=myUmapcolors,pt.size=0.7, group.by = "sample")
DimPlot(object = CAV_snRNA,cols=myUmapcolors,pt.size=0.7, group.by = "patient_group")
DimPlot(object = CAV_snRNA,cols=myUmapcolors,pt.size=0.7, group.by = "major_labl")
DimPlot(object = CAV_snRNA,cols=myUmapcolors,pt.size=0.7, group.by = "region")
dev.off()

# unsupervised clustering umap plot 
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_EC_subcluster_umap.pdf",width=12,height=10)
DimPlot(object = CAV_snRNA,cols=myUmapcolors, group.by = c("seurat_clusters","patient_group","region","major_labl"))
dev.off();

# The dot plot for annotation 
features <- c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","CD36","MGLL",# capillary Endo
    "FAM155A","TSHZ2","SLCO2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2","POSTN",# venous Endo
    "PCSK5","ARL15","LINC00639","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2","SLC45A4")# Arterial Endo)
features_label<-c(rep("capillary",10),rep("venous",9),rep("Arterial",10))

EC_atlas<- c("FABP4",# heart 
    "IFIT3","IFIT1","ISG15","STAT1","USP18",# interferon Endo
    "SPARCL1","APLNR","ADM","APLN","COL15A1","COL4A2",# angiogenic Endo
    "PRSS23","FXYD6","CP","PROX1","PDPN","LYVE1",# Lymphatic EC
    "FBLN5","HEY1","MECOM",# Artery 
    "CXCL12","RBP7",# C_A
    "KDR","ENDOU",# capillary
    "VCAM1","FMO1","FMO2",# C_V
    "MGP","CFH","BGN","VWF"# large vein
    )
EC_atlas_label<-c("heart_EC",rep("interferon",5),rep("angiogenic",6),
    rep("Lymphatic",6),rep("Artery",3),rep("C_A",2),
    rep("capillary",2),rep("C_V",3),rep("large_vein",4))

Idents(CAV_snRNA)<- CAV_snRNA$CAV
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_subcluster_annotated_dotplot.pdf",width=8,height=4)
p1<- DotPlot(CAV_snRNA, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
p1&scale_x_discrete(labels=features_label)
p2<- DotPlot(CAV_snRNA, features = EC_atlas,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p2
p2&scale_x_discrete(labels=EC_atlas_label)
dev.off()

# plot the features 
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_EC_subcluster_umap_features.pdf",width=10,height=10)
FeaturePlot(
  object = CAV_snRNA,
  features = c("nCounts_RNA","nFeaturess_RNA","nCount_RNA","nFeature_RNA"),
  order = FALSE,
  reduction='umap_harmony' ,
  cols = c("grey", "darkred"),
  ncol = 2
) & NoLegend()
dev.off()
DefaultAssay(CAV_snRNA) <- "RNA"
saveRDS(CAV_snRNA,"/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_EC_CAV_harmony.rds")

CAV_snRNA<- readRDS("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_EC_CAV_harmony.rds")
# proportion for subcluster
library(cowplot)

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_subtype_proportion.pdf",width=10,height=5)
df <- as.data.frame(CAV_snRNA@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# major_labl
df_ct <- df %>%
    group_by(major_labl, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_detail_proportion.pdf",width=10,height=5)
df <- as.data.frame(CAV_snRNA@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, CAV_detail) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=CAV_detail)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, CAV_detail) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=CAV_detail)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# major_labl
df_ct <- df %>%
    group_by(major_labl, CAV_detail) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=CAV_detail)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

dev.off()

# For ATAC subcluster 
EC_snATAC<- readRDS("./02_EC/snATAC_All_EC_harmony.rds")
Idents(EC_snATAC)<- EC_snATAC$subcluster_CAV
CAV_snATAC <- subset(EC_snATAC,idents="CAV")
DefaultAssay(CAV_snATAC)<- "ATAC"

CAV_snATAC <- RunUMAP(CAV_snATAC, 
               dims = 1:30, 
               reduction = 'harmony',
               reduction.name = "umap_harmony",
               reduction.ke = 'umapharmony_',
               verbose = FALSE,
               min.dist = 0.1)

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subcluster_reduction_umap.pdf",width=15,height=5)
p1 <- DimPlot(object = CAV_snATAC, reduction = "umap",
              group.by = "sample") +
    xlab("UMAP1") + ylab("UMAP2")
p2 <- DimPlot(object = CAV_snATAC, reduction = "umap_harmony",
              group.by = "sample", shuffle = TRUE) +
    xlab("UMAP1") + ylab("UMAP2")
p1 + p2
dev.off()
resolutions <- seq(0.2, 2, 0.2)
CAV_snATAC <- FindNeighbors(CAV_snATAC, reduction = "harmony", dims = 1:30)
CAV_snATAC <- FindClusters(CAV_snATAC,graph.name ='peaks_snn', resolution = resolutions, verbose = FALSE)
CAV_snATAC <- RunUMAP(CAV_snATAC, 
               dims = 1:20, 
               reduction = 'harmony',
               reduction.name = "umap_harmony",
               reduction.ke = 'umapharmony_',
               verbose = FALSE,metric = 'euclidean',
               min.dist = 0.05)

library(clustree)
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subcluster_resolutions_select.pdf",width=15,height=15)
clustree(CAV_snATAC, prefix = "peaks_snn_res.")
plotlist <- lapply(resolutions, function(x){
    cols <- myUmapcolors
    p <- DimPlot(CAV_snATAC, group.by = glue::glue("peaks_snn_res.{x}"), label = TRUE,
             reduction = "umap_harmony", shuffle = TRUE) +
    scale_color_manual(values = cols) +
    xlab("UMAP1") + ylab("UMAP2")
    p
})
p <- patchwork::wrap_plots(plotlist, ncol = 3)
p
dev.off()

CAV_snATAC$seurat_clusters<- CAV_snATAC$peaks_snn_res.0.8
Idents(CAV_snATAC)<- CAV_snATAC$seurat_clusters
# unsupervised clustering umap plot 
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subcluster_umap.pdf",width=5,height=5)
DimPlot(object = CAV_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "seurat_clusters",label = T)
DimPlot(object = CAV_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "patient_group",label = T)
DimPlot(object = CAV_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "patient_region",label = T)
DimPlot(object = CAV_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "region",label = T)
DimPlot(object = CAV_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "major_labl",label = T)
DimPlot(object = CAV_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "cell_type",label = T)
dev.off();

# make tree for cluster
object<- CAV_snATAC
embeddings <- Embeddings(object = object, reduction = "harmony")[,1:20]
data.dims <- lapply(X = levels(x = object), FUN = function(x) {
    cells <- WhichCells(object = object, idents = x)
    if (length(x = cells) == 1) {
        cells <- c(cells, cells)
    }
    temp <- colMeans(x = embeddings[cells, ])
})
data.dims <- do.call(what = "cbind", args = data.dims)
colnames(x = data.dims) <- levels(x = object)
library(lsa)
cosine_dist <- as.dist(1-cosine(data.dims))
data.tree <- ape::as.phylo(x = hclust(d = cosine_dist))
library(ggtree);
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subcluster_harmony_tree.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


# call peak 
DefaultAssay(CAV_snATAC)<-"ATAC"
peak<-CallPeaks(
       CAV_snATAC,
       group.by = "seurat_clusters",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 2.7e+09,
       outdir="/data/R02/nieyg/project/EC/human_MI/02_EC/CAV/peaks",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(CAV_snATAC),
     features = peak,
     cells = colnames(CAV_snATAC)
     )     
#macs2_counts<-macs2_counts[-which(rownames(macs2_counts)=="GroupUN243-311766-311965"),]
# ATAC analysis add gene annotation information

# create a new assay using the MACS2 peak set and add it to the Seurat object
CAV_snATAC[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(CAV_snATAC),
  annotation = Annotation(CAV_snATAC)
)

DefaultAssay(CAV_snATAC)<-"peaks"
gene.activities <- GeneActivity(CAV_snATAC,features=rownames(CAV_snRNA))
# add gene activities as a new assay
CAV_snATAC[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

DefaultAssay(CAV_snATAC)<-"ACTIVITY"
# The dot plot for annotation 
features <- c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","CD36","MGLL",# capillary Endo
    "FAM155A","TSHZ2","SLCO2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2","POSTN",# venous Endo
    "PCSK5","ARL15","LINC00639","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2","SLC45A4")# Arterial Endo)
features_label<-c(rep("capillary",10),rep("venous",9),rep("Arterial",10))

EC_atlas<- c("FABP4",# heart 
    "IFIT3","IFIT1","ISG15","STAT1","USP18",# interferon Endo
    "SPARCL1","APLNR","ADM","APLN","COL15A1","COL4A2",# angiogenic Endo
    "PRSS23","FXYD6","CP","PROX1","PDPN","LYVE1",# Lymphatic EC
    "FBLN5","HEY1","MECOM",# Artery 
    "CXCL12","RBP7",# C_A
    "KDR","ENDOU",# capillary
    "VCAM1","FMO1","FMO2",# C_V
    "MGP","CFH","BGN","VWF"# large vein
    )
EC_atlas_label<-c("heart_EC",rep("interferon",5),rep("angiogenic",6),
    rep("Lymphatic",6),rep("Artery",3),rep("C_A",2),
    rep("capillary",2),rep("C_V",3),rep("large_vein",4))

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subcluster_annotation_dotplot.pdf",width=8,height=6)
p1<- DotPlot(CAV_snATAC,assay="ACTIVITY" ,features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
p1&scale_x_discrete(labels=features_label)
p1<- DotPlot(CAV_snATAC,assay="RNA" ,features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
p1&scale_x_discrete(labels=features_label)
p1<- DotPlot(CAV_snATAC,assay="RNA_imputation" ,features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
p1&scale_x_discrete(labels=features_label)
dev.off()

Idents(CAV_snATAC)<-CAV_snATAC$seurat_clusters
#####further annotation########
CAV_snATAC <- RenameIdents(
  object = CAV_snATAC,
  '0' = 'C_A',
  '1' = 'Capillary',
  '2' = 'Capillary',
  '3' = 'Venous',
  '4' = 'Artery',
  '5' = 'Capillary',
  '6' = 'C_A',
  '7' = 'C_V',
  '8' = 'notsure',
  '9' = 'notsure',
  '10' = 'Capillary',
  '11' = 'notsure',
  '12' = 'Venous'
  )
CAV_snATAC@meta.data$CAV<-Idents(CAV_snATAC)
CAV_snATAC$CAV<- factor(CAV_snATAC$CAV,levels=c("Artery","C_A","Capillary","C_V","Venous","notsure"))
table(CAV_snATAC$CAV,CAV_snATAC$major_labl)

# remove the notsure
CAV_snATAC2<- subset(CAV_snATAC,idents=c("Artery","C_A","Capillary","C_V","Venous"))
table(CAV_snATAC2$CAV,CAV_snATAC2$major_labl)
CAV_snATAC<- CAV_snATAC2
CAV_snATAC$CAV<- factor(CAV_snATAC$CAV,levels=c("Artery","C_A","Capillary","C_V","Venous"))

# merge snRNA and snATAC 

# normalize gene activities
DefaultAssay(CAV_snATAC) <- "ACTIVITY"
CAV_snATAC <- CAV_snATAC %>% 
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA(verbose = FALSE)

CAV_snRNA <- FindVariableFeatures(CAV_snRNA)

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = CAV_snRNA, 
    query = CAV_snATAC,
    features = VariableFeatures(object = CAV_snRNA),
    reference.assay = "RNA",
    query.assay = "ACTIVITY", 
    reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = CAV_snRNA$CAV,
    weight.reduction = CAV_snATAC[["harmony"]], dims = 1:30)

CAV_snATAC <- AddMetaData(CAV_snATAC, metadata = celltype.predictions)
CAV_snATAC$annotation_correct <- CAV_snATAC$predicted.id == CAV_snATAC$cell_type
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/ATAC_predicted_Ground-truth_Umap.pdf",width=10,height=5)
p1 <- DimPlot(CAV_snATAC,reduction='umap_harmony' ,cols=myUmapcolors, group.by = "predicted.id", label = TRUE)  + ggtitle("Predicted annotation")
p2 <- DimPlot(CAV_snATAC,reduction='umap_harmony' ,cols=myUmapcolors,group.by = "CAV", label = TRUE)  + ggtitle("Ground-truth annotation")
p1 | p2
dev.off()

predictions <- table(CAV_snATAC$CAV, CAV_snATAC$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

order_Var1<- c(levels(Idents(CAV_snATAC)))
order_Var2<- c(levels(Idents(CAV_snATAC)))
predictions$Var1<- factor(predictions$Var1,levels=order_Var1)
predictions$Var2<- factor(predictions$Var2,levels=order_Var2)

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/ATAC_predicted_rate.pdf",width=10,height=5)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(CAV_snATAC$CAV == CAV_snATAC$predicted.id))
incorrect <- length(which(CAV_snATAC$CAV != CAV_snATAC$predicted.id))
data <- FetchData(CAV_snATAC, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
    geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2
dev.off()



Idents(CAV_snATAC)<-CAV_snATAC$seurat_clusters
table(CAV_snATAC$seurat_clusters)
#####further annotation########
CAV_snATAC <- RenameIdents(
  object = CAV_snATAC,
  '0' = 'C_A_1',
  '1' = 'Capillary_1',
  '2' = 'Capillary_2',
  '3' = 'Venous_1',
  '4' = 'Artery',
  '5' = 'Capillary_3',
  '6' = 'C_A_2',
  '7' = 'C_V',
  '10' = 'Capillary_4',
  '12' = 'Venous_2'
  )
CAV_snATAC@meta.data$CAV_detail<-Idents(CAV_snATAC)
CAV_snATAC$CAV_detail<- factor(CAV_snATAC$CAV_detail,
    levels=c("Artery","C_A_1","C_A_2",
        "Capillary_1","Capillary_2","Capillary_3","Capillary_4",
        "C_V","Venous_1","Venous_2"))

table(CAV_snRNA$CAV_detail,CAV_snRNA$major_labl)


pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subcluster_annotated_umap.pdf",width=6,height=5)
DimPlot(CAV_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "CAV")
DimPlot(CAV_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "CAV_detail")
DimPlot(CAV_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "sample")
DimPlot(CAV_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "major_labl")
DimPlot(CAV_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "patient_group")
DimPlot(CAV_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "patient_group")
DimPlot(CAV_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "development_stage")
dev.off()

pdf("./02_EC/CAV/snATAC_CAV_Unsupervised_cluster_Umap.pdf",width=6,height=5)
DimPlot(object = CAV_snATAC,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.7, group.by = "seurat_clusters")
DimPlot(object = CAV_snATAC,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.7, group.by = "sample")
DimPlot(object = CAV_snATAC,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.7, group.by = "patient_group")
DimPlot(object = CAV_snATAC,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.7, group.by = "major_labl")
DimPlot(object = CAV_snATAC,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.7, group.by = "region")
dev.off()

# unsupervised clustering umap plot 
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_EC_subcluster_umap.pdf",width=12,height=10)
DimPlot(object = CAV_snATAC,cols=myUmapcolors,reduction = "umap_harmony", group.by = c("seurat_clusters","patient_group","region","major_labl"))
dev.off();

# The dot plot for annotation 
features <- c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","CD36","MGLL",# capillary Endo
    "FAM155A","TSHZ2","SLCO2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2","POSTN",# venous Endo
    "PCSK5","ARL15","LINC00639","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2","SLC45A4")# Arterial Endo)
features_label<-c(rep("capillary",10),rep("venous",9),rep("Arterial",10))
EC_atlas<- c("FABP4",# heart 
    "IFIT3","IFIT1","ISG15","STAT1","USP18",# interferon Endo
    "SPARCL1","APLNR","ADM","APLN","COL15A1","COL4A2",# angiogenic Endo
    "PRSS23","FXYD6","CP","PROX1","PDPN","LYVE1",# Lymphatic EC
    "FBLN5","HEY1","MECOM",# Artery 
    "CXCL12","RBP7",# C_A
    "KDR","ENDOU",# capillary
    "VCAM1","FMO1","FMO2",# C_V
    "MGP","CFH","BGN","VWF"# large vein
    )
EC_atlas_label<-c("heart_EC",rep("interferon",5),rep("angiogenic",6),
    rep("Lymphatic",6),rep("Artery",3),rep("C_A",2),
    rep("capillary",2),rep("C_V",3),rep("large_vein",4))

Idents(CAV_snATAC)<- CAV_snATAC$CAV
pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subcluster_annotated_dotplot.pdf",width=8,height=4)
p1<- DotPlot(CAV_snATAC, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
p1&scale_x_discrete(labels=features_label)
p1<- DotPlot(CAV_snATAC,assay="RNA", features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
p1&scale_x_discrete(labels=features_label)
dev.off()

DefaultAssay(CAV_snATAC) <- "ATAC"
saveRDS(CAV_snATAC,"/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_harmony_annotated.rds")

# proportion for subcluster

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subtype_proportion.pdf",width=10,height=5)
df <- as.data.frame(CAV_snATAC@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# major_labl
df_ct <- df %>%
    group_by(major_labl, CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()


pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_proportion.pdf",width=10,height=5)
df <- as.data.frame(CAV_snATAC@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, CAV_detail) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=CAV_detail)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, CAV_detail) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=CAV_detail)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# major_labl
df_ct <- df %>%
    group_by(major_labl, CAV_detail) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=CAV_detail)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

dev.off()


# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(CAV_snRNA)
refdata <- GetAssayData(CAV_snRNA, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, 
    refdata = refdata, 
    weight.reduction = "cca")
CAV_snATAC[["RNA"]] <- imputation


coembed <- merge(x = CAV_snRNA, y = CAV_snATAC)
coembed<- RenameCells(coembed,
    old.names=colnames(coembed),
    new.names =paste(coembed$sample,gsub("_[1,2]*","",colnames(coembed)),sep="#"))

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- coembed %>%
    ScaleData(features = genes.use, do.scale = FALSE) %>%
    FindVariableFeatures() %>%
    RunPCA() %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/coembed_umap.pdf",width=15,height=5)
DimPlot(coembed,reduction='umap',cols=myUmapcolors,raster=FALSE, group.by = c("patient_region", "CAV"))
dev.off()

saveRDS(CAV_snATAC,"./02_EC/CAV/EC_CAV_snATAC_add_peaks_fragment.rds")
saveRDS(CAV_snRNA,"./02_EC/CAV/EC_CAV_snRNA_add_gene_symbol.rds")
saveRDS(coembed,"./02_EC/CAV/EC_CAV_snRNA_snATAC_merged.rds")


# check the snATAC annotation by trackplot of CAV marker genes
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
CAV_snRNA<- readRDS("./02_EC/CAV/EC_CAV_snRNA_add_gene_symbol.rds")
CAV_snATAC<- readRDS("./02_EC/CAV/EC_CAV_snATAC_add_peaks_fragment.rds")

features <- c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","CD36","MGLL",# capillary Endo
    "FAM155A","TSHZ2","SLCO2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2","POSTN",# venous Endo
    "PCSK5","ARL15","LINC00639","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2","SLC45A4")# Arterial Endo)
features_label<-c(rep("capillary",10),rep("venous",9),rep("Arterial",10))

Idents(CAV_snATAC)<-CAV_snATAC$CAV
library(BSgenome.Hsapiens.UCSC.hg38)
DefaultAssay(CAV_snATAC) <- "ACTIVITY"
gene<- rownames(CAV_snATAC)

DefaultAssay(CAV_snATAC) <- "peaks"
# first compute the GC content for each peak
CAV_snATAC <- RegionStats(CAV_snATAC, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
CAV_snATAC <- LinkPeaks(
  object = CAV_snATAC,
  peak.assay = "peaks",
  expression.assay = "ACTIVITY",
  genes.use = features
)
######Visulize track and RNA exp######
idents.plot <- Idents(CAV_snATAC)

pdf("./02_EC/CAV/EC_CAV_Marker_gene-peaktrack-RNAexp.pdf",height=6,width=8)
for(i in features[which(features%in% gene)]){
  print(i)
  p1 <- CoveragePlot(
  object = CAV_snATAC,
  region = i,
  features = i,
  expression.assay = "ACTIVITY",
  idents = idents.plot,
  extend.upstream = 1000,
  annotation=TRUE,
  extend.downstream = 1000
)
print(p1)}
dev.off()


# The human marker gene in https://academic.oup.com/cardiovascres/article/118/14/2960/6536962

features_CardioR <- c("CA4","INMT","NEURL1B","APLN","MEOX1","CDH13","MXRA5","IGFBP5",# capillary Endo
    "SERPINE2","PRND","SRGN","FBLN5","PIK3R3","APOA1","CXCR4","SMAD6",#A
    "INHBA","NPPC","LEPR","PLAC9","ATP1B1","CYP1B1","TGM2","EMILIN1",# Endocardial
    "KNL1","PCLAF","HIST1H1B","H2AFX","UBE2S","HIST1H1D","ARL6IP1","HIST1H4C",# proliferating 
    "ACKR1","VCAN","DUSP23","PTGIS","MMRN1","CD55",# V
    "BAMBI","SELL","C7","AKR1B1","IER3","SAT1","COL1A1",# valvular 
    "TFF3","NTS","MFAP4","CCL21","NR2F1","CLU","ARL4A","CTSL")# Arterial Endo)

# link peaks to genes
CAV_snATAC <- LinkPeaks(
  object = CAV_snATAC,
  peak.assay = "peaks",
  expression.assay = "ACTIVITY",
  genes.use = features_CardioR
)

pdf("./02_EC/CAV/EC_CAV_CardioR_Marker_gene-peaktrack-RNAexp.pdf",height=6,width=8)
for(i in features_CardioR[which(features_CardioR%in% gene)]){
  print(i)
  p1 <- CoveragePlot(
  object = CAV_snATAC,
  region = i,
  features = i,
  expression.assay = "ACTIVITY",
  idents = idents.plot,
  extend.upstream = 1000,
  annotation=TRUE,
  extend.downstream = 1000
)
print(p1)}
dev.off()

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subcluster_annotated_dotplot_CardioR_Marker.pdf",width=10,height=4)
DotPlot(CAV_snATAC,group.by="CAV",assay="ACTIVITY", features = features_CardioR,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
DotPlot(CAV_snATAC,assay="RNA",group.by="CAV", features = features_CardioR,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
DotPlot(CAV_snATAC,group.by="CAV_detail",assay="ACTIVITY", features = features_CardioR,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
DotPlot(CAV_snATAC,assay="RNA",group.by="CAV_detail", features = features_CardioR,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_subcluster_annotated_dotplot_CardioR_Marker.pdf",width=10,height=4)
DotPlot(CAV_snRNA,assay="RNA",group.by="CAV", features = features_CardioR,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
DotPlot(CAV_snRNA,assay="RNA",group.by="CAV_detail", features = features_CardioR,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()

Idents(CAV_snATAC)<- CAV_snATAC$CAV

# The marker gene in supply excel of Spatial multi-omic map of human myocardial infarction 
Findmarkes<- c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","CD36","MGLL",
"FAM155A","TSHZ2","SLCO2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2","LINC02147","POSTN",
"PCSK5","ARL15","LINC00639","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2","SLC45A4",
"LINC02388","TMEM132C","PKHD1L1","NRG3","BMPER","TMEM108","NRG1","CGNL1","CDH11","PLA2G5",
"SEMA3A","RELN","MMRN1","PIEZO2","CCL21","DOCK5","TFPI","PKHD1L1","FLT4","AKAP12")
Findmarkes<- unique(Findmarkes)
CAV_snATAC <- LinkPeaks(
  object = CAV_snATAC,
  peak.assay = "peaks",
  expression.assay = "ACTIVITY",
  genes.use = Findmarkes
)
pdf("./02_EC/CAV/EC_CAV_Findmarkes_gene-peaktrack-RNAexp.pdf",height=16,width=8)
for(i in Findmarkes[which(Findmarkes%in% gene)]){
  print(i)
  p1 <- CoveragePlot(
  object = CAV_snATAC,
  region = i,
  features = i,
  expression.assay = "ACTIVITY",
  idents = idents.plot,
  extend.upstream = 1000,
  annotation=TRUE,
  extend.downstream = 1000
)
print(p1)}
dev.off()

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snATAC_CAV_subcluster_annotated_dotplot_Findmarkes.pdf",width=10,height=4)
DotPlot(CAV_snATAC,group.by="CAV",assay="ACTIVITY", features = Findmarkes,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
DotPlot(CAV_snATAC,assay="RNA",group.by="CAV", features = Findmarkes,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
DotPlot(CAV_snATAC,group.by="CAV_detail",assay="ACTIVITY", features = Findmarkes,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
DotPlot(CAV_snATAC,assay="RNA",group.by="CAV_detail", features = Findmarkes,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()

pdf("/md01/nieyg/project/EC/human_MI/02_EC/CAV/snRNA_CAV_subcluster_annotated_dotplot_Findmarkes.pdf",width=10,height=4)
DotPlot(CAV_snRNA,assay="RNA",group.by="CAV", features = Findmarkes,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
DotPlot(CAV_snRNA,assay="RNA",group.by="CAV_detail", features = Findmarkes,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
dev.off()

DefaultAssay(CAV_snATAC)<-"ACTIVITY"
gene<- rownames(CAV_snATAC)
DefaultAssay(CAV_snATAC)<-"peaks"
known_markers<- c("VEGFA","VEGFB","WNT","NOTCH1","NOTCH3","NOTCH4","HIF1A","MECOM","KDR")
CAV_snATAC <- LinkPeaks(
  object = CAV_snATAC,
  peak.assay = "peaks",
  expression.assay = "ACTIVITY",
  genes.use = known_markers
)
Idents(CAV_snATAC)<- CAV_snATAC$CAV
idents.plot <- Idents(CAV_snATAC)

pdf("./02_EC/CAV/known_markers_EC_CAV_gene-peaktrack-RNAexp.pdf",height=10,width=8)
for(i in known_markers[which(known_markers%in% gene)]){
  print(i)
  p1 <- CoveragePlot(
  object = CAV_snATAC,
  region = i,
  features = i,
  expression.assay = "ACTIVITY",
  idents = idents.plot,
  extend.upstream = 1000,
  annotation=TRUE,
  extend.downstream = 1000
)
print(p1)}
dev.off()




# Cell cycle score of CAV in snRNA
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(CAV_snRNA)<-"RNA"
CAV_snRNA <- CellCycleScoring(CAV_snRNA, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(CAV_snRNA[[]])
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
CAV_snRNA <- ScaleData(CAV_snRNA, features = c(s.genes, g2m.genes))
CAV_snRNA <- RunPCA(CAV_snRNA, features = c(s.genes, g2m.genes))

pdf("./02_EC/CAV/CAV_snRNA_cellcycle_umap.pdf")
DimPlot(CAV_snRNA)
dev.off()

# Cell cycle score of CAV in snATAC 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
DefaultAssay(CAV_snATAC)<-"RNA"
CAV_snATAC <- CellCycleScoring(CAV_snATAC, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
# view cell cycle scores and phase assignments
head(CAV_snATAC[[]])
# Running a PCA on cell cycle genes reveals, unsurprisingly, that cells separate entirely by
# phase
CAV_snATAC <- ScaleData(CAV_snATAC, features = c(s.genes, g2m.genes))
CAV_snATAC <- RunPCA(CAV_snATAC, features = c(s.genes, g2m.genes))

pdf("./02_EC/CAV/monocle/CAV_snATAC_cellcycle_umap.pdf")
DimPlot(CAV_snATAC)
dev.off()

# ref: https://www.nature.com/articles/s41467-022-30545-8
Idents(CAV_snRNA)<- CAV_snRNA$CAV

pdf("./02_EC/CAV/snRNA_proliferative_marker_violin.pdf",width=8,height=20)
VlnPlot(object = CAV_snRNA, features = features,ncol=1)
dev.off()






features<- c(
"AURKB","MKI67",
"NR2F2","VEGFB","VEGFA","NOTCH1") 

Idents(CAV_snRNA)<- CAV_snRNA$CAV
pdf("./02_EC/CAV/snRNA_proliferative_marker_violin.pdf",width=8,height=20)
VlnPlot(object = CAV_snRNA, features = features, split.by = 'region_dai',ncol=1)
dev.off()

features<- c(
"AURKB","MKI67",# proliferative markers
"FAM174B","ADAMTS18","HEY2","ANTXR1","RASGRF2","NOTCH1","NR2F2","PSMD9","ZNRD2","DCTN6","GJA4") 
CAV_snRNA@meta.data$newgroup <- paste0(CAV_snRNA$CAV,"_",CAV_snRNA$region_dai)
Idents(CAV_snRNA) <- 'newgroup'
Idents(CAV_snRNA) <- factor(Idents(CAV_snRNA), levels = c(paste0("Artery","_",levels(CAV_snRNA$region_dai)),
                                                            paste0("C_A","_",levels(CAV_snRNA$region_dai)),
                                                            paste0("Capillary","_",levels(CAV_snRNA$region_dai)),
                                                            paste0("C_V","_",levels(CAV_snRNA$region_dai)),
                                                            paste0("Venous","_",levels(CAV_snRNA$region_dai))))        

pdf("./02_EC/CAV/snRNA_proliferative_marker_dotplot.pdf",width=12,height=12)
DotPlot(CAV_snRNA, features = features, cols = c("lightgrey", "red"))
DotPlot(CAV_snRNA, features = cc.genes$s.genes, cols = c("lightgrey", "red"))
DotPlot(CAV_snRNA, features = cc.genes$g2m.genes, cols = c("lightgrey", "red"))
dev.off()




pdf("./02_EC/CAV/CAV_splitby_region_Umap.pdf",width=20,height=5)
DimPlot(object = CAV_snRNA,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.4,split.by = "major_labl", group.by = "CAV")
DimPlot(object = CAV_snATAC,reduction = "umap_harmony",cols=myUmapcolors,pt.size=0.8,split.by = "major_labl", group.by = "CAV")
dev.off()








