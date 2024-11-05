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

snRNA<- readRDS("./01_all_celltype/all_celltype_snRNA_add_gene_symbol.rds")
snATAC<- readRDS("./01_all_celltype/all_celltype_snATAC_add_peaks_fragment.rds")
coembed<- readRDS("./01_all_celltype/all_celltype_snRNA_snATAC_merged.rds")
# For RNA subcluster 
Idents(snRNA)<- snRNA$cell_type;
CM_snRNA <- subset(snRNA,idents="cardiac muscle myoblast")

# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
CM_snRNA[["percent.mt"]] <- PercentageFeatureSet(CM_snRNA, pattern = "^MT-")

pdf("./03_CM/CM_subcluster_QC_plot.pdf",width=10,height=10)
VlnPlot(CM_snRNA, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
plot2 <- FeatureScatter(CM_snRNA, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot2
plot(density(CM_snRNA@meta.data$nCount_RNA),xlim=c(0,38000))
plot(density(CM_snRNA@meta.data$nFeature_RNA),xlim=c(0,6500))
dev.off()
CM_snRNA <- subset(CM_snRNA, subset = nFeature_RNA > 2000 & nFeature_RNA < 5000 & nCount_RNA < 20000 & nCount_RNA > 5000)

integrated_data <- SplitObject(CM_snRNA, split.by = "sample")

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

integrated_data <- NormalizeData(integrated_data)
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
                              max.iter.harmony = 20)
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

resolutions <- seq(0.2, 2, 0.2)
integrated_data <- integrated_data %>% 
  FindNeighbors(reduction = "harmony") %>% 
  FindClusters(reduction = "harmony",
    reduction.type="cca.aligned",
    k.param=15,
    resolution = resolutions) %>% 
  RunUMAP(reduction = "harmony", 
         dims = 1:20,min.dist=0.1,metric = 'euclidean',
          reduction.name = "umap_harmony") %>% 
  identity()

library(clustree)
pdf("./03_CM/resolutions_select_CM_subcluster_reintergrated.pdf",width=15,height=15)
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
DimPlot(object = integrated_data,reduction = "umap_harmony", group.by = "sample")
DimPlot(object = integrated_data,reduction = "umap_harmony", group.by = "patient_group")
DimPlot(object = integrated_data,reduction = "umap_harmony", group.by = "major_labl")
DimPlot(object = integrated_data,reduction = "umap_harmony", group.by = "region")
dev.off()

integrated_data$seurat_clusters<- integrated_data$RNA_snn_res.1
Idents(integrated_data)<- integrated_data$seurat_clusters

pdf("./03_CM/snRNA_CM_Unsupervised_cluster_Umap.pdf",width=6,height=5)
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "seurat_clusters")
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "sample")
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "patient_group")
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "major_labl")
DimPlot(object = integrated_data,cols=myUmapcolors,pt.size=0.7, group.by = "region")
dev.off()


# make tree for cluster
object<- integrated_data
embeddings <- Embeddings(object = object, reduction = "harmony")[,1:30]
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
pdf("/md01/nieyg/project/EC/human_MI/03_CM/snRNA_CM_subcluster_harmony_tree.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


CM_snRNA<- integrated_data
# The dot plot for annotation 
DefaultAssay(CM_snRNA)<- "RNA"
features <- c("LINC02208","AC022034.3","AC0107068.2","PDE1A","FHL2","GPC5","GRIK2","CORIN","FILIP1L","MGAT4C","ANKRD1",# CM1
    "LMOD2","TCAP","MYH7","ACTC1","MYL2","CRYAB","TNNI3","MYL12A","TXNIP",# CM2
    "NRXN3","GBE1","MYH9","ADAMTS9","NPPB","FLNC","SLC2A13","ATP2B4","LINC01091",# CM3
    "PCDH7","PLCG2","FN1","SPARCL1","COL6A3","TECRL","DLC1","LINC02552",# CM4
    "C1QTNF1","RRAD","CAMK1D","NAMPT","VIPR2","PDE4B","STAT3","A2M","FNBP1L","THBS1" # CM5
    )
CM_snRNA<- ScaleData(CM_snRNA)
pdf("/md01/nieyg/project/EC/human_MI/03_CM/snRNA_CM_subcluster_annotation_dotplot.pdf",width=8,height=3)
p1<- DotPlot(CM_snRNA, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
#p1&scale_x_discrete(labels=features_label)
dev.off()

Idents(CM_snRNA)<-CM_snRNA$seurat_clusters
#####further annotation########
CM_snRNA <- RenameIdents(
  object = CM_snRNA,
  '0' = 'CM2',
  '1' = 'CM1',
  '2' = 'CM1',
  '3' = 'CM3',
  '4' = 'CM4',
  '5' = 'CM1',
  '6' = 'CM3',
  '7' = 'CM2',
  '8' = 'CM1',
  '9' = 'CM5',
  '10' = 'CM2',
  '11' = 'CM3'
  )
CM_snRNA@meta.data$subcluster<-Idents(CM_snRNA)
table(CM_snRNA$subcluster,CM_snRNA$major_labl)

pdf("/md01/nieyg/project/EC/human_MI/03_CM/snRNA_CM_subcluster_annotation_umap.pdf",width=6,height=5)
DimPlot(CM_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "subcluster")
DimPlot(CM_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "major_labl")
dev.off()

DefaultAssay(CM_snRNA) <- "RNA"
saveRDS(CM_snRNA,"/md01/nieyg/project/EC/human_MI/03_CM/snRNA_CM_subtype_harmony.rds")

# proportion for subcluster
library(cowplot)
CM_snRNA$major_labl<- factor(CM_snRNA$major_labl,levels=c("CTRL","IZ","BZ","RZ","FZ"))
pdf("/md01/nieyg/project/EC/human_MI/03_CM/snRNA_CM_subtype_proportion.pdf",width=5,height=5)
df <- as.data.frame(CM_snRNA@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, subcluster) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=subcluster)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, subcluster) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=subcluster)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# major_labl
df_ct <- df %>%
    group_by(major_labl, subcluster) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=subcluster)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()

CM_snRNA$subcluster<- factor(CM_snRNA$subcluster,levels=c("CM1","CM2","CM3","CM4","CM5"))
Idents(CM_snRNA)<- CM_snRNA$subcluster
DefaultAssay(CM_snRNA)<- "RNA"
features <- c("LINC02208","AC022034.3","AC0107068.2","PDE1A","FHL2","GPC5","GRIK2","CORIN","FILIP1L","MGAT4C","ANKRD1",# CM1
    "LMOD2","TCAP","MYH7","ACTC1","MYL2","CRYAB","TNNI3","MYL12A","TXNIP",# CM2
    "NRXN3","GBE1","MYH9","ADAMTS9","NPPB","FLNC","SLC2A13","ATP2B4","LINC01091",# CM3
    "PCDH7","PLCG2","FN1","SPARCL1","COL6A3","TECRL","DLC1","LINC02552",# CM4
    "C1QTNF1","RRAD","CAMK1D","NAMPT","VIPR2","PDE4B","STAT3","A2M","FNBP1L","THBS1" # CM5
    )
CM_snRNA<- ScaleData(CM_snRNA)
pdf("/md01/nieyg/project/EC/human_MI/03_CM/snRNA_CM_subcluster_annotated_dotplot.pdf",width=10,height=4)
p1<- DotPlot(CM_snRNA, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
#p1&scale_x_discrete(labels=features_label)
dev.off()



CM_snRNA<- readRDS("/md01/nieyg/project/EC/human_MI/03_CM/snRNA_CM_subtype_harmony.rds")

# For ATAC subcluster 
Idents(snATAC)<- snATAC$cell_type
CM_snATAC <- subset(snATAC,idents="cardiac muscle myoblast")
DefaultAssay(CM_snATAC)<- "ATAC"
library(reticulate)
Python_Script_version<<-paste0("/data/R02/nieyg/ori/biosoft/conda/bin/python")
Sys.setenv(RETICULATE_PYTHON= Python_Script_version)
reticulate ::use_python(Python_Script_version, required = TRUE)
Sys.which("python")
reticulate::py_config()
sc <- import("scopen")
matDR <- sc$Main$scopen_dr(counts = CM_snATAC@assays$peaks@counts,
                           verbose = 1)
matDR <- t(as.matrix(matDR))
colnames(matDR) <- paste0("PC_", 1:ncol(matDR))
rownames(matDR) <- colnames(CM_snATAC@assays$peaks@counts)
CM_snATAC@reductions[['scopen']] <- CreateDimReducObject(embeddings = matDR,
                                                        assay = "peaks",
                                                       key = "PC_")
DepthCor(CM_snATAC, reduction = "scopen", n = 30)
CM_snATAC <- RunUMAP(object = CM_snATAC, 
                    reduction = 'scopen', 
                    dims = 1:20, 
                    verbose = FALSE)

suppressMessages(library(harmony))
CM_snATAC <- RunHarmony(
  object = CM_snATAC,
  group.by.vars = "sample",
  reduction.use='scopen',
  #reduction = 'scopen',
  project.dim = FALSE,
  assay.use = 'peaks',
  plot_convergence = FALSE,
  verbose = TRUE    
)
CM_snATAC <- RunUMAP(CM_snATAC, 
               dims = 1:20, 
               reduction = 'harmony',
               reduction.name = "umap_harmony",
               reduction.ke = 'umapharmony_',
               #metric = 'euclidean',
              verbose = FALSE)

pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_CM_subcluster_reduction_umap.pdf",width=15,height=5)
p1 <- DimPlot(object = CM_snATAC, reduction = "umap",
              group.by = "sample") +
    xlab("UMAP1") + ylab("UMAP2")
p2 <- DimPlot(object = CM_snATAC, reduction = "umap_harmony",
              group.by = "sample", shuffle = TRUE) +
    xlab("UMAP1") + ylab("UMAP2")
p1 + p2
dev.off()

resolutions <- seq(0.2, 1, 0.1)
CM_snATAC <- FindNeighbors(CM_snATAC, reduction = "harmony", dims = 1:20)
CM_snATAC <- FindClusters(CM_snATAC,graph.name ='peaks_snn', resolution = resolutions, verbose = FALSE)

library(clustree)
pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_CM_subcluster_resolutions_select.pdf",width=15,height=15)
clustree(CM_snATAC, prefix = "peaks_snn_res.")
plotlist <- lapply(resolutions, function(x){
    cols <- myUmapcolors
    p <- DimPlot(CM_snATAC, group.by = glue::glue("peaks_snn_res.{x}"), label = TRUE,
             reduction = "umap_harmony", shuffle = TRUE) +
    scale_color_manual(values = cols) +
    xlab("UMAP1") + ylab("UMAP2")
    p
})
p <- patchwork::wrap_plots(plotlist, ncol = 3)
p
dev.off()

CM_snATAC$seurat_clusters<- CM_snATAC$peaks_snn_res.0.4
Idents(CM_snATAC)<- CM_snATAC$seurat_clusters


# unsupervised clustering umap plot 
pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_CM_subcluster_umap.pdf",width=5,height=5)
DimPlot(object = CM_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "seurat_clusters",label = T)
DimPlot(object = CM_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "patient_group",label = T)
DimPlot(object = CM_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "patient_region",label = T)
DimPlot(object = CM_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "region",label = T)
DimPlot(object = CM_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "major_labl",label = T)
DimPlot(object = CM_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "cell_type",label = T)
dev.off();

# make tree for cluster
object<- CM_snATAC
Idents(CM_snATAC)<- CM_snATAC$peaks_snn_res.0.4
embeddings <- Embeddings(object = object, reduction = "harmony")[,1:30]
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
pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_CM_subcluster_harmony_tree.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


# call peak 
DefaultAssay(CM_snATAC)<-"ATAC"
peak<-CallPeaks(
       CM_snATAC,
       group.by = "seurat_clusters",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 2.7e+09,
       outdir="/data/R02/nieyg/project/EC/human_MI/03_CM/peaks",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(CM_snATAC),
     features = peak,
     cells = colnames(CM_snATAC)
     )     

# create a new assay using the MACS2 peak set and add it to the Seurat object
CM_snATAC[["peaks_CM"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(CM_snATAC),
  annotation = Annotation(CM_snATAC)
)

DefaultAssay(CM_snATAC)<-"peaks_CM"
gene.activities <- GeneActivity(CM_snATAC,features=rownames(CM_snRNA))
# add gene activities as a new assay
CM_snATAC[["ACTIVITY_CM"]] <- CreateAssayObject(counts = gene.activities)
saveRDS(CM_snATAC,"./03_CM/CM_snATAC_peakcalling.rds")
DefaultAssay(CM_snATAC)<-"RNA_imputation"
# The dot plot for annotation 
features <- c("LINC02208","AC022034.3","AC0107068.2","PDE1A","FHL2","GPC5","GRIK2","CORIN","FILIP1L","MGAT4C","ANKRD1",# CM1
    "LMOD2","TCAP","MYH7","ACTC1","MYL2","CRYAB","TNNI3","MYL12A","TXNIP",# CM2
    "NRXN3","GBE1","MYH9","ADAMTS9","NPPB","FLNC","SLC2A13","ATP2B4","LINC01091",# CM3
    "PCDH7","PLCG2","FN1","SPARCL1","COL6A3","TECRL","DLC1","LINC02552",# CM4
    "C1QTNF1","RRAD","CAMK1D","NAMPT","VIPR2","PDE4B","STAT3","A2M","FNBP1L","THBS1" # CM5
    )
CM_snATAC<- ScaleData(CM_snATAC)
pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_CM_subcluster_annotation_dotplot.pdf",width=10,height=4)
p1<- DotPlot(CM_snATAC, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
#p1&scale_x_discrete(labels=features_label)
dev.off()


Idents(CM_snATAC)<-CM_snATAC$peaks_nn_res.0.8
#####further annotation########
CM_snATAC <- RenameIdents(
  object = CM_snATAC,
  '0' = '',
  '1' = '',
  '2' = '',
  '3' = '',
  '4' = '',
  '5' = '',
  '6' = '',
  '7' = '',
  '8' = '',
  '9' = '',
  '10' = '',
  '11' = ''
  )
CM_snATAC@meta.data$subcluster<-Idents(CM_snATAC)
table(CM_snATAC$subcluster,CM_snATAC$major_labl)

pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_EC_subcluster_annotation_umap.pdf",width=6,height=5)
DimPlot(CM_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "subcluster")
DimPlot(CM_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "subcluster_CAV")
dev.off()

DefaultAssay(CM_snATAC) <- "peaks"
saveRDS(CM_snATAC,"/md01/nieyg/project/EC/human_MI/03_CM/snATAC_All_EC_harmony.rds")

CM_snATAC<- readRDS("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_All_EC_harmony.rds")


# proportion for subcluster

pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_All_EC_subtype_proportion.pdf",width=10,height=5)
df <- as.data.frame(CM_snATAC@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, subcluster_CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=subcluster_CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, subcluster) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=subcluster)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# major_labl
df_ct <- df %>%
    group_by(major_labl, subcluster) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=subcluster)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()


pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_All_EC_not_CAV_proportion.pdf",width=10,height=5)
df <- as.data.frame(CM_snRNA@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, subcluster) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=subcluster)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, subcluster_CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=subcluster_CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# major_labl
df_ct <- df %>%
    group_by(major_labl, subcluster_CAV) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=subcluster_CAV)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

dev.off()


# merge snRNA and snATAC 

# normalize gene activities
DefaultAssay(CM_snATAC) <- "ACTIVITY"
CM_snATAC <- CM_snATAC %>% 
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA(verbose = FALSE)

CM_snRNA <- FindVariableFeatures(CM_snRNA)


# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = CM_snRNA, 
    query = CM_snATAC,
    features = VariableFeatures(object = CM_snRNA),
    reference.assay = "RNA",
    query.assay = "ACTIVITY", 
    reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
	refdata = CM_snRNA$subcluster,
    weight.reduction = CM_snATAC[["scopen"]], 
    dims = 1:30)

CM_snATAC <- AddMetaData(CM_snATAC, metadata = celltype.predictions)
CM_snATAC$subcluster_annotation_correct <- CM_snATAC$predicted.id == CM_snATAC$subcluster
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = CM_snRNA$subcluster_CAV,
    weight.reduction = CM_snATAC[["scopen"]], dims = 1:30)
CM_snATAC$predicted_CAV<- CM_snATAC$predicted.id

CM_snATAC <- AddMetaData(CM_snATAC, metadata = celltype.predictions)
CM_snATAC$subcluster_CAV_annotation_correct <- CM_snATAC$predicted.id == CM_snATAC$subcluster_CAV


CM_snATAC$predicted_CAV<- factor(CM_snATAC$predicted_CAV,levels=c("Venous_Endo","C_V","Capillary_Endo","C_A","Arterial_Endo","Endocardial_Endo","Lymphatic_Endo","notsure"))
CM_snATAC$subcluster<- factor(CM_snATAC$subcluster,levels=c("Venous_Endo","C_V","Capillary_Endo","C_A","Arterial_Endo","Endocardial_Endo","Lymphatic_Endo","notsure"))

pdf("/md01/nieyg/project/EC/human_MI/03_CM/All_EC_ATAC_predicted_Ground-truth_Umap.pdf",width=10,height=5)
p1 <- DimPlot(CM_snATAC,reduction='umap_harmony' ,cols=myUmapcolors, group.by = "predicted_CAV", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(CM_snATAC,reduction='umap_harmony' ,cols=myUmapcolors,group.by = "subcluster", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2
p1 <- DimPlot(CM_snATAC,reduction='umap_harmony' ,cols=myUmapcolors, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(CM_snATAC,reduction='umap_harmony' ,cols=myUmapcolors,group.by = "subcluster_CAV", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2
dev.off()

# not split CAV 
predictions <- table(CM_snATAC$subcluster_CAV, CM_snATAC$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

order_Var1<- c(levels(Idents(CM_snATAC)))
order_Var2<- c(levels(Idents(CM_snATAC)))
predictions$Var1<- factor(predictions$Var1,levels=order_Var1)
predictions$Var2<- factor(predictions$Var2,levels=order_Var2)
library(cowplot)
pdf("/md01/nieyg/project/EC/human_MI/03_CM/All_EC_ATAC_predicted_rate_notCAV.pdf",width=10,height=5)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(CM_snATAC$subcluster_CAV == CM_snATAC$predicted.id))
incorrect <- length(which(CM_snATAC$subcluster_CAV != CM_snATAC$predicted.id))
data <- FetchData(CM_snATAC, vars = c("prediction.score.max", "subcluster_CAV_annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = subcluster_CAV_annotation_correct, colour = subcluster_CAV_annotation_correct)) +
    geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2
dev.off()

saveRDS(CM_snATAC,"/md01/nieyg/project/EC/human_MI/03_CM/snATAC_All_EC_harmony_predictions.rds")


# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(CM_snRNA)
refdata <- GetAssayData(CM_snRNA, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, 
    refdata = refdata, 
    weight.reduction = "cca")
CM_snATAC[["RNA"]] <- imputation


coembed <- merge(x = CM_snRNA, y = CM_snATAC)
coembed<- RenameCells(coembed,
    old.names=colnames(coembed),
    new.names =paste(coembed$sample,gsub("_[1,2]*","",colnames(coembed)),sep="#"))

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- coembed %>%
    ScaleData(features = genes.use, do.scale = FALSE) %>%
    RunPCA(features = genes.use,, verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

pdf("/md01/nieyg/project/EC/human_MI/03_CM/coembed_umap.pdf",width=20,height=5)
DimPlot(coembed,reduction='umap',cols=myUmapcolors,raster=FALSE, group.by = c("patient_region", "subcluster",'tech'))
DimPlot(coembed,reduction='umap',cols=myUmapcolors,raster=FALSE, group.by = c("patient_region", "subcluster_CAV",'tech'))
dev.off()
saveRDS(coembed,"./03_CM/All_CM_snRNA_snATAC_merged.rds")

CM_snRNA<- readRDS("./03_CM/snRNA_All_EC_harmony.rds")
CM_snATAC<- readRDS("./03_CM/snATAC_All_EC_harmony.rds")
coembed<- readRDS("./03_CM/All_CM_snRNA_snATAC_merged.rds")

# The features between notsure and CAV in snATAC 
    CM_snATAC <- NucleosomeSignal(CM_snATAC)
    CM_snATAC <- TSSEnrichment(CM_snATAC,fast=FALSE)

pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_features_umap.pdf",width=24,height=20)
FeaturePlot(
  object = CM_snATAC,features = c("nCount_RNA","nFeature_RNA","TSS.enrichment","nucleosome_signal",
  	"nCount_ATAC","nFeature_ATAC","nCount_peaks","nFeature_peaks",
  	"nCount_ACTIVITY","nFeature_ACTIVITY","nCounts_RNA","nFeaturess_RNA"),
  order = FALSE,reduction='umap_harmony' ,
  cols = c("grey", "darkred"),
  ncol = 4
) & NoLegend()
dev.off()
Idents(CM_snATAC)<- CM_snATAC$subcluster_CAV
pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_EC_QC_vlnplot.pdf",width=10,height=16)
VlnPlot(object = CM_snATAC,features = c("nCount_RNA","nFeature_RNA","TSS.enrichment","nucleosome_signal"), ncol = 1, pt.size = 0)
VlnPlot(object = CM_snATAC,features = c("nCount_ATAC","nFeature_ATAC","nCount_peaks","nFeature_peaks"), ncol = 1, pt.size = 0)
VlnPlot(object = CM_snATAC,features = c("nCount_ACTIVITY","nFeature_ACTIVITY","nCounts_RNA","nFeaturess_RNA"), ncol = 1, pt.size = 0)
dev.off()


# The DEP in cluster 12 and other CAV

Idents(CM_snATAC) <- CM_snATAC$subcluster_CAV

# change back to working with peaks instead of gene activities
DefaultAssay(CM_snATAC) <- 'peaks'

da_peaks <- FindMarkers(
  object = CM_snATAC,
  ident.1 = "CAV",
  ident.2 = "notsure",
  test.use = 'LR',
  latent.vars = 'nCount_peaks',
  logfc.threshold = 0,
  min.pct = 0,
  pseudocount.use = 0.01
)
head(da_peaks)

VolcanoPlot=function(dif, log2FC=log2(1.5), padj=0.05, 
                 label.symbols=NULL, label.max=30,
                 cols=c("#497aa2", "#ae3137"), title=""){
  if( all( !c("log2FoldChange", "padj", "symbol") %in% colnames(dif) )){
    stop("Colnames must include: log2FoldChange, padj, symbol")
  }
  rownames(dif)=dif$symbol
  
  # (1) define up and down
  dif$threshold="ns";
  dif[which(dif$log2FoldChange > log2FC & dif$padj <padj),]$threshold="up";
  dif[which(dif$log2FoldChange < (-log2FC) & dif$padj < padj),]$threshold="down";
  dif$threshold=factor(dif$threshold, levels=c('down','ns','up'))
  #head(dif)
  #
  tb2=table(dif$threshold); print(tb2)
  library(ggplot2)
  # (2) plot
  g1 = ggplot(data=dif, aes(x=log2FoldChange, y=-log10(padj), color=threshold)) +
    geom_point(alpha=0.8, size=0.8) +
    geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
    geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
    labs(title= ifelse(""==title, "", paste("DEG:", title)))+
    xlab(bquote(Log[2]*FoldChange))+
    ylab(bquote(-Log[10]*italic(P.adj)) )+
    theme_classic(base_size = 14) +
    theme(legend.box = "horizontal",
          legend.position="top",
          legend.spacing.x = unit(0, 'pt'),
          legend.text = element_text( margin = margin(r = 20) ),
          legend.margin=margin(b= -10, unit = "pt"),
          plot.title = element_text(hjust = 0.5, size=10)
          ) +
    scale_color_manual('',labels=c(paste0("down(",tb2[[1]],')'),'ns',
                                   paste0("up(",tb2[[3]],')' )),
                       values=c(cols[1], "grey", cols[2]) )+
    guides(color=guide_legend(override.aes = list(size=3, alpha=1))); g1;
  # (3)label genes
  if(is.null(label.symbols)){
    dif.sig=dif[which(dif$threshold != "ns" ), ]
    len=nrow(dif.sig)
    if(len<label.max){
      label.symbols=rownames(dif.sig)
    }else{
      dif.sig=dif.sig[order(dif.sig$log2FoldChange), ]
      dif.sig= rbind(dif.sig[1:(label.max/2),], dif.sig[(len-label.max/2):len,])
      label.symbols=rownames(dif.sig)
    }
  }
  dd_text = dif[label.symbols, ]
  print((dd_text))
  # add text
  library(ggrepel)
  g1 + geom_text_repel(data=dd_text,
                       aes(x=log2FoldChange, y=-log10(padj), label=row.names(dd_text)),
                         #size=2.5, 
                         colour="black",alpha=1)
}

# plot the vlocano plot 
dif=data.frame(
  symbol=rownames(da_peaks),
  log2FoldChange=da_peaks$avg_log2FC,
  padj=da_peaks$p_val_adj
)

pdf("/md01/nieyg/project/EC/human_MI/03_CM/snATAC_EC_CAV_vs_notsure_volcano.pdf",width=10,height=10)
p1 = VolcanoPlot(dif, padj=0.05, title="CAV vs notsure", label.max = 50)
p1
dev.off()






fc <- FoldChange(CM_snATAC, ident.1 = "CAV", ident.2 = "notsure")
# order by fold change
fc <- fc[order(fc$avg_log2FC, decreasing = TRUE), ]
head(fc)

open_CAV <- rownames(da_peaks[da_peaks$avg_log2FC > 3, ])
open_notsure <- rownames(da_peaks[da_peaks$avg_log2FC < -3, ])

closest_genes_cd4naive <- ClosestFeature(pbmc, regions = open_cd4naive)
closest_genes_cd14mono <- ClosestFeature(pbmc, regions = open_cd14mono)

# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
mouse_brain <- AddMotifs(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)
# test enrichment
enriched.motifs <- FindMotifs(
  object = mouse_brain,
  features = top.da.peak
)

# plot the notsure and other TSS distribution
CM_snATAC<- readRDS("./03_CM/snATAC_All_EC_harmony.rds")
CM_snATAC <- TSSEnrichment(object = CM_snATAC, fast = FALSE)
Idents(CM_snATAC) <- CM_snATAC$subcluster_CAV
CM_snATAC$TSS_group <- ifelse(CM_snATAC$subcluster_CAV == "notsure", 'notsure', 'other')
pdf("./03_CM/snATAC_notsure_TSS_distribution.pdf",width=12,height=6)
TSSPlot(CM_snATAC, group.by = 'TSS_group') + NoLegend()
TSSPlot(CM_snATAC, group.by = 'subcluster_CAV') + NoLegend()
dev.off()

CM_snATAC <- NucleosomeSignal(object = CM_snATAC)
CM_snATAC$nucleosome_group <- ifelse(CM_snATAC$nucleosome_signal == "notsure", 'notsure', 'other')
pdf("./03_CM/snATAC_notsure_NucleosomeSignal_distribution.pdf",width=12,height=6)
FragmentHistogram(object = CM_snATAC, group.by = 'nucleosome_group')
FragmentHistogram(object = CM_snATAC, group.by = 'subcluster_CAV')
dev.off()





