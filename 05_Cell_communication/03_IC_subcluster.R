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
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

# For RNA subcluster 
Idents(snRNA)<- snRNA$cell_type_original
IC_snRNA <- subset(snRNA,idents="Myeloid")
integrated_data <- SplitObject(IC_snRNA, split.by = "sample")

# Calculate HVG per sample - Here we assume that batch and patient effects aren't as strong
hvg_list <- map(integrated_data, function(x) {
  DefaultAssay(x) <- "RNA"
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
                              max.iter.harmony = 20)

# Create the UMAP with new reduction -----------
resolutions <- seq(0.2, 2, 0.2)
integrated_data <- integrated_data %>% 
  FindNeighbors(reduction = "harmony") %>% 
  FindClusters(reduction = "harmony",
    reduction.type="cca.aligned",
    k.param=15,
    resolution = resolutions) %>% 
  RunUMAP(reduction = "harmony", 
         dims = 1:10,min.dist=0.1,metric = 'euclidean',
          reduction.name = "umap_harmony") %>% 
  identity()

library(clustree)
pdf("./05_IC/resolutions_select_IC_subcluster_reintergrated.pdf",width=15,height=15)
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

integrated_data$seurat_clusters<- integrated_data$RNA_snn_res.0.6
Idents(integrated_data)<- integrated_data$seurat_clusters

# unsupervised clustering umap plot 
pdf("/md01/nieyg/project/EC/human_MI/05_IC/snRNA_IC_subcluster_umap.pdf",width=5,height=5)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "seurat_clusters",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "patient_group",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "patient_region",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "region",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "major_labl",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "cell_type",label = T)
dev.off();

# make tree for cluster
object<- integrated_data
embeddings <- Embeddings(object = object, reduction = "pca")[,1:10]
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
pdf("/md01/nieyg/project/EC/human_MI/05_IC/snRNA_IC_subcluster_harmony_tree.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


IC_snRNA<- integrated_data
# The dot plot for annotation 
features <- c("MAMDC2","PDE4D","SCN9A","RBPJ","F13A1","P2RY14","FGF13","COLEC12","NAV2","PDGFC",#LYVE_FOLR_MP
    "THBS1","PLTP","NEURL2","CR1","STAB1","FAM20A","NAMPT","TFRC","MARCO",#LYVE_PLTP_MP
    "IL1R2","CIITA","PID1","ADAM28","RTN1","FLT3","CD74","HLA-DQA1","AFF3","DOCK10",#Monocytes
    "ITGAX","MMP19","SDS","FAM20C","SPP1","TPRG1","TM4SF19","ZNF331","CD36","PLAUR",#SPP1_MP
    "KCNMA1","HS3ST2","NHSL1","CPM","GALNTL6","SNTB1","ANKRD28","ALOX15B","KCNE1","MYO1E"#CCL18_MP
    )
pdf("/md01/nieyg/project/EC/human_MI/05_IC/snRNA_IC_subcluster_annotation_dotplot.pdf",width=12,height=8)
p1<- DotPlot(IC_snRNA, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
dev.off()

Idents(IC_snRNA)<-IC_snRNA$seurat_clusters
#####further annotation########
IC_snRNA <- RenameIdents(
  object = IC_snRNA,
  '0' = 'LYVE_PLTP_MP',
  '1' = 'LYVE_FOLR_MP',
  '2' = 'LYVE_PLTP_MP',
  '3' = 'Monocytes',
  '4' = 'LYVE_FOLR_MP',
  '5' = 'SPP1_MP',
  '6' = 'SPP1_MP',
  '7' = 'CCL18_MP',
  '8' = 'LYVE_PLTP_MP',
  '9' = 'Monocytes'
  )
IC_snRNA@meta.data$subcluster<-Idents(IC_snRNA)
table(IC_snRNA$subcluster,IC_snRNA$major_labl)

pdf("/md01/nieyg/project/EC/human_MI/05_IC/snRNA_IC_subcluster_annotation_umap.pdf",width=6,height=5)
DimPlot(IC_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "subcluster")
DimPlot(IC_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "major_labl")
dev.off()

pdf("/md01/nieyg/project/EC/human_MI/05_IC/snRNA_IC_subcluster_annotated_dotplot.pdf",width=10,height=3)
p1<- DotPlot(IC_snRNA, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
dev.off()

DefaultAssay(IC_snRNA) <- "RNA"
saveRDS(IC_snRNA,"/md01/nieyg/project/EC/human_MI/05_IC/snRNA_All_MP_harmony.rds")

# proportion for subcluster

pdf("/md01/nieyg/project/EC/human_MI/05_IC/snRNA_MP_subtype_proportion.pdf",width=6,height=5)
df <- as.data.frame(IC_snRNA@meta.data)
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


IC_snRNA<- readRDS("/md01/nieyg/project/EC/human_MI/05_IC/snRNA_MP_harmony.rds")

# For ATAC subcluster 
Idents(snATAC)<- snATAC$cell_type
IC_snATAC <- subset(snATAC,idents="cardiac endothelial cell")
DefaultAssay(IC_snATAC)<- "ATAC"

library(reticulate)
Python_Script_version<<-paste0("/data/R02/nieyg/ori/biosoft/conda/bin/python")
Sys.setenv(RETICULATE_PYTHON= Python_Script_version)
reticulate ::use_python(Python_Script_version, required = TRUE)

Sys.which("python")
reticulate::py_config()

sc <- import("scopen")
matDR <- sc$Main$scopen_dr(counts = IC_snATAC@assays$peaks@counts,
                           verbose = 1)

matDR <- t(as.matrix(matDR))

colnames(matDR) <- paste0("PC_", 1:ncol(matDR))
rownames(matDR) <- colnames(IC_snATAC@assays$peaks@counts)
IC_snATAC@reductions[['scopen']] <- CreateDimReducObject(embeddings = matDR,
                                                        assay = "peaks",
                                                       key = "PC_")
DepthCor(IC_snATAC, reduction = "scopen", n = 30)
IC_snATAC <- RunUMAP(object = IC_snATAC, 
                    reduction = 'scopen', 
                    dims = 1:30, 
                    verbose = FALSE)

suppressMessages(library(harmony))

IC_snATAC <- RunHarmony(
  object = IC_snATAC,
  group.by.vars = "sample",
  reduction.use='scopen',
  #reduction = 'scopen',
  project.dim = FALSE,
  assay.use = 'peaks',
  plot_convergence = FALSE,
  verbose = TRUE    
)
IC_snATAC <- RunUMAP(IC_snATAC, 
               dims = 1:30, 
               reduction = 'harmony',
               reduction.name = "umap_harmony",
               reduction.ke = 'umapharmony_',
              verbose = FALSE,
                   min.dist = 0.4)

pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_EC_subcluster_reduction_umap.pdf",width=15,height=5)
p1 <- DimPlot(object = IC_snATAC, reduction = "umap",
              group.by = "sample") +
    xlab("UMAP1") + ylab("UMAP2")
p2 <- DimPlot(object = IC_snATAC, reduction = "umap_harmony",
              group.by = "sample", shuffle = TRUE) +
    xlab("UMAP1") + ylab("UMAP2")
p1 + p2
dev.off()

resolutions <- seq(0.2, 1, 0.1)
IC_snATAC <- FindNeighbors(IC_snATAC, reduction = "harmony", dims = 1:30)
IC_snATAC <- FindClusters(IC_snATAC,graph.name ='peaks_nn', resolution = resolutions, verbose = FALSE)

library(clustree)
pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_EC_subcluster_resolutions_select.pdf",width=15,height=15)
clustree(IC_snATAC, prefix = "peaks_snn_res.")
plotlist <- lapply(resolutions, function(x){
    cols <- myUmapcolors
    p <- DimPlot(IC_snATAC, group.by = glue::glue("peaks_snn_res.{x}"), label = TRUE,
             reduction = "umap_harmony", shuffle = TRUE) +
    scale_color_manual(values = cols) +
    xlab("UMAP1") + ylab("UMAP2")
    p
})
p <- patchwork::wrap_plots(plotlist, ncol = 3)
p
dev.off()

# unsupervised clustering umap plot 
pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_EC_subcluster_umap.pdf",width=5,height=5)
DimPlot(object = IC_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "peaks_nn_res.0.8",label = T)
DimPlot(object = IC_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "patient_group",label = T)
DimPlot(object = IC_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "patient_region",label = T)
DimPlot(object = IC_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "region",label = T)
DimPlot(object = IC_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "major_labl",label = T)
DimPlot(object = IC_snATAC,cols=myUmapcolors,reduction = "umap_harmony",  group.by = "cell_type",label = T)
dev.off();

# make tree for cluster
object<- IC_snATAC
Idents(IC_snATAC)<- IC_snATAC$peaks_nn_res.0.8
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
pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_EC_subcluster_harmony_tree.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()


# call peak 
DefaultAssay(IC_snATAC)<-"ATAC"
peak<-CallPeaks(
       IC_snATAC,
       group.by = "peaks_nn_res.0.8",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 2.7e+09,
       outdir="/data/R02/nieyg/project/EC/human_MI/02_EC/peaks",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(IC_snATAC),
     features = peak,
     cells = colnames(IC_snATAC)
     )     
#macs2_counts<-macs2_counts[-which(rownames(macs2_counts)=="GroupUN243-311766-311965"),]
# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "hg38"
annotations<- renameSeqlevels(annotations, paste0("chr",seqlevels(annotations)))

# create a new assay using the MACS2 peak set and add it to the Seurat object
IC_snATAC[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(IC_snATAC),
  annotation = Annotation(IC_snATAC)
)

DefaultAssay(IC_snATAC)<-"peaks"
gene.activities <- GeneActivity(IC_snATAC,features=rownames(IC_snRNA))
# add gene activities as a new assay
IC_snATAC[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

DefaultAssay(IC_snATAC)<-"RNA_imputation"
# The dot plot for annotation 
features <- c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","CD36","MGLL",# capillary Endo
	"FAM155A","TSHZ2","SLCO2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2","POSTN",# venous Endo
	"PCSK5","ARL15","LINC00639","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2","SLC45A4",# Arterial Endo
	"LINC02388","TMEM132C","PKHD1L1","NRG3","BMPER","TMEM108","NRG1","CGNL1","CDH11","PLA2G5",# Endocardial Endo
	"SEMA3A","RELN","MMRN1","PIEZO2","CCL21","DOCK5","TFPI","FLT4","AKAP12"# Lymphatic Endo
	)
features_label<-c(rep("capillary",10),rep("venous",9),rep("Arterial",10),rep("Endocardial",10),rep("Lymphatic",9))
features_label<- features_label[which(features%in% rownames(IC_snATAC))]
features<- features[which(features%in% rownames(IC_snATAC))]

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
 
Idents(IC_snATAC)<-IC_snATAC$peaks_nn_res.0.8

pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_EC_subcluster_annotation_dotplot_RNA_imputation.pdf",width=12,height=8)
p1<- DotPlot(IC_snATAC, features = features,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
p1&scale_x_discrete(labels=features_label)
p2<- DotPlot(IC_snATAC, features = EC_atlas,dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p2
p2&scale_x_discrete(labels=EC_atlas_label)
dev.off()

Idents(IC_snATAC)<-IC_snATAC$peaks_nn_res.0.8
#####further annotation########
IC_snATAC <- RenameIdents(
  object = IC_snATAC,
  '0' = 'Capillary_Endo',
  '1' = 'Venous_Endo',
  '2' = 'Capillary_Endo',
  '3' = 'Arterial_Endo',
  '4' = 'Capillary_Endo',
  '5' = 'notsure',
  '6' = 'notsure',
  '7' = 'Capillary_Endo',
  '8' = 'Endocardial_Endo',
  '9' = 'notsure',
  '10' = 'Lymphatic_Endo',
  '11' = 'C_A',
  '12' = 'notsure',
  '13' = 'C_V'
  )
IC_snATAC@meta.data$subcluster<-Idents(IC_snATAC)
table(IC_snATAC$subcluster,IC_snATAC$major_labl)

# not split C_A_V
Idents(IC_snATAC)<-IC_snATAC$peaks_nn_res.0.8
#####further annotation########
IC_snATAC <- RenameIdents(
  object = IC_snATAC,
  '0' = 'CAV',
  '1' = 'CAV',
  '2' = 'CAV',
  '3' = 'CAV',
  '4' = 'CAV',
  '5' = 'CAV',
  '6' = 'CAV',
  '7' = 'CAV',
  '8' = 'Endocardial_Endo',
  '9' = 'CAV',
  '10' = 'Lymphatic_Endo',
  '11' = 'CAV',
  '12' = 'notsure',
  '13' = 'CAV'
  )
IC_snATAC@meta.data$subcluster_CAV<-Idents(IC_snATAC)
table(IC_snATAC$subcluster_CAV,IC_snATAC$major_labl)

pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_EC_subcluster_annotation_umap.pdf",width=6,height=5)
DimPlot(IC_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "subcluster")
DimPlot(IC_snATAC, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "subcluster_CAV")
dev.off()

DefaultAssay(IC_snATAC) <- "peaks"
saveRDS(IC_snATAC,"/md01/nieyg/project/EC/human_MI/02_EC/snATAC_All_EC_harmony.rds")

IC_snATAC<- readRDS("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_All_EC_harmony.rds")


# proportion for subcluster

pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_All_EC_subtype_proportion.pdf",width=10,height=5)
df <- as.data.frame(IC_snATAC@meta.data)
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


pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_All_EC_not_CAV_proportion.pdf",width=10,height=5)
df <- as.data.frame(IC_snRNA@meta.data)
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
DefaultAssay(IC_snATAC) <- "ACTIVITY"
IC_snATAC <- IC_snATAC %>% 
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA(verbose = FALSE)

IC_snRNA <- FindVariableFeatures(IC_snRNA)


# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = IC_snRNA, 
    query = IC_snATAC,
    features = VariableFeatures(object = IC_snRNA),
    reference.assay = "RNA",
    query.assay = "ACTIVITY", 
    reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, 
	refdata = IC_snRNA$subcluster,
    weight.reduction = IC_snATAC[["scopen"]], 
    dims = 1:30)

IC_snATAC <- AddMetaData(IC_snATAC, metadata = celltype.predictions)
IC_snATAC$subcluster_annotation_correct <- IC_snATAC$predicted.id == IC_snATAC$subcluster
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = IC_snRNA$subcluster_CAV,
    weight.reduction = IC_snATAC[["scopen"]], dims = 1:30)
IC_snATAC$predicted_CAV<- IC_snATAC$predicted.id

IC_snATAC <- AddMetaData(IC_snATAC, metadata = celltype.predictions)
IC_snATAC$subcluster_CAV_annotation_correct <- IC_snATAC$predicted.id == IC_snATAC$subcluster_CAV


IC_snATAC$predicted_CAV<- factor(IC_snATAC$predicted_CAV,levels=c("Venous_Endo","C_V","Capillary_Endo","C_A","Arterial_Endo","Endocardial_Endo","Lymphatic_Endo","notsure"))
IC_snATAC$subcluster<- factor(IC_snATAC$subcluster,levels=c("Venous_Endo","C_V","Capillary_Endo","C_A","Arterial_Endo","Endocardial_Endo","Lymphatic_Endo","notsure"))

pdf("/md01/nieyg/project/EC/human_MI/02_EC/All_EC_ATAC_predicted_Ground-truth_Umap.pdf",width=10,height=5)
p1 <- DimPlot(IC_snATAC,reduction='umap_harmony' ,cols=myUmapcolors, group.by = "predicted_CAV", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(IC_snATAC,reduction='umap_harmony' ,cols=myUmapcolors,group.by = "subcluster", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2
p1 <- DimPlot(IC_snATAC,reduction='umap_harmony' ,cols=myUmapcolors, group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(IC_snATAC,reduction='umap_harmony' ,cols=myUmapcolors,group.by = "subcluster_CAV", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2
dev.off()

# not split CAV 
predictions <- table(IC_snATAC$subcluster_CAV, IC_snATAC$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

order_Var1<- c(levels(Idents(IC_snATAC)))
order_Var2<- c(levels(Idents(IC_snATAC)))
predictions$Var1<- factor(predictions$Var1,levels=order_Var1)
predictions$Var2<- factor(predictions$Var2,levels=order_Var2)
library(cowplot)
pdf("/md01/nieyg/project/EC/human_MI/02_EC/All_EC_ATAC_predicted_rate_notCAV.pdf",width=10,height=5)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(IC_snATAC$subcluster_CAV == IC_snATAC$predicted.id))
incorrect <- length(which(IC_snATAC$subcluster_CAV != IC_snATAC$predicted.id))
data <- FetchData(IC_snATAC, vars = c("prediction.score.max", "subcluster_CAV_annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = subcluster_CAV_annotation_correct, colour = subcluster_CAV_annotation_correct)) +
    geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2
dev.off()

saveRDS(IC_snATAC,"/md01/nieyg/project/EC/human_MI/02_EC/snATAC_All_EC_harmony_predictions.rds")


# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(IC_snRNA)
refdata <- GetAssayData(IC_snRNA, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, 
    refdata = refdata, 
    weight.reduction = "cca")
IC_snATAC[["RNA"]] <- imputation


coembed <- merge(x = IC_snRNA, y = IC_snATAC)
coembed<- RenameCells(coembed,
    old.names=colnames(coembed),
    new.names =paste(coembed$sample,gsub("_[1,2]*","",colnames(coembed)),sep="#"))

# Finally, we run PCA and UMAP on this combined object, to visualize the co-embedding of both
# datasets
coembed <- coembed %>%
    ScaleData(features = genes.use, do.scale = FALSE) %>%
    RunPCA(features = genes.use,, verbose = FALSE) %>%
    RunUMAP(dims = 1:30, verbose = FALSE)

pdf("/md01/nieyg/project/EC/human_MI/02_EC/coembed_umap.pdf",width=20,height=5)
DimPlot(coembed,reduction='umap',cols=myUmapcolors,raster=FALSE, group.by = c("patient_region", "subcluster",'tech'))
DimPlot(coembed,reduction='umap',cols=myUmapcolors,raster=FALSE, group.by = c("patient_region", "subcluster_CAV",'tech'))
dev.off()
saveRDS(coembed,"./02_EC/All_IC_snRNA_snATAC_merged.rds")

IC_snRNA<- readRDS("./02_EC/snRNA_All_EC_harmony.rds")
IC_snATAC<- readRDS("./02_EC/snATAC_All_EC_harmony.rds")
coembed<- readRDS("./02_EC/All_IC_snRNA_snATAC_merged.rds")

# The features between notsure and CAV in snATAC 
    IC_snATAC <- NucleosomeSignal(IC_snATAC)
    IC_snATAC <- TSSEnrichment(IC_snATAC,fast=FALSE)

pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_features_umap.pdf",width=24,height=20)
FeaturePlot(
  object = IC_snATAC,features = c("nCount_RNA","nFeature_RNA","TSS.enrichment","nucleosome_signal",
  	"nCount_ATAC","nFeature_ATAC","nCount_peaks","nFeature_peaks",
  	"nCount_ACTIVITY","nFeature_ACTIVITY","nCounts_RNA","nFeaturess_RNA"),
  order = FALSE,reduction='umap_harmony' ,
  cols = c("grey", "darkred"),
  ncol = 4
) & NoLegend()
dev.off()
Idents(IC_snATAC)<- IC_snATAC$subcluster_CAV
pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_EC_QC_vlnplot.pdf",width=10,height=16)
VlnPlot(object = IC_snATAC,features = c("nCount_RNA","nFeature_RNA","TSS.enrichment","nucleosome_signal"), ncol = 1, pt.size = 0)
VlnPlot(object = IC_snATAC,features = c("nCount_ATAC","nFeature_ATAC","nCount_peaks","nFeature_peaks"), ncol = 1, pt.size = 0)
VlnPlot(object = IC_snATAC,features = c("nCount_ACTIVITY","nFeature_ACTIVITY","nCounts_RNA","nFeaturess_RNA"), ncol = 1, pt.size = 0)
dev.off()


# The DEP in cluster 12 and other CAV

Idents(IC_snATAC) <- IC_snATAC$subcluster_CAV

# change back to working with peaks instead of gene activities
DefaultAssay(IC_snATAC) <- 'peaks'

da_peaks <- FindMarkers(
  object = IC_snATAC,
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

pdf("/md01/nieyg/project/EC/human_MI/02_EC/snATAC_EC_CAV_vs_notsure_volcano.pdf",width=10,height=10)
p1 = VolcanoPlot(dif, padj=0.05, title="CAV vs notsure", label.max = 50)
p1
dev.off()






fc <- FoldChange(IC_snATAC, ident.1 = "CAV", ident.2 = "notsure")
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
IC_snATAC<- readRDS("./02_EC/snATAC_All_EC_harmony.rds")
IC_snATAC <- TSSEnrichment(object = IC_snATAC, fast = FALSE)
Idents(IC_snATAC) <- IC_snATAC$subcluster_CAV
IC_snATAC$TSS_group <- ifelse(IC_snATAC$subcluster_CAV == "notsure", 'notsure', 'other')
pdf("./02_EC/snATAC_notsure_TSS_distribution.pdf",width=12,height=6)
TSSPlot(IC_snATAC, group.by = 'TSS_group') + NoLegend()
TSSPlot(IC_snATAC, group.by = 'subcluster_CAV') + NoLegend()
dev.off()

IC_snATAC <- NucleosomeSignal(object = IC_snATAC)
IC_snATAC$nucleosome_group <- ifelse(IC_snATAC$nucleosome_signal == "notsure", 'notsure', 'other')
pdf("./02_EC/snATAC_notsure_NucleosomeSignal_distribution.pdf",width=12,height=6)
FragmentHistogram(object = IC_snATAC, group.by = 'nucleosome_group')
FragmentHistogram(object = IC_snATAC, group.by = 'subcluster_CAV')
dev.off()





