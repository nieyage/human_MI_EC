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
Neuronal_snRNA <- subset(snRNA,idents="Neuronal")
integrated_data <- SplitObject(Neuronal_snRNA, split.by = "sample")

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
         dims = 1:20,min.dist=1,metric = 'euclidean',
          reduction.name = "umap_harmony") %>% 
  identity()

library(clustree)
pdf("./08_Neuron/resolutions_select_Neuron_subcluster_reintergrated.pdf",width=15,height=15)
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

integrated_data$seurat_clusters<- integrated_data$RNA_snn_res.0.8
Idents(integrated_data)<- integrated_data$seurat_clusters

# unsupervised clustering umap plot 
pdf("/md01/nieyg/project/EC/human_MI/08_Neuron/snRNA_Neuron_subcluster_umap.pdf",width=5,height=5)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "seurat_clusters",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "patient_group",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "patient_region",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "region",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "major_labl",label = T)
DimPlot(object = integrated_data,cols=myUmapcolors, group.by = "cell_type",label = T)
dev.off();

# make tree for cluster
object<- integrated_data
embeddings <- Embeddings(object = object, reduction = "pca")[,1:20]
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
pdf("/md01/nieyg/project/EC/human_MI/08_Neuron/snRNA_Neuron_subcluster_harmony_tree.pdf",width=6,height=6)
ggtree(data.tree,layout = "circular") + geom_tiplab()+ geom_treescale()
dev.off()

Neuronal_snRNA<- integrated_data

Idents(Neuronal_snRNA)<-Neuronal_snRNA$seurat_clusters
#####further annotation########
Neuronal_snRNA <- RenameIdents(
  object = Neuronal_snRNA,
  '0' = 'Neuron1',
  '1' = 'Neuron2',
  '2' = 'Neuron3',
  '3' = 'Neuron2',
  '4' = 'Neuron1',
  '5' = 'Neuron1',
  '6' = 'Neuron4',
  '7' = 'Neuron4',
  '8' = 'Neuron2',
  '9' = 'Neuron3',
  '10' = 'Neuron4',
  '11' = 'Neuron3'
  )
Neuronal_snRNA@meta.data$subcluster<-Idents(Neuronal_snRNA)
table(Neuronal_snRNA$subcluster,Neuronal_snRNA$major_labl)

pdf("/md01/nieyg/project/EC/human_MI/08_Neuron/snRNA_Neuron_subcluster_annotation_umap.pdf",width=6,height=5)
DimPlot(Neuronal_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "subcluster")
DimPlot(Neuronal_snRNA, label = T, repel = TRUE, cols=myUmapcolors, reduction = "umap_harmony",group.by = "major_labl")
dev.off()

# proportion for subcluster

pdf("/md01/nieyg/project/EC/human_MI/08_Neuron/snRNA_Neuron_subtype_proportion.pdf",width=6,height=5)
df <- as.data.frame(Neuronal_snRNA@meta.data)
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

# Find subtype markers and GO enrichment 
combined<- SCTransform(Neuronal_snRNA)
DefaultAssay(combined)<- "SCT"
combined<- ScaleData(combined)

markers <- FindAllMarkers(combined, only.pos = TRUE,logfc.threshold = 0.25)
#markers<- markers[markers$p_val_adj<0.05&markers$avg_log2FC>1,]

write.csv(markers,"./08_Neuron/Neuron_subtype_DEG_markers-log2FC1adjP005.csv")

top10 <- markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC);

blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")
pdf("./08_Neuron/Neuron_subtype_markers-log2FC1adjP005_top_markers-DEG_heatmap.pdf",width=10,height=10)
DoHeatmap(object = combined,features=markers$gene,label=T,size = 2,group.by = "subcluster") 
DoHeatmap(object = combined,features=top10$gene,label=T, group.colors =c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494","#B3B3B3" ),
  disp.min = -2,disp.max = 2,size = 2,group.by = "subcluster") 
DoHeatmap(object = combined,features=markers$gene,label=T, group.colors =c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494","#B3B3B3" ),
  disp.min = -2,disp.max = 2,size = 2,group.by = "subcluster") 
dev.off();


# GO and KEGG for subcluster
Neuron1<-markers[markers$cluster=="Neuron1",7]
Neuron2<-markers[markers$cluster=="Neuron2",7]
Neuron3<-markers[markers$cluster=="Neuron3",7]
Neuron4<-markers[markers$cluster=="Neuron4",7]

library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
Neuron_list<- list(Neuron1,Neuron2,Neuron3,Neuron4)
subtype<- c("Neuron1","Neuron2","Neuron3","Neuron4")
all_ego<- data.frame()
pdf("./08_Neuron/Neuron_subtype_DEG_GO_BP.pdf")
for(i in 1:length(subtype)){
    gene<- Neuron_list[[i]];
    gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.5,
                readable = TRUE)
    p<- barplot(ego, showCategory=20)
    print(p)
    write.csv(ego,paste0(subtype[i],"-GO-BP.csv",sep=""))
    ego<- as.data.frame(ego)
    ego$celltype<- subtype[i];
    all_ego<- rbind(all_ego,ego);

}
dev.off()

write.csv(all_ego,"./08_Neuron/Allsubtype_GO.csv")

#top_term <- read.csv("./08_Neuron/Allsubtype_GO-selected.csv") 
library(ggplot2)
#top_term$Description<- factor(top_term$Description,levels=rev(top_term$Description))

top_term <- all_ego %>% group_by(celltype) %>% top_n(n = -5, wt = p.adjust);
top_term<- top_term[!duplicated(top_term$Description),]
top_term$Description<- factor(top_term$Description,levels=rev(top_term$Description))
top_term$celltype<- factor(top_term$celltype,levels=levels(combined))
pdf("./08_Neuron/All_Neuron_subtype_GO.pdf",width=14,height=12)
p <- ggplot(top_term,aes(y=Count,x=Description,fill=pvalue)) + 
      geom_bar(stat="identity",position = "dodge") +
      facet_grid(celltype~.,scales = "free",space = "free") + 
      coord_flip() + 
      theme_bw() +
      scale_fill_gradient(low = '#FF0000', high = '#1202FF')+
      theme(plot.title = element_text(hjust = 0.5),
            strip.text.y = element_text(size = 14),
            legend.position="right",
            legend.title = element_text(size=18),
            legend.text = element_text(size=14),
            axis.text.x = element_text(size=14),
            axis.text.y = element_text(size=18),
            axis.title.x = element_text(size=14),
            axis.title.y = element_text(size=14))
p
dev.off()





pdf("/md01/nieyg/project/EC/human_MI/08_Neuron/snRNA_Neuron_subcluster_annotated_dotplot.pdf",width=10,height=3)
p1<- DotPlot(Neuronal_snRNA, features = unique(top10$gene),dot.scale = 3) +theme(axis.text.x=element_text(size=8,angle=90,hjust=1,vjust=0.5))
p1
dev.off()
DefaultAssay(Neuronal_snRNA) <- "RNA"
saveRDS(Neuronal_snRNA,"/md01/nieyg/project/EC/human_MI/08_Neuron/snRNA_Neuron_harmony.rds")
