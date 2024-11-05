library(monocle3)
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
library(monocle)
library(SeuratWrappers)
CAV_snRNA<- readRDS("./02_EC/CAV/EC_CAV_snRNA_add_gene_symbol.rds")
CAV_snATAC<- readRDS("./02_EC/CAV/EC_CAV_snATAC_add_peaks_fragment.rds")

# Building trajectories with Monocle 2 for snRNA of CAV 
DefaultAssay(CAV_snRNA)<-"RNA"
CAV_snRNA.list <- SplitObject(CAV_snRNA,split.by = "sample")
for (i in 1:length(CAV_snRNA.list)) {
  DefaultAssay(CAV_snRNA.list[[i]]) <- "RNA"
  CAV_snRNA.list[[i]] <- SCTransform(CAV_snRNA.list[[i]], verbose = FALSE)
}
for (i in 1:length(CAV_snRNA.list)) {
  DefaultAssay(CAV_snRNA.list[[i]]) <- "SCT"
}

CAV_snRNA.anchors <- FindIntegrationAnchors(object.list = CAV_snRNA.list,anchor.features = 2000,dims = 1:30)
CAV_snRNA.integrated <- IntegrateData(anchorset = CAV_snRNA.anchors, dims = 1:30,normalization.method = "SCT",new.assay.name = "integrated_CAV")
head(CAV_snRNA.integrated@meta.data)
DefaultAssay(CAV_snRNA.integrated) <- "integrated_CAV"

# scale and center features in the dataset
CAV_snRNA.integrated <- ScaleData(CAV_snRNA.integrated, features =rownames(CAV_snRNA.integrated),verbose = FALSE)

# only keep the Artery_1
# DEG among subtype
Idents(CAV_snRNA.integrated)<- CAV_snRNA.integrated$CAV_detail
CAV_snRNA_injury.integrated<- subset(CAV_snRNA.integrated,idents=levels(CAV_snRNA.integrated)[-9])
DefaultAssay(CAV_snRNA_injury.integrated)<-"RNA"
Idents(CAV_snRNA_injury.integrated)<- CAV_snRNA_injury.integrated$CAV
CAV_snRNA_injury.integrated <- ScaleData(CAV_snRNA_injury.integrated,features=rownames(CAV_snRNA_injury.integrated))
markers <- FindAllMarkers(CAV_snRNA_injury.integrated, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
table(markers$cluster)
CAV_subtype_DEG<- markers$gene;
# DEG among cluster
Idents(CAV_snRNA_injury.integrated)<- CAV_snRNA_injury.integrated$CAV_detail
markers <- FindAllMarkers(CAV_snRNA_injury.integrated, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
table(markers$cluster)
CAV_cluster_DEG<- markers$gene;

# Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(CAV_snRNA_injury.integrated@assays$integrated_CAV@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = CAV_snRNA_injury.integrated@meta.data)
fData <- data.frame(gene_short_name = row.names(CAV_snRNA_injury.integrated), row.names = row.names(CAV_snRNA_injury.integrated))
fd <- new('AnnotatedDataFrame', data = fData)
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FMPCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFMP3", "#BEBADA", "#MP8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

monocle_cds <- monocle::newCellDataSet(data,
	phenoData = pd,
	featureData = fd,
	lowerDetectionLimit = 0.5,
	expressionFamily = uninormal()) 
monocle_cds<- estimateSizeFactors(monocle_cds)
#monocle_cds<- estimateDispersions(monocle_cds)
monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
expressed_genes <- row.names(subset(featureData(monocle_cds)@data,
                                     num_cells_expressed >= 10))
# 1:seurat DEG
## 1) subtype DEG 
monocle_cds_tmp<- setOrderingFilter(monocle_cds,CAV_subtype_DEG)
monocle_cds_tmp <- reduceDimension(monocle_cds_tmp,max_components=2,norm_method  ="none" ,method="DDRTree") 
monocle_cds_tmp <- orderCells(monocle_cds_tmp)
pdf("./02_EC/CAV/05_trajectory/02_monocle2/02_injury_CAV/seurat_subtypeDEG/injury_seurat_DEG_subtype_pseudotime.pdf",width=12,height=12)
plot_cell_trajectory(monocle_cds_tmp, color_by = "seurat_clusters",size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",,size=1,show_backbone=T)+scale_color_manual(values=myUmapcolors)
plot_complex_cell_trajectory(monocle_cds_tmp,x=1,y=1,color_by="CAV")+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "State",,size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "Pseudotime",,size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",size=0.01) +facet_wrap(~major_labl, nrow =4)+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV_detail",size=0.01) +facet_wrap(~major_labl, nrow =4)+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV_detail",size=0.5,show_backbone=T)+scale_color_manual(values=myUmapcolors)
dev.off()
df<- monocle_cds_tmp@phenoData@data
df<- as.data.frame(df)
pdf("./02_EC/CAV/05_trajectory/02_monocle2/02_injury_CAV/seurat_subtypeDEG/injury_seurat_DEG_CAV_subtype_cell_density_by_pseudotime.pdf",width=10,height=5)
ggplot(df,aes(Pseudotime,colour=CAV,fill=CAV))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()+scale_color_manual(values=myUmapcolors)+scale_fill_manual(values=myUmapcolors)
dev.off()
saveRDS(monocle_cds_tmp,"./02_EC/CAV/05_trajectory/02_monocle2/02_injury_CAV/seurat_subtypeDEG/injury_integrated_seurat_DEG_CAV_subtype_pseudotime_monocle.rds")

## 2) cluster DEG 
monocle_cds_tmp<- setOrderingFilter(monocle_cds,CAV_cluster_DEG)
monocle_cds_tmp <- reduceDimension(monocle_cds_tmp,max_components=2,norm_method  ="none" ,method="DDRTree") 
monocle_cds_tmp <- orderCells(monocle_cds_tmp)
pdf("./02_EC/CAV/05_trajectory/02_monocle2/02_injury_CAV/seurat_clusterDEG/injury_seurat_DEG_cluster_pseudotime.pdf",width=12,height=12)
plot_cell_trajectory(monocle_cds_tmp, color_by = "seurat_clusters",size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",,size=1,show_backbone=T)+scale_color_manual(values=myUmapcolors)
plot_complex_cell_trajectory(monocle_cds_tmp,x=1,y=1,color_by="CAV")+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "State",,size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "Pseudotime",,size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",size=0.01) +facet_wrap(~major_labl, nrow =4)+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV_detail",size=0.01) +facet_wrap(~major_labl, nrow =4)+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV_detail",size=0.5,show_backbone=T)+scale_color_manual(values=myUmapcolors)
dev.off()
df<- monocle_cds_tmp@phenoData@data
df<- as.data.frame(df)
pdf("./02_EC/CAV/05_trajectory/02_monocle2/02_injury_CAV/seurat_clusterDEG/injury_seurat_DEG_CAV_cluster_cell_density_by_pseudotime.pdf",width=10,height=5)
ggplot(df,aes(Pseudotime,colour=CAV,fill=CAV))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()+scale_color_manual(values=myUmapcolors)+scale_fill_manual(values=myUmapcolors)
dev.off()
saveRDS(monocle_cds_tmp,"./02_EC/CAV/05_trajectory/02_monocle2/02_injury_CAV/seurat_clusterDEG/injury_integrated_seurat_DEG_CAV_cluster_pseudotime_monocle.rds")

# monocle_cds_tmp<- readRDS("./02_EC/CAV/monocle2/integrated_seurat_DEG_CAV_subtype_pseudotime_monocle.rds")






















































# Select the DEG for trajectory 
library(monocle3)
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
library(monocle)
library(SeuratWrappers)
monocle_cds_tmp<- readRDS("./02_EC/CAV/monocle2/integrated_seurat_DEG_CAV_subtype_pseudotime_monocle.rds")
df<- monocle_cds_tmp@phenoData@data
df<- as.data.frame(df)

# pdf("./02_EC/CAV/monocle2/check_seurat_DEG_CAV_subtype_pseudotime.pdf",width=12,height=12)
# plot_cell_trajectory(monocle_cds_tmp, color_by = "seurat_clusters",size=1)
# plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",,size=1,show_backbone=T)+scale_color_manual(values=myUmapcolors)
# plot_complex_cell_trajectory(monocle_cds_tmp,x=1,y=1,color_by="CAV")+scale_color_manual(values=myUmapcolors)
# plot_cell_trajectory(monocle_cds_tmp, color_by = "State",,size=1)
# plot_cell_trajectory(monocle_cds_tmp, color_by = "Pseudotime",,size=1)
# plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",size=0.01) +facet_wrap(~major_labl, nrow =4)+scale_color_manual(values=myUmapcolors)
# plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV_detail",size=0.01) +facet_wrap(~major_labl, nrow =4)+scale_color_manual(values=myUmapcolors)
# plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV_detail",size=0.5,show_backbone=T)+scale_color_manual(values=myUmapcolors)
# dev.off()


pdf("./02_EC/CAV/monocle2/check_seurat_DEG_CAV_subtype_pseudotime.pdf",width=12,height=12)
plot_cell_trajectory(monocle_cds_tmp, color_by = "State",,size=1)
plot_complex_cell_trajectory(monocle_cds_tmp,x=1,y=1,color_by="CAV",size=0.01)+scale_color_manual(values=myUmapcolors)+facet_wrap(~major_labl, nrow =4)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV_detail",size=0.01) +facet_wrap(~major_labl, nrow =4)+scale_color_manual(values=myUmapcolors)
dev.off()

# for C_A end 

df<- monocle_cds_tmp@phenoData@data
df<- as.data.frame(df)
valid_cells <- row.names(subset(df,
            State %in%c("1","3","4","5")))
C_A_end_cds <- monocle_cds_tmp[,valid_cells]
# split by condition 
df<- C_A_end_cds@phenoData@data
df<- as.data.frame(df)
# stacked density
CTRL<- ggplot(df[df$major_labl=="CTRL",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("CTRL")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)
BZ<- ggplot(df[df$major_labl=="BZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("BZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)
FZ<- ggplot(df[df$major_labl=="FZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("FZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+xlim(0,11)
IZ_df<- df[df$major_labl=="IZ",]
IZ<- ggplot(IZ_df,aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("IZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)
RZ<- ggplot(df[df$major_labl=="RZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("RZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)
global<- ggplot(df,aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+theme_classic()+ggtitle("global")+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)
pdf("./02_EC/CAV/monocle2/C_A_end_stacked_density_by_pseudotime.pdf",width=20,height=2)
global|CTRL|IZ|BZ|RZ|FZ
dev.off()

CTRL<- ggplot(df[df$major_labl=="CTRL",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("CTRL")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)+ylim(0,0.6)
BZ<- ggplot(df[df$major_labl=="BZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("BZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)+ylim(0,0.6)
FZ<- ggplot(df[df$major_labl=="FZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("FZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+xlim(0,11)+ylim(0,0.6)
IZ<- ggplot(df[df$major_labl=="IZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("IZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)+ylim(0,0.6)
RZ<- ggplot(df[df$major_labl=="RZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("RZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)+ylim(0,0.6)
global<- ggplot(df,aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+theme_classic()+ggtitle("global")+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,11)+ylim(0,0.6)
pdf("./02_EC/CAV/monocle2/C_A_end_density_by_pseudotime.pdf",width=20,height=2)
global|CTRL|IZ|BZ|RZ|FZ
dev.off()


pdf("./02_EC/CAV/monocle2/C_A_end_ECDF_test.pdf",width=20,height=2)
global<- ggplot(df,aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("global")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,11)
CTRL<- ggplot(df[df$major_labl=="CTRL",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("CTRL")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,11)
IZ<- ggplot(df[df$major_labl=="IZ",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("IZ")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,11)
BZ<- ggplot(df[df$major_labl=="BZ",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("BZ")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,11)
RZ<- ggplot(df[df$major_labl=="RZ",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("RZ")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,11)
FZ<- ggplot(df[df$major_labl=="FZ",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("FZ")+theme_classic()+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,11)
global|CTRL|IZ|BZ|RZ|FZ
dev.off()


df$Pseudotime <- as.integer(df$Pseudotime)
pdf("./02_EC/CAV/monocle2/C_A_end_barplot.pdf",width=20,height=2)
global<- ggplot(df, aes(x = as.character(Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("global") +theme_bw() +theme( legend.position =" none ")
CTRL<- ggplot(df[df$major_labl=="CTRL",], aes(x = as.character(Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("CTRL") +theme_bw() +theme( legend.position =" none ")
IZ<- ggplot(df[df$major_labl=="IZ",], aes(x = as.character(Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("IZ") +theme_bw() +theme( legend.position =" none ")
BZ<- ggplot(df[df$major_labl=="BZ",], aes(x = as.character(Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("BZ") +theme_bw() +theme( legend.position =" none ")
RZ<- ggplot(df[df$major_labl=="RZ",], aes(x = as.character(Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("RZ") +theme_bw() +theme( legend.position =" none ")
FZ<- ggplot(df[df$major_labl=="FZ",], aes(x = as.character(Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("FZ") +theme_bw()
global|CTRL|IZ|BZ|RZ|FZ
dev.off()




df<- monocle_cds_tmp@phenoData@data
df<- as.data.frame(df)

# for Artery end 
valid_cells <- row.names(subset(df,
            State %in%c("1","3","7")))
Artery_end_cds <- monocle_cds_tmp[,valid_cells]
# split by condition 
df<- Artery_end_cds@phenoData@data
df<- as.data.frame(df)
CTRL<- ggplot(df[df$major_labl=="CTRL",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("CTRL")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
BZ<- ggplot(df[df$major_labl=="BZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("BZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
FZ<- ggplot(df[df$major_labl=="FZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("FZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+xlim(0,13)
IZ_df<- df[df$major_labl=="IZ",]
IZ<- ggplot(IZ_df,aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("IZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
RZ<- ggplot(df[df$major_labl=="RZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+ggtitle("RZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
global<- ggplot(df,aes(Pseudotime,group=CAV,fill=CAV))+geom_density(adjust=1.5,position="fill")+theme_classic()+ggtitle("global")+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
pdf("./02_EC/CAV/monocle2/Artery_end_stacked_density_by_pseudotime.pdf",width=20,height=2)
global|CTRL|IZ|BZ|RZ|FZ
dev.off()

CTRL<- ggplot(df[df$major_labl=="CTRL",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("CTRL")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
BZ<- ggplot(df[df$major_labl=="BZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("BZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
FZ<- ggplot(df[df$major_labl=="FZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("FZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+xlim(0,13)
IZ<- ggplot(df[df$major_labl=="IZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("IZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
RZ<- ggplot(df[df$major_labl=="RZ",],aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+ggtitle("RZ")+theme_classic()+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
global<- ggplot(df,aes(Pseudotime,group=CAV,fill=CAV))+geom_density(bw=0.5,alpha = 0.5)+theme_classic()+ggtitle("global")+scale_fill_manual(values=myUmapcolors)+theme( legend.position =" none ")+xlim(0,13)
pdf("./02_EC/CAV/monocle2/Artery_end_density_by_pseudotime.pdf",width=20,height=2)
global|CTRL|IZ|BZ|RZ|FZ
dev.off()
pdf("./02_EC/CAV/monocle2/Artery_end_ECDF_test.pdf",width=20,height=2)
global<- ggplot(df,aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("global")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,13)
CTRL<- ggplot(df[df$major_labl=="CTRL",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("CTRL")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,13)
IZ<- ggplot(df[df$major_labl=="IZ",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("IZ")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,13)
BZ<- ggplot(df[df$major_labl=="BZ",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("BZ")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,13)
RZ<- ggplot(df[df$major_labl=="RZ",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("RZ")+theme_classic()+theme( legend.position =" none ")+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,13)
FZ<- ggplot(df[df$major_labl=="FZ",],aes(Pseudotime,col=CAV))+stat_ecdf(geom="smooth", se=F, size=1.2)+ggtitle("FZ")+theme_classic()+  labs(x="Pseudotime",y="ECDF",col="") +scale_fill_manual(values=myUmapcolors)+scale_color_manual(values=myUmapcolors)+xlim(0,13)
global|CTRL|IZ|BZ|RZ|FZ
dev.off()

df<- Artery_end_cds@phenoData@data
df<- as.data.frame(df)
df$Pseudotime <- as.integer(df$Pseudotime)
df$Pseudotime <- as.factor(df$Pseudotime)
df$Pseudotime <- factor(df$Pseudotime,levels=as.character(0:12))
pdf("./02_EC/CAV/monocle2/Artery_end_barplot.pdf",width=20,height=2)
global<- ggplot(df, aes(x = Pseudotime, fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("global") +theme_bw() +theme( legend.position =" none ")
CTRL<- ggplot(df[df$major_labl=="CTRL",], aes(x = (Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("CTRL") +theme_bw() +theme( legend.position =" none ")
IZ<- ggplot(df[df$major_labl=="IZ",], aes(x = (Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("IZ") +theme_bw() +theme( legend.position =" none ")
BZ<- ggplot(df[df$major_labl=="BZ",], aes(x = (Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("BZ") +theme_bw() +theme( legend.position =" none ")
RZ<- ggplot(df[df$major_labl=="RZ",], aes(x = (Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("RZ") +theme_bw() +theme( legend.position =" none ")
FZ<- ggplot(df[df$major_labl=="FZ",], aes(x = (Pseudotime), fill = as.factor(CAV))) +geom_bar(stat = "count",width = 0.8, position = position_dodge(width = 0.8)) +scale_fill_manual(values = c(myUmapcolors)) +ggtitle("FZ") +theme_bw()
global|CTRL|IZ|BZ|RZ|FZ
dev.off()


# find the Pseudotime related DEG 
valid_cells <- row.names(subset(df,
            State %in%c("1","3","4","5","7")))

BEAM_cds <- monocle_cds_tmp[,valid_cells]


BEAM_res <- BEAM(BEAM_cds, branch_point = 2, cores = 1,progenitor_method = 'duplicate')
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
  pdf("./02_EC/CAV/monocle2/branched_heatmap.pdf",width=5,height=8)
  plot_genes_branched_heatmap(monocle_cds_tmp[row.names(subset(BEAM_res,
                                            pval < 1e-1)),],
                                            branch_point = 1,
                                            num_clusters = 4,
                                            cores = 1,
                                            use_gene_short_name = T,
                                            show_rownames = T)
  dev.off()

show_genes <- row.names(subset(monocle_cds_tmp@featureData@data,
          gene_short_name %in% c("CDK6", "YY1AP1", "MYO1E","COX3")))



pdf("./02_EC/CAV/monocle2/branched_DEG_pseudotime.pdf",width=8,height=8)
plot_genes_branched_pseudotime(monocle_cds_tmp[show_genes,],
                       branch_point = 2,
                       color_by = "State",
                       ncol = 2)
plot_genes_branched_pseudotime(monocle_cds_tmp[show_genes,],
                       branch_point = 2,
                       color_by = "CAV",
                       ncol = 2)
dev.off()




# Timediff_test_res <- differentialGeneTest(monocle_cds[ordering_genes,],
#               fullModelFormulaStr = "~sm.ns(Pseudotime)")
# write.csv(Timediff_test_res,"./06_monocle/Pseudotime_related_DEG.csv")
# 
# sig_gene_names <- row.names(subset(Timediff_test_res, qval < 0.1))
# pdf("./06_monocle/Pseudotime_related_DEG_heatmap.pdf",width=5,height=10)
# plot_pseudotime_heatmap(monocle_cds[sig_gene_names,],
#                 num_clusters = 3,
#                 cores = 1,
#                 show_rownames = T)
# dev.off()





















