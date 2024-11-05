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


# Extract data, phenotype data, and feature data from the SeuratObject
data <- as(as.matrix(CAV_snRNA.integrated@assays$integrated_CAV@data), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = CAV_snRNA.integrated@meta.data)
fData <- data.frame(gene_short_name = row.names(CAV_snRNA.integrated), row.names = row.names(CAV_snRNA.integrated))
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

# 1: diff_gene_monocle 
diff_test_res <- differentialGeneTest(monocle_cds[expressed_genes,],fullModelFormulaStr = "~CAV",cores=1)
diff_test_res_batch <- differentialGeneTest(monocle_cds[expressed_genes,], fullModelFormulaStr = " ~ CAV + sample",reducedModelFormulaStr = " ~ sample", relative_expr=TRUE, cores=4)

head(diff_test_res)
deg <- subset(diff_test_res, qval < 0.001)
deg <- deg[order(deg$qval,decreasing=F),]
write.csv(deg,"./02_EC/CAV/monocle2/integrated_train_monocle_DEG.csv")

deg<- read.csv("./02_EC/CAV/monocle2/integrated_train_monocle_DEG.csv")
head(deg)

#disp_table <- dispersionTable(monocle_cds)
#disp.genes <- subset(disp_table, mean_expression >= 0.1 &dispersion_empirical>= 1*dispersion_fit)$gene_id

ordering_genes<- rownames(deg)
#ordering_genes<- rownames(deg)[which(rownames(deg)%in% disp.genes)]

monocle_cds_tmp<- setOrderingFilter(monocle_cds,ordering_genes)
# pdf("./02_EC/CAV/monocle2/integrated_monocle_plot_ordering_monocle_DEG_genes.pdf")
# plot_ordering_genes(monocle_cds_tmp)
# #plot_pc_variance_explained(monocle_cds, return_all = F)
# dev.off()
# reduceDimension（cds ， max_components  =  2 ， reduction_method  =  c（“DDRTree” ，
  “ICA” ， “tSNE” ， “SimplePPT” ， “L1-graph” ， “SGL-tree” ）， norm_method  =  c（“log” ，
  “vstExprs” ， “none” ）， residualModelFormulaStr  =  NULL， pseudo_expr  =  1 ，
  relative_expr  =  TRUE， auto_param_selection  =  TRUE， verbose  =  FALSE，
  scaling  =  TRUE， ...）



monocle_cds_tmp <- reduceDimension(monocle_cds_tmp,max_components=2,norm_method  ="none" ,method="DDRTree",residualModelFormulaStr="~donor_id") 

monocle_cds_tmp <- orderCells(monocle_cds_tmp)
pdf("./02_EC/CAV/monocle2/integrated_DEG_monocle2_CAV_subtype_pseudotime.pdf",width=10,height=10)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",size=0.6,show_backbone=T)+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV_detail",size=0.6,show_backbone=T)+scale_color_manual(values=myUmapcolors)
plot_complex_cell_trajectory(monocle_cds_tmp,x=1,y=1,size=0.6,color_by="CAV")+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "State",size=0.6)
plot_cell_trajectory(monocle_cds_tmp, color_by = "Pseudotime",size=0.6)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",size=0.1) +facet_wrap(~major_labl, nrow =3)+scale_color_manual(values=myUmapcolors)
dev.off()

library(ggpubr)
df<- monocle_cds_tmp@phenoData@data
df<- as.data.frame(df)
pdf("./02_EC/CAV/monocle2/integrated_DEG_monocle2_cell_density_by_pseudotime.pdf",width=10,height=5)
ggplot(df,aes(Pseudotime,colour=CAV,fill=CAV))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()+scale_color_manual(values=myUmapcolors)+scale_fill_manual(values=myUmapcolors)
dev.off()

# 2: top variable features in seurat 
variable_feature_seurat<- VariableFeatures(CAV_snRNA)
monocle_cds_tmp<- setOrderingFilter(monocle_cds,variable_feature_seurat[which(variable_feature_seurat%in%expressed_genes)])

monocle_cds_tmp <- reduceDimension(monocle_cds_tmp,max_components=2,norm_method  ="none" ,method="DDRTree") 
monocle_cds_tmp <- orderCells(monocle_cds_tmp)

pdf("./02_EC/CAV/monocle2/integrated_seurat_top_vari_CAV_subtype_pseudotime.pdf",width=12,height=12)
plot_cell_trajectory(monocle_cds_tmp, color_by = "seurat_clusters",size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",,size=1,show_backbone=T)+scale_color_manual(values=myUmapcolors)
plot_complex_cell_trajectory(monocle_cds_tmp,x=1,y=1,color_by="CAV")+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "State",,size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "Pseudotime",,size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",,size=0.5) +facet_wrap(~major_labl, nrow =4)+scale_color_manual(values=myUmapcolors)
dev.off()
df<- monocle_cds_tmp@phenoData@data
df<- as.data.frame(df)
pdf("./02_EC/CAV/monocle2/seurat_top_vari_CAV_subtype_cell_density_by_pseudotime.pdf",width=10,height=5)
ggplot(df,aes(Pseudotime,colour=CAV,fill=CAV))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()+scale_color_manual(values=myUmapcolors)+scale_fill_manual(values=myUmapcolors)
dev.off()

saveRDS(monocle_cds_tmp,"./02_EC/CAV/monocle2/integrated_seurat_top_vari_CAV_subtype_pseudotime_monocle.rds")


# 3:seurat DEG
cluster_deg<- read.csv("./02_EC/CAV/DEG_DEP/FindAllMarkers_gene.csv")
#cluster_deg<-  cluster_deg[cluster_deg$avg_log2FC>0.3,]
#cluster_deg<-  cluster_deg[cluster_deg$p_val_adj<0.05,]$gene
monocle_cds_tmp<- setOrderingFilter(monocle_cds,cluster_deg)

monocle_cds_tmp <- reduceDimension(monocle_cds_tmp,max_components=2,norm_method  ="none" ,method="DDRTree") 
monocle_cds_tmp <- orderCells(monocle_cds_tmp)
pdf("./02_EC/CAV/monocle2/seurat_DEG_CAV_subtype_pseudotime.pdf",width=12,height=12)
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
pdf("./02_EC/CAV/monocle2/seurat_DEG_CAV_subtype_cell_density_by_pseudotime.pdf",width=10,height=5)
ggplot(df,aes(Pseudotime,colour=CAV,fill=CAV))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()+scale_color_manual(values=myUmapcolors)+scale_fill_manual(values=myUmapcolors)
dev.off()
saveRDS(monocle_cds_tmp,"./02_EC/CAV/monocle2/integrated_seurat_DEG_CAV_subtype_pseudotime_monocle.rds")

monocle_cds_tmp<- readRDS("./02_EC/CAV/monocle2/integrated_seurat_DEG_CAV_subtype_pseudotime_monocle.rds")

# 5: CAV marker gene
marker_gene<- c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","CD36","MGLL",# capillary Endo
	"FAM155A","TSHZ2","SLCO2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2","POSTN",# venous Endo
	"PCSK5","ARL15","LINC00639","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2","SLC45A4",
	"FBLN5","HEY1","MECOM",# Artery 
	"CXCL12","RBP7",# C_A
	"KDR","ENDOU",# capillary
	"VCAM1","FMO1","FMO2",# C_V
	"MGP","CFH","BGN","VWF"# large vein
	)# Arterial Endo

monocle_cds_tmp<- setOrderingFilter(monocle_cds,marker_gene)

monocle_cds_tmp <- reduceDimension(monocle_cds_tmp,max_components=2,norm_method  ="none" ,method="DDRTree") 
monocle_cds_tmp <- orderCells(monocle_cds_tmp)
pdf("./02_EC/CAV/monocle2/CAV_marker_CAV_subtype_pseudotime.pdf",width=6,height=6)
plot_cell_trajectory(monocle_cds_tmp, color_by = "seurat_clusters",size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",,size=1,show_backbone=T)+scale_color_manual(values=myUmapcolors)
plot_complex_cell_trajectory(monocle_cds_tmp,x=1,y=1,color_by="CAV")+scale_color_manual(values=myUmapcolors)
plot_cell_trajectory(monocle_cds_tmp, color_by = "State",,size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "Pseudotime",,size=1)
plot_cell_trajectory(monocle_cds_tmp, color_by = "CAV",,size=0.5) +facet_wrap(~major_labl, nrow =4)+scale_color_manual(values=myUmapcolors)
dev.off()
df<- monocle_cds_tmp@phenoData@data
df<- as.data.frame(df)
pdf("./02_EC/CAV/monocle2/CAV_marker_CAV_subtype_cell_density_by_pseudotime.pdf",width=10,height=5)
ggplot(df,aes(Pseudotime,colour=CAV,fill=CAV))+geom_density(bw=0.5,size=1,alpha=0.5)+theme_classic()+scale_color_manual(values=myUmapcolors)+scale_fill_manual(values=myUmapcolors)
dev.off()













# find the Pseudotime related DEG 
Timediff_test_res <- differentialGeneTest(monocle_cds[ordering_genes,],
              fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(Timediff_test_res,"./06_monocle/Pseudotime_related_DEG.csv")

sig_gene_names <- row.names(subset(Timediff_test_res, qval < 0.1))
pdf("./06_monocle/Pseudotime_related_DEG_heatmap.pdf",width=5,height=10)
plot_pseudotime_heatmap(monocle_cds[sig_gene_names,],
                num_clusters = 3,
                cores = 1,
                show_rownames = T)
dev.off()





















