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
library(SeuratWrappers)
CAV_snRNA<- readRDS("./02_EC/CAV/EC_CAV_snRNA_add_gene_symbol.rds")
CAV_snATAC<- readRDS("./02_EC/CAV/EC_CAV_snATAC_add_peaks_fragment.rds")
# Building trajectories with Monocle 3 for snRNA of CAV 
Idents(CAV_snRNA)<- CAV_snRNA$CAV
CAV_snRNA@reductions$UMAP<- CAV_snRNA@reductions$umap_harmony
CAV_snRNA.cds <- as.cell_data_set(CAV_snRNA)
CAV_snRNA.cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(CAV_snRNA.cds)
CAV_snRNA.cds <- cluster_cells(cds = CAV_snRNA.cds,
  cluster_method =  "louvain",num_iter=5,
 reduction_method = "UMAP")# partition
CAV_snRNA.cds <- learn_graph(CAV_snRNA.cds, use_partition = T)
# plot trajectories colored by pseudotime
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, days="Venous"){
  cell_ids <- which(colData(cds)[, "CAV"] == days)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

 CAV_snRNA.cds <- order_cells(CAV_snRNA.cds, root_pr_nodes=get_earliest_principal_node(CAV_snRNA.cds))

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

#CAV_snRNA.cds <- order_cells(CAV_snRNA.cds, reduction_method = "UMAP", root_cells = Capillary)
pdf("./02_EC/CAV/monocle/snRNA_pseudotime_umap.pdf",width=6,height=5)
plot_cells(cds = CAV_snRNA.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
plot_cells(cds = CAV_snRNA.cds,label_leaves=FALSE,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
dev.off()
pdf("./02_EC/CAV/monocle/snRNA_trajectory_umap.pdf",width=5,height=5)
plot_cells(cds = CAV_snRNA.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = CAV_snRNA.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV_detail",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = CAV_snRNA.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "days.after.infarction",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
dev.off()
# add pseudotime info in seurat object 

CAV_snRNA <- AddMetaData(
  object = CAV_snRNA,
  metadata = CAV_snRNA.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "CAV_pseudotime"
)

saveRDS(CAV_snRNA.cds,"./02_EC/CAV/monocle/CAV_snRNA_add_pseudotime_info_monocle.rds")

# split the C to A and C to V
# PART1: C to A 
Idents(CAV_snRNA)<- CAV_snRNA$CAV
C2A_snRNA<- subset(CAV_snRNA,idents=c("Capillary","C_A","Artery"))
#C2A_snRNA@reductions$UMAP<- C2A_snRNA@reductions$umap_harmony
C2A_snRNA.cds <- as.cell_data_set(C2A_snRNA)
C2A_snRNA.cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(C2A_snRNA.cds)
C2A_snRNA.cds <- cluster_cells(cds = C2A_snRNA.cds, reduction_method = "UMAP")# partition
C2A_snRNA.cds <- learn_graph(C2A_snRNA.cds, use_partition = T)
# plot trajectories colored by pseudotime
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, days="control"){
  cell_ids <- which(colData(cds)[, "days.after.infarction"] == days)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

C2A_snRNA.cds <- order_cells(C2A_snRNA.cds, root_pr_nodes=get_earliest_principal_node(C2A_snRNA.cds))

pdf("./02_EC/CAV/monocle/C_to_A_snRNA_pseudotime_umap.pdf",width=6,height=5)
plot_cells(cds = C2A_snRNA.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
plot_cells(cds = C2A_snRNA.cds,label_leaves=FALSE,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
dev.off()
pdf("./02_EC/CAV/monocle/C_to_A_snRNA_trajectory_umap.pdf",width=5,height=5)
plot_cells(cds = C2A_snRNA.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = C2A_snRNA.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV_detail",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = C2A_snRNA.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "days.after.infarction",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
dev.off()

# PART2: C to V 
Idents(CAV_snRNA)<- CAV_snRNA$CAV
C2V_snRNA<- subset(CAV_snRNA,idents=c("Capillary","C_V","Venous"))
#C2V_snRNA@reductions$UMAP<- C2V_snRNA@reductions$umap_harmony
C2V_snRNA.cds <- as.cell_data_set(C2V_snRNA)
C2V_snRNA.cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(C2V_snRNA.cds)
C2V_snRNA.cds <- cluster_cells(cds = C2V_snRNA.cds, reduction_method = "UMAP")# partition
C2V_snRNA.cds <- learn_graph(C2V_snRNA.cds, use_partition = T)
# plot trajectories colored by pseudotime
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, days="control"){
  cell_ids <- which(colData(cds)[, "days.after.infarction"] == days)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

C2V_snRNA.cds <- order_cells(C2V_snRNA.cds, root_pr_nodes=get_earliest_principal_node(C2V_snRNA.cds))

pdf("./02_EC/CAV/monocle/C_to_V_snRNA_pseudotime_umap.pdf",width=6,height=5)
plot_cells(cds = C2V_snRNA.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
plot_cells(cds = C2V_snRNA.cds,label_leaves=FALSE,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
dev.off()
pdf("./02_EC/CAV/monocle/C_to_V_snRNA_trajectory_umap.pdf",width=5,height=5)
plot_cells(cds = C2V_snRNA.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = C2V_snRNA.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV_detail",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = C2V_snRNA.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "days.after.infarction",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
dev.off()
# add pseudotime info in seurat object 

CAV_snRNA <- AddMetaData(
  object = CAV_snRNA,
  metadata = C2A_snRNA.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "C_to_A"
)

CAV_snRNA <- AddMetaData(
  object = CAV_snRNA,
  metadata = C2V_snRNA.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "C_to_V"
)

pdf("./02_EC/CAV/monocle/C_to_V_or_A_snRNA_trajectory_umap.pdf",width=10,height=5)
FeaturePlot(CAV_snRNA, c("C_to_A", "C_to_V"), pt.size = 0.1) & scale_color_viridis_c()
dev.off()

saveRDS(CAV_snRNA,"./02_EC/CAV/monocle/CAV_snRNA_add_pseudotime_info.rds")
saveRDS(C2A_snRNA.cds,"./02_EC/CAV/monocle/C_to_A_CAV_snRNA_add_pseudotime_info_monocle.rds")
saveRDS(C2V_snRNA.cds,"./02_EC/CAV/monocle/C_to_V_CAV_snRNA_add_pseudotime_info_monocle.rds")


# Differential expression analysis

gene_fits <- fit_models(CAV_snRNA.cds, model_formula_str = "~days.after.infarction")
fit_coefs <- coefficient_table(gene_fits)
# Finding genes that change as a function of pseudotime
Track_genes <- graph_test(CAV_snRNA.cds, neighbor_graph="principal_graph", cores=1)
Track_genes <- Track_genes[,c(5,2,3,4,1,6)] %>% filter(q_value<1e-3)
Track_genes_sig<- Track_genes%>% top_n(n=10,morans_I) %>% pull(gene_ID) %>% as.character()

features <- c("CA4","PKD1L1","BTNL9","CLIC5","RGCC","SLC9C1","LNX1","F8","CD36","MGLL",# capillary Endo
	"FAM155A","TSHZ2","SLCO2A1","ZBTB7C","KCNIP4","LYST","IGFBP5","SNTG2","POSTN",# venous Endo
	"PCSK5","ARL15","LINC00639","SMAD6","NEBL","MECOM","FUT8","PRDM16","PDZD2","SLC45A4")# Arterial Endo)
cds_subset <- CAV_snRNA.cds[row.names(subset(rowData(CAV_snRNA.cds), gene_short_name %in% features)),]
pdf("./02_EC/CAV/monocle/all_CAV/snRNA_CAV_markers_trajectory_umap.pdf",width=20,height=10)
plot_genes_violin(cds_subset, group_cells_by="days.after.infarction", ncol=4) + theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()
pdf("./02_EC/CAV/monocle/all_CAV/snRNA_CAV_markers_plot_genes_in_pseudotime.pdf",width=20,height=10)
plot_genes_in_pseudotime(CAV_snRNA.cds[features,],color_cells_by="CAV",min_expr=0.5,ncol=4)
dev.off()
pdf("./02_EC/CAV/monocle/all_CAV/snRNA_CAV_markers_plot_genes_in_umap.pdf",width=10,height=10)
plot_cells(CAV_snRNA.cds, genes=features[1:4],show_trajectory_graph=FALSE,label_cell_groups=FALSE,label_leaves=FALSE)
dev.off()

gene_module_df <- find_gene_modules(CAV_snRNA.cds[features,], resolution=c(10^seq(-6,-1)))
cell_group_df <- tibble::tibble(cell=row.names(colData(CAV_snRNA.cds)), cell_group=colData(CAV_snRNA.cds)$CAV)
agg_mat <- aggregate_gene_expression(CAV_snRNA.cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
pheatmap::pheatmap(agg_mat,scale="column", clustering_method="ward.D2")



# Building trajectories with Monocle 3 for snATAC of CAV 
Idents(CAV_snATAC)<- CAV_snATAC$CAV
CAV_snATAC@reductions$UMAP<- CAV_snATAC@reductions$umap_harmony
CAV_snATAC.cds <- as.cell_data_set(CAV_snATAC)
CAV_snATAC.cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(CAV_snATAC.cds)
CAV_snATAC.cds <- cluster_cells(cds = CAV_snATAC.cds, reduction_method = "UMAP")# partition
CAV_snATAC.cds <- learn_graph(CAV_snATAC.cds, use_partition = T)
# plot trajectories colored by pseudotime
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, days="control"){
  cell_ids <- which(colData(cds)[, "days.after.infarction"] == days)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}
CAV_snATAC.cds <- order_cells(CAV_snATAC.cds, root_pr_nodes=get_earliest_principal_node(CAV_snATAC.cds))
pdf("./02_EC/CAV/monocle/snATAC_pseudotime_umap.pdf",width=6,height=5)
plot_cells(cds = CAV_snATAC.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
plot_cells(cds = CAV_snATAC.cds,label_leaves=FALSE,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
dev.off()
pdf("./02_EC/CAV/monocle/snATAC_trajectory_umap.pdf",width=5,height=5)
plot_cells(cds = CAV_snATAC.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = CAV_snATAC.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV_detail",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = CAV_snATAC.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "days.after.infarction",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
dev.off()

# add pseudotime info in seurat object 
CAV_snATAC <- AddMetaData(
  object = CAV_snATAC,
  metadata = CAV_snATAC.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "CAV_pseudotime"
)

saveRDS(CAV_snATAC.cds,"./02_EC/CAV/monocle/CAV_snATAC_add_pseudotime_info_monocle.rds")

# split the C to A and C to V
# PART1: C to A 
Idents(CAV_snATAC)<- CAV_snATAC$CAV
C2A_snATAC<- subset(CAV_snATAC,idents=c("Capillary","C_A","Artery"))
#C2A_snATAC@reductions$UMAP<- C2A_snATAC@reductions$umap_harmony
C2A_snATAC.cds <- as.cell_data_set(C2A_snATAC)
C2A_snATAC.cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(C2A_snATAC.cds)
C2A_snATAC.cds <- cluster_cells(cds = C2A_snATAC.cds, reduction_method = "UMAP")# partition
C2A_snATAC.cds <- learn_graph(C2A_snATAC.cds, use_partition = T)
# plot trajectories colored by pseudotime
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, days="control"){
  cell_ids <- which(colData(cds)[, "days.after.infarction"] == days)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

C2A_snATAC.cds <- order_cells(C2A_snATAC.cds, root_pr_nodes=get_earliest_principal_node(C2A_snATAC.cds))

pdf("./02_EC/CAV/monocle/C_to_A_snATAC_pseudotime_umap.pdf",width=6,height=5)
plot_cells(cds = C2A_snATAC.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
plot_cells(cds = C2A_snATAC.cds,label_leaves=FALSE,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
dev.off()
pdf("./02_EC/CAV/monocle/C_to_A_snATAC_trajectory_umap.pdf",width=5,height=5)
plot_cells(cds = C2A_snATAC.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = C2A_snATAC.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV_detail",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = C2A_snATAC.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "days.after.infarction",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
dev.off()

# PART2: C to V 
Idents(CAV_snATAC)<- CAV_snATAC$CAV
C2V_snATAC<- subset(CAV_snATAC,idents=c("Capillary","C_V","Venous"))
#C2V_snATAC@reductions$UMAP<- C2V_snATAC@reductions$umap_harmony
C2V_snATAC.cds <- as.cell_data_set(C2V_snATAC)
C2V_snATAC.cds@rowRanges@elementMetadata@listData$gene_short_name <- rownames(C2V_snATAC.cds)
C2V_snATAC.cds <- cluster_cells(cds = C2V_snATAC.cds, reduction_method = "UMAP")# partition
C2V_snATAC.cds <- learn_graph(C2V_snATAC.cds, use_partition = T)
# plot trajectories colored by pseudotime
# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, days="control"){
  cell_ids <- which(colData(cds)[, "days.after.infarction"] == days)
  closest_vertex <-
  cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
  igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
  (which.max(table(closest_vertex[cell_ids,]))))]
  root_pr_nodes
}

C2V_snATAC.cds <- order_cells(C2V_snATAC.cds, root_pr_nodes=get_earliest_principal_node(C2V_snATAC.cds))

pdf("./02_EC/CAV/monocle/C_to_V_snATAC_pseudotime_umap.pdf",width=6,height=5)
plot_cells(cds = C2V_snATAC.cds,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
plot_cells(cds = C2V_snATAC.cds,label_leaves=FALSE,color_cells_by = "pseudotime",show_trajectory_graph = TRUE)
dev.off()
pdf("./02_EC/CAV/monocle/C_to_V_snATAC_trajectory_umap.pdf",width=5,height=5)
plot_cells(cds = C2V_snATAC.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = C2V_snATAC.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "CAV_detail",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
plot_cells(cds = C2V_snATAC.cds,group_label_size = 4.5,label_leaves=FALSE,color_cells_by = "days.after.infarction",show_trajectory_graph = TRUE)+ scale_color_manual(values = myUmapcolors)
dev.off()
# add pseudotime info in seurat object 

CAV_snATAC <- AddMetaData(
  object = CAV_snATAC,
  metadata = C2A_snATAC.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "C_to_A"
)

CAV_snATAC <- AddMetaData(
  object = CAV_snATAC,
  metadata = C2V_snATAC.cds@principal_graph_aux@listData$UMAP$pseudotime,
  col.name = "C_to_V"
)

pdf("./02_EC/CAV/monocle/C_to_V_or_A_snATAC_trajectory_umap.pdf",width=10,height=5)
FeaturePlot(CAV_snATAC, c("C_to_A", "C_to_V"), pt.size = 0.1) & scale_color_viridis_c()
dev.off()

saveRDS(CAV_snATAC,"./02_EC/CAV/monocle/CAV_snATAC_add_pseudotime_info.rds")
saveRDS(C2A_snATAC.cds,"./02_EC/CAV/monocle/C_to_A_CAV_snATAC_add_pseudotime_info_monocle.rds")
saveRDS(C2V_snATAC.cds,"./02_EC/CAV/monocle/C_to_V_CAV_snATAC_add_pseudotime_info_monocle.rds")


