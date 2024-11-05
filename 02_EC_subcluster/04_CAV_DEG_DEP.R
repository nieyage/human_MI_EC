library(monocle)
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
# DEG for each CAV 
DefaultAssay(CAV_snRNA)<-"RNA"
Idents(CAV_snRNA)<- CAV_snRNA$CAV
CAV_snRNA <- ScaleData(CAV_snRNA,features=rownames(CAV_snRNA))

markers <- FindAllMarkers(CAV_snRNA, only.pos = TRUE, min.pct = 0.05, logfc.threshold = 0.1)
table(markers$cluster)
write.csv(markers,"./02_EC/CAV/DEG_DEP/FindAllMarkers_gene.csv")

# verification
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
signif_markers <- markers[markers$p_val_adj<0.05,] 
top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

DefaultAssay(CAV_snRNA)<- "RNA"
pdf("./02_EC/CAV/DEG_CAV_markers_top_markers-DEG_heatmap.pdf",width=10,height=10)
DoHeatmap(object = CAV_snRNA,features=top10$gene,group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = CAV_snRNA,features=markers$gene, group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
dev.off();


# GO and KEGG for CM subtype
Artery<-markers[markers$cluster=="Artery",7]
C_A<-markers[markers$cluster=="C_A",7]
Capillary<-markers[markers$cluster=="Capillary",7]
C_V<-markers[markers$cluster=="C_V",7]
Venous<-markers[markers$cluster=="Venous",7]

library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

CAV_list<- list(Artery,C_A,Capillary,C_V,Venous)
subtype<- c("Artery","C_A","Capillary","C_V","Venous")
all_ego<- data.frame()
pdf("./02_EC/CAV/EC_CAV_subtype_DEG_GO_BP.pdf")
for(i in 1:length(subtype)){
  gene<- CAV_list[[i]];
  gene.df <- bitr(gene, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)
    p<- barplot(ego, showCategory=20)
    print(p)
    write.csv(ego,paste0("./02_EC/CAV/",subtype[i],"-GO-BP.csv",sep=""))
    ego<- as.data.frame(ego)
    ego$celltype<- subtype[i];
    all_ego<- rbind(all_ego,ego);
}
dev.off()

write.csv(all_ego,"./02_EC/CAV/Allsubtype_GO.csv")

top_term <- all_ego %>% group_by(celltype) %>% top_n(n = -5, wt = pvalue);
top_term<- top_term[!duplicated(top_term$Description),]
top_term$Description<- factor(top_term$Description,levels=rev(top_term$Description))
top_term$celltype<- factor(top_term$celltype,levels=levels(CAV_snRNA))

top_term <- read.csv("/data/R02/nieyg/project/EC/human_MI/02_EC/CAV/Allsubtype_GO_select.csv") 
#top_term<- top_term[order(top_term$pvalue),]
library(ggplot2)
top_term$celltype<- factor(top_term$celltype,levels=levels(CAV_snRNA))
top_term$Description<- factor(top_term$Description,levels=rev(unique(top_term$Description)))
top_term<- top_term[order(top_term$Description),]
top_term<- top_term[order(top_term$celltype),]

pdf("./02_EC/CAV/EC_CAV_subtype_GO_top5_term.pdf",width=14,height=12)
p <- ggplot(top_term,aes(y=Count,x=Description,fill=-log10(pvalue))) + 
      geom_bar(stat="identity",position = "dodge") +
      facet_grid(celltype~.,scales = "free",space = "free") + 
      coord_flip() + 
      theme_bw() +
      scale_fill_gradient(high = '#FF0000', low = '#1202FF')+
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


DefaultAssay(CAV_snATAC)<-"peaks"
da_peaks <- FindAllMarkers(
  object = CAV_snATAC,
  test.use = 'LR',
  #logfc.threshold = 0.1,
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)
write.csv(da_peaks,"./02_EC/CAV/snATAC_CAV_DEP_FindAllMarkers.csv")
#markers <- FindAllMarkers(ORN, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
da_peaks<-da_peaks[da_peaks$p_val_adj<0.05,]
#da_peaks<- da_peaks[da_peaks$avg_log2FC>1,]
table(da_peaks$cluster)
#verification
peak2show<- rownames(da_peaks)


# plot the avg heatmap 
peak_Avg <- AverageExpression(CAV_snATAC,features=peak2show,assays = "peaks")
count=t(scale(t(peak_Avg$peaks),scale = T,center = F))
pdf("./02_EC/CAV/snATAC_CAV_DEP_FindAllMarkers_heatmap_avg.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
dev.off();

venous<- da_peaks[da_peaks$cluster=="Venous",]
# test enrichment
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
CAV_snATAC <- AddMotifs(
  object = CAV_snATAC,
  genome = library(BSgenome.Hsapiens.UCSC.hg38),
  pfm = pfm
)
enriched.motifs <- FindMotifs(
  object = CAV_snATAC,
  features = venous
)

MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(enriched.motifs))
)

mouse_brain <- RunChromVAR(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(mouse_brain) <- 'chromvar'

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = mouse_brain,
  features = "MA0497.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

differential.activity <- FindMarkers(
  object = mouse_brain,
  ident.1 = 'Pvalb',
  ident.2 = 'Sst',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)



# The DEGs between CA and A 


markers <- FindMarkers(CAV_snRNA,ident.1="C_A",ident.2="Artery", min.pct = 0.05, logfc.threshold = 0.1)
table(markers$cluster)
write.csv(markers,"./02_EC/CAV/DEG_DEP/FindAllMarkers_gene.csv")

# verification
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
signif_markers <- markers[markers$p_val_adj<0.05,] 
#top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

#自定义阈值
log2FC = 0.1
padj = 0.05 

markers$threshold="ns";
markers[which(markers$avg_log2FC  > log2FC & markers$p_val_adj <padj),]$threshold="up";
markers[which(markers$avg_log2FC  < (-log2FC) & markers$p_val_adj < padj),]$threshold="down";
markers$threshold=factor(markers$threshold, levels=c('down','up'))

pdf("./02_EC/CAV/CAvsA_DEG_volcano.pdf",width=5,height=5)
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  #labs(title= ifelse(""==title, "", paste("DEG:", title)))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)) )+
  theme_classic(base_size = 14) +
  scale_color_manual('',labels=c(paste0("down(",table(markers$threshold)[[1]],')'),
                                 paste0("up(",table(markers$threshold)[[2]],')' )),
                     values=c("blue","red" ) )+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))
dev.off();


DefaultAssay(CAV_snRNA)<- "RNA"
markers<- markers[order(markers$avg_log2FC),]
pdf("./02_EC/CAV/CAvsA_DEG_heatmap.pdf",width=10,height=10)
DoHeatmap(object = CAV_snRNA,features=rownames(markers),group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = CAV_snRNA,cells=rownames(CAV_snRNA@meta.data[CAV_snRNA$CAV%in%c("C_A","Artery"),]),features=rownames(markers), group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
dev.off();

# GO and KEGG for CM subtype
Artery<- rownames(markers[markers$avg_log2FC<0,])
#C_A<-markers[markers$cluster=="C_A",7]

library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

pdf("./02_EC/CAV/Artery_CAvsA_DEG_GO_BP.pdf")

  gene.df <- bitr(Artery, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)
    p<- barplot(ego, showCategory=20)
    print(p)
    write.csv(ego,"Artery_CAvsA_DEG_GO_BP")
dev.off()


DefaultAssay(CAV_snATAC)<-"peaks"
da_peaks <- FindMarkers(
  object = CAV_snATAC,
  test.use = 'LR',ident.1="C_A",ident.2="Artery",
  #logfc.threshold = 0.1,
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)


write.csv(da_peaks,"./02_EC/CAV/CAvsA_DEP_FindAllMarkers.csv")
da_peaks<- read.csv("./02_EC/CAV/CAvsA_DEP_FindAllMarkers.csv",row.names=1)
#markers <- FindAllMarkers(ORN, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
da_peaks<-da_peaks[da_peaks$p_val_adj<0.05,]
da_peaks<- da_peaks[order(da_peaks$avg_log2FC),]
#verification
peak2show<- rownames(da_peaks)
CAV_snATAC<- ScaleData(CAV_snATAC)
pdf("./02_EC/CAV/CAvsA_DEP_heatmap.pdf",width=10,height=10)
#DoHeatmap(object = CAV_snATAC,features=peak2show,group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
#DoHeatmap(object = CAV_snATAC,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),cells=rownames(CAV_snATAC@meta.data[CAV_snRNA$CAV%in%c("C_A","Artery"),]),features=peak2show, group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = CAV_snATAC,cells=rownames(CAV_snATAC@meta.data[CAV_snATAC$CAV%in%c("C_A","Artery"),]),features=peak2show, group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = "#377EB8", "white", "#E41A1C")
DoHeatmap(object = CAV_snATAC,features=peak2show, group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = "#377EB8", "white", "#E41A1C")
dev.off();

log2FC = 0.1
padj = 0.05 

da_peaks$threshold="ns";
da_peaks[which(da_peaks$avg_log2FC  > log2FC & da_peaks$p_val_adj <padj),]$threshold="up";
da_peaks[which(da_peaks$avg_log2FC  < (-log2FC) & da_peaks$p_val_adj < padj),]$threshold="down";
da_peaks$threshold=factor(da_peaks$threshold, levels=c('down','up'))

pdf("./02_EC/CAV/CAvsA_DEP_volcano.pdf",width=5,height=5)
ggplot(data=da_peaks, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  #labs(title= ifelse(""==title, "", paste("DEG:", title)))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)) )+
  theme_classic(base_size = 14) +
  scale_color_manual('',labels=c(paste0("down(",table(da_peaks$threshold)[[1]],')'),
                                 paste0("up(",table(da_peaks$threshold)[[2]],')' )),
                     values=c("blue","red" ) )+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))
dev.off();



# plot the avg heatmap 
peak_Avg <- AverageExpression(CAV_snATAC,features=peak2show,assays = "peaks")
count=t(scale(t(peak_Avg$peaks),scale = T,center = F))
pdf("./02_EC/CAV/CAvsA_snATAC_CAV_DEP_FindAllMarkers_heatmap_avg.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
dev.off();

C_A<- rownames(da_peaks[da_peaks$avg_log2FC>0.1,])
Artery<- rownames(da_peaks[da_peaks$avg_log2FC< -0.1,])

head(C_A)
C_A_df <- data.frame(
  chromosome = sub("-.*", "", C_A),
  start = as.numeric(sub(".*-(\\d+)-.*", "\\1", C_A)),
  end = as.numeric(sub(".*-\\d+-(\\d+)", "\\1", C_A))
)
write.table(C_A_df, file = "./02_EC/CAV/CAvsA_C_A.bed", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

Artery_df <- data.frame(
  chromosome = sub("-.*", "", Artery),
  start = as.numeric(sub(".*-(\\d+)-.*", "\\1", Artery)),
  end = as.numeric(sub(".*-\\d+-(\\d+)", "\\1", Artery))
)
write.table(Artery_df, file = "./02_EC/CAV/CAvsA_Artery.bed", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

findMotifsGenome.pl CAvsA_Artery.bed  hg19 CAvsA_Artery_Homer -len 6,8,10,12  
findMotifsGenome.pl CAvsA_C_A.bed  hg19 CAvsA_C_A_Homer -len 6,8,10,12  


# test enrichment
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

# add motif information
CAV_snATAC <- AddMotifs(
  object = CAV_snATAC,
  genome = library(BSgenome.Hsapiens.UCSC.hg38),
  pfm = pfm
)
enriched.motifs <- FindMotifs(
  object = CAV_snATAC,
  features = C_A
)

MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(enriched.motifs))
)

mouse_brain <- RunChromVAR(
  object = mouse_brain,
  genome = BSgenome.Mmusculus.UCSC.mm10
)

DefaultAssay(mouse_brain) <- 'chromvar'

# look at the activity of Mef2c
p2 <- FeaturePlot(
  object = mouse_brain,
  features = "MA0497.1",
  min.cutoff = 'q10',
  max.cutoff = 'q90',
  pt.size = 0.1
)
p1 + p2

differential.activity <- FindMarkers(
  object = mouse_brain,
  ident.1 = 'Pvalb',
  ident.2 = 'Sst',
  only.pos = TRUE,
  mean.fxn = rowMeans,
  fc.name = "avg_diff"
)

MotifPlot(
  object = mouse_brain,
  motifs = head(rownames(differential.activity)),
  assay = 'peaks'
)



# The DEGs between Artery_1 and Artery_2 
Idents(CAV_snRNA)<- CAV_snRNA$CAV_detail
CAV_snRNA <- ScaleData(CAV_snRNA,features=rownames(CAV_snRNA))

CAV_snRNA<- SCTransform(CAV_snRNA, verbose = FALSE)
CAV_snRNA <- ScaleData(CAV_snRNA,features=rownames(CAV_snRNA))
markers <- FindMarkers(CAV_snRNA,slot = "scale.data",assays="SCT",ident.1="Artery_1",ident.2="Artery_2", min.pct = 0, logfc.threshold = 0)
table(markers$cluster)



DefaultAssay(CAV_snRNA)<- "RNA"
CAV_snRNA<- NormalizeData(CAV_snRNA, verbose = FALSE)
markers <- FindMarkers(CAV_snRNA,slot = "data",assays="RNA",ident.1="Artery_1",ident.2="Artery_2", min.pct = 0, logfc.threshold = 0, test.use = "negbinom")
table(markers$cluster)

write.csv(markers,"./02_EC/CAV/DEG_DEP/Artery_1_vs_Artery_2_FindAllMarkers_gene.csv")

#Artery_markers<- c("GABBR2","GRIA2","SSUH2","JAG1",#embryonic
#  "SEMA3G","EFNB2","DLL4"#adult
#  )

Artery_markers<- c("EFNB2","DLL4")
Idents(CAV_snRNA)<- factor(Idents(CAV_snRNA),levels=rev(levels(CAV_snRNA)))
pdf("./02_EC/CAV/Artery_1_vs_Artery_2_makrer_volcano.pdf",width=4,height=4)
VlnPlot(CAV_snRNA, features = Artery_markers,idents=c("Artery_1","Artery_2"),
 stack = TRUE, flip = TRUE) +
  theme(legend.position = "none") +geom_boxplot(width=0.06,fill="white",outlier.size=0)
dev.off()

# verification
myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")
signif_markers <- markers[markers$p_val_adj<0.05,] 
#top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC);
blueYellow = c("1"="#352A86","2"="#343DAE","3"="#0262E0","4"="#1389D2","5"="#2DB7A3","6"="#A5BE6A","7"="#F8BA43","8"="#F6DA23","9"="#F8FA0D")
solarExtra<- c("#3361A5", "#248AF3", "#14B3FF" ,"#88CEEF" ,"#C1D5DC", "#EAD397" ,"#FDB31A" ,"#E42A2A","#A31D1D")

#自定义阈值
log2FC = 0.1
padj = 0.05 

markers$threshold="ns";
markers[which(markers$avg_log2FC  > log2FC & markers$p_val_adj <padj),]$threshold="up";
markers[which(markers$avg_log2FC  < (-log2FC) & markers$p_val_adj < padj),]$threshold="down";
markers$threshold=factor(markers$threshold, levels=c('down','up'))

pdf("./02_EC/CAV/CAvsA_DEG_volcano.pdf",width=5,height=5)
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  #labs(title= ifelse(""==title, "", paste("DEG:", title)))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)) )+
  theme_classic(base_size = 14) +
  scale_color_manual('',labels=c(paste0("down(",table(markers$threshold)[[1]],')'),
                                 paste0("up(",table(markers$threshold)[[2]],')' )),
                     values=c("blue","red" ) )+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))
dev.off();


DefaultAssay(CAV_snRNA)<- "RNA"
markers<- markers[order(markers$avg_log2FC),]
pdf("./02_EC/CAV/CAvsA_DEG_heatmap.pdf",width=10,height=10)
DoHeatmap(object = CAV_snRNA,features=rownames(markers),group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = CAV_snRNA,cells=rownames(CAV_snRNA@meta.data[CAV_snRNA$CAV%in%c("C_A","Artery"),]),features=rownames(markers), group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
dev.off();

# GO and KEGG for CM subtype
Artery<- rownames(markers[markers$avg_log2FC<0,])
#C_A<-markers[markers$cluster=="C_A",7]

library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

pdf("./02_EC/CAV/Artery_CAvsA_DEG_GO_BP.pdf")

  gene.df <- bitr(Artery, fromType = "SYMBOL",
                toType = c("ENSEMBL", "ENTREZID"),
                OrgDb = org.Hs.eg.db)
    ego <- enrichGO(gene.df$ENTREZID,
                keyType = 'ENTREZID',
                OrgDb = org.Hs.eg.db,
                ont = "BP", ###BP,MF,CC
                pAdjustMethod  = "BH",
                pvalueCutoff = 0.5,
                qvalueCutoff = 0.5,
                readable = TRUE)
    p<- barplot(ego, showCategory=20)
    print(p)
    write.csv(ego,"Artery_CAvsA_DEG_GO_BP")
dev.off()


DefaultAssay(CAV_snATAC)<-"peaks"
da_peaks <- FindMarkers(
  object = CAV_snATAC,
  test.use = 'LR',ident.1="C_A",ident.2="Artery",
  #logfc.threshold = 0.1,
  min.pct = 0.05,
  latent.vars = 'nCount_peaks'
)


write.csv(da_peaks,"./02_EC/CAV/CAvsA_DEP_FindAllMarkers.csv")
da_peaks<- read.csv("./02_EC/CAV/CAvsA_DEP_FindAllMarkers.csv",row.names=1)
#markers <- FindAllMarkers(ORN, only.pos = TRUE, min.pct = 0.2, logfc.threshold = 0.2)
da_peaks<-da_peaks[da_peaks$p_val_adj<0.05,]
da_peaks<- da_peaks[order(da_peaks$avg_log2FC),]
#verification
peak2show<- rownames(da_peaks)
CAV_snATAC<- ScaleData(CAV_snATAC)
pdf("./02_EC/CAV/CAvsA_DEP_heatmap.pdf",width=10,height=10)
#DoHeatmap(object = CAV_snATAC,features=peak2show,group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
#DoHeatmap(object = CAV_snATAC,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),cells=rownames(CAV_snATAC@meta.data[CAV_snRNA$CAV%in%c("C_A","Artery"),]),features=peak2show, group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = solarExtra[3:8])
DoHeatmap(object = CAV_snATAC,cells=rownames(CAV_snATAC@meta.data[CAV_snATAC$CAV%in%c("C_A","Artery"),]),features=peak2show, group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = "#377EB8", "white", "#E41A1C")
DoHeatmap(object = CAV_snATAC,features=peak2show, group.colors =myUmapcolors,size = 2,group.by = "CAV") #+scale_fill_gradientn(colors = "#377EB8", "white", "#E41A1C")
dev.off();

log2FC = 0.1
padj = 0.05 

da_peaks$threshold="ns";
da_peaks[which(da_peaks$avg_log2FC  > log2FC & da_peaks$p_val_adj <padj),]$threshold="up";
da_peaks[which(da_peaks$avg_log2FC  < (-log2FC) & da_peaks$p_val_adj < padj),]$threshold="down";
da_peaks$threshold=factor(da_peaks$threshold, levels=c('down','up'))

pdf("./02_EC/CAV/CAvsA_DEP_volcano.pdf",width=5,height=5)
ggplot(data=da_peaks, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  #labs(title= ifelse(""==title, "", paste("DEG:", title)))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)) )+
  theme_classic(base_size = 14) +
  scale_color_manual('',labels=c(paste0("down(",table(da_peaks$threshold)[[1]],')'),
                                 paste0("up(",table(da_peaks$threshold)[[2]],')' )),
                     values=c("blue","red" ) )+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))
dev.off();



# plot the avg heatmap 
peak_Avg <- AverageExpression(CAV_snATAC,features=peak2show,assays = "peaks")
count=t(scale(t(peak_Avg$peaks),scale = T,center = F))
pdf("./02_EC/CAV/CAvsA_snATAC_CAV_DEP_FindAllMarkers_heatmap_avg.pdf",width=15,height=15)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "firebrick3"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
pheatmap(count,cluster_cols = F,cluster_rows = F,color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),show_rownames=F,show_colnames=T)
dev.off();

C_A<- rownames(da_peaks[da_peaks$avg_log2FC>0.1,])
Artery<- rownames(da_peaks[da_peaks$avg_log2FC< -0.1,])

head(C_A)
C_A_df <- data.frame(
  chromosome = sub("-.*", "", C_A),
  start = as.numeric(sub(".*-(\\d+)-.*", "\\1", C_A)),
  end = as.numeric(sub(".*-\\d+-(\\d+)", "\\1", C_A))
)
write.table(C_A_df, file = "./02_EC/CAV/CAvsA_C_A.bed", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

Artery_df <- data.frame(
  chromosome = sub("-.*", "", Artery),
  start = as.numeric(sub(".*-(\\d+)-.*", "\\1", Artery)),
  end = as.numeric(sub(".*-\\d+-(\\d+)", "\\1", Artery))
)
write.table(Artery_df, file = "./02_EC/CAV/CAvsA_Artery.bed", row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)

findMotifsGenome.pl CAvsA_Artery.bed  hg19 CAvsA_Artery_Homer -len 6,8,10,12  
findMotifsGenome.pl CAvsA_C_A.bed  hg19 CAvsA_C_A_Homer -len 6,8,10,12  


# test enrichment
library(JASPAR2020)
library(TFBSTools)
library(patchwork)
library(BSgenome.Hsapiens.UCSC.hg38)
# Get a list of motif position frequency matrices from the JASPAR database
pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

