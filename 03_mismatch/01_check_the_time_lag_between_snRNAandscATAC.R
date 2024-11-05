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
DefaultAssay(CAV_snRNA)<-"RNA"
Idents(CAV_snRNA)<- CAV_snRNA$CAV
CV_V_snRNA<- subset(CAV_snRNA,idents=c("C_V","Venous"))
Idents(CV_V_snRNA)<- CV_V_snRNA$major_labl
IZvsCTRL_CVvsV_snRNA<- subset(CV_V_snRNA,idents=c("IZ","CTRL"))

DefaultAssay(CAV_snATAC)<-"RNA"
Idents(CAV_snATAC)<- CAV_snATAC$CAV
CV_V_snATAC<- subset(CAV_snATAC,idents=c("C_V","Venous"))
Idents(CV_V_snATAC)<- CV_V_snATAC$major_labl
IZvsCTRL_CVvsV_snATAC<- subset(CV_V_snATAC,idents=c("IZ","CTRL"))

## snRNA:
# The DEG between IZ and CTRL;
Idents(IZvsCTRL_CVvsV_snRNA)<- IZvsCTRL_CVvsV_snRNA$major_labl
IZvsCTRL_snRNA_markers <- FindMarkers(IZvsCTRL_CVvsV_snRNA,ident.1="IZ",ident.2="CTRL", min.pct = 0.05, logfc.threshold = 0)
log2FC = 0.1
padj = 0.05 
markers<- IZvsCTRL_snRNA_markers
markers$threshold="ns";
markers[which(markers$avg_log2FC  > log2FC & markers$p_val_adj <padj),]$threshold="up";
markers[which(markers$avg_log2FC  < (-log2FC) & markers$p_val_adj < padj),]$threshold="down";
markers$threshold=factor(markers$threshold, levels=c('down','up'))

pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/IZvsCTRL_DEG_volcano.pdf",width=5,height=5)
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +xlim(-0.5,0.5)+ylim(0,100)+ggtitle("C_V vs Venous")+
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

# The DEG between CV and V;
Idents(IZvsCTRL_CVvsV_snRNA)<- IZvsCTRL_CVvsV_snRNA$CAV
CVvsV_snRNA_markers <- FindMarkers(IZvsCTRL_CVvsV_snRNA,ident.1="C_V",ident.2="Venous", min.pct = 0.05, logfc.threshold = 0)
log2FC = 0.1
padj = 0.05 
markers<- CVvsV_snRNA_markers
markers$threshold="ns";
markers[which(markers$avg_log2FC  > log2FC & markers$p_val_adj <padj),]$threshold="up";
markers[which(markers$avg_log2FC  < (-log2FC) & markers$p_val_adj < padj),]$threshold="down";
markers$threshold=factor(markers$threshold, levels=c('down','up'))

pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/CVvsV_DEG_volcano.pdf",width=5,height=5)
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +xlim(-0.5,0.5)+ylim(0,100)+ggtitle("C_V vs Venous")+
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

# quad plot 
library(tidyverse)
library(Seurat)
library(ggrepel)

# parpare input data

IZvsCTRL_snRNA_gene <- rownames(IZvsCTRL_snRNA_markers)
CVvsV_snRNA_gene <- rownames(CVvsV_snRNA_markers)
geneset <- c(unique(intersect(IZvsCTRL_snRNA_gene, CVvsV_snRNA_gene)))
IZvsCTRL_snRNA_markers <- IZvsCTRL_snRNA_markers[rownames(IZvsCTRL_snRNA_markers) %in% geneset, ]
IZvsCTRL_snRNA_markers$gene_name <- rownames(IZvsCTRL_snRNA_markers)
colnames(IZvsCTRL_snRNA_markers)[1:5] <- c("p_val_IZvsCTRL", "avg_log2FC_IZvsCTRL", "pct.1_IZvsCTRL", "pct.2_IZvsCTRL", "p_val_adj_IZvsCTRL")
IZvsCTRL_snRNA_markers <- IZvsCTRL_snRNA_markers %>%arrange(gene_name)
CVvsV_snRNA_markers <- CVvsV_snRNA_markers[rownames(CVvsV_snRNA_markers) %in% geneset, ]
CVvsV_snRNA_markers$gene_name <- rownames(CVvsV_snRNA_markers)
colnames(CVvsV_snRNA_markers)[1:5] <- c("p_val_CVvsV", "avg_log2FC_CVvsV", "pct.1_CVvsV", "pct.2_CVvsV", "p_val_adj_CVvsV")
CVvsV_snRNA_markers <- CVvsV_snRNA_markers %>%arrange(gene_name)
quad_df <- data.frame(cbind(IZvsCTRL_snRNA_markers, CVvsV_snRNA_markers))
head(quad_df)
quad_df_sign<- quad_df[quad_df$p_val_IZvsCTRL<0.05&quad_df$p_val_CVvsV<0.05,]
quad_df_sign$type="nosig";
log2FC=0.1
quad_df_sign[which(quad_df_sign$avg_log2FC_IZvsCTRL  > log2FC & quad_df_sign$avg_log2FC_CVvsV >log2FC),]$type="IZ_CV";
quad_df_sign[which(quad_df_sign$avg_log2FC_IZvsCTRL  > log2FC & quad_df_sign$avg_log2FC_CVvsV < -log2FC),]$type="IZ_V";
quad_df_sign[which(quad_df_sign$avg_log2FC_IZvsCTRL  < -log2FC & quad_df_sign$avg_log2FC_CVvsV >log2FC),]$type="CTRL_CV";
quad_df_sign[which(quad_df_sign$avg_log2FC_IZvsCTRL  < -log2FC & quad_df_sign$avg_log2FC_CVvsV < -log2FC),]$type="CTRL_V";
quad_df_sign$type<- factor(quad_df_sign$type,levels=c("IZ_CV","IZ_V","CTRL_CV","CTRL_V","nosig"))
quad_plot <- quad_df_sign %>%
  ggplot(aes(x = avg_log2FC_CVvsV, y = avg_log2FC_IZvsCTRL,color=type)) + #, label = label
  #annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "#457B9D", alpha = 0.65) +
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= "#A8DADC", alpha = 0.65) +
  #annotate("rect", xmin = Inf, xmax = 0, ymin = -Inf, ymax = 0, fill= "#F1FAEE", alpha = 0.65) +
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = 0 , fill= "#E63946", alpha = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.25, color = "black", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.25, color = "black", alpha = 0.85) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual('',labels=c(paste0("IZ_CV(",table(quad_df_sign$type)[[1]],')'),
                                 paste0("IZ_V(",table(quad_df_sign$type)[[2]],')' ),
                                 paste0("CTRL_CV(",table(quad_df_sign$type)[[3]],')'),
                                 paste0("CTRL_V(",table(quad_df_sign$type)[[4]],')' ),
                                 paste0("not_sig(",table(quad_df_sign$type)[[5]],')' )),
                     values=c("#88B04B","#92A8D1","#6B5B95","#D64161","grey") )+
  scale_x_continuous(limits = c(-0.8, 0.8)) +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  xlab("Log2FC(C_V / Venous)") +
  ylab("Log2FC(IZ / CTRL)") +
  ggtitle("snRNA:C_V vs Venous & IZ vs CTRL (avg_log2FC>0.1)")+
  theme_linedraw() +
  theme(axis.text = element_text(size = 15), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        aspect.ratio = 5/5)
pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/IZvsCTRL_CVvsV_DEG_quadplot.pdf",width=7,height=7)
quad_plot
dev.off()
quad_df_sign$CVvsV_IZvsCTRL<- quad_df_sign$avg_log2FC_IZvsCTRL*quad_df_sign$avg_log2FC_CVvsV
IZ_CV_data <- quad_df_sign[quad_df_sign$type=="IZ_CV",]
IZ_CV_data<- IZ_CV_data[order(IZ_CV_data$CVvsV_IZvsCTRL,decreasing=T),]
IZ_CV <- rownames(IZ_CV_data)

IZ_V_data <- quad_df_sign[quad_df_sign$type=="IZ_V",]
IZ_V_data<- IZ_V_data[order(IZ_V_data$CVvsV_IZvsCTRL,decreasing=F),]
IZ_V <- rownames(IZ_V_data)

CTRL_CV_data <- quad_df_sign[quad_df_sign$type=="CTRL_CV",]
CTRL_CV_data<- CTRL_CV_data[order(CTRL_CV_data$CVvsV_IZvsCTRL,decreasing=F),]
CTRL_CV <- rownames(CTRL_CV_data)

CTRL_V_data <- quad_df_sign[quad_df_sign$type=="CTRL_V",]
CTRL_V_data<- CTRL_V_data[order(CTRL_V_data$CVvsV_IZvsCTRL,decreasing=T),]
CTRL_V <- rownames(CTRL_V_data)

IZvsCTRL_CVvsV_snRNA <- SCTransform(IZvsCTRL_CVvsV_snRNA, verbose = FALSE)
IZvsCTRL_CVvsV_snRNA <- ScaleData(IZvsCTRL_CVvsV_snRNA)
library(DoMultiBarHeatmap)
# DefaultAssay(IZvsCTRL_CVvsV_snRNA)<- "RNA"
# IZvsCTRL_CVvsV_snRNA <- ScaleData(IZvsCTRL_CVvsV_snRNA)

IZvsCTRL_CVvsV_snRNA$CAV <- as.character(IZvsCTRL_CVvsV_snRNA$CAV)
IZvsCTRL_CVvsV_snRNA$major_labl <- as.character(IZvsCTRL_CVvsV_snRNA$major_labl)

IZvsCTRL_CVvsV_snRNA$CAV <- factor(IZvsCTRL_CVvsV_snRNA$CAV,levels=c("C_V","Venous"))
IZvsCTRL_CVvsV_snRNA$major_labl <- factor(IZvsCTRL_CVvsV_snRNA$major_labl,levels=c("IZ","CTRL"))

if (!require(devtools)) {
  install.packages("devtools")
}
# devtools::install_github("xmc811/Scillus", ref = "development")
library(Scillus)
DefaultAssay(IZvsCTRL_CVvsV_snRNA) <- "RNA"
IZvsCTRL_CVvsV_snRNA <- SCTransform(IZvsCTRL_CVvsV_snRNA, verbose = FALSE)
IZvsCTRL_CVvsV_snRNA <- ScaleData(IZvsCTRL_CVvsV_snRNA,features=rownames(IZvsCTRL_CVvsV_snRNA))

pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/IZvsCTRL_CVvsV_DEG_heatmap.pdf",width=10,height=10)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snRNA, 
              markers = IZ_CV,
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$blueYellow)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snRNA, 
              markers = IZ_V[1:200],
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$blueYellow)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snRNA, 
              markers = CTRL_CV[1:200],
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$blueYellow)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snRNA, 
              markers = CTRL_V,
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$blueYellow)
dev.off()

# cols.use <- list(CAV=c(C_V='#DAEBA7', Venous='#60CAF8'),major_labl=c(IZ='#BE017C', CTRL='#00A741'))
# pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/IZvsCTRL_CVvsV_DEG_heatmap.pdf",width=10,height=10)
# DoMultiBarHeatmap(IZvsCTRL_CVvsV_snRNA,disp.max=1,disp.min=-1,label=TRUE,cols.use= cols.use,features=IZ_CV, group.by='CAV', additional.group.by = 'major_labl')
# DoMultiBarHeatmap(IZvsCTRL_CVvsV_snRNA, features=IZ_V, group.by='CAV', additional.group.by = 'major_labl')
# DoMultiBarHeatmap(IZvsCTRL_CVvsV_snRNA, features=CTRL_CV, group.by='CAV', additional.group.by = 'major_labl')
# DoMultiBarHeatmap(IZvsCTRL_CVvsV_snRNA, features=CTRL_V, group.by='CAV', additional.group.by = 'major_labl')
# dev.off();

# gene to peak and peak heatmap 
DEG<- c(IZ_CV,IZ_V,CTRL_CV,CTRL_V)
DefaultAssay(IZvsCTRL_CVvsV_snATAC) <- "peaks"
# first compute the GC content for each peak
IZvsCTRL_CVvsV_snATAC <- RegionStats(IZvsCTRL_CVvsV_snATAC, genome = BSgenome.Hsapiens.UCSC.hg38)
IZvsCTRL_CVvsV_snATAC <- LinkPeaks(
  object = IZvsCTRL_CVvsV_snATAC,
  peak.assay = "peaks",
  expression.assay = "RNA",pvalue_cutoff = 1,min.cells = 5,
  genes.use = DEG,score_cutoff = 0
)
lnk <- Links(object = IZvsCTRL_CVvsV_snATAC[["peaks"]])
highest_score_link <- as.data.frame(lnk) %>%
  group_by(gene) %>%
  filter(score == max(score))
highest_score_link<- as.data.frame(highest_score_link)

rownames(highest_score_link)<- highest_score_link$gene;
IZ_CV_highest_score_peak<- na.omit(highest_score_link[IZ_CV,]$peak)
IZ_V_highest_score_peak<- na.omit(highest_score_link[IZ_V,]$peak)
CTRL_CV_highest_score_peak<- na.omit(highest_score_link[CTRL_CV,]$peak)
CTRL_V_highest_score_peak<- na.omit(highest_score_link[CTRL_V,]$peak)

DefaultAssay(IZvsCTRL_CVvsV_snATAC) <- "peaks"

IZvsCTRL_CVvsV_snATAC$CAV <- as.character(IZvsCTRL_CVvsV_snATAC$CAV)
IZvsCTRL_CVvsV_snATAC$major_labl <- as.character(IZvsCTRL_CVvsV_snATAC$major_labl)
IZvsCTRL_CVvsV_snATAC$CAV <- factor(IZvsCTRL_CVvsV_snATAC$CAV,levels=c("C_V","Venous"))
IZvsCTRL_CVvsV_snATAC$major_labl <- factor(IZvsCTRL_CVvsV_snATAC$major_labl,levels=c("IZ","CTRL"))

#IZvsCTRL_CVvsV_snATAC <- SCTransform(IZvsCTRL_CVvsV_snATAC, verbose = FALSE)
IZvsCTRL_CVvsV_snATAC <- ScaleData(IZvsCTRL_CVvsV_snATAC,features=rownames(IZvsCTRL_CVvsV_snATAC))

pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/IZvsCTRL_CVvsV_DEG_linkedPeaks_heatmap.pdf",width=10,height=10)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snATAC, 
              markers = IZ_CV_highest_score_peak,
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$solarExtra)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snATAC, 
              markers = IZ_V_highest_score_peak,
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$solarExtra)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snATAC, 
              markers = CTRL_CV_highest_score_peak[1:200],
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$solarExtra)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snATAC, 
              markers = CTRL_V_highest_score_peak,
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$solarExtra)
dev.off()

# all_peak<- GRangesToString(granges(IZvsCTRL_CVvsV_snATAC), sep = c("-", "-"))
# closest_genes_peak <- ClosestFeature(IZvsCTRL_CVvsV_snATAC, regions = all_peak)
# IZ_CV_peak<- closest_genes_peak[closest_genes_peak$gene_name%in%IZ_CV, ]
# IZ_V_peak<- closest_genes_peak[closest_genes_peak$gene_name%in%IZ_V, ]
# CTRL_CV_peak<- closest_genes_peak[closest_genes_peak$gene_name%in%CTRL_CV, ]
# CTRL_V_peak<- closest_genes_peak[closest_genes_peak$gene_name%in%CTRL_V, ]


## snATAC:
# The DEP between IZ and CTRL;
Idents(IZvsCTRL_CVvsV_snATAC)<- IZvsCTRL_CVvsV_snATAC$major_labl
IZvsCTRL_snATAC_markers <- FindMarkers(IZvsCTRL_CVvsV_snATAC,ident.1="IZ",ident.2="CTRL", min.pct = 0.05, test.use = 'wilcox',logfc.threshold = 0)
log2FC = 0.1
padj = 0.05 
markers<- IZvsCTRL_snATAC_markers
markers$threshold="ns";
markers[which(markers$avg_log2FC  > log2FC & markers$p_val_adj <padj),]$threshold="up";
markers[which(markers$avg_log2FC  < (-log2FC) & markers$p_val_adj < padj),]$threshold="down";
markers$threshold=factor(markers$threshold, levels=c('down','up'))

pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/IZvsCTRL_DEP_volcano.pdf",width=5,height=5)
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +
  xlim(-12,12)+ylim(0,10)+
  ggtitle("C_V vs Venous")+
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  #labs(title= ifelse(""==title, "", paste("DEP:", title)))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)) )+
  theme_classic(base_size = 14) +
  scale_color_manual('',labels=c(paste0("down(",table(markers$threshold)[[1]],')'),
                                 paste0("up(",table(markers$threshold)[[2]],')' )),
                     values=c("blue","red" ) )+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))
dev.off();

# The DEP between CV and V;
Idents(IZvsCTRL_CVvsV_snATAC)<- IZvsCTRL_CVvsV_snATAC$CAV
CVvsV_snATAC_markers <- FindMarkers(IZvsCTRL_CVvsV_snATAC,ident.1="C_V",ident.2="Venous", min.pct = 0.05, test.use = 'wilcox', logfc.threshold = 0)
log2FC = 0.1
padj = 0.05 
markers<- CVvsV_snATAC_markers
markers$threshold="ns";
markers[which(markers$avg_log2FC  > log2FC & markers$p_val_adj <padj),]$threshold="up";
markers[which(markers$avg_log2FC  < (-log2FC) & markers$p_val_adj < padj),]$threshold="down";
markers$threshold=factor(markers$threshold, levels=c('down','up'))

pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/CVvsV_DEP_volcano.pdf",width=5,height=5)
ggplot(data=markers, aes(x=avg_log2FC, y=-log10(p_val_adj), color=threshold)) +
  geom_point(alpha=0.8, size=0.8) +ggtitle("C_V vs Venous")+
  #xlim(-0.5,0.5)+ylim(0,100)+
  geom_vline(xintercept = c(-log2FC, log2FC), linetype=2, color="grey")+
  geom_hline(yintercept = -log10(padj), linetype=2, color="grey")+
  #labs(title= ifelse(""==title, "", paste("DEP:", title)))+
  xlab(bquote(Log[2]*FoldChange))+
  ylab(bquote(-Log[10]*italic(P.adj)) )+
  theme_classic(base_size = 14) +
  scale_color_manual('',labels=c(paste0("down(",table(markers$threshold)[[1]],')'),
                                 paste0("up(",table(markers$threshold)[[2]],')' )),
                     values=c("blue","red" ) )+
  guides(color=guide_legend(override.aes = list(size=3, alpha=1)))
dev.off();

# quad plot 
library(tidyverse)
library(Seurat)
library(ggrepel)
# parpare input data
IZvsCTRL_peak <- rownames(IZvsCTRL_snATAC_markers)
CVvsV_peak <- rownames(CVvsV_snATAC_markers)
peakset <- c(unique(intersect(IZvsCTRL_peak, CVvsV_peak)))

IZvsCTRL_snATAC_markers <- IZvsCTRL_snATAC_markers[rownames(IZvsCTRL_snATAC_markers) %in% peakset, ]
IZvsCTRL_snATAC_markers$peak_name <- rownames(IZvsCTRL_snATAC_markers)
colnames(IZvsCTRL_snATAC_markers)[1:5] <- c("p_val_IZvsCTRL", "avg_log2FC_IZvsCTRL", "pct.1_IZvsCTRL", "pct.2_IZvsCTRL", "p_val_adj_IZvsCTRL")
IZvsCTRL_snATAC_markers <- IZvsCTRL_snATAC_markers %>%arrange(peak_name)
CVvsV_snATAC_markers <- CVvsV_snATAC_markers[rownames(CVvsV_snATAC_markers) %in% peakset, ]
CVvsV_snATAC_markers$peak_name <- rownames(CVvsV_snATAC_markers)
colnames(CVvsV_snATAC_markers)[1:5] <- c("p_val_CVvsV", "avg_log2FC_CVvsV", "pct.1_CVvsV", "pct.2_CVvsV", "p_val_adj_CVvsV")
CVvsV_snATAC_markers <- CVvsV_snATAC_markers %>%arrange(peak_name)
quad_df <- data.frame(cbind(IZvsCTRL_snATAC_markers, CVvsV_snATAC_markers))
head(quad_df)
quad_df_sign<- quad_df[quad_df$p_val_IZvsCTRL<0.05&quad_df$p_val_CVvsV<0.05,]
quad_df_sign$type="nosig";
log2FC=0.1
quad_df_sign[which(quad_df_sign$avg_log2FC_IZvsCTRL  > log2FC & quad_df_sign$avg_log2FC_CVvsV >log2FC),]$type="IZ_CV";
quad_df_sign[which(quad_df_sign$avg_log2FC_IZvsCTRL  > log2FC & quad_df_sign$avg_log2FC_CVvsV < -log2FC),]$type="IZ_V";
quad_df_sign[which(quad_df_sign$avg_log2FC_IZvsCTRL  < -log2FC & quad_df_sign$avg_log2FC_CVvsV >log2FC),]$type="CTRL_CV";
quad_df_sign[which(quad_df_sign$avg_log2FC_IZvsCTRL  < -log2FC & quad_df_sign$avg_log2FC_CVvsV < -log2FC),]$type="CTRL_V";
quad_df_sign$type<- factor(quad_df_sign$type,levels=c("IZ_CV","IZ_V","CTRL_CV","CTRL_V","nosig"))
quad_plot <- quad_df_sign %>%
  ggplot(aes(x = avg_log2FC_CVvsV, y = avg_log2FC_IZvsCTRL,color=type)) + #, label = label
  #annotate("rect", xmin = Inf, xmax = 0, ymin = Inf, ymax = 0, fill= "#457B9D", alpha = 0.65) +
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = -Inf, ymax = 0 , fill= "#A8DADC", alpha = 0.65) +
  #annotate("rect", xmin = Inf, xmax = 0, ymin = -Inf, ymax = 0, fill= "#F1FAEE", alpha = 0.65) +
  #annotate("rect", xmin = -Inf, xmax = 0, ymin = Inf, ymax = 0 , fill= "#E63946", alpha = 0.65) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1.25, color = "black", alpha = 0.85) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 1.25, color = "black", alpha = 0.85) +
  geom_point(size = 1, alpha = 0.8) +
  scale_color_manual('',labels=c(paste0("IZ_CV(",table(quad_df_sign$type)[[1]],')'),
                                 paste0("IZ_V(",table(quad_df_sign$type)[[2]],')' ),
                                 paste0("CTRL_CV(",table(quad_df_sign$type)[[3]],')'),
                                 paste0("CTRL_V(",table(quad_df_sign$type)[[4]],')' ),
                                 paste0("not_sig(",table(quad_df_sign$type)[[5]],')' )),
                     values=c("#88B04B","#92A8D1","#6B5B95","#D64161","grey") )+
  scale_x_continuous(limits = c(-0.8, 0.8)) +
  scale_y_continuous(limits = c(-1.5, 1.5)) +
  xlab("Log2FC(C_V / Venous)") +
  ylab("Log2FC(IZ / CTRL)") +
  ggtitle("snATAC:C_V vs Venous & IZ vs CTRL (avg_log2FC>0.1)")+
  theme_linedraw() +
  theme(axis.text = element_text(size = 15), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(linetype = "solid", colour = "black", size = 1),
        axis.title = element_text(size = 15),
        legend.text = element_text(size = 15),
        aspect.ratio = 5/5)
pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/IZvsCTRL_CVvsV_DEP_quadplot.pdf",width=7,height=7)
quad_plot
dev.off()

IZ_CV <- rownames(quad_df_sign[quad_df_sign$type=="IZ_CV",])
IZ_V <- rownames(quad_df_sign[quad_df_sign$type=="IZ_V",])
CTRL_CV <- rownames(quad_df_sign[quad_df_sign$type=="CTRL_CV",])
CTRL_V <- rownames(quad_df_sign[quad_df_sign$type=="CTRL_V",])

quad_df_sign$CVvsV_IZvsCTRL<- quad_df_sign$avg_log2FC_IZvsCTRL*quad_df_sign$avg_log2FC_CVvsV
IZ_CV_data <- quad_df_sign[quad_df_sign$type=="IZ_CV",]
IZ_CV_data<- IZ_CV_data[order(IZ_CV_data$CVvsV_IZvsCTRL,decreasing=T),]
IZ_CV_peak <- rownames(IZ_CV_data)
IZ_V_data <- quad_df_sign[quad_df_sign$type=="IZ_V",]
IZ_V_data<- IZ_V_data[order(IZ_V_data$CVvsV_IZvsCTRL,decreasing=F),]
IZ_V_peak <- rownames(IZ_V_data)
CTRL_CV_data <- quad_df_sign[quad_df_sign$type=="CTRL_CV",]
CTRL_CV_data<- CTRL_CV_data[order(CTRL_CV_data$CVvsV_IZvsCTRL,decreasing=F),]
CTRL_CV_peak <- rownames(CTRL_CV_data)
CTRL_V_data <- quad_df_sign[quad_df_sign$type=="CTRL_V",]
CTRL_V_data<- CTRL_V_data[order(CTRL_V_data$CVvsV_IZvsCTRL,decreasing=T),]
CTRL_V_peak <- rownames(CTRL_V_data)

DefaultAssay(IZvsCTRL_CVvsV_snATAC) <- "peaks"
IZvsCTRL_CVvsV_snATAC <- ScaleData(IZvsCTRL_CVvsV_snATAC,features=rownames(IZvsCTRL_CVvsV_snATAC))

pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/IZvsCTRL_CVvsV_DEP_heatmap.pdf",width=10,height=10)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snATAC, 
              markers = IZ_CV_peak[1:200],
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$solarExtra)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snATAC, 
              markers = IZ_V_peak,
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$solarExtra)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snATAC, 
              markers = CTRL_CV_peak,
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$solarExtra)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snATAC, 
              markers = CTRL_V_peak[1:200],
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$solarExtra)
dev.off()

closest_IZ_CV <- ClosestFeature(IZvsCTRL_CVvsV_snATAC, IZ_CV_peak)
closest_IZ_V <- ClosestFeature(IZvsCTRL_CVvsV_snATAC, IZ_V_peak)
closest_CTRL_CV <- ClosestFeature(IZvsCTRL_CVvsV_snATAC, CTRL_CV_peak)
closest_CTRL_V <- ClosestFeature(IZvsCTRL_CVvsV_snATAC, CTRL_V_peak)
closest_IZ_CV<- closest_IZ_CV$gene_name
closest_IZ_V<- closest_IZ_V$gene_name
closest_CTRL_CV<- closest_CTRL_CV$gene_name
closest_CTRL_V<- closest_CTRL_V$gene_name

DefaultAssay(IZvsCTRL_CVvsV_snRNA) <- "RNA"
IZvsCTRL_CVvsV_snRNA <- SCTransform(IZvsCTRL_CVvsV_snRNA, verbose = FALSE)
IZvsCTRL_CVvsV_snRNA <- ScaleData(IZvsCTRL_CVvsV_snRNA,features=rownames(IZvsCTRL_CVvsV_snRNA))

pdf("./02_EC/CAV/03_DEG_DEP/05_IZvsCTRL_CVvsV/IZvsCTRL_CVvsV_DEP_linkedGene_heatmap.pdf",width=10,height=10)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snRNA, 
              markers = closest_IZ_CV,
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$blueYellow)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snRNA, 
              markers = closest_IZ_V[1:200],
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$blueYellow)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snRNA, 
              markers = closest_CTRL_CV[1:200],
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$blueYellow)
plot_heatmap(dataset = IZvsCTRL_CVvsV_snRNA, 
              markers = closest_CTRL_V,
              sort_var = c("CAV","major_labl"),
              anno_var = c("CAV","major_labl"),
              anno_colors = list(c('#DAEBA7', '#60CAF8'),
                                 c('#BE017C', '#00A741')),
              hm_limit = seq(-1.5, 1.5, by = 0.375),hm_colors =ArchRPalettes$blueYellow)
dev.off()

