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
CAV_snRNA <- ScaleData(CAV_snRNA,features=rownames(CAV_snRNA))

CAV_snRNA$days.after.infarction<- factor(CAV_snRNA$days.after.infarction,levels=c("control","2","3","4","5","6","11","31","40","45","62","88","101","153","166"))
CAV_snRNA$region_dai<- paste(CAV_snRNA$major_labl,CAV_snRNA$days.after.infarction,sep="_dai")
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai=="CTRL_daicontrol")]<-"CTRL"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("BZ_dai2","BZ_dai5","BZ_dai31"))]<-"BZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("RZ_dai2" ,"RZ_dai3" ,"RZ_dai5","RZ_dai31"))]<-"RZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("FZ_dai40","FZ_dai62","FZ_dai88","FZ_dai101","FZ_dai153","FZ_dai166"))]<-"FZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("IZ_dai2"))]<-"IZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("IZ_dai11","IZ_dai4","IZ_dai5","IZ_dai6"))]<-"IZ"
CAV_snRNA$region_dai[which(CAV_snRNA$region_dai%in%c("IZ_dai45"))]<-"IZ"
CAV_snRNA$region_dai<- factor(CAV_snRNA$region_dai,levels=c("CTRL","IZ","BZ","RZ","FZ"))


CAV_snRNA@meta.data$newgroup <- paste0(CAV_snRNA$CAV,"_",CAV_snRNA$region_dai)
Idents(CAV_snRNA) <- 'newgroup'
Idents(CAV_snRNA) <- factor(Idents(CAV_snRNA), levels = c(paste0("Artery","_",levels(CAV_snRNA$region_dai)),
                                                            paste0("C_A","_",levels(CAV_snRNA$region_dai)),
                                                            paste0("Capillary","_",levels(CAV_snRNA$region_dai)),
                                                            paste0("C_V","_",levels(CAV_snRNA$region_dai)),
                                                            paste0("Venous","_",levels(CAV_snRNA$region_dai))))        
# RNA expression part:
# Ki67 + cellcycle gene 
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
gene<- gene[!duplicated(c(s.genes,g2m.genes))]

Idents(CAV_snRNA) <- CAV_snRNA$CAV
pdf("./02_EC/CAV/snRNA_cellcycle_marker_dotplot.pdf",width=20,height=6)
DotPlot(CAV_snRNA, features = gene, cols = c("lightgrey", "red"))
DotPlot(CAV_snRNA, features = gene, cols = c('#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87','#E95C59'),split.by = "region_dai")& theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

markers <- c("MCM6","HELLS","RFC2","SLBP","CASP8AP2","USP1","SMC4","NCAPD2","ECT2","CKAP5","CTCF","GAS2L3","CBX5","")
pdf("./02_EC/CAV/snRNA_cellcycle_marker_dotplot.pdf",width=12,height=10)
DotPlot(CAV_snRNA, features = markers, cols = c("lightgrey", "red"))& theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(CAV_snRNA, features = markers, cols = c("lightgrey", "red"),group.by="newgroup")& theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
selcet_markers <- c("RFC2","SLBP","CASP8AP2","SMC4","CKAP5")


library(DOSE)
library(GOSemSim)
library(clusterProfiler)
library(org.Hs.eg.db)
# library(org.Mm.eg.db)
# library(org.Rn.eg.db)
library(dplyr)
library(GO.db)
get_GO_data <- function(OrgDb, ont, keytype) {
  GO_Env <- get_GO_Env()
  use_cached <- FALSE
  
  if (exists("organism", envir=GO_Env, inherits=FALSE) &&
      exists("keytype", envir=GO_Env, inherits=FALSE)) {
    
    org <- get("organism", envir=GO_Env)
    kt <- get("keytype", envir=GO_Env)
    
    if (org == DOSE:::get_organism(OrgDb) &&
        keytype == kt &&
        exists("goAnno", envir=GO_Env, inherits=FALSE)) {
      ## https://github.com/GuangchuangYu/clusterProfiler/issues/182
      ## && exists("GO2TERM", envir=GO_Env, inherits=FALSE)){
      
      use_cached <- TRUE
    }
  }
  
  if (use_cached) {
    goAnno <- get("goAnno", envir=GO_Env)
  } else {
    OrgDb <- GOSemSim:::load_OrgDb(OrgDb)
    kt <- keytypes(OrgDb)
    if (! keytype %in% kt) {
      stop("keytype is not supported...")
    }
    
    kk <- keys(OrgDb, keytype=keytype)
    goAnno <- suppressMessages(
      AnnotationDbi::select(OrgDb, keys=kk, keytype=keytype,
             columns=c("GOALL", "ONTOLOGYALL")))
    
    goAnno <- unique(goAnno[!is.na(goAnno$GOALL), ])
    
    assign("goAnno", goAnno, envir=GO_Env)
    assign("keytype", keytype, envir=GO_Env)
    assign("organism", DOSE:::get_organism(OrgDb), envir=GO_Env)
  }
  
  if (ont == "ALL") {
    GO2GENE <- unique(goAnno[, c(2,1)])
  } else {
    GO2GENE <- unique(goAnno[goAnno$ONTOLOGYALL == ont, c(2,1)])
  }
  
  GO_DATA <- DOSE:::build_Anno(GO2GENE, get_GO2TERM_table())
  
  goOnt.df <- goAnno[, c("GOALL", "ONTOLOGYALL")] %>% unique
  goOnt <- goOnt.df[,2]
  names(goOnt) <- goOnt.df[,1]
  assign("GO2ONT", goOnt, envir=GO_DATA)
  return(GO_DATA)
}
get_GO_Env <- function () {
  if (!exists(".GO_clusterProfiler_Env", envir = .GlobalEnv)) {
    pos <- 1
    envir <- as.environment(pos)
    assign(".GO_clusterProfiler_Env", new.env(), envir=envir)
  }
  get(".GO_clusterProfiler_Env", envir = .GlobalEnv)
}

get_GO2TERM_table <- function() {
  GOTERM.df <- get_GOTERM()
  GOTERM.df[, c("go_id", "Term")] %>% unique
}

get_GOTERM <- function() {
  pos <- 1
  envir <- as.environment(pos)
  if (!exists(".GOTERM_Env", envir=envir)) {
    assign(".GOTERM_Env", new.env(), envir)
  }
  GOTERM_Env <- get(".GOTERM_Env", envir = envir)
  if (exists("GOTERM.df", envir = GOTERM_Env)) {
    GOTERM.df <- get("GOTERM.df", envir=GOTERM_Env)
  } else {
    GOTERM.df <- toTable(GOTERM)
    assign("GOTERM.df", GOTERM.df, envir = GOTERM_Env)
  }
  return(GOTERM.df)
}

GO_DATA <- get_GO_data("org.Hs.eg.db", "ALL", "SYMBOL")


findGO <- function(pattern, method = "key"){
    if(!exists("GO_DATA"))
        load("GO_DATA.RData")
    if(method == "key"){
        pathways = cbind(GO_DATA$PATHID2NAME[grep(pattern, GO_DATA$PATHID2NAME)])
    } else if(method == "gene"){
        pathways = cbind(GO_DATA$PATHID2NAME[GO_DATA$EXTID2PATHID[[pattern]]])
    }

    colnames(pathways) = "pathway"

    if(length(pathways) == 0){
        cat("No results!\n")
    } else{
        return(pathways)
    }
}

getGO <- function(ID){

    if(!exists("GO_DATA"))
        load("GO_DATA.RData")
    allNAME = names(GO_DATA$PATHID2EXTID)
    if(ID %in% allNAME){
        geneSet = GO_DATA$PATHID2EXTID[ID]
        names(geneSet) = GO_DATA$PATHID2NAME[ID]
        return(geneSet)     
    } else{
        cat("No results!\n")
    }
}

findGO("endothelial cell proliferation")

#load("GO_DATA.RData") # 载入数据 GO_DATA
#findGO("endothelial cell proliferation") # 寻找 含有指定关键字的 pathway name 的 pathway
#findGO("INS", method = "gene") # 寻找含有指定基因名的 pathway
#getGO("GO:0008152") # 获取指定 GO ID 的 gene set
#$`insulin metabolic process`
#[1] "CEACAM1" "CPE"     "ERN1"    "IDE"     "PCSK2"   "ERO1B"
#gene_related_endothelial_cell_proliferation<-unlist(getGO("GO:0001935"))
posi_gene_related_endothelial_cell_proliferation<-unlist(getGO("GO:0001938"))
nega_gene_related_endothelial_cell_proliferation<-unlist(getGO("GO:0001937"))

Idents(CAV_snRNA) <- CAV_snRNA$CAV
pdf("./02_EC/CAV/snRNA_endothelial_cell_proliferation_dotplot.pdf",width=30,height=4)
DotPlot(CAV_snRNA, features = posi_gene_related_endothelial_cell_proliferation, cols = c("lightgrey", "red"))& theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(CAV_snRNA, features = nega_gene_related_endothelial_cell_proliferation, cols = c("lightgrey", "red"))& theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

posi_markers <- c("BMP6","HIF1A","NF1","DYSF","NRP2","NRP1")
pdf("./02_EC/CAV/snRNA_posi_endothelial_cell_proliferation_dotplot.pdf",width=12,height=10)
DotPlot(CAV_snRNA, features = posi_markers, cols = c("lightgrey", "red"))& theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(CAV_snRNA, features = posi_markers, cols = c("lightgrey", "red"),group.by="newgroup")& theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()


# GO:migration
findGO("endothelial cell migration")

# #GO:0003348 "cardiac endothelial cell differentiation"
# # GO:0045603 "positive regulation of endothelial cell differentiation"
# 
# cardiac_EC_differentiation<-unlist(getGO("GO:0003348"))
# EC_differentiation<-unlist(getGO("GO:0045603"))
# 
# pdf("./02_EC/CAV/snRNA_EC_differentiation_dotplot.pdf",width=12,height=10)
# DotPlot(CAV_snRNA, features = cardiac_EC_differentiation, cols = c("lightgrey", "red"))& theme(axis.text.x = element_text(angle = 90, hjust = 1))
# DotPlot(CAV_snRNA, features = EC_differentiation, cols = c("lightgrey", "red"),group.by="newgroup")& theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()

#GO:1904989 "positive regulation of endothelial cell activation"
positive_EC_activation<- unlist(getGO("GO:1904989"))

#GO:0010595 "positive regulation of endothelial cell migration"
positive_EC_migration<-unlist(getGO("GO:0010595"))

pdf("./02_EC/CAV/snRNA_EC_migration_dotplot.pdf",width=30,height=8)
DotPlot(CAV_snRNA, features = positive_EC_migration, cols = c("lightgrey", "red"))& theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(CAV_snRNA, features = positive_EC_migration, cols = c("lightgrey", "red"),group.by="newgroup")& theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
select_positive_EC_migration<-c("ABL1","ETS1","PIK3C2A","SASH1","FOXP1","HDAC7","RHOJ")
pdf("./02_EC/CAV/snRNA_EC_migration_dotplot.pdf",width=12,height=10)
DotPlot(CAV_snRNA, features = select_positive_EC_migration, cols = c("lightgrey", "red"))& theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(CAV_snRNA, features = select_positive_EC_migration, cols = c("lightgrey", "red"),group.by="newgroup")& theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
# # GO:0060842 "arterial endothelial cell differentiation"
# # GO:0060844 "arterial endothelial cell fate commitment"
# # GO:0060853 "Notch signaling pathway involved in arterial endothelial cell fate commitment"
# arterial<- unlist(getGO("GO:0060842"))
# # GO:0060843 "venous endothelial cell differentiation"
# venous<- unlist(getGO("GO:0060843"))
# # GO:0072577 "endothelial cell apoptotic process"
# # GO:2000353 "positive regulation of endothelial cell apoptotic process"
# #apoptotic<- unlist(getGO("GO:0072577"))
# apoptotic<- unlist(getGO("GO:2000353"))
# pdf("./02_EC/CAV/snRNA_EC_arterial_venous_apoptotic_dotplot.pdf",width=30,height=8)
# DotPlot(CAV_snRNA, features = c(arterial,apoptotic), cols = c("lightgrey", "red"))& theme(axis.text.x = element_text(angle = 90, hjust = 1))
# DotPlot(CAV_snRNA, features = c(arterial,apoptotic), cols = c("lightgrey", "red"),group.by="newgroup")& theme(axis.text.x = element_text(angle = 90, hjust = 1))
# dev.off()

last<- c("NR2F2",posi_markers,select_positive_EC_migration)

pdf("./02_EC/CAV/snRNA_cellcycle_proli_migration_dotplot.pdf",width=8,height=10)
DotPlot(CAV_snRNA, features = last, cols = c("lightgrey", "red"))& theme(axis.text.x = element_text(angle = 90, hjust = 1))
DotPlot(CAV_snRNA, features = last, cols = c("lightgrey", "red"),group.by="newgroup")& theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()



CAV_snRNA$CAV<- factor(CAV_snRNA$CAV,levels=rev(levels(CAV_snRNA$CAV)))
#CAV_snRNA$CAV<- factor(CAV_snRNA$CAV,levels=rev(levels(CAV_snRNA$CAV)))

Idents(CAV_snRNA)<- CAV_snRNA$CAV
set<- myUmapcolors[c(5,4,3,2,1)]

pdf("./02_EC/CAV/snRNA_cellcycle_proli_migration_violin.pdf",width=20,height=40)
VlnPlot(object = CAV_snRNA, features = last, add.noise =F,split.by = 'region_dai',log=TRUE,stack =TRUE)+geom_boxplot(width=0.03,fill="white",outlier.size=0)+ scale_x_continuous(limits = c(1, 5))
VlnPlot(object = CAV_snRNA, features = last, add.noise =TRUE,split.by = 'region_dai',log=TRUE,ncol=1)+scale_fill_manual(values=set)+geom_boxplot(width=0.03,fill="white",outlier.size=0)+ scale_y_continuous(limits = c(1, 5))
dev.off()


Idents(CAV_snRNA)<- CAV_snRNA$CAV
pdf("./02_EC/CAV/snRNA_proliferative_marker_violin.pdf",width=12,height=20)
VlnPlot(object = CAV_snRNA, features = features, split.by = 'region_dai',ncol=1,log=TRUE)
VlnPlot(object = CAV_snRNA, features = features, split.by = 'region_dai',ncol=1)
VlnPlot(object = CAV_snRNA, features = features,ncol=1)
dev.off()


# ATAC accessibility part:

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
  genes.use = last
)
######Visulize track and RNA exp######
Idents(CAV_snATAC)<- factor(Idents(CAV_snATAC),levels=c("Venous","C_V","Capillary","C_A","Artery"))
idents.plot <-Idents(CAV_snATAC)

pdf("./02_EC/CAV/snATAC_proli_migration-peaktrack-ACTIVITY.pdf",height=6,width=8)
for(i in last[which(last%in% gene)]){
  print(i)
  p1 <- CoveragePlot(
  object = CAV_snATAC,
  region = i,
  features = i,
  expression.assay = "ACTIVITY",
  idents = idents.plot,
  extend.upstream = 3000,
  annotation=TRUE,
  extend.downstream = -3000
)
print(p1)}
dev.off()

DefaultAssay(CAV_snATAC) <- "RNA_imputation"
gene<- rownames(CAV_snATAC)

DefaultAssay(CAV_snATAC) <- "peaks"
# first compute the GC content for each peak
CAV_snATAC <- RegionStats(CAV_snATAC, genome = BSgenome.Hsapiens.UCSC.hg38)

# link peaks to genes
CAV_snATAC <- LinkPeaks(
  object = CAV_snATAC,
  peak.assay = "peaks",
  expression.assay = "RNA_imputation",
  genes.use = last
)
######Visulize track and RNA exp######
idents.plot <- Idents(CAV_snATAC)

pdf("./02_EC/CAV/snATAC_proli_migration-peaktrack-RNA_imputation.pdf",height=6,width=8)
for(i in last[which(last%in% gene)]){
  print(i)
  p1 <- CoveragePlot(
  object = CAV_snATAC,
  region = i,
  features = i,
  expression.assay = "RNA_imputation",
  idents = idents.plot,
  extend.upstream = 1000,
  annotation=TRUE,
  extend.downstream = 1000
)
print(p1)}
dev.off()

p1 <- CoveragePlot(
  object = CAV_snATAC,
  region = "chr15-96326500-96328000",
  #features = i,
  #expression.assay = "RNA_imputation",
  idents = idents.plot,
  extend.upstream = 100,
  annotation=TRUE,
  extend.downstream = 100
)
p2 <- CoveragePlot(
  object = CAV_snATAC,
  region = "chr2-205675000-205700000",
  #features = i,
  #expression.assay = "RNA_imputation",
  idents = idents.plot,
  extend.upstream = 100,
  annotation=TRUE,
  extend.downstream = 100
)
p3 <- CoveragePlot(
  object = CAV_snATAC,
  region = "chr12-47805000-47815000",
  #features = i,
  #expression.assay = "RNA_imputation",
  idents = idents.plot,
  extend.upstream = 100,
  annotation=TRUE,
  extend.downstream = 100
)
p4 <- CoveragePlot(
  object = CAV_snATAC,
  region = "chr14-63225000-63235000",
  #features = i,
  #expression.assay = "RNA_imputation",
  idents = idents.plot,
  extend.upstream = 100,
  annotation=TRUE,
  extend.downstream = 100
)
set<- myUmapcolors[c(5,4,3,2,1)]
p1<-p1& scale_fill_manual(values=set)
p2<-p2& scale_fill_manual(values=set)& theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p3<-p3& scale_fill_manual(values=set)& theme(strip.text.y.left = element_blank(),strip.background = element_blank())
p4<-p4& scale_fill_manual(values=set)& theme(strip.text.y.left = element_blank(),strip.background = element_blank())
pdf("./02_EC/CAV/snATAC_show_migration-peaktrack-RNA_imputation.pdf",height=6,width=6)

p1|p2|p3|p4
dev.off()









