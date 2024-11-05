conda activate r43
cd /md01/nieyg/project/EC/human_MI
R

.libPaths("/data/R02/nieyg/ori/biosoft/conda/envs/r43/lib/R/library")

library(Seurat)
library(Signac)
library(SeuratDisk)

snRNA <- readRDS("/md01/nieyg/project/EC/human_MI/00_data/All-snRNA-seq.rds")
snATAC<- readRDS("/md01/nieyg/project/EC/human_MI/00_data/All-snATAC-seq.rds")
# add the ATAC matrix for snATAC obj
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK168.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK171.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK173.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK336.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK339.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK340.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK342.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK344.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK345.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK347.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK348.tsv.gz
wget https://zenodo.org/records/6578553/files/10X_ATAC_CK380.tsv.gz
wget https://zenodo.org/records/6578617/files/10X_ATAC_CK385.tsv.gz

wget https://zenodo.org/records/6578047/files/Visium_control_P1.h5ad
wget https://zenodo.org/records/6578047/files/Visium_control_P17.h5ad
wget https://zenodo.org/records/6578047/files/Visium_control_P7.h5ad
wget https://zenodo.org/records/6578047/files/Visium_control_P8.h5ad
wget https://zenodo.org/records/6578047/files/Visium_FZ_GT_P19.h5ad
wget https://zenodo.org/records/6578047/files/Visium_FZ_GT_P4.h5ad
wget https://zenodo.org/records/6578047/files/Visium_FZ_P14.h5ad
wget https://zenodo.org/records/6578047/files/Visium_FZ_P18.h5ad
wget https://zenodo.org/records/6578047/files/Visium_FZ_P20.h5ad
wget https://zenodo.org/records/6578047/files/Visium_GT_IZ_P13.h5ad
wget https://zenodo.org/records/6578047/files/Visium_GT_IZ_P15.h5ad
wget https://zenodo.org/records/6578047/files/Visium_GT_IZ_P9.h5ad
wget https://zenodo.org/records/6578047/files/Visium_GT_IZ_P9_rep2.h5ad
wget https://zenodo.org/records/6578047/files/Visium_IZ_BZ_P2.h5ad
wget https://zenodo.org/records/6578047/files/Visium_IZ_P10.h5ad
wget https://zenodo.org/records/6578047/files/Visium_IZ_P15.h5ad
wget https://zenodo.org/records/6578047/files/Visium_IZ_P16.h5ad
wget https://zenodo.org/records/6578047/files/Visium_IZ_P3.h5ad
wget https://zenodo.org/records/6578047/files/Visium_RZ_BZ_P12.h5ad
wget https://zenodo.org/records/6578047/files/Visium_RZ_BZ_P2.h5ad
wget https://zenodo.org/records/6578047/files/Visium_RZ_BZ_P3.h5ad
wget https://zenodo.org/records/6578047/files/Visium_RZ_FZ_P5.h5ad
wget https://zenodo.org/records/6578047/files/Visium_RZ_GT_P2.h5ad
wget https://zenodo.org/records/6578047/files/Visium_RZ_P11.h5ad
wget https://zenodo.org/records/6578047/files/Visium_RZ_P3.h5ad
wget https://zenodo.org/records/6578047/files/Visium_RZ_P6.h5ad
wget https://zenodo.org/records/6578047/files/Visium_RZ_P9.h5ad


# visium rawdata 

wget https://zenodo.org/records/6580069/files/10X_Visium_10X001.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_10X0017.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_10X0018.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_10X0020.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_10X0025.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_10X0026.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_10X0027.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_10X009.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0010.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0011.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0012.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0013.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0014.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0015.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0016.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0019.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH002.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0021.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0022.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0023.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0024.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH0028.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH003.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH004.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH005.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH006.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH007.tar.gz
wget https://zenodo.org/records/6580069/files/10X_Visium_ACH008.tar.gz
wget https://zenodo.org/records/6580069/files/metadata-Visium.csv

# snRNA rawdata 
wget https://zenodo.org/records/6578617/files/10X_RNA_CK158.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK159.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK160.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK161.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK162.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK163.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK164.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK165.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK356.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK357.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK358.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK359.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK360.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK361.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK362.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK363.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK364.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK365.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK366.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK367.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK368.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK369.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK370.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK371.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK372.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK373.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK374.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK375.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_CK376.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_KL001.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_KL002.tar.gz
wget https://zenodo.org/records/6578617/files/10X_RNA_KL003.tar.gz

gunzip -d *

sed -i '/#/d;' *

awk '{$4 = "CK166#" $4; print}' 10X_ATAC_CK166.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK166.tsv

 awk '{$4 = "CK167#" $4; print}' 10X_ATAC_CK167.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK167.tsv
 awk '{$4 = "CK168#" $4; print}' 10X_ATAC_CK168.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK168.tsv
 awk '{$4 = "CK169#" $4; print}' 10X_ATAC_CK169.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK169.tsv
 awk '{$4 = "CK170#" $4; print}' 10X_ATAC_CK170.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK170.tsv
 awk '{$4 = "CK171#" $4; print}' 10X_ATAC_CK171.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK171.tsv
 awk '{$4 = "CK173#" $4; print}' 10X_ATAC_CK173.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK173.tsv
 awk '{$4 = "CK174#" $4; print}' 10X_ATAC_CK174.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK174.tsv
 awk '{$4 = "CK336#" $4; print}' 10X_ATAC_CK336.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK336.tsv
 awk '{$4 = "CK337#" $4; print}' 10X_ATAC_CK337.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK337.tsv
 awk '{$4 = "CK338#" $4; print}' 10X_ATAC_CK338.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK338.tsv
 awk '{$4 = "CK339#" $4; print}' 10X_ATAC_CK339.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK339.tsv
 awk '{$4 = "CK340#" $4; print}' 10X_ATAC_CK340.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK340.tsv
 awk '{$4 = "CK341#" $4; print}' 10X_ATAC_CK341.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK341.tsv
 awk '{$4 = "CK342#" $4; print}' 10X_ATAC_CK342.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK342.tsv
 awk '{$4 = "CK343#" $4; print}' 10X_ATAC_CK343.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK343.tsv
 awk '{$4 = "CK344#" $4; print}' 10X_ATAC_CK344.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK344.tsv
 awk '{$4 = "CK345#" $4; print}' 10X_ATAC_CK345.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK345.tsv
 awk '{$4 = "CK346#" $4; print}' 10X_ATAC_CK346.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK346.tsv
 awk '{$4 = "CK347#" $4; print}' 10X_ATAC_CK347.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK347.tsv
 awk '{$4 = "CK348#" $4; print}' 10X_ATAC_CK348.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK348.tsv
 awk '{$4 = "CK349#" $4; print}' 10X_ATAC_CK349.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK349.tsv
 awk '{$4 = "CK350#" $4; print}' 10X_ATAC_CK350.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK350.tsv
 awk '{$4 = "CK351#" $4; print}' 10X_ATAC_CK351.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK351.tsv
 awk '{$4 = "CK352#" $4; print}' 10X_ATAC_CK352.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK352.tsv
 awk '{$4 = "CK353#" $4; print}' 10X_ATAC_CK353.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK353.tsv
 awk '{$4 = "CK354#" $4; print}' 10X_ATAC_CK354.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK354.tsv
 awk '{$4 = "CK355#" $4; print}' 10X_ATAC_CK355.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK355.tsv
 awk '{$4 = "CK380#" $4; print}' 10X_ATAC_CK380.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK380.tsv
 awk '{$4 = "CK381#" $4; print}' 10X_ATAC_CK381.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK381.tsv
 awk '{$4 = "CK382#" $4; print}' 10X_ATAC_CK382.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK382.tsv
 awk '{$4 = "CK383#" $4; print}' 10X_ATAC_CK383.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK383.tsv
 awk '{$4 = "CK385#" $4; print}' 10X_ATAC_CK385.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK385.tsv
 awk '{$4 = "CK386#" $4; print}' 10X_ATAC_CK386.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK386.tsv
 awk '{$4 = "CK387#" $4; print}' 10X_ATAC_CK387.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK387.tsv
 awk '{$4 = "CK388#" $4; print}' 10X_ATAC_CK388.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK388.tsv
 awk '{$4 = "CK389#" $4; print}' 10X_ATAC_CK389.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK389.tsv
 awk '{$4 = "CK390#" $4; print}' 10X_ATAC_CK390.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK390.tsv
 awk '{$4 = "CK391#" $4; print}' 10X_ATAC_CK391.tsv > tmp.txt && mv tmp.txt 10X_ATAC_CK391.tsv

cat *tsv > fragments.tsv
awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$5}' fragments.tsv >fragments_bed.tsv
cat fragments.tsv |sort -k1,1 -k2,2n | bgzip > fragments.tsv.gz 
tabix -s 1 -b 2 -e 3 -p bed fragments.tsv.gz

awk '{print $1,$2,$3}' fragments.tsv | sort | uniq > peaks.bed

frags.snATAC <- CreateFragmentObject(
  path = "/md01/nieyg/project/EC/human_MI/00_data/fragments/fragments.tsv.gz",
  cells = colnames(snATAC)
)

peaks <- read.table(
  file = "/md01/nieyg/project/EC/human_MI/00_data/fragments/peaks.bed",
  col.names = c("chr", "start", "end")
)
# convert to genomic ranges
gr <- makeGRangesFromDataFrame(peaks)

# Create a unified set of peaks to quantify in each dataset
combined.peaks <- reduce(x = gr)

# Filter out bad peaks based on length
peakwidths <- width(combined.peaks)
combined.peaks <- combined.peaks[peakwidths  < 10000 & peakwidths > 20]
combined.peaks

snATAC.counts <- FeatureMatrix(
  fragments = frags.snATAC,
  features = combined.peaks,
  cells = colnames(snATAC)
)
snATAC_assay <- CreateChromatinAssay(snATAC.counts, fragments = frags.snATAC)
snATAC[["ATAC"]] <- snATAC_assay

# callpeaks for snATAC 

# call peak 

DefaultAssay(snATAC)<-"ATAC"
peak<-CallPeaks(
       snATAC,
       group.by = "cell_type",
       macs2.path = "/public/home/nieyg/biosoft/conda/bin/macs2",
       broad = FALSE,
       format = "BED",
       fragment.tempdir = tempdir(),
       effective.genome.size = 2.7e+09,
       outdir="/data/R02/nieyg/project/EC/human_MI/01_all_celltype/peaks",
       combine.peaks=TRUE
)
macs2_counts <- FeatureMatrix(
     fragments = Fragments(snATAC),
     features = peak,
     cells = colnames(snATAC)
     )     
#macs2_counts<-macs2_counts[-which(rownames(macs2_counts)=="GroupUN243-311766-311965"),]
# ATAC analysis add gene annotation information
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
genome(annotations) <- "hg38"
annotations<- renameSeqlevels(annotations, paste0("chr",seqlevels(annotations)))
#seqlevelsStyle(annotations) <- "UCSC"


#Integrating scRNA-seq and scATAC-seq data
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
 snRNA$cell_type<- factor(snRNA$cell_type,levels=c(levels(Idents(snATAC)),"mast cell",
    "adipocyte of epicardial fat of left ventricle","unknown"))
Idents(snRNA)<- snRNA$cell_type
Idents(snATAC)<- snATAC$cell_type

myUmapcolors <- c(  '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87',
         '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658',
         '#9FA3A8', '#E0D4CA', '#5F3D69', '#58A4C3', '#AA9A59', '#E63863', '#E39A35', 
         '#C1E6F3', '#6778AE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#625D9E', 
         '#68A180', '#3A6963', '#968175', '#161853', '#FF9999', '#344CB7', '#FFCC1D', 
         '#116530', '#678983', '#A19882', '#FFBCBC', '#24A19C', '#FF9A76', "#8DD3C7",
         "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
         "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F", "#E41A1C", "#377EB8", "#4DAF4A", 
         "#FF7F00", "#FFFF33", "#A65628", "#F781BF")

p1 <- DimPlot(snRNA,  cols=myUmapcolors,reduction="umap" ,raster=FALSE,label = TRUE) + NoLegend() + ggtitle("RNA")
p2 <- DimPlot(snATAC,  cols=myUmapcolors,group.by = "cell_type", label = TRUE) + NoLegend() + ggtitle("ATAC")
pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/RNA_ATAC_split_Umap.pdf",width=10,height=5)
p1 + p2
plot <- (p1 + p2) & xlab("UMAP 1") & ylab("UMAP 2") & theme(axis.title = element_text(size = 18))
plot
dev.off()

# rename the gene symbol in snRNA
library(org.Hs.eg.db)
head(rownames(snRNA))
ids=select(org.Hs.eg.db,keys = rownames(snRNA),
           columns = c('ENSEMBL','SYMBOL'),
           keytype = 'ENSEMBL')
head(ids)
dim(ids) 
ids=na.omit(ids)
dim(ids) 
length(unique(ids$SYMBOL)) # [1] 15494 
ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
pos=match(ids$ENSEMBL,rownames(snRNA) )
snRNA=snRNA[pos,]
snRNA 
RenameGenesSeurat <- function(obj , 
                              newnames ) { 
print("Run this before integration. It only changes obj@assays$RNA@counts, @data and @scale.data.")
RNA <- obj@assays$RNA
if (nrow(RNA) == length(newnames)) {
    if (length(RNA@counts)) RNA@counts@Dimnames[[1]]            <- newnames
    if (length(RNA@data)) RNA@data@Dimnames[[1]]                <- newnames
    if (length(RNA@scale.data)) RNA@scale.data@Dimnames[[1]]    <- newnames
  } else {"Unequal gene sets: nrow(RNA) != nrow(newnames)"}
  obj@assays$RNA <- RNA
  return(obj)
}
snRNA_tmp=RenameGenesSeurat(obj = snRNA, 
                  newnames = ids$SYMBOL)
snRNA<- snRNA_tmp
# quantify gene activity
rownames(snRNA[["RNA"]]@meta.features) <- snRNA[["RNA"]]@meta.features$feature_name
snRNA <- FindVariableFeatures(snRNA)
#data.frame(row.names = rownames(snRNA[["RNA"]]))

# rename the gene symbol in snATAC
DefaultAssay(snATAC)<-"RNA"
head(rownames(snATAC))
ids=select(org.Hs.eg.db,keys = rownames(snATAC),
           columns = c('ENSEMBL','SYMBOL'),
           keytype = 'ENSEMBL')
head(ids)
dim(ids) 
ids=na.omit(ids)
dim(ids) 
length(unique(ids$SYMBOL)) # [1] 15494 
ids=ids[!duplicated(ids$SYMBOL),]
ids=ids[!duplicated(ids$ENSEMBL),]
pos=match(ids$ENSEMBL,rownames(snATAC) )
snATAC=snATAC[pos,]
snATAC 
snATAC_tmp=RenameGenesSeurat(obj = snATAC, 
                  newnames = ids$SYMBOL)
snATAC<- snATAC_tmp
# quantify gene activity
rownames(snATAC[["RNA"]]@meta.features) <- snATAC[["RNA"]]@meta.features$feature_name
snATAC <- FindVariableFeatures(snATAC)

# create a new assay using the MACS2 peak set and add it to the Seurat object
snATAC[["peaks"]] <- CreateChromatinAssay(
  counts = macs2_counts,
  #fragments = fragpath[[i]],
  fragments = Fragments(snATAC),
  annotation = Annotation(snATAC)
)

DefaultAssay(snATAC)<-"peaks"
gene.activities <- GeneActivity(snATAC,features=rownames(snRNA))
# add gene activities as a new assay
snATAC[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

# normalize gene activities
DefaultAssay(snATAC) <- "ACTIVITY"
snATAC <- snATAC %>% 
        NormalizeData() %>%
        FindVariableFeatures() %>%
        ScaleData() %>%
        RunPCA(verbose = FALSE)

# Identify anchors
transfer.anchors <- FindTransferAnchors(reference = snRNA, 
    query = snATAC,
    features = VariableFeatures(object = snRNA),
    reference.assay = "RNA",
    query.assay = "ACTIVITY", 
    reduction = "cca")

celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = snRNA$cell_type,
    weight.reduction = snRNA[["harmony"]], dims = 1:30)

snATAC <- AddMetaData(snATAC, metadata = celltype.predictions)
snATAC$annotation_correct <- snATAC$predicted.id == snATAC$cell_type
pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/ATAC_predicted_Ground-truth_Umap.pdf",width=10,height=5)
p1 <- DimPlot(snATAC,reduction='umap' , group.by = "predicted.id", label = TRUE) + NoLegend() + ggtitle("Predicted annotation")
p2 <- DimPlot(snATAC,reduction='umap' ,group.by = "cell_type", label = TRUE) + NoLegend() + ggtitle("Ground-truth annotation")
p1 | p2
dev.off()

predictions <- table(snATAC$cell_type, snATAC$predicted.id)
predictions <- predictions/rowSums(predictions)  # normalize for number of cells in each cell type
predictions <- as.data.frame(predictions)

order_Var1<- c(levels(Idents(snATAC)),"mast cell",
    "adipocyte of epicardial fat of left ventricle","unknown")
order_Var2<- c(levels(Idents(snATAC)))
predictions$Var1<- factor(predictions$Var1,levels=order_Var1)
predictions$Var2<- factor(predictions$Var2,levels=order_Var2)

pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/ATAC_predicted_rate.pdf",width=10,height=5)
p1 <- ggplot(predictions, aes(Var1, Var2, fill = Freq)) + geom_tile() + scale_fill_gradient(name = "Fraction of cells",
    low = "#ffffc8", high = "#7d0025") + xlab("Cell type annotation (RNA)") + ylab("Predicted cell type label (ATAC)") +
    theme_cowplot() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

correct <- length(which(snATAC$cell_type == snATAC$predicted.id))
incorrect <- length(which(snATAC$cell_type != snATAC$predicted.id))
data <- FetchData(snATAC, vars = c("prediction.score.max", "annotation_correct"))
p2 <- ggplot(data, aes(prediction.score.max, fill = annotation_correct, colour = annotation_correct)) +
    geom_density(alpha = 0.5) + theme_cowplot() + scale_fill_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + scale_color_discrete(name = "Annotation Correct",
    labels = c(paste0("FALSE (n = ", incorrect, ")"), paste0("TRUE (n = ", correct, ")"))) + xlab("Prediction Score")
p1 + p2
dev.off()

# note that we restrict the imputation to variable genes from scRNA-seq, but could impute the
# full transcriptome if we wanted to
genes.use <- VariableFeatures(snRNA)
refdata <- GetAssayData(snRNA, assay = "RNA", slot = "data")[genes.use, ]

# refdata (input) contains a scRNA-seq expression matrix for the scRNA-seq cells.  imputation
# (output) will contain an imputed scRNA-seq matrix for each of the ATAC cells
imputation <- TransferData(anchorset = transfer.anchors, 
    refdata = refdata, 
    weight.reduction = "cca")
snATAC[["RNA"]] <- imputation

# change the metadata column in snRNA and snATAC to keep same 
overlap_col<- intersect(colnames(snRNA@meta.data),colnames(snATAC@meta.data))
setdiff(colnames(snRNA@meta.data),overlap_col)
setdiff(colnames(snATAC@meta.data),overlap_col)

colnames(snRNA@meta.data)[8]<- "patient_region"
## update the same cell name
snATAC$major_labl<- snATAC$region 
snATAC$major_labl <- stringr::str_replace_all(snATAC$major_labl,
                                              c("control" = "CTRL",
                                                "RZ/BZ" = "BZ",
                                                "RZ/FZ" = "FZ",
                                                "RZ/GT" = "RZ",
                                                "FZ/GT" = "FZ",
                                                "GT/IZ" = "IZ",
                                                "IZ/BZ" = "IZ"))
snRNA$region <- gsub("_P.*","",snRNA$patient_region)
snRNA$region <- factor(snRNA$region,levels=levels(snATAC$region))
snATAC$major_labl <- factor(snATAC$major_labl,levels=levels(snRNA$major_labl))
snRNA$tech <- "RNA"
snATAC$tech <- "ATAC"
overlap_col<- intersect(colnames(snRNA@meta.data),colnames(snATAC@meta.data))
snATAC@meta.data<- snATAC@meta.data[,overlap_col]
snRNA@meta.data<- snRNA@meta.data[,overlap_col]

coembed <- merge(x = snRNA, y = snATAC)
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

pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/coembed_umap.pdf",width=15,height=5)
DimPlot(coembed,reduction='umap',cols=myUmapcolors,raster=FALSE, group.by = c("patient_region", "cell_type"))
dev.off()

saveRDS(snATAC,"./01_all_celltype/all_celltype_snATAC_add_peaks_fragment.rds")
saveRDS(snRNA,"./01_all_celltype/all_celltype_snRNA_add_gene_symbol.rds")
saveRDS(coembed,"./01_all_celltype/all_celltype_snRNA_snATAC_merged.rds")

# The celltype proportion in different samples
# in RNA 
pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/RNA_celltype_proportion.pdf",width=10,height=5)
df <- as.data.frame(snRNA@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

# region
df_ct <- df %>%
    group_by(region, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(region, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# major_labl
df_ct <- df %>%
    group_by(major_labl, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

# donor_id
df_ct <- df %>%
    group_by(donor_id, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(donor_id, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

# patient_region
df_ct <- df %>%
    group_by(patient_region, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_region, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()

# in ATAC 

pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/ATAC_celltype_proportion.pdf",width=10,height=5)
df <- as.data.frame(snATAC@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

# region
df_ct <- df %>%
    group_by(region, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(region, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# major_labl
df_ct <- df %>%
    group_by(major_labl, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# donor_id
df_ct <- df %>%
    group_by(donor_id, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(donor_id, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_region
df_ct <- df %>%
    group_by(patient_region, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_region, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()

# in coembed


pdf("/md01/nieyg/project/EC/human_MI/01_all_celltype/coembed_celltype_proportion.pdf",width=10,height=5)
df <- as.data.frame(coembed@meta.data)
# development_stage
df_ct <- df %>%
    group_by(development_stage, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(development_stage, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# patient_group
df_ct <- df %>%
    group_by(patient_group, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_group, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

# region
df_ct <- df %>%
    group_by(region, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(region, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

# major_labl
df_ct <- df %>%
    group_by(major_labl, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(major_labl, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
# donor_id
df_ct <- df %>%
    group_by(donor_id, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(donor_id, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p

# patient_region
df_ct <- df %>%
    group_by(patient_region, cell_type) %>%
    summarise(counts = n()) %>%
    mutate(cell_proportion = counts / sum(counts))
cols <- myUmapcolors
p <- ggplot(df_ct, aes(patient_region, cell_proportion, fill=cell_type)) + 
    geom_bar(stat="identity", position = position_stack(reverse = TRUE)) +
    scale_fill_manual(values = cols) +
    theme_cowplot() +
    xlab("") + ylab("") +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
p
dev.off()





















