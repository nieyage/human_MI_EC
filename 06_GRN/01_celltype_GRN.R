conda activate r4-base
# Step1: Seurat 2 loom and scanpy (scRNA)
library(loomR)
# Suppress loading messages when building the HTML
suppressPackageStartupMessages({
  library(SCopeLoomR)
  library(BiocParallel)
  library(ggplot2)
  library(data.table)
  library(grid)
  library(Signac)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(JASPAR2020)
  library(TFBSTools)
  library(patchwork)
  set.seed(1234)
})
library(SeuratDisk)
snRNA<- readRDS("../01_all_celltype/all_celltype_snRNA_add_gene_symbol.rds")

library(SeuratDisk)
sdata.loom <- as.loom(x = snRNA, filename = "/md01/nieyg/project/EC/human_MI/07_scenicplus/01_all_celltype/All_celltype_seurat2loom.loom", verbose = FALSE)
# Always remember to close loom files when done
sdata.loom$close_all()

conda activate scenicplus
# in python
import scanpy as sc
adata = sc.read_loom("/md01/nieyg/project/EC/human_MI/07_scenicplus/01_all_celltype/All_celltype_seurat2loom.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
adata.write('/md01/nieyg/project/EC/human_MI/07_scenicplus/01_all_celltype/All_celltype_snRNA_adata.h5ad', compression='gzip')













# Step2: Signac 2 pycisTopic 
conda activate scenicplus
#supress warnings
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import sys
import os
_stderr = sys.stderr
null = open(os.devnull,'wb')
work_dir = '/md01/nieyg/project/honeybee/honebee-latest-Version/09_GRN/02_SCENIC/04_scenicplus/'
import pycisTopic
#make a directory for to store the processed scRNA-seq data.
#if not os.path.exists(os.path.join(work_dir, '02_scATAC')):
#    os.makedirs(os.path.join(work_dir, '02_scATAC'))
tmp_dir = '/md01/nieyg/tmp'

# in shell
# merge the fragment file 
cp /md01/nieyg/project/honeybee/data/cellranger/Forager/Forager-NCBI-manually/outs/fragments.tsv .
mv fragments.tsv Forager_fragments.tsv
cp /md01/nieyg/project/honeybee/data/cellranger/NE/NE-NCBI-manually/outs/fragments.tsv .
mv fragments.tsv NE_fragments.tsv
cp /md01/nieyg/project/honeybee/data/cellranger/Nurse/Nurse-NCBI-manually/outs/fragments.tsv .
mv fragments.tsv Nurse_fragments.tsv
# add sample info in barcode 
awk '$4="NE_"$4' NE_fragments.tsv > NE_fragments_2.tsv
awk '$4="Nurse_"$4' Nurse_fragments.tsv > Nurse_fragments_2.tsv
awk '$4="Forager_"$4' Forager_fragments.tsv > Forager_fragments_2.tsv
cat NE_fragments_2.tsv Nurse_fragments_2.tsv Forager_fragments_2.tsv > last_fragments.tsv
awk '{print$1"\t"$2"\t"$3"\t"$4"\t"$5}' last_fragments.tsv > last_fragments_2.tsv
cat last_fragments_2.tsv |sort -k1,1 -k2,2n | bgzip > last_fragments.tsv.gz
tabix -b 2 -e 3 -p bed last_fragments.tsv.gz

# in python 
fragments_dict = {'ORN': os.path.join(work_dir, '02_scATAC/fragments_file/last_fragments.tsv.gz')}

# 1.generate pseudobulk ATAC-seq profiles per cell type and call peaks

# load the cell type annotation we generated in the scRNA-seq analysis above.
import scanpy as sc
adata = sc.read_h5ad(os.path.join(work_dir, '01_scRNA/adata.h5ad'))
cell_data = adata.obs
cell_data['sample_id'] = 'ORN'
cell_data['subcluster'] = cell_data['subcluster'].astype(str) # set data type of the celltype column to str, otherwise the export_pseudobulk function will complain.


