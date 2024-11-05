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
CAV_snRNA<- readRDS("../02_EC/CAV/EC_CAV_snRNA_add_gene_symbol.rds")
library(SeuratDisk)
sdata.loom <- as.loom(x = CAV_snRNA, filename = "/md01/nieyg/project/EC/human_MI/07_scenicplus/02_CAV_subtype/CAV_subtype_seurat2loom.loom", verbose = FALSE)
# Always remember to close loom files when done
sdata.loom$close_all()

conda activate scenicplus
# in python
import scanpy as sc
adata = sc.read_loom("/md01/nieyg/project/EC/human_MI/07_scenicplus/02_CAV_subtype/CAV_subtype_seurat2loom.loom", sparse=True, cleanup=False, X_name='spliced', obs_names='CellID', var_names='Gene', dtype='float32')
adata.write('/md01/nieyg/project/EC/human_MI/07_scenicplus/02_CAV_subtype/CAV_subtype_snRNA_adata.h5ad', compression='gzip')

# Step2: Signac 2 pycisTopic (snATAC)
library(cisTopic)
library(Signac)
library(Seurat)
library(tidyverse)
library(arrow)
## read in processed signac data
CAV_snATAC<- readRDS("../../02_EC/CAV/EC_CAV_snATAC_add_peaks_fragment.rds")
CAV_snATAC$CAV <- CAV_snATAC$CAV %>% as.character()
CAV_snATAC$CAV_detail <- CAV_snATAC$CAV_detail %>% as.character()
## create cisTopic from matrix
counts <- CAV_snATAC@assays$ATAC@counts
rownames(counts) <- rownames(counts) %>% sub("-", ":", .)
cisTopicObject <- createcisTopicObject(counts, 
                                       project.name='CAV')
## add metadata
cellData_md <- CAV_snATAC[[]]
cisTopicObject <- addCellMetadata(cisTopicObject, 
                                  cell.data = cellData_md)
## build models
cisTopicObject <- runWarpLDAModels(cisTopicObject, topic=c(16:24), 
                                   seed=123, nCores=5, 
                                   addModels=FALSE, 
                                   tmp = "tmp_models")
saveRDS(cisTopicObject, "CAV_snATAC_cisTopicObject.rds")
cisTopicObject <- selectModel(cisTopicObject, type='derivative')
library(arrow)
path <- "feathers/"
cisTopic_obj <- readRDS("CAV_snATAC_cisTopicObject.rds")
modelMat <- modelMatSelection(cisTopicObject, 'cell', 'Probability')
modelMat <- as.data.frame(modelMat)
write_feather(modelMat, sink=paste0(path, 'cell_topic.feather'))
modelMat <- modelMatSelection(cisTopicObject, 'region', 'Probability', all.regions=TRUE)
modelMat <- as.data.frame(modelMat)
write_feather(modelMat, sink=paste0(path, 'topic_region.feather'))
## cell data
celldata <- cisTopicObject@cell.data
celldata <- as.data.frame(celldata)
write_feather(celldata, sink=paste0(path, 'celldata.feather'))
ct <- cisTopicObject@count.matrix
regions <- rownames(ct)
ct <- as.data.frame(ct)
write_feather(ct, sink=paste0(path, 'count_matrix.feather'))
#write_feather(regions, sink=paste0(path, 'regions.feather'))
write.csv(regions, "regions.txt", quote = FALSE, row.names = FALSE)

conda deactivate
conda activate scenicplus


import scenicplus
import scanpy
import pycisTopic
## 1. Initialize cisTopic object
from pycisTopic.cistopic_class import *
from pycisTopic.utils import *
from pycisTopic.lda_models import CistopicLDAModel
work_dir="/data/R02/nieyg/project/EC/human_MI/07_scenicplus/02_CAV_subtype"
# Load count matrix
matrix_path="feathers/count_matrix.feather"
fragment_matrix = pd.read_feather(matrix_path)
regions=pd.read_csv("regions.txt")
region_names = regions["x"].tolist()
fragment_matrix.index = region_names
cisTopic_obj = pycisTopic.cistopic_class.create_cistopic_object(fragment_matrix)
# Also add the cell annotation, cell_data should be a pandas df with cells as rows (cell names as index) and variables as columns
cell_data = pd.read_feather("feathers/celldata.feather")
cell_data.index = fragment_matrix.columns
cisTopic_obj.cell_data.index = fragment_matrix.columns
cisTopic_obj.add_cell_data(cell_data)
cisTopic_obj.cell_names = fragment_matrix.columns
pickle.dump(cisTopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))


import pycisTopic
## 1. Initialize cisTopic object
from pycisTopic.cistopic_class import *
# Load count matrix
matrix_path=PATH_TO_FRAGMENTS_MATRIX
fragment_matrix = pd.read_feather(matrix_path)
cisTopic_obj = create_cistopic_object(fragment_matrix)
# Also add the cell annotation, cell_data should be a pandas df with cells as rows (cell names as index) and variables as columns
cisTopic_obj.add_cell_data(cell_data)
## 2. Add model
from pycisTopic.utils import *
model = load_cisTopic_model(PATH_TO_THE_FOLDER_WHERE_YOU_SAVED_DATA_FROM_R)
cistopic_obj.add_LDA_model(model)






