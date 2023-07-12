# **_HELLO - MY NAME IS LE PHUOC QUYEN_**

### This script has been developed to perform a comprehensive large-scale copy number variation analysis on the epithelial cell population from 21 tumor samples of pancreatic ductal adenocarcinoma (PDAC) using single-cell RNA sequencing (scRNA-seq) data.

```R
# Path into cellranger output (aggregation 21 tumors sample)
setwd("/media/bio03/DATA/quyen/THESIS/DATASETS/PRJCA001063/seurat/aggr_21tumors/saveRDS/epithelial_21tumor_re0.5_16cluster.rds")
```
I.    Subset the cluster "2", "5", "13", "14", "0", "1", "15"                        from epithelial_21tumor    ---------  seurat.combined_1
II.   Subset the cluster "3", "4", "6", "7", "8", "9", "10", "11", "12"              from epithelial_21tumor    ---------  seurat.combined_2

# **I. QUALITY CONTROL + NORMALIZATION + DIMENSIONALITY REDUCTION + CLUSTERING**
## 1. Library packages
```R
library(ggplot2)
library(SingleR)
library(celldex)	
library(RColorBrewer)
library(SingleCellExperiment)
library(dplyr)
library(Seurat)
library(patchwork)
library(scater)
library(ExperimentHub)
library(clustifyr)
library(infercnv)
library(copykat)
library(Seurat)
library(Matrix)
library(stringr)
main_color = c("#7fcdbb", "#c994c7", "#3182bd", "#fa9fb5", "#41b6c4", "#31a354", "#fb6a4a", "#54278f", "#9ecae1", "#c51b8a", "#636363", "#fec44f", "#66c2a4", "#1c9099", "#6baed6", "#f03b20", "#980043", "#fee391", "#d95f0e", "#4eb3d3", "#41ab5d", "#a6bddb", "#016450", "#0570b0", "#78c679", "#8c6bb1", "#6e016b", "#fc9272", "#cb181d", "#ef6548", "#990000", "#525252", "#993404", "#756bb1", "#addd8e")
# Color_epi ( only for epithelial objective)
color_epi = c("#c994c7", "#7fcdbb",  "#41b6c4",  "#fb6a4a", "#9ecae1", "#fec44f",  "#addd8e", "#c51b8a", "#636363",  "#fee391", "#66c2a4", "#1c9099", "#6baed6", "#fa9fb5", "#41ab5d", "#f03b20", "#980043", "#d95f0e", "#4eb3d3", "#a6bddb", "#016450", "#0570b0", "#78c679", "#8c6bb1", "#6e016b", "#fc9272", "#cb181d", "#ef6548", "#990000", "#525252", "#993404", "#756bb1", "#addd8e")
color_srat_1 = c("#c994c7", "#7fcdbb",  "#41b6c4",  "#fec44f", "#fa9fb5", "#41ab5d", "#f03b20", "#980043", "#d95f0e", "#4eb3d3", "#a6bddb", "#016450", "#0570b0", "#78c679", "#8c6bb1", "#6e016b", "#fc9272", "#cb181d", "#ef6548", "#990000", "#525252", "#993404", "#756bb1", "#addd8e")
color_srat_2 = c("#fb6a4a", "#9ecae1",  "#addd8e", "#c51b8a", "#636363",  "#fee391", "#66c2a4", "#1c9099", "#6baed6", "#980043", "#d95f0e", "#4eb3d3", "#a6bddb", "#016450", "#0570b0", "#78c679", "#8c6bb1", "#6e016b", "#fc9272", "#cb181d", "#ef6548", "#990000", "#525252", "#993404", "#756bb1", "#addd8e")

```
## 2. ReadRDS epithelial_21tumor ( epithelial cluster from aggr_21tumor_norm)
```r
epithelial_21tumor <- readRDS("/media/bio03/DATA/quyen/THESIS/DATASETS/PRJCA001063/seurat/aggr_21tumors/saveRDS/epithelial_21tumor_re0.5_16cluster.rds")
DimPlot(epithelial_21tumor, label.size = 4, cols = color_epi, repel = T, label = T, reduction = "umap")
```
## 3. Subset the cluster "2", "5", "13", "14", "0", "1", "15" from epithelial_21tumor
```R
Idents(object = epithelial_21tumor) <- "seurat_clusters"
srat_1_epi_21tumor <- subset(x = epithelial_21tumor, idents = c("2", "5", "13", "14", "0", "1", "15"))
srat_1_epi_21tumor
DimPlot(srat_1_epi_21tumor, label.size = 3, cols = color_srat_1, repel = T, label = T, reduction = "umap")
```

# **II. Subsetting T and Endocrine cell from aggr_7normal_norm**
### The reference group cell for inferCNV and Copykat analysis
```R
normal_T_endo <- readRDS("/media/bio03/DATA/quyen/THESIS/DATASETS/PRJCA001063/seurat/aggr_7normal/saveRDS/normal_T_endo.rds")
DimPlot(normal_T_endo, label.size = 3, cols = main_color, repel = T, label = T, reduction = "umap")
```

# **III.  Seurat combined 2 objectives**
```r
# Merge srat_1 object with reference object by Seurat
seurat.combined_1 <- merge(srat_1_epi_21tumor, y = normal_T_endo, add.cell.ids = c("P", "N"), project = "merge")
seurat.combined_1
head(colnames(seurat.combined_1))
table(seurat.combined_1$orig.ident)
# Normalization, runPCA
seurat.combined_1 <- NormalizeData(seurat.combined_1)
seurat.combined_1 <- FindVariableFeatures(seurat.combined_1, selection.method = "vst", nfeatures = 3000)
all.genes <- rownames(seurat.combined_1)
seurat.combined_1 <- ScaleData(seurat.combined_1, features = all.genes)
seurat.combined_1 <- RunPCA(seurat.combined_1, features = VariableFeatures(object = seurat.combined_1))
ElbowPlot(seurat.combined_1)
# Clustering
seurat.combined_1 <- FindNeighbors(seurat.combined_1, dims = 1:10)
seurat.combined_1 <- FindClusters(seurat.combined_1, resolution = 0.5)
seurat.combined_1 <- RunUMAP(seurat.combined_1, dims = 1:10)
DimPlot(seurat.combined_1, label.size = 3, cols = main_color, repel = T, label = T, reduction = "umap")
# Identify T & endocrine cell cluster in combined object by level expression of marker genes
VlnPlot(seurat.combined_1, features = c("INS"), raster = FALSE)  # Endocrine cluster 11
VlnPlot(seurat.combined_1, features = c("CD8B", "CD8A", "CD3E"), raster = FALSE)   # T cell cluster 4, 13

Idents(object = seurat.combined_1) <- "seurat_clusters"
cluster12.markers <- FindMarkers(seurat.combined_1, ident.1 = 12, min.pct = 0.25)
head(cluster12.markers, n = 15)

# Add the cell_types feature 
seurat.combined_1[["cell_types"]] <- ""
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '0','cluster0',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '1','cluster1',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '2','cluster2',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '3','cluster3',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '4','T_cell_1',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '5','cluster5',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '6','cluster6',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '7','cluster7',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '8','cluster8',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '9','cluster9',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '10','cluster10',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '11','endocrine_normal',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '12','cluster12',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '13','T_cell_2',seurat.combined_1@meta.data$cell_types)
seurat.combined_1[["cell_types"]] <- ifelse(seurat.combined_1@meta.data$seurat_clusters == '14','cluster14',seurat.combined_1@meta.data$cell_types)
table(seurat.combined_1$cell_types)
# Save seurat.combined_1 object by RDS file
saveRDS(seurat.combined_1, file="/media/bio03/DATA/quyen/THESIS/DATASETS/PRJCA001063/CNVs_analysis/inferCNV/aggr_21tumor/seurat.combined_1/seurat.combined_1.rds")
seurat.combined_1 <- readRDS("/media/bio03/DATA/quyen/THESIS/DATASETS/PRJCA001063/CNVs_analysis/inferCNV/aggr_21tumor/seurat.combined_1/seurat.combined_1.rds")
```

# **IV.  Running inferCNV**
### seurat.combined_1
```r
library(infercnv)
# 1. Create raw counts matrix (cell-genes)
seurat.combined_1_matrix <- data.frame(seurat.combined_1@assays$RNA@counts)
write.table(seurat.combined_1_matrix, file="/media/bio03/DATA/quyen/THESIS/DATASETS/PRJCA001063/CNVs_analysis/inferCNV/aggr_21tumor/seurat.combined_1/run_infercnv/seurat.combined_1_matrix.txt")

# 2. Create annotations file which indicates which cells types
celltypes <- data.frame(colnames(seurat.combined_1), seurat.combined_1[["cell_types"]])   #2 column: 1 for cellbarcore, 1 for cell type
write.table(celltypes, file="/media/bio03/DATA/quyen/THESIS/DATASETS/PRJCA001063/CNVs_analysis/inferCNV/aggr_21tumor/seurat.combined_1/run_infercnv/seurat.combined_1_celltypes.txt", sep = "\t", row.names = FALSE, col.names = FALSE)

# 3. Gene/chromosome positions file    #https://data.broadinstitute.org/Trinity/CTAT/cnv/
curl -JLO  https://data.broadinstitute.org/Trinity/CTAT/cnv/gencode_v19_gene_pos.txt

# 4. inferCNV analysis
# Create the infercnv object
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=read.table("seurat.combined_1_matrix.txt"),
                                    annotations_file="seurat.combined_1_celltypes.txt",
                                    delim="\t",
                                    gene_order_file="gencode_v19_gene_pos.txt",
                                    ref_group_names= c("T_cell_1", "T_cell_2", "endocrine_normal")
)

infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir="/media/bio03/DATA/quyen/THESIS/DATASETS/PRJCA001063/CNVs_analysis/inferCNV/aggr_21tumor/seurat.combined_1/output",
                             cluster_by_groups=T,
                             denoise=T,
                             HMM=T,
                             analysis_mode="samples",
                             num_threads=12)

