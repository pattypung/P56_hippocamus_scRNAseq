library(Seurat)
library(dplyr)
library(ggplot2)
library(tidyverse)

options(Seurat.object.assay.version = "v3")
#### import scRNA-seq data from the output files of 10X Cellranger count ####
### sample from P56HippoM1
input.data <- Read10X(data.dir = "/Users/boom7/Documents/YHH_lab/RNA-seq_data/P56_M1_HIP") 
P56HippoM1 <- CreateSeuratObject(counts = input.data, project = "P56HippoM1")
P56HippoM1
P56HippoM1$source <- "P56HippoM1"
P56HippoM1[["percent.mt"]] <- PercentageFeatureSet(P56HippoM1, pattern = "^mt-")
VlnPlot(P56HippoM1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P56HippoM1 <- subset(P56HippoM1, subset = nFeature_RNA > 200 & percent.mt < 25) ###QC check: percentage of mtRNA<25% or15/10
P56HippoM1 <- subset(P56HippoM1, subset = nFeature_RNA > 200)  ###QC check: >200gene per cell
P56HippoM1 <- NormalizeData(P56HippoM1, verbose = FALSE)
P56HippoM1 <- FindVariableFeatures(P56HippoM1, selection.method = "vst")

### sample from P56HippoM2
input.data <- Read10X(data.dir = "/Users/boom7/Documents/YHH_lab/RNA-seq_data/P56_M2_HIP")
P56HippoM2 <- CreateSeuratObject(counts = input.data, project = "P56HippoM2")
P56HippoM2
P56HippoM2$source <- "P56HippoM2"
P56HippoM2[["percent.mt"]] <- PercentageFeatureSet(P56HippoM2, pattern = "^mt-")
VlnPlot(P56HippoM2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
P56HippoM2 <- subset(P56HippoM2, subset = nFeature_RNA > 200 & percent.mt < 25)
P56HippoM2 <- subset(P56HippoM2, subset = nFeature_RNA > 200)
P56HippoM2 <- NormalizeData(P56HippoM2, verbose = FALSE)
P56HippoM2 <- FindVariableFeatures(P56HippoM2, selection.method = "vst")


#### Integration of samples ####
features <- SelectIntegrationFeatures(object.list = list(P56HippoM1,P56HippoM2)) ###Integrate Gene Read From Samples
P56HippoM1M2.anchors <- FindIntegrationAnchors(object.list = list(P56HippoM1,P56HippoM2), anchor.features = features) ###Gene For Cluster
P56HippoM1M2.combined <- IntegrateData(anchorset = P56HippoM1M2.anchors) ###Combined different samples

saveRDS(P56HippoM1M2.combined, "20240424_P56HippoM1M2_combined_CCA.rds")
P56HippoM1M2.combined <- readRDS("/Users/boom7/Documents/YHH_lab/scRNA-seq/NeMo_RRID_SCR_016152_P56_hippo/20240424_P56HippoM1M2_combined_CCA.rds")

DefaultAssay(P56HippoM1M2.combined) <- "integrated"
P56HippoM1M2.combined <- ScaleData(P56HippoM1M2.combined, verbose = FALSE)
P56HippoM1M2.combined <- RunPCA(P56HippoM1M2.combined, npcs = 30, verbose = FALSE)
P56HippoM1M2.combined <- RunUMAP(P56HippoM1M2.combined, reduction = "pca", dims = 1:30)
P56HippoM1M2.combined <- FindNeighbors(P56HippoM1M2.combined, reduction = "pca", dims = 1:30)
P56HippoM1M2.combined <- FindClusters(P56HippoM1M2.combined, resolution = 0.01)

P56HippoM1M2.combined$orig.ident <- factor(x = P56HippoM1M2.combined$orig.ident, levels = c("P56HippoM1", "P56HippoM2"))

### Draw UMAP
DimPlot(P56HippoM1M2.combined, reduction = "umap", label = T)
DimPlot(P56HippoM1M2.combined, reduction = "umap", label = T) + NoLegend()
DimPlot(P56HippoM1M2.combined, reduction = "umap", group.by = "orig.ident")
DimPlot(P56HippoM1M2.combined, split.by = "orig.ident", ncol = 2, label = F)
DimPlot(P56HippoM1M2.combined, split.by = "orig.ident", ncol = 2, label = T) + NoLegend()
DimPlot(P56HippoM1M2.combined, split.by = "orig.ident") + NoLegend()
DimPlot(P56HippoM1M2.combined, reduction = "umap")
DimPlot(P56HippoM1M2.combined, reduction = "umap", label = TRUE, label.size = 6) + NoLegend()


DefaultAssay(P56HippoM1M2.combined) <- "RNA"
P56HippoM1M2.combined <- ScaleData(P56HippoM1M2.combined, verbose = FALSE)

### check excitatory neuron marker
FeaturePlot(P56HippoM1M2.combined, features = c("Grin1"), pt.size = 0.4, cols = c("lightgrey","red"))
FeaturePlot(P56HippoM1M2.combined, features = c("Slc17a7"), pt.size = 0.4, cols = c("lightgrey","red"))
FeaturePlot(P56HippoM1M2.combined, features = c("Camk2a"), pt.size = 0.4, cols = c("lightgrey","red"))


saveRDS(P56HippoM1M2.combined, "20240424_P56Hippo_combined.rds")

#### Get the pseudobulk from a specific cluster ####
specific_cluster <- WhichCells(P56HippoM1M2.combined, idents = "0") # Replace "1" with your cluster number
counts_in_cluster <- GetAssayData(P56HippoM1M2.combined, slot = "counts")[, specific_cluster]
bulk_counts <- rowSums(counts_in_cluster)
bulk_data <- data.frame(gene = rownames(counts_in_cluster), bulk_count = bulk_counts)
write.csv(bulk_data, "bulk_rna_seq_data.csv", row.names = FALSE)


