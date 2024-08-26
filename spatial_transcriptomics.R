library(Seurat)
library(SeuratData)
library(ggplot2)
library(patchwork)
library(dplyr)
library(glmGamPoi)
library(future)

options(timeout = 12000)
InstallData("stxBrain")
#Loading the Data
brain <- LoadData("stxBrain", type = "anterior1")

#Data Preprocessing
plot1 <- VlnPlot(brain, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(brain, features = "nCount_Spatial") + theme(legend.position = "right")
png(filename ="nCount_Spatial.png")
dev.off()

options(future.globals.maxSize = 3 * 1024^3)
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE) #Normalisation

#Gene expression visualization
png(filename ="Spatial_Feature_Plot_Hpca_Ttr.png")
SpatialFeaturePlot(brain, features = c("Hpca", "Ttr"))
dev.off()

#Improving the visualization of the histology image
png(filename ="Spatial_Feature_Plot_Hpca.png")
p1 <- SpatialFeaturePlot(brain, features = "Hpca", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Hpca", alpha = c(0.1, 1))
p1 + p2
dev.off()

png(filename ="Spatial_Feature_Plot_Ttr.png")
p1 <- SpatialFeaturePlot(brain, features = "Ttr", pt.size.factor = 1)
p2 <- SpatialFeaturePlot(brain, features = "Ttr", alpha = c(0.1, 1))
p1 + p2
dev.off()

#Dimensionality reduction, clustering, and visualization
##UMAP
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

png(filename ="SpatialDimPlot_DimPlot_UMAP.png")
p1 <- DimPlot(brain, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
dev.off()
png(filename ="SpatialDimPlot_DimPlot_UMAP_Highlighted_Cells.png")
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,5, 8)), facet.highlight = TRUE, ncol = 3)
dev.off()

##TSNE

brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunTSNE(brain, reduction = "pca", dims = 1:30)

png(filename ="SpatialDimPlot_DimPlot_TSNE.png")
p1 <- DimPlot(brain, reduction = "tsne", label = TRUE)
p2 <- SpatialDimPlot(brain, label = TRUE, label.size = 3)
p1 + p2
dev.off()
png(filename ="SpatialDimPlot_DimPlot_TSNE_Highlighted_Cells.png")
SpatialDimPlot(brain, cells.highlight = CellsByIdentities(object = brain, idents = c(2, 1, 4, 3,5, 8)), facet.highlight = TRUE, ncol = 3)
dev.off()

#Identification of Spatially Variable Features
de_markers <- FindMarkers(brain, ident.1 = 5, ident.2 = 6)
png(filename ="SpatialFeaturePlot_Markers.png")
SpatialFeaturePlot(object = brain, features = rownames(de_markers)[1:3], alpha = c(0.1, 1), ncol = 3)
dev.off()

brain <- FindSpatiallyVariableFeatures(brain, assay = "SCT", features = VariableFeatures(brain)[1:1000], selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(brain, selection.method = "moransi"), 5)
png(filename ="SpatialVariableFeaturePlot_Top5_features.png")
SpatialFeaturePlot(brain, features = top.features, ncol = 3, alpha = c(0.1, 1))
dev.off()

#Subsetting out anatomical features
cortex <- subset(brain, idents = c(1, 2, 3, 4, 6, 7))
cortex <- subset(cortex, anterior1_imagerow > 400 | anterior1_imagecol < 150, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 275 & anterior1_imagecol > 370, invert = TRUE)
cortex <- subset(cortex, anterior1_imagerow > 250 & anterior1_imagecol > 440, invert = TRUE)
png(filename ="SpatialDimPlot_Anatomical_Regions.png")
p1 <- SpatialDimPlot(cortex, crop = TRUE, label = TRUE)
p2 <- SpatialDimPlot(cortex, crop = FALSE, label = TRUE, pt.size.factor = 1, label.size = 3)
p1 + p2
dev.off()

#Integration with scRNA data
allen_reference <- readRDS("allen_cortex.rds")
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>% RunPCA(verbose = FALSE) %>% RunUMAP(dims = 1:30)
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>% RunPCA(verbose = FALSE) #normalization after subsetting
png(filename ="SubClass_Annotations.png")
DimPlot(allen_reference, group.by = "subclass", label = TRUE)
dev.off()

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$subclass, prediction.assay = TRUE, weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay
DefaultAssay(cortex) <- "predictions" 
png(filename ="SpatialFeaturePlot_Prediction_Scores.png")
SpatialFeaturePlot(cortex, features = c("L2/3 IT", "L4"), pt.size.factor = 1.6, ncol = 2, crop = TRUE)
dev.off()

cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "moransi", features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex, selection.method = "moransi"), 4)
png(filename ="SpatialFeaturePlot_Prediction_Scores_TopClusters.png")
SpatialPlot(object = cortex, features = top.clusters, ncol = 2)
dev.off()

png(filename ="SpatialFeaturePlot_Localized_Patterns_All_Subsets.png")
SpatialFeaturePlot(cortex, features = c("Astro", "L2/3 IT", "L4", "L5 PT", "L5 IT", "L6 CT", "L6 IT", "L6b", "Oligo"), pt.size.factor = 1, ncol = 2, crop = FALSE, alpha = c(0.1, 1))
dev.off()

#Working with multiple slices in Seurat
brain2 <- LoadData("stxBrain", type = "posterior1")
brain2 <- SCTransform(brain2, assay = "Spatial", verbose = FALSE)
brain.merge <- merge(brain, brain2)

DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(brain), VariableFeatures(brain2))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)

png(filename ="SpatialDimPlot_Multiple_Slices.png")
DimPlot(brain.merge, reduction = "umap", group.by = c("ident", "orig.ident"))
dev.off()

png(filename ="SpatialDimPlot_Merge.png")
SpatialDimPlot(brain.merge)
dev.off()

png(filename ="SpatialFeaturePlot_Merge.png")
SpatialFeaturePlot(brain.merge, features = c("Hpca", "Plp1"))
dev.off()



