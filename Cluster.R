install.packages("BiocManager")
BiocManager::install("GenomicRanges")
if (!requireNamespace("GenomicRanges", quietly = TRUE)) {
  BiocManager::install("GenomicRanges")
}

library(GenomicRanges)
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(patchwork)

obj@assays$ATAC@fragments # check obj path

# Adjust path to where your .rds file actually is
obj <- readRDS("/Volumes/TOSHIBA/Human_scATAC_macula_RGC_coembed.rds")

Idents(obj) # cluster IDs for all cells
table(Idents(obj)) # counts of cells per cluster

predicted_id_counts <- table(obj$predicted.id)

# Identify the predicted.id values that have more than 20 cells
major_predicted_ids <- names(predicted_id_counts[predicted_id_counts > 20])
obj <- obj[, obj$predicted.id %in% major_predicted_ids]

# change cell identities to the per-cell predicted labels
Idents(obj) <- obj$predicted.id

# change to working with peaks instead of gene activities
DefaultAssay(obj) <- 'peaks'
obj <- RunTFIDF(obj)

clusters_to_refine <- major_predicted_ids

head(clusters_to_refine)

# refine clusters

# Run SVD (LSI) for dimensionality reduction
obj <- FindTopFeatures(obj, min.cutoff = 'q0')  # select variable peaks
obj <- RunSVD(obj)

# Construct a shared nearest neighbor graph
obj <- FindNeighbors(obj, reduction = "lsi", dims = 2:30)

# Cluster cells (adjust resolution as needed)
obj <- FindClusters(obj, resolution = 0.5)

obj <- RunUMAP(obj, reduction = "lsi", dims = 2:30)

p1 <- DimPlot(obj, label = TRUE) + ggtitle("Refined Clusters UMAP")
DimPlot(obj, group.by = "predicted.id")

ggsave(
  filename = "/Volumes/TOSHIBA/Refined_Clusters_UMAP.png",
  plot = p1,
  width = 8,
  height = 6,
  dpi = 300
)

# refine the first cluster ()

# 1. Subset the cells
RGC_cells <- subset(obj, subset = predicted.id == "ON_MGC")

# 2. Set the assay to peaks for ATAC data
DefaultAssay(RGC_cells) <- "peaks"

# 3. Normalize and reduce dimensions
RGC_cells <- RunTFIDF(RGC_cells)
RGC_cells <- FindTopFeatures(RGC_cells, min.cutoff = 10)
RGC_cells <- RunSVD(RGC_cells)

# 4. Build neighbor graph and recluster
RGC_cells <- FindNeighbors(RGC_cells, reduction = "lsi", dims = 2:30)
RGC_cells <- FindClusters(RGC_cells, resolution = 0.5)  # tweak resolution as needed

# 5. Run UMAP and visualize refined clusters
RGC_cells <- RunUMAP(RGC_cells, reduction = "lsi", dims = 2:30)
DimPlot(RGC_cells, label = TRUE, pt.size = 0.5)

# wilcox is the default option for test.use
cell_types <- c("ipRGCs", "mRGC5", "mRGC6", "OFF_MGC", "ON_MGC", "ON_PGC", "OFF_PGC", "SBCs")

# Create a folder to save the results
if(!dir.exists("Top_Peaks")) dir.create("Top_Peaks")

for (cluster_name in clusters_to_refine) {
  da_peaks <- FindMarkers(
    object = obj,
    ident.1 = cluster_name,
    test.use = 'wilcox',
    min.pct = 0.1
  )
  top_peaks <- head(da_peaks[order(da_peaks$p_val_adj), ], 20)
  write.csv(top_peaks, paste0("Top_Peaks/TopPeaks_", cluster_name, ".csv"))
}

# wilcox is the default option for test.use, use differential accessibility (differential peaks)
da_peaks <- FindMarkers(
  object = obj,
  ident.1 = "OFF_MGC",
  test.use = 'wilcox',
  min.pct = 0.1
)

head(da_peaks)

# sort plotting order according to similarities across annnotated cell types
obj <- SortIdents(obj)

DefaultAssay(obj) <- "peaks"

Idents(obj) <- obj$predicted.id



DefaultAssay(obj) = "ATAC"
Idents(obj) <- obj$predicted.id
CoveragePlot(
  object = obj,
  region=c("EOMES"),
  extend.upstream = 2000,
  extend.downstream = 11000,
  ncol=1
)

# calculate phastcons and phylop of the peaks i found in primates mostly?
# can try using a different graph if i think it's better

Fragments(obj)[[1]]@path <- "/Volumes/TOSHIBA/A7/fragments.tsv.gz"