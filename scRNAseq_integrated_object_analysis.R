library("Seurat")
library("dplyr")
library("tidyr")
library("ggplot2")
library("ggplotify")
library("viridis")
library("RColorBrewer")
library("tibble")
library("pheatmap")
library("stringr")
library("tidyverse")
library("patchwork")
library("reshape2")
library("data.table")
library("Matrix")
library("cowplot")
library("ggrepel")
library("GGally")
library("scales")
library("ggnewscale")
library("forcats")
library("philentropy")
library("ggbreak")
library("ggpubr")
library("ggpattern")
library("UCell")


renv::init()
renv::snapshot()


packages <- c(
  "Seurat", "dplyr", "tidyr", "ggplot2", "ggplotify", "viridis",
  "RColorBrewer", "tibble", "pheatmap", "stringr", "tidyverse",
  "patchwork", "reshape2", "data.table", "Matrix", "cowplot",
  "ggrepel", "GGally", "scales", "ggnewscale", "forcats",
  "philentropy", "ggbreak", "ggpubr", "ggpattern", "UCell"
)

# Install from CRAN (and Bioconductor for UCell)
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

for (pkg in packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    tryCatch({
      install.packages(pkg)
    }, error = function(e) {
      message(paste("Trying Bioconductor for:", pkg))
      BiocManager::install(pkg, ask = FALSE)
    })
  } else {
    message(paste("Already installed:", pkg))
  }
}


setwd("~/Lab-files/Seurat/integrated")
setwd("~/")

# Functions from Long's util code
  # Exclude TCR genes from Highly variable genes
  quietTCRgenes <- function(matrix) {
    unwanted_genes <- c()
    for(file in list.files("long_excluded_gene_lists/TCellReceptors/")){
      csv <- read.csv(paste0("long_excluded_gene_lists/TCellReceptors/", file))
      unwanted_genes <- c(unwanted_genes, csv[["Approved.symbol"]])
    }
    matrix <- matrix[!rownames(matrix) %in% unwanted_genes, ]
    return(matrix)
  }
  
  # Exclude sex genes from highly variable genes
  quietSexgenes <- function(matrix){
    unwanted_genes <- "XIST|RPS4Y1|RPS4Y"
    unwanted_genes <- grep(pattern = unwanted_genes, x = rownames(matrix), value = TRUE)
    matrix <- matrix[!rownames(matrix) %in% unwanted_genes, ]
    return(matrix)
  }
  
  # Exclude BCR genes from highly variable genes
  quietBCRgenes <- function(matrix){
    unwanted_genes <- c()
    for(file in list.files("long_excluded_gene_lists/BCellImunoglobulin/")){
      csv <- read.csv(paste0("long_excluded_gene_lists/BCellImunoglobulin/", file))
      unwanted_genes <- c(unwanted_genes, csv[["Approved.symbol"]])
    }
    matrix <- matrix[!rownames(matrix) %in% unwanted_genes, ]
    return(matrix)
  }
  
  # Exclude TCR and BCR, as well as sex related genes
  ExcludeUnwantedGenes <- function(matrix, TCR=TRUE, BCR=TRUE, Sex=TRUE){
    if(TCR) 
    {matrix <- quietTCRgenes(matrix)}
    if(BCR) 
    {matrix <- quietBCRgenes(matrix)}
    if(Sex) 
    {matrix <- quietSexgenes(matrix)}
    return(matrix)
  }

#Updated TCR quieting function
  quietTCRgenes_new <- function(seurat_obj) {
    unwanted_genes <- c()
    for (file in list.files("long_excluded_gene_lists/TCellReceptors/", full.names = TRUE)) {
      csv <- read.csv(file)
      unwanted_genes <- c(unwanted_genes, csv[["Approved.symbol"]])
    }
    
    DefaultAssay(seurat_obj) <- "RNA"
    
    # Get current RNA features
    current_genes <- Features(seurat_obj[["RNA"]])
    genes_to_keep <- setdiff(current_genes, unwanted_genes)
    
    # Subset all layers safely
    for (layer in Layers(seurat_obj[["RNA"]])) {
      layer_data <- LayerData(seurat_obj[["RNA"]], layer = layer)
      genes_in_layer <- intersect(rownames(layer_data), genes_to_keep)
      
      # Only subset if genes exist in this layer
      if (length(genes_in_layer) > 0) {
        new_layer_data <- layer_data[genes_in_layer, , drop = FALSE]
        seurat_obj[["RNA"]] <- SetAssayData(
          seurat_obj[["RNA"]],
          layer = layer,
          new.data = new_layer_data
        )
      }
    }
    
    return(seurat_obj)
  }

#Updated sex genes quieting function
  quietSexgenes_new <- function(seurat_obj) {
    unwanted_genes <- c("XIST", "RPS4Y1", "RPS4Y")
    
    DefaultAssay(seurat_obj) <- "RNA"
    
    # Get current RNA features
    current_genes <- Features(seurat_obj[["RNA"]])
    genes_to_keep <- setdiff(current_genes, unwanted_genes)
    
    # Subset all layers safely
    for (layer in Layers(seurat_obj[["RNA"]])) {
      layer_data <- LayerData(seurat_obj[["RNA"]], layer = layer)
      genes_in_layer <- intersect(rownames(layer_data), genes_to_keep)
      
      # Only subset if genes exist in this layer
      if (length(genes_in_layer) > 0) {
        new_layer_data <- layer_data[genes_in_layer, , drop = FALSE]
        seurat_obj[["RNA"]] <- SetAssayData(
          seurat_obj[["RNA"]],
          layer = layer,
          new.data = new_layer_data
        )
      }
    }
    
    return(seurat_obj)
  }
  




# Load the Seurat object from the .rds file
mdata <- readRDS(file = "integratedP565_TCRGraph_v7.Rds")
ldata <- readRDS(file = "seurat_integrated.rds")

mdata <- readRDS(file = "final_mdata.rds")
mdata_pp65 <- readRDS(file = "final_mdata_pp65.rds")
mdata_CMV <- readRDS(file = "final_mdata_CMV.rds")
mdata_old <- readRDS (file = "integratedP565_TCRGraph_v7.Rds")

saveRDS(mdata, file = "final_mdata.rds")
saveRDS(mdata_pp65, file = "final_mdata_pp65.rds")
saveRDS(mdata_CMV, file = "final_mdata_CMV.rds")
saveRDS(mdata_CMV, file = "/box/Campbell/Sam K/CMV_project/####_RNAseq/# Cross-experiment analyses/### Seurat analyses/integrated/final_mdata_CMV.rds")
saveRDS(mdata_15, file = "mdata_15.rds")
saveRDS(mdata_20, file = "mdata_20.rds")
saveRDS(pseudobulk, file = "pseudobulk.rds")

# Exclude unwanted genes using Long's functions
mdata_pp65@assays$RNA@data <- ExcludeUnwantedGenes(mdata@assays$RNA@data)
mdata@assays$RNA <- ExcludeUnwantedGenes(mdata@assays$RNA)

mdata_old <- quietTCRgenes_new(mdata)
mdata_old <- quietSexgenes_new(mdata)

mdata <- subset(mdata, subset = percent_ribo >= 15)
mdata_20 <- subset(mdata, subset = percent_ribo >= 20)

mdata_D4_clone <- subset(mdata_CMV, subset = cloneId == "clone17donor4")
mdata_D5_clone <- subset(mdata_CMV, subset = cloneId == "clone48donor5")


# Create or update the 'specificity' column with multiple conditions
mdata$specificity <- ifelse(mdata$sampleBarcodeThresholdLogCpm %in% 
                              c("HTO-ET.pp65.specific.CD8s", "HTO-LT.pp65.specific.CD8s",
                                "HTO-Donor.4.Early.pp65", "HTO-Donor.4.Late.pp65",
                                "HTO-Donor.5.Early.pp65", "HTO-Donor.5.Late.pp65"), 
                            "pp65",
                            ifelse(mdata$sampleBarcodeThresholdLogCpm %in% 
                                  c("HTO-ET.IE1.specific.CD8s", "HTO-LT.IE1.specific.CD8s",
                                   "HTO-Donor.4.Early.IE1", "HTO-Donor.4.Late.IE1",
                                   "HTO-Donor.5.Early.IE1", "HTO-Donor.5.Late.IE1"), 
                                   "IE1", 
                                   ifelse(mdata$sampleBarcodeThresholdLogCpm %in% 
                                            c("HTO-ET.CD8s", "HTO-LT.CD8s",
                                              "HTO-Donor.4.Early.Bulk", "HTO-Donor.4.Late.Bulk",
                                              "HTO-Donor.5.Early.Bulk", "HTO-Donor.5.Late.Bulk"), 
                                              "Bulk", 
                                          NA)))

# Create or update the 'donor_identity' column with multiple conditions
mdata$donor_identity <- ifelse(grepl("donor1", mdata@meta.data[["cloneId"]]), "Donor 1",
                               ifelse(grepl("donor2", mdata@meta.data[["cloneId"]]), "Donor 2",
                                      ifelse(grepl("donor3", mdata@meta.data[["cloneId"]]), "Donor 3",
                                             ifelse(grepl("Donor.4", mdata@meta.data[["sampleBarcodeThresholdLogCpm"]]), "Donor 4",
                                                    ifelse(grepl("Donor.5", mdata@meta.data[["sampleBarcodeThresholdLogCpm"]]), "Donor 5",
                                          "NA")))))

# Create or update the 'timepoint' column with multiple conditions
mdata$timepoint <- ifelse(mdata$sampleBarcodeThresholdLogCpm %in% 
                              c("HTO-ET.pp65.specific.CD8s", "HTO-ET.IE1.specific.CD8s", "HTO-ET.CD8s",
                                "HTO-Donor.4.Early.pp65", "HTO-Donor.4.Early.IE1", "HTO-Donor.4.Early.Bulk",
                                "HTO-Donor.5.Early.pp65", "HTO-Donor.5.Early.IE1", "HTO-Donor.5.Early.Bulk"), 
                            "Early",
                            ifelse(mdata$sampleBarcodeThresholdLogCpm %in% 
                                     c("HTO-LT.pp65.specific.CD8s", "HTO-LT.IE1.specific.CD8s", "HTO-LT.CD8s",
                                       "HTO-Donor.4.Late.pp65", "HTO-Donor.4.Late.IE1", "HTO-Donor.4.Late.Bulk",
                                       "HTO-Donor.5.Late.pp65", "HTO-Donor.5.Late.IE1", "HTO-Donor.5.Late.Bulk"), 
                                   "Late", 
                                          NA))



#
# Reintegration for mdata
#
# Split and pre-process
seurat_list <- SplitObject(mdata, split.by = "donor_identity")
shared_genes <- Reduce(intersect, lapply(seurat_list, rownames))

# Process each subset
seurat_list <- lapply(seurat_list, function(x) {
  x <- subset(x, features = shared_genes)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})

# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_list)
# Scale and PCA
seurat_list <- lapply(seurat_list, function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
  return(x)
})
# Save other assays before integration
original_assays <- mdata@assays
non_rna_assays <- original_assays[!names(original_assays) %in% c("RNA")]
# Integration
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
mdata_integrated <- IntegrateData(anchorset = anchors)
# Reattach original non-RNA assays
for (assay_name in names(non_rna_assays)) {
  mdata_integrated[[assay_name]] <- non_rna_assays[[assay_name]]
}
# Replace mdata with the enriched object
mdata <- mdata_integrated


# Visualization for mdata
DefaultAssay(mdata) <- "FB"
mdata <- NormalizeData(mdata)
mdata <- FindVariableFeatures(mdata, selection.method = "vst")
mdata <- ScaleData(mdata, verbose = FALSE)
mdata <- RunPCA(mdata, npcs = 30, verbose = FALSE)
mdata <- RunUMAP(mdata, reduction = "pca", dims = 1:30, reduction.name = "umap")
mdata <- FindNeighbors(mdata, reduction = "pca", dims = 1:30)
mdata <- FindClusters(mdata, resolution = 1.0)

mdata <- RunUMAP(mdata, reduction = "integrated.cca", dims = 1:30)
mdata <- FindNeighbors(mdata, reduction = "integrated.cca", dims = 1:30)
mdata <- FindClusters(mdata, resolution = 0.3)

DefaultAssay(mdata) <- "integrated"
DimPlot(mdata, reduction = "umap", group.by = "donor_identity")
DimPlot(mdata, reduction = "umap", group.by = "specificity")
DimPlot(mdata, reduction = "umap", group.by = "timepoint")
DimPlot(mdata, reduction = "umap", group.by = "seurat_clusters")
DoHeatmap(mdata, features = c("LT-pp65-specific-CD8s", "LT-CD8s"))

FeaturePlot(mdata_CMV, features = "PD-1", reduction = "umap",
            cols = c("lightgrey", "darkred"))

DimPlot(mdata, reduction = "umap", 
        cells.highlight = WhichCells(mdata, expression = cloneId == 'clone48donor5'), 
          cols.highlight = "red", cols = "lightgray")
DimPlot(mdata, reduction = "umap", 
        cells.highlight = WhichCells(mdata, expression = cloneId == 'clone382donor5'), 
        cols.highlight = "red", cols = "lightgray")


# Subset the Seurat object to retain only the cells with the specified values in 'sampleBarcodeThresholdLogCpm'
mdata_pp65 <- subset(mdata, subset = specificity == 'pp65')
mdata_CMV <- subset(mdata, subset = specificity %in% c('pp65', 'IE1'))


mdata_CMV <- FindNeighbors(mdata, reduction = "integrated.cca", dims = 1:30)
mdata_CMV <- FindClusters(mdata, resolution = 0.5)

#
# Reintegration for mdata_CMV
#
seurat_list <- SplitObject(mdata_CMV, split.by = "donor_identity")
shared_genes <- Reduce(intersect, lapply(seurat_list, rownames))
seurat_list <- lapply(seurat_list, function(x) {
  x <- subset(x, features = shared_genes)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})
features <- SelectIntegrationFeatures(object.list = seurat_list)
seurat_list <- lapply(seurat_list, function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
  return(x)
})
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
mdata_CMV <- IntegrateData(anchorset = anchors)


# Visualization for mdata_CMV
DefaultAssay(mdata_CMV) <- "FB"
mdata_CMV <- NormalizeData(mdata_CMV)
mdata_CMV <- FindVariableFeatures(mdata_CMV, selection.method = "vst")
mdata_CMV <- ScaleData(mdata_CMV, verbose = FALSE)
mdata_CMV <- RunPCA(mdata_CMV, npcs = 30, verbose = FALSE)
mdata_CMV <- RunUMAP(mdata_CMV, reduction = "pca", dims = 1:30, reduction.name = "umap")
mdata_CMV <- FindNeighbors(mdata_CMV, reduction = "pca", dims = 1:30)
mdata_CMV <- FindClusters(mdata_CMV, graph.name = "integrated_snn", resolution = 0.3)


DefaultAssay(mdata_CMV) <- "RNA"
DimPlot(mdata_CMV, reduction = "umap", group.by = "donor_identity")
DimPlot(mdata_CMV, reduction = "umap", group.by = "specificity")
DimPlot(mdata_CMV, reduction = "umap", group.by = "seurat_clusters")
FeaturePlot(mdata_CMV, features = "CMC1", reduction = "umap",
            cols = c("lightgrey", "darkred"))

DimPlot(mdata_CMV, reduction = "umap", cells.highlight = WhichCells(mdata_CMV, idents = "4"), cols.highlight = "red", cols = "lightgray")
mdata_CMV <- subset(mdata_CMV, subset = seurat_clusters != 9)


# Exclude unwanted genes using Long's functions
mdata_pp65@assays$integrated@data <- ExcludeUnwantedGenes(mdata_pp65@assays$integrated@data)
mdata_pp65@assays$integrated@meta.features <- ExcludeUnwantedGenes(mdata_pp65@assays$integrated@meta.features)


#
# Reintegration for mdata_pp65
#
seurat_list <- SplitObject(mdata_pp65, split.by = "donor_identity")
shared_genes <- Reduce(intersect, lapply(seurat_list, rownames))
seurat_list <- lapply(seurat_list, function(x) {
  x <- subset(x, features = shared_genes)
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
  return(x)
})
features <- SelectIntegrationFeatures(object.list = seurat_list)
seurat_list <- lapply(seurat_list, function(x) {
  x <- ScaleData(x, features = features)
  x <- RunPCA(x, features = features)
  return(x)
})
anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
mdata_pp65 <- IntegrateData(anchorset = anchors)


# Visualization for mdata_pp65
DefaultAssay(mdata_pp65) <- "integrated"
#mdata_pp65 <- NormalizeData(mdata_pp65@assays$RNA)
#mdata_pp65 <- FindVariableFeatures(mdata_pp65, selection.method = "vst")
mdata_pp65 <- ScaleData(mdata_pp65, verbose = FALSE)
mdata_pp65 <- RunPCA(mdata_pp65, npcs = 30, verbose = FALSE)
mdata_pp65 <- RunUMAP(mdata_pp65, reduction = "pca", dims = 1:30, reduction.name = "umap")
mdata_pp65 <- FindNeighbors(mdata_pp65, reduction = "pca", dims = 1:30)
mdata_pp65 <- FindClusters(mdata_pp65, resolution = 0.5)

mdata_pp65 <- RunUMAP(mdata_pp65, reduction = "integrated.cca", dims = 1:30)
mdata_pp65 <- FindNeighbors(mdata_pp65, reduction = "pca", dims = 1:30)

mdata_pp65 <- FindClusters(mdata_pp65, resolution = 0.325) # adjust resolution to reduce clusters


DefaultAssay(mdata_pp65) <- "RNA"
DimPlot(mdata_pp65, reduction = "umap", group.by = "seurat_clusters")
FeaturePlot(mdata_pp65, features = "CMC1", reduction = "umap",
            cols = c("lightgrey", "darkred"))

#
# Plot the relative amounts of each clone in each cluster
#
# Extract metadata
meta_data <- mdata@meta.data

# Define clone IDs
clone_ids <- paste0("D2_pp65-clone-", 1:5)

# Count number of cells per cluster for each Global_Clone_ID
heatmap_data <- meta_data %>%
  filter(Global_Clone_ID %in% clone_ids) %>%
  group_by(seurat_clusters, Global_Clone_ID) %>%
  summarise(Cell_Count = n(), .groups = "drop")

# Extract numeric part of Global_Clone_ID for correct ordering
heatmap_data <- heatmap_data %>%
  mutate(Global_Clone_ID_Num = as.numeric(gsub("D2_pp65-clone-", "", Global_Clone_ID))) %>%
  arrange(Global_Clone_ID_Num)  # Order by numeric clone ID

# Create wide-format data (clone IDs as rows, clusters as columns)
heatmap_matrix <- heatmap_data %>%
  select(seurat_clusters, Global_Clone_ID, Cell_Count) %>%
  pivot_wider(names_from = seurat_clusters, values_from = Cell_Count, values_fill = 0)

# Save row names and remove them from the matrix
rownames(heatmap_matrix) <- heatmap_matrix$Global_Clone_ID
heatmap_matrix <- as.matrix(heatmap_matrix[, -1])  # Remove Global_Clone_ID column

# Sort columns numerically (clusters)
colnames(heatmap_matrix) <- as.character(sort(as.numeric(colnames(heatmap_matrix))))

# Plot heatmap with correctly ordered rows
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,  # Keep rows in specified order
         cluster_cols = FALSE,  # Keep clusters in numeric order
         scale = "none",        # Keep absolute cell counts
         main = "Cell Counts per Cluster for Selected Global_Clone_IDs",
         fontsize_row = 10, 
         fontsize_col = 10)

# Save
# Set output file
png("exports/integrated_v5/stacked_bar_graphs/D2_pp65_clones_per_cluster.png", width = 1200, height = 1000, res = 150)  # Adjust size and resolution

# Generate heatmap
pheatmap(heatmap_matrix,
         cluster_rows = FALSE,  
         cluster_cols = FALSE,  
         scale = "none",        
         main = "Cell Counts per Cluster for Selected Global_Clone_IDs",
         fontsize_row = 10, 
         fontsize_col = 10)

# Save and close the file
dev.off()





#
#Plot cells based on antigen specificity
#
DimPlot(mdata, group.by = 'seurat_clusters', 
        reduction = 'umap',
        cells.highlight = WhichCells(mdata, 
                                     expression = #sampleBarcodeThresholdLogCpm == 'HTO-ET.IE1.specific.CD8s' | 
                                                  sampleBarcodeThresholdLogCpm == 'HTO-ET.IE1.specific.CD8s'),
        sizes.highlight = 1.25)  # Set size for highlighted cells (default is usually 1)



FeaturePlot(mdata_pp65, features = "CD45RA", reduction = "umap",
            cols = c("lightgrey", "darkred"))




#
# Find markers for all clusters
#
Idents(mdata) <- "metacluster"
Idents(mdata) <- "seurat_clusters"
Idents(mdata_CMV) <- "seurat_clusters"

cluster_markers <- FindAllMarkers(mdata, only.pos = TRUE, logfc.threshold = 0.25)
filtered_markers <- cluster_markers %>% filter(pct.1 > 0.1)

mdata_CMV <- JoinLayers(mdata_CMV)
cluster_markers_CMV <- FindAllMarkers(mdata_CMV, only.pos = TRUE, logfc.threshold = 0.1, min.pct = 0.1)
filtered_markers <- cluster_markers_CMV %>% filter(pct.1 > 0.1)

cluster_markers_CMV <- FindAllMarkers(mdata_CMV, only.pos = TRUE, logfc.threshold = 0.1, min.pct = 0.1)
filtered_markers <- cluster_markers_CMV

mdata_pp65 <- JoinLayers(mdata_pp65)
cluster_markers_pp65 <- FindAllMarkers(mdata_pp65, only.pos = TRUE, logfc.threshold = 0.1, min.pct = 0.1)
filtered_markers_pp65 <- cluster_markers_pp65 %>% filter(pct.1 > 0.1)

# Plot a heatmap of the top 10 markers per cluster
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- filtered_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- filtered_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(mdata, features = top10$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

top10_CMV <- cluster_markers_CMV %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_CMV <- filtered_markers_CMV %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(mdata_CMV, features = top10_CMV$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

top10_pp65 <- cluster_markers_pp65 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_pp65 <- filtered_markers_pp65 %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
DoHeatmap(mdata_pp65, features = top10_pp65$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

D1_subset <- subset(mdata_pp65, idents = paste0("D1_global-clone-", 1:14))
DoHeatmap(D1_subset, features = top10_pp65$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))

D2_subset <- subset(mdata_pp65, idents = paste0("D2_global-clone-", 1:9))
DoHeatmap(D2_subset, features = top10_pp65$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))


Idents(mdata_CMV) <- "specificity"
cluster_CMV_markers <- FindMarkers(mdata_CMV, ident.1 = "pp65", ident.2 = "IE1")
top_genes <- cluster_CMV_markers %>%
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(10) %>%
  pull(gene)
DoHeatmap(mdata_CMV, features = top_genes) + 
  scale_fill_gradientn(colors = c("blue", "white", "red"))

top_genes <- filtered_markers %>%
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(10) %>%
  pull(gene)



  # Set the identity class of the Seurat object to the 't_clonotype_id' column
  Idents(mdata_pp65) <- "t_clonotype_id"
  # Get unique clonotype IDs, removing NA values
  clonotype_ids <- unique(na.omit(mdata_pp65$t_clonotype_id))
  # Initialize a list to store the results
  all_markers <- list()
  # Loop through each clonotype and find markers compared to all other cells
  for (clonotype in clonotype_ids) {
    # Ensure the clonotype is valid and not NA
    if (!is.na(clonotype) && clonotype != "") {
      # Perform differential expression analysis for clonotype vs. all other cells
      markers <- FindMarkers(mdata_pp65, ident.1 = clonotype)
      
      # Store the markers in a list
      all_markers[[clonotype]] <- markers
      
      # Print the top markers for the current clonotype
      print(paste("Top markers for", clonotype))
      print(head(markers))
    }
  }
  # To access markers for a specific clonotype, you can reference the list:
  # e.g., markers for clonotype1:
  head(all_markers[['clonotype1']])

# Plot the top marker for cluster n
top_gene <- cluster_markers %>% top_n(n = 1, wt = avg_log2FC) %>% pull(gene)
VlnPlot(mdata, features = top_gene)





#
# R script to ADD GLOBAL CLONES and generate clonotype plots
#

# Load the CSV file for adding global clone names to the object
csv_path <- "exports/Final_Global_Clone_IDs.csv"
global_clones <- read.csv(csv_path, stringsAsFactors = FALSE)

# Extract metadata from Seurat object
metadata <- mdata_CMV@meta.data

# Initialize new metadata column with NA
metadata$Global_Clone_ID <- NA

# Create a named lookup vector for fast matching
clone_lookup <- setNames(global_clones$GClone_name, global_clones$Clone_ID)

# Update 'Global_Clone_ID' based on matches in 'cloneId'
metadata$Global_Clone_ID <- clone_lookup[metadata$cloneId]

# Assign updated metadata back to Seurat object
mdata_CMV@meta.data <- metadata


# Open PDF file
pdf(file = "exports/integrated_v5/Clonotype_plots.pdf", width = 8, height = 6)

# Loop through donors
for(x in 1:5) {
  
  # Loop through IE1 clones
  for(i in 1:20) {
    clonotype_id <- paste0("D", x, "_IE1-clone-", i)
    
    # Check if there are any cells with this clonotype
    cells_to_highlight <- WhichCells(mdata_CMV, expression = Global_Clone_ID == clonotype_id)
    
    # Only plot if there are cells to highlight
    if (length(cells_to_highlight) > 0) {
      # Generate the DimPlot for each clonotype
      plot <- DimPlot(mdata_CMV, reduction = 'umap', 
                      cells.highlight = cells_to_highlight, 
                      sizes.highlight = 1.25) + ggtitle(clonotype_id)
      
      # Print each plot to the PDF file
      print(plot)
    }
  }
  
  #Loop through pp65 clones
  for(i in 1:20) {
    clonotype_id <- paste0("D", x, "_pp65-clone-", i)
    
    # Check if there are any cells with this clonotype
    cells_to_highlight <- WhichCells(mdata_CMV, expression = Global_Clone_ID == clonotype_id)
    
    # Only plot if there are cells to highlight
    if (length(cells_to_highlight) > 0) {
      # Generate the DimPlot for each clonotype
      plot <- DimPlot(mdata_CMV, reduction = 'umap', 
                      cells.highlight = cells_to_highlight, 
                      sizes.highlight = 1.25) + ggtitle(clonotype_id)
      
      # Print each plot to the PDF file
      print(plot)
    }
  }
}


# Close the PDF file
dev.off()


#
# Investigating clones with mislabeled specificity
#

DimPlot(mdata_CMV, reduction = 'umap', 
        cells.highlight = WhichCells(mdata_CMV, expression = cloneId == "clone17donor4"), 
        sizes.highlight = 1.25)

DimPlot(mdata_CMV, reduction = 'umap', 
        cells.highlight = WhichCells(mdata_CMV, expression = cloneId == "clone48donor5"), 
        sizes.highlight = 1.25)

DimPlot(mdata_D4_clone, reduction = "umap", group.by = "specificity")
DimPlot(mdata_D5_clone, reduction = "umap", group.by = "specificity")


# Extract metadata
meta <- mdata@meta.data

# Count cells by cloneId and specificity
counts_table <- meta %>%
  group_by(cloneId, specificity) %>%
  summarise(cell_count = n(), .groups = "drop") %>%
  pivot_wider(names_from = specificity, values_from = cell_count, values_fill = 0) %>%
  arrange(cloneId)

# Add a total cell count column
counts_table$total_cells <- rowSums(counts_table[ , !(names(counts_table) %in% c("cloneId"))])

# Remove rows with 'notCalled' in cloneId
counts_table <- counts_table %>%
  filter(!grepl("notCalled", cloneId))

# List of donor identifiers
donors <- paste0("donor", 1:5)
specificities <- c("Bulk", "pp65", "IE1")

# Loop through each donor
for (donor in donors) {
  # Filter for cloneIds containing the donor string
  donor_table <- counts_table %>%
    filter(grepl(donor, cloneId)) %>%
    arrange(desc(total_cells))
  
  # Determine which specificity is highest for each cloneId
  donor_table$max_specificity <- apply(
    donor_table[, specificities], 1,
    function(row) {
      if (all(is.na(row))) return(NA)
      spec <- names(row)[which.max(row)]
      return(spec)
    }
  )
  
  # Convert "donor1" to "Donor_1" for filenames
  donor_label <- gsub("donor(\\d+)", "Donor_\\1", donor)
  
  # Loop through each specificity to create subgroups
  for (spec in specificities) {
    spec_table <- donor_table %>%
      filter(max_specificity == spec)
    
    # Export the subgroup
    filename <- paste0("exports/TCRgraph_clonotype_counts/Clone_specificity_counts__", donor_label, "_", spec, ".csv")
    write.csv(spec_table, filename, row.names = FALSE)
  }
}


#
#
# Feature plot saving
#
#

FeaturePlot(mdata_pp65, features = "SERTAD1", reduction = "umap",
            cols = c("blue", "white", "red"))

#Create and save the FeaturePlot without storing it as a variable
plot_feature <- function(gene) {
  p <- FeaturePlot(
    mdata,
    features = gene,
    reduction = "umap",
    cols = c("lightgrey", "darkred"),
    min.cutoff = "0"
  )
  
  ggsave(
    filename = paste0("exports/integrated_v5/UMAP_", gene, ".png"),
    plot = p,
    width = 6, height = 5, dpi = 300
  )
}

# Example use:
plot_feature("KLRG1")


ggsave("exports/integrated_v5/UMAP_ZBTB16.png", plot = 
         FeaturePlot(mdata, features = "ZBTB16", reduction = "umap",
                     cols = c("lightgrey", "darkred")), 
       width = 6, height = 5, dpi = 300)

# Set desired order
mdata_CMV$metacluster <- factor(
  mdata_CMV$metacluster,
  levels = c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")
)

ggsave("exports/integrated_v5/CMV_UMAP_by_metacluster.png", plot = 
         DimPlot(mdata_CMV, group.by = "metacluster", reduction = "umap"), 
       width = 6, height = 5, dpi = 300)


ggsave(
  "exports/integrated_v5/CMV_DEGs_per_specificity.png",
  plot = DoHeatmap(
    mdata_CMV,
    features = top10$gene,
  ) +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +  # expression values
    theme(
      axis.text.y = element_text(size = 12)  # gene names
      # legend.text and legend.title can be adjusted here if needed
    ),
  width = 10, height = 20, dpi = 300
)

ggsave(
  "exports/integrated_v5/CMV_DEGs_per_specificity.png",
  plot = DoHeatmap(
    mdata_CMV,
    features = top10$gene,
  ) +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +  # expression values
    theme(
      axis.text.y = element_text(size = 20)  # gene names
      # legend.text and legend.title can be adjusted here if needed
    ),
  width = 10, height = 10, dpi = 300
)

ggsave("exports/integrated_v2/IE1_specific_cells.png", plot = 
         DimPlot(mdata, group.by = 'seurat_clusters', 
                 reduction = 'umap',
                 cells.highlight = WhichCells(mdata, 
                                              expression = sampleBarcodeThresholdLogCpm == 'HTO-ET.IE1.specific.CD8s' | 
                                                           sampleBarcodeThresholdLogCpm == 'HTO-LT.IE1.specific.CD8s'),
                 sizes.highlight = 1),  # Set size for highlighted cells (default is usually 1) 
       width = 6, height = 5, dpi = 300)


# Order the metaclusters
metacluster_order <- c("Naive", "Central Memory", "Effector Memory", 
                       "NK-like", "ISG+", "MAIT")
mdata$metacluster <- factor(mdata$metacluster,
                                     levels = metacluster_order)

metacluster_order <- c("Tcm", "Tem", "T-CMC1", 
                       "NK-like", "ISG+")
mdata_CMV$metacluster <- factor(mdata_CMV$metacluster,
                            levels = metacluster_order)

# Example DimPlot to extract colors
dp <- DimPlot(mdata_CMV, group.by = "specificity", reduction = "umap")

# Extract the mapping of colors from the ggplot object
# In Seurat v5, the color scale is stored in dp$layers[[1]]$mapping$colour
# But safer is to use the scales via the ggplot_build object

ggb <- ggplot_build(dp)  # build the plot

# The first layer contains the fill colors actually used
# Extract unique colors assigned to each metacluster
layer_data <- ggb$data[[1]]
color_mapping <- unique(layer_data$colour)  # colors actually used in the plot

# Make sure colors are in order of metacluster levels
metacluster_levels <- levels(factor(mdata_CMV$specificity))
names(color_mapping) <- metacluster_levels  # assign names

# Reverse colors
reversed_colors <- rev(color_mapping)


# Make sure Idents are set to specificity
Idents(mdata_CMV) <- mdata_CMV$specificity

# Get the current Idents as a factor
current_levels <- levels(Idents(mdata_CMV))

# Reverse the factor levels
Idents(mdata_CMV) <- factor(Idents(mdata_CMV), levels = rev(current_levels))




# Set the specificity factor levels in desired order
mdata_CMV$specificity <- factor(mdata_CMV$specificity, levels = c("pp65", "IE1"))

# Generate individual plots for each specificity in order
plots <- lapply(levels(mdata_CMV$specificity), function(sp) {
  DimPlot(
    mdata_CMV,
    group.by = "specificity",
    cells = WhichCells(
      mdata_CMV,
      expression = specificity == sp
    )
  ) +
    scale_color_manual(values = c(
      "pp65" = "#0047FF",
      "IE1"  = "#FF1A1A"
    )) +
    ggtitle(paste0(sp, "-Specific Cells")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
})

# Combine plots side-by-side
combined_plot <- wrap_plots(plots, nrow = 1)

# Display the plot
print(combined_plot)

# Save to file (optional)
ggsave("exports/integrated_v5/CMV_UMAPs_split_by_specificity.png", combined_plot, width = 8, height = 4)




#
# Export total UMAP plots subsetted on spec or donor
#
export_dir <- "exports/integrated_v5/total_UMAP_plots_split"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# List of specificities
specs <- c("pp65", "IE1", "Bulk")

# Loop over each specificity and save a UMAP plot
for (spec in specs) {
  # Subset the object
  mdata_sub <- subset(mdata, subset = specificity == spec)
  
  # Generate UMAP plot
  p <- DimPlot(mdata_sub, group.by = "seurat_clusters", reduction = "umap") +
    ggtitle(paste("UMAP -", spec))
  
  # Save
  ggsave(
    filename = file.path(export_dir, paste0("Total_UMAP_split_", spec, ".png")),
    plot = p,
    width = 6,
    height = 5,
    dpi = 300
  )
}

#
# Explicitly define all combos
donors <- unique(mdata$donor_identity)
timepoints <- unique(mdata$timepoint)

# Get global UMAP coordinates from the full object
umap_coords <- Embeddings(mdata, reduction = "umap")
x_range <- range(umap_coords[, 1])
y_range <- range(umap_coords[, 2])

for (donor in donors) {
  for (tp in timepoints) {
    message("Processing: ", donor, " - ", tp)
    
    # Get cells for this donor-timepoint
    cells_use <- colnames(mdata)[mdata$donor_identity == donor & mdata$timepoint == tp]
    
    # Subset using cells, not metadata expression
    mdata_sub <- subset(mdata, cells = cells_use)
    
    # Plot
    p <- DimPlot(mdata_sub, group.by = "seurat_clusters", reduction = "umap") +
      ggtitle(paste("UMAP -", donor, "-", tp)) +
      xlim(x_range) +
      ylim(y_range)
    
    # Filename
    fname <- paste0("Total_UMAP_split_", gsub(" ", "_", donor), "_", tp, ".png")
    
    # Save
    ggsave(
      filename = file.path(export_dir, fname),
      plot = p,
      width = 6,
      height = 5,
      dpi = 300
    )
  }
}



# Subset the Seurat object for cells belonging to the two groups
subset_obj <- subset(mdata, subset = specificity %in% c("pp65", "IE1"))

# Set the identity class based on the 'specificity' metadata
Idents(subset_obj) <- "specificity"

# Perform differential expression analysis
de_results <- FindMarkers(subset_obj, ident.1 = "pp65", ident.2 = "IE1")

# View the top differentially expressed genes
head(de_results)

# Filter for top markers based on adjusted p-value and log fold change
top_genes <- de_results %>%
  rownames_to_column("gene") %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(avg_log2FC)) %>%
  head(10)  # Change the number to get more top genes if needed

# Print the top differentially expressed genes
print(top_genes)

# Save the results to a CSV file
write.csv(de_results, "differential_expression_results.csv", row.names = TRUE)




# CLONOTYPE ANALYSIS
# Extract numbers for each TCRgraph clonotype
#
# Step 1: Extract metadata from Seurat object
meta <- mdata@meta.data

# Step 2: Filter to only cells with specificity == "pp65"
meta_pp65 <- meta %>% filter(specificity == "IE1")

# Step 3: Remove rows where cloneId contains "notCalled"
meta_pp65 <- meta_pp65 %>% filter(!grepl("notCalled", cloneId))

# Step 4: Count Early and Late timepoint cells per cloneId
clone_id_counts <- meta_pp65 %>%
  group_by(cloneId) %>%
  summarize(
    Early_count = sum(timepoint == "Early"),
    Late_count  = sum(timepoint == "Late"),
    .groups = "drop"
  )

# Step 5: Split into donor-specific tables and sort
Donor_1 <- clone_id_counts %>%
  filter(grepl("donor1", cloneId)) %>%
  arrange(desc(Early_count + Late_count))

Donor_2 <- clone_id_counts %>%
  filter(grepl("donor2", cloneId)) %>%
  arrange(desc(Early_count + Late_count))

Donor_3 <- clone_id_counts %>%
  filter(grepl("donor3", cloneId)) %>%
  arrange(desc(Early_count + Late_count))

Donor_4 <- clone_id_counts %>%
  filter(grepl("donor4", cloneId)) %>%
  arrange(desc(Early_count + Late_count))

Donor_5 <- clone_id_counts %>%
  filter(grepl("donor5", cloneId)) %>%
  arrange(desc(Early_count + Late_count))

# Step 6: Pad shorter donor tables to equal length
max_rows <- max(nrow(Donor_1), nrow(Donor_2), nrow(Donor_3), nrow(Donor_4), nrow(Donor_5))

pad_df <- function(df) {
  rows_to_add <- max_rows - nrow(df)
  if (rows_to_add > 0) {
    empty <- data.frame(
      cloneId = rep(NA, rows_to_add),
      Early_count = rep(NA, rows_to_add),
      Late_count = rep(NA, rows_to_add)
    )
    df <- rbind(df, empty)
  }
  return(df)
}

Donor_1 <- pad_df(Donor_1)
Donor_2 <- pad_df(Donor_2)
Donor_3 <- pad_df(Donor_3)
Donor_4 <- pad_df(Donor_4)
Donor_5 <- pad_df(Donor_5)

# Step 7: Combine all donor tables into a wide format
master_results <- data.frame(
  "Donor 1 cloneId"    = Donor_1$cloneId,
  "Donor 1 Early"      = Donor_1$Early_count,
  "Donor 1 Late"       = Donor_1$Late_count,
  
  "Donor 2 cloneId"    = Donor_2$cloneId,
  "Donor 2 Early"      = Donor_2$Early_count,
  "Donor 2 Late"       = Donor_2$Late_count,
  
  "Donor 3 cloneId"    = Donor_3$cloneId,
  "Donor 3 Early"      = Donor_3$Early_count,
  "Donor 3 Late"       = Donor_3$Late_count,
  
  "Donor 4 cloneId"    = Donor_4$cloneId,
  "Donor 4 Early"      = Donor_4$Early_count,
  "Donor 4 Late"       = Donor_4$Late_count,
  
  "Donor 5 cloneId"    = Donor_5$cloneId,
  "Donor 5 Early"      = Donor_5$Early_count,
  "Donor 5 Late"       = Donor_5$Late_count
)

write.csv(master_results, "exports/TCRgraph_clonotype_counts/IE1_TCRgraph_clonotype_counts.csv", row.names = FALSE)




# CLONOTYPE ANALYSIS
# Extract TCR sequences for each clonotype

# Extract metadata
metadata <- mdata_CMV@meta.data

# Initialize empty data frames for AA and NT results
results_aa <- data.frame(Specificity = character(), Clone_ID = character(),
                         Global_Clone_ID = character(), Donor = character(),
                         stringsAsFactors = FALSE)
results_nt <- data.frame(Specificity = character(), Clone_ID = character(),
                         Global_Clone_ID = character(), Donor = character(),
                         stringsAsFactors = FALSE)

# Function to update or add a row
update_results <- function(results, clone_id, sequence, specificity, global_clone_id, donor, seq_prefix) {
  
  row_index <- which(results$Clone_ID == clone_id)
  
  if (length(row_index) == 0) {
    # Clone_ID not found — add new row
    new_row <- data.frame(Specificity = specificity, Clone_ID = clone_id,
                          Global_Clone_ID = global_clone_id, Donor = donor,
                          stringsAsFactors = FALSE)
    new_row[paste0(seq_prefix, "_sequence_1")] <- sequence
    new_row["Count_1"]                         <- 1
    results <- bind_rows(results, new_row)
    
  } else {
    row_index <- row_index[1]
    
    existing_seq_cols <- grep(paste0("^", seq_prefix, "_sequence_"), names(results), value = TRUE)
    
    match_found <- FALSE
    for (seq_col in existing_seq_cols) {
      count_col <- gsub(paste0(seq_prefix, "_sequence_"), "Count_", seq_col)
      if (!is.na(results[row_index, seq_col]) && results[row_index, seq_col] == sequence) {
        results[row_index, count_col] <- results[row_index, count_col] + 1
        match_found <- TRUE
        break
      }
    }
    
    if (!match_found) {
      # Use max existing index + 1 to avoid overwrite bug
      existing_indices <- as.integer(gsub(paste0(seq_prefix, "_sequence_"), "", existing_seq_cols))
      filled_indices   <- existing_indices[!is.na(results[row_index, existing_seq_cols])]
      next_index       <- if (length(filled_indices) == 0) 1 else max(filled_indices) + 1
      
      new_seq_col   <- paste0(seq_prefix, "_sequence_", next_index)
      new_count_col <- paste0("Count_", next_index)
      
      if (!(new_seq_col %in% names(results))) {
        results[[new_seq_col]]   <- NA
        results[[new_count_col]] <- NA
      }
      
      results[row_index, new_seq_col]   <- sequence
      results[row_index, new_count_col] <- 1
    }
  }
  
  return(results)
}

# Iterate through Seurat metadata
for (i in 1:nrow(metadata)) {
  clone_id        <- metadata$cloneId[i]
  cdr3_aa         <- metadata$t_cdr3s_aa[i]
  cdr3_nt         <- metadata$t_cdr3s_nt[i]
  specificity     <- metadata$specificity[i]
  global_clone_id <- metadata$Global_Clone_ID[i]
  donor           <- metadata$donor_identity[i]
  
  if (!is.na(clone_id) && !is.na(cdr3_aa)) {
    results_aa <- update_results(results_aa, clone_id, cdr3_aa, specificity, global_clone_id, donor, "CDR3_aa")
  }
  if (!is.na(clone_id) && !is.na(cdr3_nt)) {
    results_nt <- update_results(results_nt, clone_id, cdr3_nt, specificity, global_clone_id, donor, "CDR3_nt")
  }
}

# Shared post-processing function
process_results <- function(results) {
  
  donor_order       <- c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5")
  specificity_order <- c("pp65", "IE1")
  
  results %>%
    mutate(Total_count = rowSums(select(., starts_with("Count_")), na.rm = TRUE)) %>%
    relocate(Total_count, .after = Global_Clone_ID) %>%
    filter(Total_count >= 5) %>%
    filter(!grepl("notCalled", Clone_ID)) %>%
    select(where(~ !all(is.na(.)))) %>%
    # Derive specificity for sorting from Global_Clone_ID (reliable)
    # rather than the Specificity metadata column (can be mislabeled)
    # e.g. "D3_pp65-clone-7" -> "pp65", "D4_IE1-clone-2" -> "IE1"
    mutate(
      sort_specificity = str_extract(Global_Clone_ID, "pp65|IE1"),
      Donor            = factor(Donor,            levels = donor_order),
      sort_specificity = factor(sort_specificity, levels = specificity_order)
    ) %>%
    arrange(Donor, sort_specificity, desc(Total_count)) %>%
    # Drop the temporary sort column and convert factors back to character
    select(-sort_specificity) %>%
    mutate(Donor = as.character(Donor))
}
filtered_aa <- process_results(results_aa)
filtered_nt <- process_results(results_nt)

# Export
write.csv(filtered_aa, "exports/5-D_clonotypes_TCR_aa_sequences_filtered.csv", row.names = FALSE)
write.csv(filtered_nt, "exports/5-D_clonotypes_TCR_nt_sequences_filtered.csv", row.names = FALSE)



#
# Matching aa sequence values??
#
find_matching_values <- function(file_path) {
  
  # Read the CSV file
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Normalize column names in case R has converted "_" to "." (e.g., "CDR3_aa-sequence_1" → "CDR3.aa.sequence.1")
  colnames(df) <- gsub("\\.", "-", colnames(df))  # Replace all dots with underscores
  
  # Identify columns that match the pattern "CDR3_aa-sequence_x" (x = 1 to 100)
  scan_cols <- grep("^CDR3_aa-sequence_[1-9][0-9]?$|^CDR3_aa-sequence_100$", colnames(df), value = TRUE)
  
  if (length(scan_cols) == 0) {
    stop("Error: No 'CDR3_aa-sequence_x' columns found in the dataset.")
  }
  
  # Initialize results dataframe
  results <- data.frame(Match_1_gclone = character(),
                        Match_1_cloneID = character(),
                        Match_2_gclone = character(),
                        Match_2_cloneID = character(),
                        Matching_aa_sequence = character(),
                        stringsAsFactors = FALSE)
  
  # Create a list of all (row, column) pairs in the scanned range
  scanned_cells <- expand.grid(Row = 1:nrow(df), Column = scan_cols, stringsAsFactors = FALSE)
  
  # Iterate over all unique pairs of scanned cells
  for (idx1 in 1:(nrow(scanned_cells) - 1)) {
    for (idx2 in (idx1 + 1):nrow(scanned_cells)) {
      row1 <- scanned_cells$Row[idx1]
      col1 <- scanned_cells$Column[idx1]
      row2 <- scanned_cells$Row[idx2]
      col2 <- scanned_cells$Column[idx2]
      
      val1 <- df[row1, col1]
      val2 <- df[row2, col2]
      
      # Check if both values are non-NA before proceeding
      if (!is.na(val1) && !is.na(val2)) {
        # Check if both values are non-empty and match
        if (val1 != "" && val2 != "" && val1 == val2) {
          results <- rbind(results, data.frame(
            Match_1_gclone = df$GClone_name[row1],
            Match_1_cloneID = df$Clone_ID[row1],
            Match_2_gclone = df$GClone_name[row2],
            Match_2_cloneID = df$Clone_ID[row2],
            Matching_aa_sequence = val1,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Create output filename
  file_name <- tools::file_path_sans_ext(basename(file_path))
  output_file <- paste0(dirname(file_path), "/", file_name, "_AA-SEQUENCE-MATCHES.csv")
  
  # Export results if there are matches
  if (nrow(results) > 0) {
    write.csv(results, output_file, row.names = FALSE)
    message("Results saved to: ", output_file)
  } else {
    message("No matches found. No output file created.")
  }
  
  return(results)
}


find_matching_values("exports/D3_pp65_global_clones_master_list.csv")



#
# Code to count clones between timepoints
#
# Define the clone IDs to include
clone_ids <- paste0("D1_global-clone-", 1:14)

# Remove NA values before counting
df_filtered <- mdata_pp65@meta.data[!is.na(mdata_pp65@meta.data$Global_Clone_ID) & 
                                      !is.na(mdata_pp65@meta.data$sampleBarcodeThresholdLogCpm), ]

# Initialize a data frame to store the counts
D1_pp65_clones_early_vs_late_counts <- data.frame(
  Global_Clone_ID = clone_ids,
  HTO_ET_Counts = numeric(length(clone_ids)),
  HTO_LT_Counts = numeric(length(clone_ids)),
  stringsAsFactors = FALSE
)

# Loop through each clone ID and count cells based on the conditions
for (i in seq_along(clone_ids)) {
  clone <- clone_ids[i]
  
  # Count cells for HTO-ET.pp65.specific.CD8s
  D1_pp65_clones_early_vs_late_counts$HTO_ET_Counts[i] <- sum(
    df_filtered$Global_Clone_ID == clone &
      df_filtered$sampleBarcodeThresholdLogCpm == "HTO-ET.pp65.specific.CD8s"
  )
  
  # Count cells for HTO-LT.pp65.specific.CD8s
  D1_pp65_clones_early_vs_late_counts$HTO_LT_Counts[i] <- sum(
    df_filtered$Global_Clone_ID == clone &
      df_filtered$sampleBarcodeThresholdLogCpm == "HTO-LT.pp65.specific.CD8s"
  )
}

# Export the table to a CSV file
write.csv(D1_pp65_clones_early_vs_late_counts, "exports/integrated_v3/D1_pp65_clones_early-vs-late_counts.csv", row.names = FALSE)

# Print the table to check the output
print(D1_pp65_clones_early_vs_late_counts)



# Create or update the 'timepoint' column with multiple conditions
mdata$timepoint <- ifelse(mdata$sampleBarcodeThresholdLogCpm %in% 
                      c("HTO-ET.pp65.specific.CD8s", "HTO-ET.IE1.specific.CD8s", "HTO-ET.CD8s"), 
                    "Early",
                    ifelse(mdata$sampleBarcodeThresholdLogCpm %in% 
                            c("HTO-LT.pp65.specific.CD8s", "HTO-LT.IE1.specific.CD8s", "HTO-LT.CD8s"), 
                           "Late", 
                      NA))


#
# Stacked bar graphs for clones in each cluster (by cluster)
#
# Define the Global_Clone_IDs to keep
clone_ids <- paste0("D3_global-clone-", 1:17)

# Subset the Seurat object
D3_clones <- subset(mdata_pp65, subset = Global_Clone_ID %in% clone_ids)

# Create a dataframe for plotting
df <- D3_clones@meta.data %>%
  dplyr::count(seurat_clusters, Global_Clone_ID) %>%
  dplyr::mutate(seurat_clusters = as.factor(seurat_clusters))

# Generate the stacked bar plot
plot <- ggplot(df, aes(x = seurat_clusters, y = n, fill = Global_Clone_ID)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic() +
  labs(title = "Distribution of Global Clones Across Clusters - Donor 3",
       x = "Seurat Clusters",
       y = "Cell Count",
       fill = "Global Clone ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Save the plot as a .png file
ggsave("exports/integrated_v3/D3_clones_stacked_bar.png", plot = plot, width = 8, height = 6, dpi = 300)


#
# Stacked bar graphs for clones in each cluster (by clone)
#
# Define the Global_Clone_IDs to keep
clone_ids <- paste0("D2_pp65-clone-", 1:5)

# Subset the Seurat object
D2_clones <- subset(mdata, subset = Global_Clone_ID %in% clone_ids)

# Create a dataframe for plotting
df <- D2_clones@meta.data %>%
  dplyr::count(Global_Clone_ID, seurat_clusters) %>%
  dplyr::group_by(Global_Clone_ID) %>%
  dplyr::mutate(percent = n / sum(n) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(Global_Clone_ID = factor(Global_Clone_ID, levels = clone_ids)) # Ensure order

# Generate the stacked bar plot
plot <- ggplot(df, aes(x = Global_Clone_ID, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic(base_size = 14) +
  labs(title = "Cluster Composition of Global Clones - Donor 2 - pp65",
       x = "Global Clone ID",
       y = "Percentage of Cells",
       fill = "Seurat Cluster") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA))

# Save the plot as a .png file
ggsave("exports/integrated_v5/stacked_bar_graphs/D2_stacked_bar_by_clone.png", plot = plot, width = 10, height = 6, dpi = 300, bg = "white")



#
# Stacked bar graphs for clones in each cluster (split by timepoint)
#
# Define the Global_Clone_IDs to keep
clone_ids <- paste0("D5_pp65-clone-", 1:10)

# Subset the Seurat object
D5_pp65_clones <- subset(mdata_CMV, subset = Global_Clone_ID %in% clone_ids)

# Assign "Early" and "Late" labels based on 'timepoint'
D5_pp65_clones@meta.data <- D5_pp65_clones@meta.data %>%
  dplyr::mutate(Stage = case_when(
    timepoint == "Early" ~ "Early",
    timepoint == "Late" ~ "Late",
    TRUE ~ NA_character_  # Exclude other values
  )) %>%
  dplyr::filter(!is.na(Stage))  # Remove rows without Early/Late classification

# Create a dataframe for plotting
df <- D5_pp65_clones@meta.data %>%
  dplyr::count(Global_Clone_ID, Stage, seurat_clusters) %>%
  dplyr::group_by(Global_Clone_ID, Stage) %>%
  dplyr::mutate(percent = n / sum(n) * 100) %>%
  dplyr::ungroup()

# Create a formatted x-axis label
df <- df %>%
  dplyr::mutate(Clone_Number = as.numeric(str_extract(Global_Clone_ID, "\\d+$"))) %>%
  dplyr::arrange(Clone_Number, Stage) %>%
  dplyr::mutate(Formatted_Clone = factor(paste(Global_Clone_ID, Stage, sep = " - "), 
                                         levels = unique(paste(Global_Clone_ID, Stage, sep = " - "))))

# Generate the stacked bar plot
plot <- ggplot(df, aes(x = Formatted_Clone, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_classic(base_size = 14) +
  labs(title = "Early vs Late Cluster Composition of pp65 Clones - Donor 5",
       x = "Clone ID - Timepoint",
       y = "Percentage of Cells",
       fill = "Seurat Cluster") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.background = element_rect(fill = "white", color = NA),
        plot.background = element_rect(fill = "white", color = NA),
        axis.ticks.x = element_blank()) +
  scale_x_discrete(expand = expansion(mult = c(0.05, 0.05))) +  # Add space for breaks
  geom_vline(xintercept = seq(2.5, length(clone_ids) * 2, by = 2), linetype = "dashed", color = "black")

# Save the plot as a .png file
ggsave("exports/integrated_v5/stacked_bar_graphs/D5_pp65_clones_early_late_stacked_bar.png", plot = plot, width = 12, height = 6, dpi = 300, bg = "white")



#
# Code to create pie charts for clonal cluster frequencies across timepoints,
# with a stability metric comparing Early and Late timepoints,
# and exporting all metrics to a single CSV
#

# Extract metadata
meta <- mdata_CMV@meta.data

# Parse Global_Clone_ID into donor, specificity, clone_number
meta <- meta %>%
  mutate(
    donor = str_extract(Global_Clone_ID, "D\\d") %>% str_remove("D") %>% as.integer(),
    specificity = str_extract(Global_Clone_ID, "(pp65|IE1)"),
    clone_number = str_extract(Global_Clone_ID, "clone-\\d+") %>% str_remove("clone-") %>% as.integer(),
    cluster = metacluster
  )

# Initialize master table
all_stability_metrics <- list()

# Define consistent colors for all clusters
all_clusters <- sort(unique(meta$cluster))
cluster_colors <- setNames(
  RColorBrewer::brewer.pal(n = max(3, min(12, length(all_clusters))), "Set3")[1:length(all_clusters)],
  all_clusters
)

# Loop through donors and specificities
for (d in sort(unique(meta$donor))) {
  for (spec in unique(meta$specificity)) {
    
    df_sub <- meta %>%
      filter(donor == d, specificity == spec)
    
    if (nrow(df_sub) == 0) next
    
    plots <- list()
    
    for (cl in sort(unique(df_sub$clone_number))) {
      pie_early <- pie_late <- NULL
      freq_early <- freq_late <- NULL
      
      for (tp in c("Early", "Late")) {
        df_plot <- df_sub %>%
          filter(clone_number == cl, timepoint == tp) %>%
          count(cluster)
        
        if (nrow(df_plot) == 0 || sum(df_plot$n) == 0) {
          df_plot <- tibble(cluster = NA, freq = 1)
          n_cells <- 0
          p <- ggplot() + theme_void() +
            ggtitle(paste("Clone", cl, "-", tp, "\nNo data"))
        } else {
          df_plot <- df_plot %>%
            mutate(freq = n / sum(n))
          
          n_cells <- sum(df_plot$n)
          
          p <- ggplot(df_plot, aes(x = "", y = freq, fill = factor(cluster, levels = all_clusters))) +
            geom_bar(stat = "identity", width = 1) +
            coord_polar("y") +
            scale_fill_manual(values = cluster_colors, drop = FALSE) +  # <- consistent colors
            theme_void() +
            theme(legend.position = "none") +
            ggtitle(paste("Clone", cl, "-", tp, "\n(n =", n_cells, ")")) +
            theme(plot.title = element_text(size = 10))
        }
        
        p <- p + theme(
          plot.margin = margin(0, 0, 0, 0),
          panel.spacing = unit(0, "lines")
        )
        
        if (tp == "Early") {
          pie_early <- p
          freq_early <- df_plot %>% select(cluster, freq)
        } else {
          pie_late <- p
          freq_late <- df_plot %>% select(cluster, freq)
        }
      }
      
      # Compute stability info
      freq_all <- full_join(freq_early, freq_late, by = "cluster", suffix = c("_early", "_late")) %>%
        replace_na(list(freq_early = 0, freq_late = 0)) %>%
        mutate(diff = abs(freq_early - freq_late))
      
      stability <- sum(freq_all$diff)
      
      n_cells_early <- df_sub %>% filter(clone_number == cl, timepoint == "Early") %>% nrow()
      n_cells_late  <- df_sub %>% filter(clone_number == cl, timepoint == "Late") %>% nrow()
      
      label <- paste("Donor", d, "-", spec, "- clone", cl)
      
      # Save summary to master list
      all_stability_metrics[[label]] <- list(
        donor = d,
        specificity = spec,
        clone_number = cl,
        stability_metric = stability,
        n_cells_early = n_cells_early,
        n_cells_late = n_cells_late
      )
      
      # Write per-clone CSV
      breakdown <- freq_all %>%
        arrange(cluster) %>%
        mutate(
          stability_metric = stability,
          n_cells_early = n_cells_early,
          n_cells_late = n_cells_late
        ) %>%
        select(
          cluster,
          freq_early,
          freq_late,
          diff,
          stability_metric,
          n_cells_early,
          n_cells_late
        )
      
      write_csv(
        breakdown,
        file = paste0(
          "exports/integrated_v5/clone_cluster_pie_charts/cluster_stability_breakdowns/",
          "Donor", d, "_", spec, "_clone", cl, "_stability_breakdown.csv"
        )
      )
      
      # Add stability metric plot
      metric_plot <- ggplot() +
        annotate("text", x = 0.5, y = 0.5, label = sprintf("Stability:\n%.2f", stability), size = 4) +
        theme_void()
      
      # Combine plots
      combined <- pie_early + pie_late + metric_plot + plot_layout(ncol = 3, widths = c(1, 1, 0.7))
      plots[[paste0("clone", cl)]] <- combined
    }
    
    # Stack all clone plots vertically
    final_plot <- wrap_plots(plots, ncol = 1)
    
    final_plot <- final_plot + plot_annotation(
      title = paste("Donor", d, "-", spec),
      theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
    )
    
    ggsave(
      filename = paste0("exports/integrated_v5/clone_cluster_pie_charts/Donor", d, "_", spec, "_clone_pies.png"),
      plot = final_plot,
      width = 10,
      height = 2 * length(plots)
    )
    
    print(paste("Saved plot for Donor", d, "Specificity", spec))
  }
}

# ----- Export legend -----

# Generate a dummy plot to extract the legend
legend_df <- tibble(
  cluster = factor(all_clusters, levels = all_clusters),
  freq = 1  # dummy values
)

# White background theme applied
legend_plot <- ggplot(legend_df, aes(x = "", y = freq, fill = cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  guides(fill = guide_legend(title = "Cluster")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),  # <- white background
    plot.background = element_rect(fill = "white", color = NA)     # <- white outer area
  )

# Function to extract legend from ggplot object
extract_legend <- function(p) {
  g <- ggplotGrob(p)
  legends <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[legends]]
}

legend_only <- extract_legend(legend_plot)

# Save legend as .png
ggsave(
  filename = "exports/integrated_v5/clone_cluster_pie_charts/cluster_color_legend.png",
  plot = as.ggplot(legend_only),
  width = 3,
  height = length(all_clusters) * 0.3 + 0.5
)

# ----- Export combined stability metrics -----

stability_df <- tibble(
  label = names(all_stability_metrics),
  donor = sapply(all_stability_metrics, function(x) x$donor),
  specificity = sapply(all_stability_metrics, function(x) x$specificity),
  clone_number = sapply(all_stability_metrics, function(x) x$clone_number),
  stability_metric = sapply(all_stability_metrics, function(x) x$stability_metric),
  n_cells_early = sapply(all_stability_metrics, function(x) x$n_cells_early),
  n_cells_late = sapply(all_stability_metrics, function(x) x$n_cells_late)
)

write_csv(
  stability_df,
  file = "exports/integrated_v5/clone_cluster_pie_charts/all_stability_metrics.csv"
)

#
# Pie charts for both timpoints combined
#

# Define export directory
export_dir <- "exports/integrated_v5/clone_cluster_pie_charts/both_timepoints_combined"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# Create subdirectory for per-clone CSVs
csv_dir <- file.path(export_dir, "cluster_breakdowns_combined_timepoints")
dir.create(csv_dir, showWarnings = FALSE)

# Extract metadata
meta <- mdata_CMV@meta.data

# Parse Global_Clone_ID into donor, specificity, clone_number
meta <- meta %>%
  mutate(
    donor = str_extract(Global_Clone_ID, "D\\d") %>% str_remove("D") %>% as.integer(),
    specificity = str_extract(Global_Clone_ID, "(pp65|IE1)"),
    clone_number = str_extract(Global_Clone_ID, "clone-\\d+") %>% str_remove("clone-") %>% as.integer(),
    cluster = metacluster
  )

# Initialize export table
all_clone_metrics <- list()

# Define consistent colors for all clusters
all_clusters <- sort(unique(meta$cluster))
cluster_colors <- setNames(
  RColorBrewer::brewer.pal(n = max(3, min(12, length(all_clusters))), "Set3")[1:length(all_clusters)],
  all_clusters
)

# Loop through donors and specificities
for (d in sort(unique(meta$donor))) {
  for (spec in unique(meta$specificity)) {
    
    df_sub <- meta %>%
      filter(donor == d, specificity == spec)
    
    if (nrow(df_sub) == 0) next
    
    plots <- list()
    
    for (cl in sort(unique(df_sub$clone_number))) {
      
      df_plot <- df_sub %>%
        filter(clone_number == cl) %>%
        count(cluster)
      
      if (nrow(df_plot) == 0 || sum(df_plot$n) == 0) {
        df_plot <- tibble(cluster = NA, freq = 1)
        n_cells <- 0
        p <- ggplot() + theme_void() +
          ggtitle(paste("Clone", cl, "\nNo data"))
      } else {
        df_plot <- df_plot %>%
          mutate(freq = n / sum(n))
        
        n_cells <- sum(df_plot$n)
        
        p <- ggplot(df_plot, aes(x = "", y = freq, fill = factor(cluster, levels = all_clusters))) +
          geom_bar(stat = "identity", width = 1) +
          coord_polar("y") +
          scale_fill_manual(values = cluster_colors, drop = FALSE) +
          theme_void() +
          theme(legend.position = "none") +
          ggtitle(paste("Clone", cl, "\n(n =", n_cells, ")")) +
          theme(plot.title = element_text(size = 10))
      }
      
      p <- p + theme(
        plot.margin = margin(0, 0, 0, 0),
        panel.spacing = unit(0, "lines")
      )
      
      plots[[paste0("clone", cl)]] <- p
      
      # Save summary for export
      summary_tbl <- df_plot %>%
        arrange(cluster) %>%
        mutate(
          donor = d,
          specificity = spec,
          clone_number = cl,
          n_cells_total = n_cells
        )
      
      all_clone_metrics[[paste0("Donor", d, "_", spec, "_clone", cl)]] <- summary_tbl
      
      # Write per-clone breakdown
      write_csv(
        summary_tbl,
        file = file.path(csv_dir, paste0("Donor", d, "_", spec, "_clone", cl, "_combined_breakdown.csv"))
      )
    }
    
    # Stack all clone plots vertically
    final_plot <- wrap_plots(plots, ncol = 1)
    
    final_plot <- final_plot + plot_annotation(
      title = paste("Donor", d, "-", spec, "(Combined Timepoints)"),
      theme = theme(plot.title = element_text(size = 14, hjust = 0.5))
    )
    
    ggsave(
      filename = file.path(export_dir, paste0("Donor", d, "_", spec, "_combined_pies.png")),
      plot = final_plot,
      width = 6,
      height = 2 * length(plots)
    )
    
    print(paste("Saved combined pie chart for Donor", d, "Specificity", spec))
  }
}

# ----- Export legend -----

legend_df <- tibble(
  cluster = factor(all_clusters, levels = all_clusters),
  freq = 1
)

legend_plot <- ggplot(legend_df, aes(x = "", y = freq, fill = cluster)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cluster_colors) +
  guides(fill = guide_legend(title = "Cluster")) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

extract_legend <- function(p) {
  g <- ggplotGrob(p)
  legends <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[legends]]
}

legend_only <- extract_legend(legend_plot)

ggsave(
  filename = file.path(export_dir, "cluster_color_legend_combined_timepoints.png"),
  plot = as.ggplot(legend_only),
  width = 3,
  height = length(all_clusters) * 0.3 + 0.5
)

# ----- Export combined metrics -----

all_metrics_df <- bind_rows(all_clone_metrics, .id = "label")

write_csv(
  all_metrics_df,
  file = file.path(export_dir, "all_clone_metrics_combined_timepoints.csv")
)


#
# Stacked bar graphs for metacluster breakdown (per clone per donor)
#

# Define export directory
export_dir <- "exports/integrated_v5/clone_cluster_stacked_bars/corrected_colors/stacked_bars_thin"

# Extract metadata and parse fields
meta <- mdata_CMV@meta.data %>%
  mutate(
    donor = str_extract(Global_Clone_ID, "D\\d") %>% str_remove("D") %>% as.integer(),
    specificity = str_extract(Global_Clone_ID, "(pp65|IE1)"),
    clone_number = str_extract(Global_Clone_ID, "clone-\\d+") %>% str_remove("clone-") %>% as.integer(),
    cluster = metacluster
  )

# -------------------------------------------------------------------------
# Manual cluster colors (sampled from uploaded UMAP image legend)
# -------------------------------------------------------------------------
all_clusters <- sort(unique(meta$cluster))

# Manually assign the desired colors (hex values sampled from your image)
# Make sure the names on the left match the metacluster labels in your data.
manual_colors <- c(
  "Tcm"     = "#F8766D",
  "Tem"     = "#A3A500",
  "T-CMC1"  = "#00BF7D",
  "NK-like" = "#00B0F6",
  "ISG+"    = "#E76BF3"
)

# Keep only colors matching clusters present in this dataset and preserve order
cluster_colors <- manual_colors[intersect(names(manual_colors), all_clusters)]

# Reorder to match all_clusters exactly
cluster_colors <- cluster_colors[all_clusters]
# -------------------------------------------------------------------------

# Loop through donors
for (d in sort(unique(meta$donor))) {
  donor_data <- meta %>% filter(donor == d)
  
  # Summarize clone sizes
  clone_summary <- donor_data %>%
    count(specificity, clone_number, name = "n_cells_clone") %>%
    group_by(specificity) %>%
    mutate(clone_size_percent = n_cells_clone / sum(n_cells_clone) * 100) %>%
    ungroup()
  
  # Cluster composition per clone
  comp_data <- donor_data %>%
    count(specificity, clone_number, cluster, name = "n_cells") %>%
    group_by(specificity, clone_number) %>%
    mutate(freq = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    left_join(clone_summary, by = c("specificity", "clone_number")) %>%
    mutate(
      # Use formatted clone labels for plotting
      clone_label_raw = paste0(specificity, "_clone", clone_number),
      clone_label = paste0(specificity, " Clone ", clone_number),
      label_color = dplyr::case_when(
        specificity == "IE1"  ~ "red",
        specificity == "pp65" ~ "blue",
        TRUE                  ~ "black"
      )
    )
  
  # Sort clone labels
  ie1_labels <- comp_data %>%
    filter(specificity == "IE1") %>%
    distinct(clone_number, clone_label) %>%
    arrange(clone_number) %>%
    pull(clone_label)
  
  pp65_labels <- comp_data %>%
    filter(specificity == "pp65") %>%
    distinct(clone_number, clone_label) %>%
    arrange(clone_number) %>%
    pull(clone_label)
  
  # Create factor levels with gap if both groups exist
  if (length(ie1_labels) > 0 && length(pp65_labels) > 0) {
    levels_vec <- c(rev(pp65_labels), "gap_row", rev(ie1_labels))
    gap_row <- tibble(
      specificity = NA_character_,
      clone_number = NA_integer_,
      cluster = NA_character_,
      n_cells = 0,
      n_cells_clone = 0,
      clone_size_percent = NA_real_,
      freq = 0,
      clone_label = "gap_row",
      label_color = NA_character_
    )
    comp_data <- bind_rows(comp_data, gap_row)
    draw_gap_line <- TRUE
  } else if (length(ie1_labels) > 0) {
    levels_vec <- rev(ie1_labels)
    draw_gap_line <- FALSE
  } else if (length(pp65_labels) > 0) {
    levels_vec <- rev(pp65_labels)
    draw_gap_line <- FALSE
  } else {
    message("No IE1 or pp65 clones for donor ", d, " — skipping.")
    next
  }
  
  # Apply factor levels
  comp_data <- comp_data %>%
    mutate(clone_label = factor(as.character(clone_label), levels = levels_vec))
  
  # Prepare labels for plotting
  y_labels_df <- tibble(clone_label = levels_vec) %>%
    mutate(
      label_text = ifelse(clone_label == "gap_row", "", clone_label),
      label_color = sapply(clone_label, function(l) {
        v <- comp_data$label_color[as.character(comp_data$clone_label) == l][1]
        if (is.null(v) || is.na(v)) "black" else v
      })
    )
  
  # Percent labels per clone
  percent_text_df <- comp_data %>%
    filter(!is.na(specificity)) %>%
    group_by(clone_label) %>%
    slice(1) %>%
    ungroup()
  
  # Plot
  bar_plot <- ggplot(comp_data, aes(x = freq, y = clone_label, fill = factor(cluster, levels = all_clusters))) +
    geom_bar(stat = "identity", width = 0.8) +
    
    # Percent labels to the right
    geom_text(
      data = percent_text_df,
      aes(x = 1.02, y = clone_label, label = sprintf("%.1f%%", clone_size_percent), color = label_color),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3
    ) +
    
    # Clone labels at start of bar, shifted left
    geom_text(
      data = y_labels_df,
      aes(x = 0, y = clone_label, label = label_text, color = label_color),
      inherit.aes = FALSE,
      hjust = 1.12,
      size = 3
    ) +
    
    scale_fill_manual(values = cluster_colors, drop = FALSE) +
    scale_color_identity() +
    scale_x_continuous(
      expand = expansion(mult = c(0.02, 0.15)),
      labels = scales::percent_format(accuracy = 1)
    ) +
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_line(color = "black"),
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.margin = margin(10, 30, 10, 70)
    ) +
    guides(fill = guide_legend(title = "Metacluster")) +
    ggtitle(paste("Donor", d, "- Clone Composition by Cluster (Combined Timepoints)"))
  
  # Dashed line in gap row
  if (draw_gap_line) {
    gap_index <- which(levels_vec == "gap_row")
    bar_plot <- bar_plot +
      geom_hline(yintercept = gap_index, linetype = "dashed", color = "gray50", size = 0.4)
  }
  
  # Save plot
  ggsave(
    filename = file.path(export_dir, paste0("Donor", d, "_stacked_bars.png")),
    plot = bar_plot,
    width = 4.75,
    height = 0.3 * length(levels_vec) + 2
  )
  
  message("Saved bar plot for Donor ", d)
}





# Export directory
export_dir <- "exports/integrated_v5/clone_cluster_stacked_bars/corrected_colors/ET_vs_LT_single_graph_thin"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# Extract metadata
meta <- mdata_CMV@meta.data %>%
  mutate(
    donor = str_extract(Global_Clone_ID, "D\\d") %>% str_remove("D") %>% as.integer(),
    specificity = str_extract(Global_Clone_ID, "(pp65|IE1)"),
    clone_number = str_extract(Global_Clone_ID, "clone-\\d+") %>% str_remove("clone-") %>% as.integer(),
    cluster = metacluster
  )

# -------------------------------------------------------------------------
# Manual cluster colors (sampled from uploaded UMAP image legend)
# -------------------------------------------------------------------------
all_clusters <- sort(unique(meta$cluster))

# Manually assign the desired colors (hex values sampled from your image)
# Make sure the names on the left match the metacluster labels in your data.
manual_colors <- c(
  "Tcm"     = "#F8766D",
  "Tem"     = "#A3A500",
  "T-CMC1"  = "#00BF7D",
  "NK-like" = "#00B0F6",
  "ISG+"    = "#E76BF3"
)

# Keep only colors matching clusters present in this dataset and preserve order
cluster_colors <- manual_colors[intersect(names(manual_colors), all_clusters)]

# Reorder to match all_clusters exactly
cluster_colors <- cluster_colors[all_clusters]
# -------------------------------------------------------------------------

# Function to generate combined Early/Late plot per donor
generate_combined_plot <- function(donor_meta, d) {
  if (nrow(donor_meta) == 0) return(NULL)
  
  # Summarize cluster composition per clone/timepoint
  comp_data <- donor_meta %>%
    count(specificity, clone_number, timepoint, cluster, name = "n_cells") %>%
    group_by(specificity, clone_number, timepoint) %>%
    mutate(freq = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    mutate(
      clone_label_base = paste0(specificity, " Clone ", clone_number),
      y_label_text = paste0(clone_label_base, " - ", timepoint),
      label_color = case_when(
        specificity == "IE1" ~ "red",
        specificity == "pp65" ~ "blue",
        TRUE ~ "black"
      )
    )
  
  # Compute percent of each clone out of all cells of same specificity/timepoint
  clone_percent <- donor_meta %>%
    count(specificity, clone_number, timepoint, name = "clone_cells") %>%
    group_by(specificity, timepoint) %>%
    mutate(clone_percent = clone_cells / sum(clone_cells) * 100) %>%
    ungroup() %>%
    mutate(y_label_text = paste0(specificity, " Clone ", clone_number, " - ", timepoint))
  
  # Build factor levels with Early/Late + small gaps + full gap between IE1/pp65
  y_levels <- c()
  dummy_rows <- tibble()
  
  for (spec in c("IE1", "pp65")) {
    clones <- sort(unique(comp_data$clone_number[comp_data$specificity == spec]))
    for (i in seq_along(clones)) {
      cl <- clones[i]
      # Early/Late bars
      y_levels <- c(
        y_levels,
        paste0(spec, " Clone ", cl, " - Early"),
        paste0(spec, " Clone ", cl, " - Late")
      )
      # Small gap after each clone except last in group
      if (i < length(clones)) {
        gap_label <- paste0("gap_", spec, "_", cl)
        y_levels <- c(y_levels, gap_label)
        dummy_rows <- bind_rows(dummy_rows,
                                tibble(
                                  specificity = NA_character_,
                                  clone_number = NA_integer_,
                                  timepoint = NA_character_,
                                  cluster = NA_character_,
                                  n_cells = 0,
                                  freq = 0,
                                  clone_label_base = NA_character_,
                                  y_label_text = gap_label,
                                  label_color = NA_character_
                                ))
      }
    }
    # Full gap between IE1 and pp65 (two dummy rows)
    if (spec == "IE1") {
      full_gap_label1 <- "gap_IE1_pp65_1"
      full_gap_label2 <- "gap_IE1_pp65_2"
      y_levels <- c(y_levels, full_gap_label1, full_gap_label2)
      
      dummy_rows <- bind_rows(dummy_rows,
                              tibble(
                                specificity = NA_character_,
                                clone_number = NA_integer_,
                                timepoint = NA_character_,
                                cluster = NA_character_,
                                n_cells = 0,
                                freq = 0,
                                clone_label_base = NA_character_,
                                y_label_text = full_gap_label1,
                                label_color = NA_character_
                              ),
                              tibble(
                                specificity = NA_character_,
                                clone_number = NA_integer_,
                                timepoint = NA_character_,
                                cluster = NA_character_,
                                n_cells = 0,
                                freq = 0,
                                clone_label_base = NA_character_,
                                y_label_text = full_gap_label2,
                                label_color = NA_character_
                              ))
    }
  }
  
  # Combine actual data + dummy rows
  full_plot_data <- bind_rows(comp_data, dummy_rows) %>%
    mutate(y_label_text = factor(y_label_text, levels = rev(y_levels)))
  
  # Dashed line: middle of full gap
  full_gap_pos <- which(levels(full_plot_data$y_label_text) %in% c("gap_IE1_pp65_1", "gap_IE1_pp65_2"))
  dashed_line_position <- mean(full_gap_pos)
  
  # Percent labels (corrected per specificity/timepoint)
  percent_text_df <- clone_percent %>%
    select(y_label_text, clone_percent, specificity) %>%
    distinct()
  
  # Create plot
  ggplot(full_plot_data, aes(x = freq, y = y_label_text, fill = factor(cluster, levels = all_clusters))) +
    geom_bar(stat = "identity", width = 0.8) +
    
    # Percent labels on right
    geom_text(
      data = percent_text_df,
      aes(x = 1.02, y = y_label_text, label = sprintf("%.1f%%", clone_percent), color = case_when(specificity=="IE1"~"red", specificity=="pp65"~"blue", TRUE~"black")),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3
    ) +
    
    # Clone labels on left
    geom_text(
      data = comp_data %>% distinct(y_label_text, label_color),
      aes(x = 0, y = y_label_text, label = y_label_text, color = label_color),
      inherit.aes = FALSE,
      hjust = 1.1,  # reduced from 1.25
      size = 3
    ) +
    
    # Fill + color scales
    scale_fill_manual(values = cluster_colors, drop = FALSE) +
    scale_color_identity() +
    
    # X-axis expansion for percent labels
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.15)),
                       labels = scales::percent_format(accuracy = 1)) +
    
    # Remove y-axis ticks and text
    coord_cartesian(clip = "off") +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_line(color = "black"),
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.margin = margin(10, 30, 10, 100)
    ) +
    
    guides(fill = guide_legend(title = "Metacluster")) +
    
    geom_hline(yintercept = dashed_line_position, linetype = "dashed", color = "gray50", size = 0.4) +
    
    ggtitle(paste("Donor", d, "- Clone Composition by Cluster (Early vs Late)"))
}

# Loop through donors
for (d in sort(unique(meta$donor))) {
  donor_meta <- meta %>% filter(donor == d)
  
  combined_plot <- generate_combined_plot(donor_meta, d)
  
  if (!is.null(combined_plot)) {
    ggsave(
      filename = file.path(export_dir, paste0("Donor", d, "_ET_vs_LT_single_graph.png")),
      plot = combined_plot,
      width = 4.75, # normal = 6
      height = 0.25 * length(unique(donor_meta$Global_Clone_ID)) + 2
    )
    message("Saved combined Early/Late plot for Donor ", d)
  }
}




#
# Gene set scores by cluster
#
# Extract scores and cluster identity
df <- mdata_CMV@meta.data %>%
  select(seurat_clusters, Tcm_Score1, Tem_Score1, NKlike_Score1, ISG_Score1) %>%
  rename(Tcm = Tcm_Score1, Tem = Tem_Score1, `NK-like` = NKlike_Score1, ISG = ISG_Score1) %>%
  melt(id.vars = "seurat_clusters", variable.name = "Score_Type", value.name = "Score")

# Ensure factor order for proper facet stacking
df$Score_Type <- factor(df$Score_Type, levels = c("Tcm", "Tem", "NK-like", "ISG"))

# Plot
p <- ggplot(df, aes(x = factor(seurat_clusters), y = Score, fill = factor(seurat_clusters))) +
  geom_violin(scale = "width", trim = FALSE, color = "black") +  # black outline
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.1) +             # individual points
  facet_wrap(~ Score_Type, ncol = 1, scales = "free_y", strip.position = "left") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  ) +
  labs(
    title = "Stacked Violin Plots of Gene Set Scores by Seurat Cluster",
    x = "Seurat Cluster",
    y = "Score"
  )

# Save the plot
ggsave("exports/integrated_v5/gene_set_scores_by_seurat_cluster.png", plot = p, width = 10, height = 8)


#
# Stacked violin plots for gene expression
#
# List of genes
genes_of_interest <- c("TCF7", "GZMB", "CMC1", "KLRC2", "IFIT1")

# Extract expression values for these genes
expr_data <- FetchData(mdata_CMV, vars = genes_of_interest)
expr_data$seurat_clusters <- mdata_CMV@meta.data$seurat_clusters

# Reshape to long format
df_long <- reshape2::melt(expr_data, id.vars = "seurat_clusters", variable.name = "Gene", value.name = "Expression")

# Ensure gene order
df_long$Gene <- factor(df_long$Gene, levels = genes_of_interest)

# Optional: convert clusters to factor for consistent ordering
df_long$seurat_clusters <- factor(df_long$seurat_clusters)

# Plot
p <- ggplot(df_long, aes(x = seurat_clusters, y = Expression, fill = seurat_clusters)) +
  geom_violin(scale = "width", trim = FALSE, color = "black") +
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.1) +  # smaller points
  facet_wrap(~ Gene, ncol = 1, scales = "free_y", strip.position = "left") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  ) +
  labs(
    title = "Stacked Violin Plots of Gene Expression by Seurat Cluster",
    x = "Seurat Cluster",
    y = "Expression Level"
  )

# Save the plot
ggsave("exports/integrated_v5/metacluster_gene_expression_by_seurat_cluster.png", plot = p, width = 10, height = 8)


#
# Stacked violin plots for KLRC2/CD160 gene expression, per donor, per clone
#
# List of genes
genes_of_interest <- c("KLRC2", "CD160")

# Loop over donors D1 through D5
for (i in 1:5) {
  donor_label <- paste0("Donor ", i)
  donor_id <- paste0("D", i)
  
  # Subset to cells from this donor
  donor_subset <- subset(mdata_CMV, subset = donor_identity == donor_label)
  
  # Skip if no relevant cells
  if (ncol(donor_subset) == 0) {
    message(paste("Skipping", donor_label, "- no matching cells"))
    next
  }
  
  # Extract expression values and Global_Clone_ID for donor cells
  expr_data <- FetchData(donor_subset, vars = genes_of_interest)
  expr_data$Global_Clone_ID <- donor_subset@meta.data$Global_Clone_ID
  
  # Remove NA Global_Clone_IDs
  expr_data <- expr_data %>% filter(!is.na(Global_Clone_ID))
  
  # Skip if no clones with expression data
  if (nrow(expr_data) == 0) {
    message(paste("Skipping", donor_label, "- no clone data after NA filter"))
    next
  }
  
  # Melt for ggplot
  df <- melt(expr_data, id.vars = "Global_Clone_ID", variable.name = "Gene", value.name = "Expression")
  
  # Order clones: IE1 first, then pp65
  unique_ids <- sort(unique(df$Global_Clone_ID))
  ie1_ids <- unique_ids[grepl("IE1", unique_ids)]
  pp65_ids <- unique_ids[grepl("pp65", unique_ids)]
  ordered_ids <- c(ie1_ids, pp65_ids)
  df$Global_Clone_ID <- factor(df$Global_Clone_ID, levels = ordered_ids)
  
  # Determine separator line location
  separator_x <- length(ie1_ids) + 0.5
  
  # Plot
  p <- ggplot(df, aes(x = Global_Clone_ID, y = Expression, fill = Global_Clone_ID)) +
    geom_violin(scale = "width", trim = FALSE, color = "black") +
    geom_jitter(width = 0.2, size = 0.2, alpha = 0.5) +
    geom_vline(xintercept = separator_x, linetype = "dashed", color = "black", size = 0.5) +
    facet_wrap(~ Gene, ncol = 1, scales = "free_y", strip.position = "left") +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.placement = "outside",
      strip.background = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = "none"
    ) +
    labs(title = paste("Expression of KLRC2/CD160 -", donor_label),
         x = "Global Clone ID", y = "Expression")
  
  # Export
  output_path <- paste0("exports/integrated_v5/KLRC2_CD160_expression_by_clone/", donor_id, "_clone_expression_violin.png")
  ggsave(output_path, plot = p, width = 10, height = 8)
  
  message(paste("Saved expression plot for", donor_label, "to", output_path))
}



#
# D5 bulk clones - KLRC2 expression?
# 

# Subset
clone_subset <- subset(mdata, subset = cloneId %in% c('clone48donor5', 'clone382donor5'))

# List of genes
genes_of_interest <- c("GZMB", "GZMK", "CMC1", "KLRC2")

# Extract expression values for these genes
expr_data <- FetchData(clone_subset, vars = genes_of_interest)
expr_data$cloneId <- clone_subset@meta.data$cloneId

# Reshape to long format
df_long <- melt(expr_data, id.vars = "cloneId", variable.name = "Gene", value.name = "Expression")

# Ensure gene order
df_long$Gene <- factor(df_long$Gene, levels = genes_of_interest)

# Optional: convert clusters to factor for consistent ordering
df_long$cloneId <- factor(df_long$cloneId)

# Plot
p <- ggplot(df_long, aes(x = cloneId, y = Expression, fill = cloneId)) +
  geom_violin(scale = "width", trim = FALSE, color = "black") +
  geom_jitter(width = 0.2, size = 0.2, alpha = 0.5) +  # smaller points
  facet_wrap(~ Gene, ncol = 1, scales = "free_y", strip.position = "left") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.placement = "outside",
    strip.background = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    legend.position = "none"
  ) +
  labs(
    title = "Stacked Violin Plots of Gene Expression for D5 'Bulk' clones",
    x = "cloneId",
    y = "Expression Level"
  )

# Save the plot
ggsave("exports/integrated_v5/D5_Bulk_clones_gene_expression.png", plot = p, width = 10, height = 8)




#
# Metrics for each clone (including percent dominance)
#

# Initialize empty list to collect all data
all_metrics_list <- list()

# Loop over donors D1 through D5
for (i in 1:5) {
  donor_label <- paste0("Donor ", i)
  donor_id <- paste0("D", i)
  
  # Subset to cells from this donor
  donor_subset <- subset(mdata_CMV, subset = donor_identity == donor_label)
  
  # Skip if no relevant cells
  if (ncol(donor_subset) == 0) {
    message(paste("Skipping", donor_label, "- no matching cells"))
    next
  }
  
  # Prepare metadata
  df <- donor_subset@meta.data %>%
    select(Global_Clone_ID, Tcm_Score1, Tem_Score1, NKlike_Score1) %>%
    filter(!is.na(Global_Clone_ID)) %>%
    rename(Tcm = Tcm_Score1, Tem = Tem_Score1, `NK-like` = NKlike_Score1)
  
  for (specificity in c("pp65", "IE1")) {
    # Filter by specificity
    df_spec <- df %>%
      filter(grepl(specificity, Global_Clone_ID, ignore.case = TRUE))
    
    # Skip if none
    if (nrow(df_spec) == 0) {
      message(paste("Skipping", donor_id, specificity, "- no matching Global_Clone_IDs"))
      next
    }
    
    # Total number of relevant cells
    total_cells <- nrow(df_spec)
    
    # Calculate metrics per clone
    clone_metrics <- df_spec %>%
      group_by(Global_Clone_ID) %>%
      summarise(
        Clone_Size = n(),
        Clone_Percentage = Clone_Size / total_cells,
        Mean_Tcm = mean(Tcm, na.rm = TRUE),
        Mean_Tem = mean(Tem, na.rm = TRUE),
        Mean_NKlike = mean(`NK-like`, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Donor = donor_id, Specificity = specificity)
    
    # Append to the full list
    all_metrics_list[[length(all_metrics_list) + 1]] <- clone_metrics
  }
}

# Combine all donor/specificity data into one data frame
all_metrics_df <- bind_rows(all_metrics_list)

# Export to one CSV
final_output_path <- "exports/integrated_v5/gene_set_score_analyses/All_donors_clone_percent_metrics.csv"
fwrite(all_metrics_df, file = final_output_path)
message(paste("Final combined table saved to:", final_output_path))



#
# Add 5 metacluster categories
#

mdata_CMV$metacluster <- case_when(
  mdata_CMV$seurat_clusters %in% c(3, 4, 7)     ~ "Tcm",
  mdata_CMV$seurat_clusters %in% c(0, 2, 8)     ~ "Tem",
  mdata_CMV$seurat_clusters == 1                ~ "T-CMC1",
  mdata_CMV$seurat_clusters == 5                ~ "NK-like",
  mdata_CMV$seurat_clusters == 6                ~ "ISG+",
  TRUE                                  ~ NA_character_  # Catch any other clusters
)

DimPlot(mdata_CMV, reduction = "umap", group.by = "metacluster")

FeaturePlot(mdata_CMV, features = "CMC1", reduction = "umap",
            cols = c("lightgrey", "darkred"))
Idents(mdata_CMV) <- "seurat_clusters"
DimPlot(mdata_CMV, reduction = "umap", 
        cells.highlight = WhichCells(mdata_CMV, idents = 7), 
        cols.highlight = "red", cols = "lightgray")



mdata$UMAP_2 <- Embeddings(mdata, "umap")[,2]

mdata$metacluster <- case_when(
  mdata$seurat_clusters %in% c(0, 3, 10)     ~ "Effector Memory",
  mdata$seurat_clusters %in% c(8, 9)     ~ "Central Memory",
  mdata$seurat_clusters %in% c(1, 5)     ~ "Naive",
  mdata$seurat_clusters == 6                ~ "MAIT",
  mdata$seurat_clusters == 4                ~ "NK-like",
  mdata$seurat_clusters == 7                ~ "ISG+",
  
  mdata$seurat_clusters %in% c(2) & mdata$UMAP_2 > -5   ~ "Central Memory",
  mdata$seurat_clusters %in% c(2) & mdata$UMAP_2 <= -5  ~ "Naive",
  
  TRUE                                                        ~ NA_character_  # Catch any other clusters
)

Idents(mdata) <- "seurat_clusters"
DimPlot(mdata, reduction = "umap", 
        cells.highlight = WhichCells(mdata, idents = 10), 
        cols.highlight = "red", cols = "lightgray")

DimPlot(mdata, reduction = "umap", group.by = "metacluster")



mdata_CMV$UMAP_2 <- Embeddings(mdata_CMV, "umap")[,2]

mdata_CMV$boot_cluster <- case_when(
  mdata_CMV$UMAP_2 < -7.5   ~ "Y",
  mdata_CMV$UMAP_2 >= -7.5   ~ "N",
  
  TRUE                                                        ~ NA_character_  # Catch any other clusters
)

Idents(mdata_CMV) <- "seurat_clusters"
DimPlot(mdata_CMV, reduction = "umap", 
        cells.highlight = WhichCells(mdata_CMV, idents = 10), 
        cols.highlight = "red", cols = "lightgray")

DimPlot(mdata_CMV, reduction = "umap", group.by = "boot_cluster")


#
# Pseudobulk analysis
#

# Step 1: Extract expression matrix and metadata
expr_matrix <- GetAssayData(mdata_CMV, slot = "data")  # or "counts"
meta_data <- mdata_CMV@meta.data

# Step 2: Filter out NA clone IDs and NA timepoints
meta_data <- meta_data %>% filter(!is.na(Global_Clone_ID), !is.na(timepoint))

# Step 3: Identify all (clone, timepoint) combinations
clone_time_combos <- meta_data %>%
  distinct(Global_Clone_ID, timepoint) %>%
  mutate(pseudo_label = paste0(Global_Clone_ID, "_", timepoint))

# Initialize lists
pseudobulk_expr_list <- list()
pseudobulk_meta_list <- list()

for (i in seq_len(nrow(clone_time_combos))) {
  clone_id <- clone_time_combos$Global_Clone_ID[i]
  timepoint <- clone_time_combos$timepoint[i]
  label <- clone_time_combos$pseudo_label[i]
  
  # Get cells matching this clone and timepoint
  cells <- rownames(meta_data)[meta_data$Global_Clone_ID == clone_id & meta_data$timepoint == timepoint]
  
  if (length(cells) == 0) next  # skip if none
  
  # Expression averaging
  expr_subset <- expr_matrix[, cells, drop = FALSE]
  avg_expr <- Matrix::rowMeans(expr_subset)
  
  # Metadata averaging
  meta_subset <- meta_data[cells, c("Tcm_Score1", "Tem_Score1", "NKlike_Score1", "ISG_Score1"), drop = FALSE]
  avg_meta <- colMeans(meta_subset, na.rm = TRUE)
  
  # Add timepoint and original clone info
  avg_meta <- c(avg_meta,
                Global_Clone_ID = clone_id,
                timepoint = timepoint,
                pseudo_label = label,
                cell_count = length(cells))
  
  # Store
  pseudobulk_expr_list[[label]] <- avg_expr
  pseudobulk_meta_list[[label]] <- avg_meta
}

# Combine expression matrix and metadata
pseudobulk_expr_matrix <- do.call(cbind, pseudobulk_expr_list)
pseudobulk_expr_matrix <- as.matrix(pseudobulk_expr_matrix)

pseudobulk_metadata <- do.call(rbind, pseudobulk_meta_list)
pseudobulk_metadata <- as.data.frame(pseudobulk_metadata)
rownames(pseudobulk_metadata) <- colnames(pseudobulk_expr_matrix)

# Create Seurat object
pseudobulk <- CreateSeuratObject(counts = pseudobulk_expr_matrix, meta.data = pseudobulk_metadata)

# Step 6: Average UMAP coordinates per (Global_Clone_ID, timepoint)
umap_coords <- Embeddings(mdata_CMV, reduction = "umap") %>% as.data.frame()
umap_coords$Global_Clone_ID <- mdata_CMV$Global_Clone_ID
umap_coords$timepoint <- mdata_CMV$timepoint
umap_coords$pseudo_label <- paste0(umap_coords$Global_Clone_ID, "_", umap_coords$timepoint)

# Remove NA labels
umap_coords <- umap_coords[complete.cases(umap_coords$pseudo_label), ]

# Compute average UMAP position for each (clone, timepoint)
avg_umap <- umap_coords %>%
  group_by(pseudo_label) %>%
  summarise(across(starts_with("UMAP"), mean, na.rm = TRUE), .groups = "drop")

# Format as matrix with rownames = pseudo_label (e.g., clone123_Early)
avg_umap_mat <- as.matrix(avg_umap[, -1])
rownames(avg_umap_mat) <- avg_umap$pseudo_label

# Ensure it matches pseudobulk object
avg_umap_mat <- avg_umap_mat[colnames(pseudobulk), , drop = FALSE]

# Step 7: Add averaged UMAP to the new Seurat object
pseudobulk[["umap"]] <- CreateDimReducObject(
  embeddings = avg_umap_mat,
  key = "UMAP_",
  assay = DefaultAssay(pseudobulk)
)

# Step 8: Add cell counts per (clone, timepoint)
clone_counts <- meta_data %>%
  filter(!is.na(Global_Clone_ID), !is.na(timepoint)) %>%
  mutate(pseudo_label = paste0(Global_Clone_ID, "_", timepoint)) %>%
  group_by(pseudo_label) %>%
  summarise(cell_count = n(), .groups = "drop")

# Add to metadata of pseudobulk
pseudobulk$cell_count <- clone_counts$cell_count[match(colnames(pseudobulk), clone_counts$pseudo_label)]

# Pseudobulk visualization
pseudobulk <- NormalizeData(pseudobulk)  # Applies to default assay

# Extract axis limits from original UMAP (mdata_CMV)
umap_coords_original <- Embeddings(mdata_CMV, reduction = "umap")
x_limits <- range(umap_coords_original[, 1])
y_limits <- range(umap_coords_original[, 2])

# Plot pseudobulk UMAP with matching axis limits
DimPlot(pseudobulk, reduction = "umap") +
  xlim(x_limits) +
  ylim(y_limits)

DimPlot(mdata_CMV, reduction = "umap")

DimPlot(pseudobulk, reduction = "umap")

FeaturePlot(pseudobulk, features = "TCF7", reduction = "umap",
            cols = c("lightgrey", "darkred"))


# Per-donor plots

# Set output directory
output_dir <- "exports/integrated_v5/pseudobulk_analysis"

# Full UMAP from mdata_CMV for gray background
umap_bg_df <- Embeddings(mdata_CMV, reduction = "umap") %>%
  as.data.frame() %>%
  rename(UMAP_1 = 1, UMAP_2 = 2)

# Get UMAP coordinates
umap_df <- Embeddings(pseudobulk, reduction = "umap") %>% 
  as.data.frame() %>%
  tibble::rownames_to_column("pseudo_label")

# Combine with metadata
plot_df <- pseudobulk_filtered@meta.data %>%
  left_join(umap_df, by = "pseudo_label") %>%
  mutate(
    timepoint_shape = timepoint,
    log10_cell_count = log10(as.numeric(cell_count)),
    donor_id = factor(orig.ident, levels = paste0("D", 1:5))
  )

# Add simplified clone labels
plot_df <- plot_df %>%
  mutate(clone_label = Global_Clone_ID %>%
           str_remove("^D\\d+_") %>%          # remove donor info
           str_replace("-clone-", " clone "))       # format "IE1-clone-1" → "IE1 clone 1"

# Extract clone number and specificity for sorting
plot_df <- plot_df %>%
  mutate(
    specificity = str_extract(clone_label, "^[^ ]+"),
    clone_number = as.numeric(str_extract(clone_label, "\\d+$"))
  ) %>%
  arrange(specificity, clone_number) %>%
  mutate(clone_label = factor(clone_label, levels = unique(clone_label)))  # set factor with correct order

# Get size scale guide
size_guide <- data.frame(
  cell_count = c(10, 100, 1000),
  log10_cell_count = log10(c(10, 100, 1000)),
  UMAP_1 = min(plot_df$UMAP_1) - 2,
  UMAP_2 = seq(min(plot_df$UMAP_2), length.out = 3, by = 1)
)

# For connecting lines between timepoints per clone
line_df <- plot_df %>%
  select(Global_Clone_ID, timepoint, UMAP_1, UMAP_2, orig.ident, clone_label, donor_id) %>%
  filter(!is.na(Global_Clone_ID)) %>%
  group_by(orig.ident, Global_Clone_ID, clone_label, donor_id) %>%
  filter(n() == 2) %>%
  summarise(
    x = UMAP_1[timepoint == "Early"],
    y = UMAP_2[timepoint == "Early"],
    xend = UMAP_1[timepoint == "Late"],
    yend = UMAP_2[timepoint == "Late"],
    donor_id = first(donor_id),  # Keep donor_id for plotting
    .groups = "drop"
  )

# Use full UMAP limits from mdata_CMV for axis consistency
umap_limits <- umap_bg_df %>%
  summarise(
    xmin = min(UMAP_1) - 0.5,
    xmax = max(UMAP_1) + 0.5,
    ymin = min(UMAP_2) - 0.5,
    ymax = max(UMAP_2) + 0.5
  )

# Loop through each donor
for (donor in levels(plot_df$donor_id)) {
  donor_df <- plot_df %>% filter(donor_id == donor)
  donor_lines <- line_df %>% filter(orig.ident == donor)
  
  # Dummy points for consistent size legend
  dummy_sizes <- data.frame(
    UMAP_1 = rep(min(umap_limits$xmin) - 10, 3),  # Place off UMAP plot
    UMAP_2 = rep(min(umap_limits$ymin) - 10, 3),
    log10_cell_count = log10(c(10, 100, 1000)),
    timepoint_shape = "Late",  # Arbitrary, doesn’t affect size scale
    clone_label = NA,          # So they don't appear in color legend
    orig.ident = donor
  )
  
  # Append to donor_df (only used for legend)
  donor_df_plot <- bind_rows(donor_df, dummy_sizes)
  
  p <- ggplot(donor_df_plot, aes(x = UMAP_1, y = UMAP_2)) +
    geom_point(
      data = umap_bg_df,
      aes(x = UMAP_1, y = UMAP_2),
      color = "lightgray",
      size = 0.5,
      alpha = 0.3
    ) +
    geom_segment(data = donor_lines, aes(x = x, y = y, xend = xend, yend = yend, color = clone_label), size = 0.3) +
    geom_point(aes(color = clone_label, shape = timepoint_shape, size = sqrt(cell_count)), stroke = 0.5) +
    scale_shape_manual(
      name = "Timepoint",
      values = c("Early" = 1, "Late" = 16)
    ) +
    scale_size_continuous(
      name = "Clone Size",
      range = c(2, 10),
      breaks = sqrt(c(10, 100, 1000)),
      labels = c("10", "100", "1000"),
      limits = sqrt(c(10, 5000))
    ) +
    scale_color_discrete(name = "Global_Clone_ID", na.translate = FALSE) +
    guides(
      shape = guide_legend(override.aes = list(size = 4)),
      size = guide_legend(override.aes = list(shape = 16)),
      color = guide_legend(override.aes = list(size = 3))
    ) +
    guides(
      color = guide_legend(order = 1, override.aes = list(size = 3)),
      shape = guide_legend(order = 2, override.aes = list(size = 4)),
      size = guide_legend(order = 3, override.aes = list(shape = 16))
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "right") +
    ggtitle(paste("Donor", gsub("D", "", donor), "- Clones (Pseudobulked)")) +
    coord_cartesian(
      xlim = c(umap_limits$xmin, umap_limits$xmax),
      ylim = c(umap_limits$ymin, umap_limits$ymax),
      expand = FALSE
    )
  
  # Save
  ggsave(
    filename = file.path(output_dir, paste0("Donor_", gsub("D", "", donor), "_clones_pseudobulked.png")),
    plot = p,
    width = 9,
    height = 7
  )
}

# Unified plot with all clones colored by donor
p_all <- ggplot(plot_df, aes(x = UMAP_1, y = UMAP_2)) +
  geom_point(
    data = umap_bg_df,
    aes(x = UMAP_1, y = UMAP_2),
    color = "lightgray",
    size = 0.5,
    alpha = 0.3
  ) +
  geom_segment(
    data = line_df,
    aes(x = x, y = y, xend = xend, yend = yend, color = donor_id),
    size = 0.3
  ) +
  geom_point(
    aes(color = donor_id, shape = timepoint_shape, size = sqrt(cell_count)),
    stroke = 0.5
  ) +
  scale_color_brewer(name = "Donor", palette = "Set1") +  # Or use scale_color_manual() for custom colors
  scale_shape_manual(
    name = "Timepoint",
    values = c("Early" = 1, "Late" = 16)
  ) +
  scale_size_continuous(
    name = "Clone Size",
    range = c(2, 10),
    breaks = sqrt(c(10, 100, 1000)),
    labels = c("10", "100", "1000"),
    limits = sqrt(c(10, 5000))
  ) +
  guides(
    color = guide_legend(order = 1, override.aes = list(size = 3)),
    shape = guide_legend(order = 2, override.aes = list(size = 4)),
    size = guide_legend(order = 3, override.aes = list(shape = 16))
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "right") +
  ggtitle("All Donors – Clones (Pseudobulked)") +
  coord_cartesian(
    xlim = c(umap_limits$xmin, umap_limits$xmax),
    ylim = c(umap_limits$ymin, umap_limits$ymax),
    expand = FALSE
  )

# Save
ggsave(
  filename = file.path(output_dir, "All_donors_clones_pseudobulked.png"),
  plot = p_all,
  width = 9,
  height = 7
)


# Modify pseudobulk to remove small clones
# Step 1: Extract metadata
meta <- pseudobulk@meta.data %>%
  select(-matches("^pseudo_label$")) %>%  # Remove all columns named exactly 'pseudo_label'
  tibble::rownames_to_column("pseudo_label")

# Step 2: Parse Global_Clone_ID and timepoint
meta <- meta %>%
  mutate(
    clone_base = str_replace(pseudo_label, "_(Early|Late)$", ""),
    timepoint = ifelse(str_detect(pseudo_label, "_Late$"), "Late", "Early")
  )

# Step 3: Count cells per clone per timepoint
clone_counts <- meta %>%
  group_by(clone_base, timepoint) %>%
  summarise(cell_count = sum(as.numeric(cell_count), na.rm = TRUE), .groups = "drop")

# Step 4: Keep only clones with ≥10 cells at BOTH timepoints
valid_clones <- clone_counts %>%
  group_by(clone_base) %>%
  filter(all(c("Early", "Late") %in% timepoint)) %>%
  filter(all(cell_count > 15)) %>%
  distinct(clone_base)

# Step 5: Subset pseudobulk to only valid clones
valid_labels <- meta %>%
  filter(clone_base %in% valid_clones$clone_base) %>%
  pull(pseudo_label)

pseudobulk_filtered <- subset(pseudobulk, cells = valid_labels)



#
#
# Metacluster statistical analysis
#
#

# Jensen-Shannon divergence function (base 2)
js_divergence <- function(p, q) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  
  # Add small pseudocount to avoid log(0)
  epsilon <- 1e-10
  p <- p + epsilon
  q <- q + epsilon
  
  m <- 0.5 * (p + q)
  js <- 0.5 * sum(p * log2(p / m)) + 0.5 * sum(q * log2(q / m))
  return(js)
}

# File paths
input_dir <- "exports/integrated_v5/clone_cluster_pie_charts/cluster_stability_breakdowns"
output_file <- "exports/integrated_v5/timepoint_stability_plots/metacluster_statistics.csv"

# List files
csv_files <- list.files(input_dir, pattern = "Donor\\d+_(pp65|IE1)_clone\\d+_stability_breakdown.csv", full.names = TRUE)

results <- list()

for (file in csv_files) {
  
  file_name <- basename(file)
  matches <- str_match(file_name, "Donor(\\d+)_(pp65|IE1)_clone(\\d+)_stability_breakdown.csv")
  donor <- as.integer(matches[2])
  specificity <- matches[3]
  clone <- as.integer(matches[4])
  
  # Read the file
  df <- read.csv(file)
  
  # Get total cell counts
  n_cells_early <- df$n_cells_early[1]
  n_cells_late  <- df$n_cells_late[1]
  
  # Determine dominant cluster
  dominant_cluster <- "Heterogeneous"
  for (i in 1:nrow(df)) {
    if (df$freq_early[i] > 0.5 && df$freq_late[i] > 0.5) {
      dominant_cluster <- df$cluster[i]
      break
    }
  }
  
  # Compute JS divergence
  jsd <- js_divergence(df$freq_early, df$freq_late)
  
  # Store result
  results[[length(results) + 1]] <- data.frame(
    donor = donor,
    specificity = specificity,
    clone = clone,
    n_cells_early = n_cells_early,
    n_cells_late = n_cells_late,
    dominant_cluster = dominant_cluster,
    js_divergence = jsd
  )
}

# Combine and export
output_df <- bind_rows(results)
write.csv(output_df, output_file, row.names = FALSE)










#
# DONOR 0 ground-up analysis
#

# 1. Read raw 10x data
Input_10x_data <- Seurat::Read10X_h5("Donor_0/P565_4_1_raw_data/P565_4_1_raw_data.h5")

# 2. Rename cell barcodes in the raw matrices by prepending prefix
prefix <- "P565_4_1"

# Gene Expression matrix
gene_expr <- Input_10x_data$`Gene Expression`
colnames(gene_expr) <- paste0(prefix, "_", colnames(gene_expr))

# Antibody Capture matrix
ab_capture <- Input_10x_data$`Antibody Capture`
colnames(ab_capture) <- paste0(prefix, "_", colnames(ab_capture))

# 3. Create Seurat object with renamed cells
P565_4_1_raw <- CreateSeuratObject(counts = gene_expr)

# 4. Add assays, making sure you keep renamed columns
P565_4_1_raw[["FB"]] <- CreateAssayObject(counts = ab_capture[1:8, , drop = FALSE])
P565_4_1_raw[["Hashes"]] <- CreateAssayObject(counts = ab_capture[9:14, , drop = FALSE])

P565_4_1_raw <- 
  PercentageFeatureSet(P565_4_1_raw, "^MT-", col.name = "percent_mito") # percent mitochondrial reads
P565_4_1_raw <-
  PercentageFeatureSet(P565_4_1_raw, "^RP[SL]", col.name = "percent_ribo") # percent ribosomal protein reads
P565_4_1_raw <-
  PercentageFeatureSet(P565_4_1_raw, "^HB[^(P)]", col.name = "percent_hb") # percent hemoglobin reads

qcMetricsThresholds <-
  c(
    "min_nFeature_RNA" = 250,
    "max_nFeature_RNA" = 6000,
    "min_nCount_RNA" = 450,
    "max_nCount_RNA" = 35000,
    # "nFeature_ADT" = 50,
    # "nCount_ADT" = 500,
    "percent_mito" = 12.5,
    "min_percent_ribo" = 7.5,
    "max_percent_ribo" = 50,
    "percent_hb" = 10
  )

qcMetricsDfForPlot <-
  data.frame(
    metric = str_replace(names(qcMetricsThresholds), "^(min|max)_", ""),
    threshold = unname(qcMetricsThresholds)
  )

qcMetrics.tmp <-
  c("nFeature_RNA", "nCount_RNA",
    # "nFeature_ADT", "nCount_ADT",
    "percent_mito", "percent_ribo", "percent_hb")

rasterResolutionDpi <- 300
plot.tmp <-
  P565_4_1_raw@meta.data %>%
  dplyr::select(all_of(qcMetrics.tmp)) %>%
  pivot_longer(
    cols = one_of(qcMetrics.tmp),
    names_to = "metric", values_to = "value") %>%
  mutate(metric = metric %>% factor(levels = qcMetrics.tmp)) %>%
  ggplot(mapping = aes(x = "sample", y = value)) +
  geom_violin() +
  ggrastr::rasterise(
    geom_jitter(size = 0.5), dpi = rasterResolutionDpi) +
  facet_wrap(facets = vars(metric), ncol = 2, scales = "free_y") +
  geom_hline(
    data =
      qcMetricsDfForPlot %>%
      dplyr::filter(metric %in% qcMetrics.tmp) %>%
      mutate(metric = metric %>% factor(levels = qcMetrics.tmp)),
    mapping = aes(yintercept = threshold),
    color = "red", linetype = "dashed") +
  labs(x = NULL, y = NULL)

print(plot.tmp)

cellsKeepPool1 <-
  WhichCells(
    P565_4_1_raw,
    expression = nFeature_RNA >= qcMetricsThresholds["min_nFeature_RNA"]) %>% 
  intersect(
    WhichCells(
      P565_4_1_raw,
      expression = nFeature_RNA <= qcMetricsThresholds["max_nFeature_RNA"])) %>%
  intersect(
    WhichCells(
      P565_4_1_raw,
      expression = nCount_RNA >= qcMetricsThresholds["min_nCount_RNA"])) %>% 
  intersect(
    WhichCells(
      P565_4_1_raw,
      expression = nCount_RNA <= qcMetricsThresholds["max_nCount_RNA"])) %>% 
  intersect(
    WhichCells(
      P565_4_1_raw, expression = percent_ribo >= qcMetricsThresholds["min_percent_ribo"])) %>%
  intersect(
    WhichCells(
      P565_4_1_raw, expression = percent_ribo <= qcMetricsThresholds["max_percent_ribo"])) %>%
  intersect(
    WhichCells(
      P565_4_1_raw, expression = percent_mito <= qcMetricsThresholds["percent_mito"])) %>%
  intersect(
    WhichCells(
      P565_4_1_raw, expression = percent_hb <= qcMetricsThresholds["percent_hb"]))

cellsDropPool1 <-
  setdiff(colnames(P565_4_1_raw), cellsKeepPool1)

P565_4_1_includingMultiplets  <-
  subset(P565_4_1_raw, cells = cellsKeepPool1)

offset.tmp <- 1
P565_4_1_includingMultiplets <- P565_4_1_includingMultiplets %>%
  AddMetaData(
    metadata = log(t(as.matrix(P565_4_1_includingMultiplets[["Hashes"]]@counts + offset.tmp)) /
                     P565_4_1_includingMultiplets$nCount_Hashes * 1e6) %>%
      as.data.frame() %>%
      { colnames(.) <- make.names(paste0("logCpm", colnames(.))); . }
  )

thresholdsHtoPool1 <-
  list(logCounts = c("posAbove" = 5.5, "negBelow" = 4.5),
       logCpm = c("posAbove" = 11.5, "negBelow" = 11, "minGap" = 1.5))

P565_4_1_includingMultiplets@meta.data %>%
  pivot_longer(
    cols = matches("logCpm"),
    names_to = "hashtag", names_prefix = "log",
    values_to = "count") %>%
  mutate(hashtag = factor(hashtag, levels = str_sort(unique(hashtag), numeric=TRUE))) %>%
  ggplot(mapping = aes(x = count)) +
  geom_density(size=1) +
  facet_wrap(vars(hashtag), nrow = 2) +
  labs(x = "Hashtag counts (log)") +
  geom_vline(
    xintercept = thresholdsHtoPool1[["logCpm"]]["negBelow"],
    linetype = "dashed", color = "red") +
  geom_vline(
    xintercept = thresholdsHtoPool1[["logCpm"]]["posAbove"],
    linetype = "dashed", color = "blue")

logCpmHTO.tmp <-
  P565_4_1_includingMultiplets@meta.data %>%
  dplyr::select(matches("logCpm"))

thresholdLogCpmHTOStatus.tmp <-
  matrix(
    "ambiguous",
    nrow = nrow(logCpmHTO.tmp), ncol = ncol(logCpmHTO.tmp),
    dimnames = dimnames(logCpmHTO.tmp))

thresholdLogCpmHTOStatus.tmp[logCpmHTO.tmp > thresholdsHtoPool1[["logCpm"]]["posAbove"]] <- "positive"
thresholdLogCpmHTOStatus.tmp[logCpmHTO.tmp < thresholdsHtoPool1[["logCpm"]]["negBelow"]] <- "negative"

thresholdLogCpmHTOStatus.tmp <-
  matrix(
    "ambiguous",
    nrow = nrow(logCpmHTO.tmp), ncol = ncol(logCpmHTO.tmp),
    dimnames = dimnames(logCpmHTO.tmp))
thresholdLogCpmHTOStatus.tmp[logCpmHTO.tmp > thresholdsHtoPool1[["logCpm"]]["posAbove"]] <- "positive"
thresholdLogCpmHTOStatus.tmp[
  t(apply(logCpmHTO.tmp, MARGIN = 1, function(x) {max(x) - x})) >
    thresholdsHtoPool1[["logCpm"]]["minGap"]] <-
  "negative"

sampleBarcode.tmp <- rep("unassigned", ncol(P565_4_1_includingMultiplets))

cells.tmp <-
  which(apply(thresholdLogCpmHTOStatus.tmp, 
              MARGIN = 1, FUN = function(x) {sum(x %in% c("ambiguous", "positive")) > 1}))
sampleBarcode.tmp[cells.tmp] <- "multiplePositive"

cells.tmp <-
  which(apply(thresholdLogCpmHTOStatus.tmp,
              MARGIN = 1, 
              FUN = function(x) {(sum(x %in% "positive") == 1) & (sum(x %in% "ambiguous") == 0)}))
sampleBarcode.tmp[cells.tmp] <-
  colnames(thresholdLogCpmHTOStatus.tmp)[
    apply(thresholdLogCpmHTOStatus.tmp[cells.tmp,],
          MARGIN = 1, function(x) which(x %in% "positive"))] %>%
  str_replace_all("logCpm", "HTO-")

sampleBarcode.tmp <-
  sampleBarcode.tmp %>%
  factor(levels = str_sort(unique(.), numeric=TRUE))

P565_4_1_includingMultiplets <- P565_4_1_includingMultiplets %>%
  AddMetaData(
    metadata = sampleBarcode.tmp,
    col.name = "sampleBarcodeThresholdLogCpm")

Idents(P565_4_1_includingMultiplets) <- "sampleBarcodeThresholdLogCpm"
P565_4_1 <- subset(x = P565_4_1_includingMultiplets, idents = "multiplePositive", invert = TRUE)

add_clonotype <- function(tcr_prefix, seurat_obj, prefix = NULL, type = "t") {
  # Read filtered contig annotations
  tcr <- read.csv(paste0(tcr_prefix, "filtered_contig_annotations.csv"))
  
  # Remove duplicated barcodes, keep first occurrence
  tcr <- tcr[!duplicated(tcr$barcode), ]
  
  # Add prefix to barcodes if provided (to match Seurat cell names)
  if (!is.null(prefix)) {
    tcr$barcode <- paste0(prefix, "_", tcr$barcode)
  }
  
  # Keep only barcode and clonotype columns
  tcr <- tcr[, c("barcode", "raw_clonotype_id")]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  
  # Read clonotypes data
  clono <- read.csv(paste0(tcr_prefix, "clonotypes.csv"))
  
  # Merge clonotype info (cdr3s_aa) onto tcr by clonotype_id
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")], by = "clonotype_id", all.x = TRUE)
  
  # Rearrange columns: barcode first, set as rownames
  tcr <- tcr[, c("barcode", "clonotype_id", "cdr3s_aa")]
  rownames(tcr) <- tcr$barcode
  tcr$barcode <- NULL
  
  # Rename columns with type prefix (e.g., t_clonotype_id)
  colnames(tcr) <- paste0(type, "_", colnames(tcr))
  
  # Add metadata to Seurat object
  seurat_obj <- AddMetaData(object = seurat_obj, metadata = tcr)
  
  return(seurat_obj)
}

P565_4_1 <- add_clonotype("Donor_0/P565_4_1_raw_data/", P565_4_1, prefix = "P565_4_1", type = "t")

#TCRgraph
if (!require("devtools")) install.packages("devtools")
devtools::install_github("BenaroyaResearch/tcrGraph", build_vignettes = TRUE)
library("tcrGraph")

P565_4_1@meta.data$barcode <- rownames(P565_4_1@meta.data)

P565_4_1_Contigs <- read.table("Donor_0/P565_4_1_raw_data/airr_rearrangement.tsv", sep = '\t', header = TRUE)
prefix <- "P565_4_1_"
P565_4_1_Contigs$cell_id <- paste0(prefix, P565_4_1_Contigs$cell_id)

P565_4_1_TCRgraph <- P565_4_1@meta.data %>% dplyr::select(barcode) %>% rename(cell_id = barcode)

P565_4_1_TCRgraph <- inner_join(P565_4_1_Contigs, P565_4_1_TCRgraph, by = "cell_id")

P565_4_1_TCRgraph <- P565_4_1_TCRgraph %>% select(cell_id, v_call, j_call, sequence) %>% rename(libid = cell_id, v_gene = v_call, j_gene = j_call)

P565_4_1_TCRGraphOutput <- makeTcrGraph(P565_4_1_TCRgraph, link = "sequence")

plot(P565_4_1_TCRGraphOutput)

clones <- getClonesFromTcrGraph(P565_4_1_TCRGraphOutput, maxB = 1, format = "full")

clones <- clones %>% select(libid, cloneId) %>% rename(barcode = libid)

clones <- clones[!duplicated(clones$barcode),]

P565_4_1@meta.data <- left_join(P565_4_1@meta.data, clones, by = "barcode")

P565_4_1@meta.data$t_clonotype_id[is.na(P565_4_1@meta.data$t_clonotype_id)] <- "notCalled"

P565_4_1@meta.data$cloneId[is.na(P565_4_1@meta.data$cloneId)] <- "notCalled"

P565_4_1$cloneId <- paste0("run1_", P565_4_1$cloneId)


#Merging

rownames(P565_4_1@meta.data) <- Cells(P565_4_1)
rownames(P565_4_2@meta.data) <- Cells(P565_4_2)

P565_4_merged <- merge(
  x = P565_4_1,
  y = P565_4_2,
  merge.data = TRUE
)


# Extract the separate count layers
counts_1 <- GetAssayData(P565_4_merged[["RNA"]], layer = "counts.1")
counts_2 <- GetAssayData(P565_4_merged[["RNA"]], layer = "counts.2")

# Combine counts by column-binding (combining cells side-by-side)
combined_counts <- Matrix::cbind2(counts_1, counts_2)

# Replace the RNA assay counts slot with the combined counts
P565_4_merged[["RNA"]] <- SetAssayData(
  P565_4_merged[["RNA"]],
  layer = "counts",
  new.data = combined_counts
)

# Remove the extra count layers if needed
P565_4_merged[["RNA"]]@layers$counts.1 <- NULL
P565_4_merged[["RNA"]]@layers$counts.2 <- NULL

# TCR gene quieting

P565_4_merged <- quietTCRgenes_new(P565_4_merged)

# Load the TCR gene list
get_unwanted_tcr_genes <- function() {
  unwanted_genes <- c()
  for (file in list.files("long_excluded_gene_lists/TCellReceptors/", full.names = TRUE)) {
    csv <- read.csv(file)
    unwanted_genes <- c(unwanted_genes, csv[["Approved.symbol"]])
  }
  unique(unwanted_genes)
}

# Get unwanted gene list
unwanted_genes <- get_unwanted_tcr_genes()

# Get RNA data and remove unwanted genes
rna_counts <- GetAssayData(P565_4_merged, assay = "RNA", slot = "counts")
filtered_counts <- rna_counts[!rownames(rna_counts) %in% unwanted_genes, ]

# Recreate RNA assay with filtered genes
P565_4_merged[["RNA"]] <- CreateAssayObject(counts = filtered_counts)

# Create or update the 'specificity' column with multiple conditions
P565_4_merged$specificity <- ifelse(grepl("pp65", P565_4_merged@meta.data[["sampleBarcodeThresholdLogCpm"]]), "pp65",
                               ifelse(grepl("IE1", P565_4_merged@meta.data[["sampleBarcodeThresholdLogCpm"]]), "IE1",
                                      ifelse(grepl("Bulk", P565_4_merged@meta.data[["sampleBarcodeThresholdLogCpm"]]), "Bulk",
                                            "NA")))

# Create or update the 'timepoint' column with multiple conditions
P565_4_merged$timepoint <- ifelse(grepl("visit.1", P565_4_merged@meta.data[["sampleBarcodeThresholdLogCpm"]]), "Visit 1",
                                    ifelse(grepl("visit.2", P565_4_merged@meta.data[["sampleBarcodeThresholdLogCpm"]]), "Visit 2",
                                           ifelse(grepl("visit.3", P565_4_merged@meta.data[["sampleBarcodeThresholdLogCpm"]]), "Visit 3",
                                                  ifelse(grepl("visit.4", P565_4_merged@meta.data[["sampleBarcodeThresholdLogCpm"]]), "Visit 4",
                                                        "NA"))))

#DefaultAssay(P565_4_merged) <- "RNA"
P565_4_merged <- NormalizeData(P565_4_merged)
P565_4_merged <- FindVariableFeatures(P565_4_merged, selection.method = "vst")
P565_4_merged <- ScaleData(P565_4_merged, verbose = FALSE)
P565_4_merged <- RunPCA(P565_4_merged, npcs = 30, verbose = FALSE)
P565_4_merged <- RunUMAP(P565_4_merged, reduction = "pca", dims = 1:30, reduction.name = "umap")
P565_4_merged <- FindNeighbors(P565_4_merged, reduction = "pca", dims = 1:30)
P565_4_merged <- FindClusters(P565_4_merged, resolution = 0.5)

#DefaultAssay(P565_4_merged) <- "RNA"
DimPlot(P565_4_merged, reduction = "umap", group.by = "specificity")
DimPlot(P565_4_merged, reduction = "umap", group.by = "timepoint")
DimPlot(P565_4_merged, reduction = "umap", group.by = "seurat_clusters")


# #
# # Integration for P565_4_merged
# #
# # Split and pre-process
# seurat_list <- SplitObject(P565_4_merged, split.by = "timepoint")
# shared_genes <- Reduce(intersect, lapply(seurat_list, rownames))
# 
# # Process each subset
# seurat_list <- lapply(seurat_list, function(x) {
#   x <- subset(x, features = shared_genes)
#   x <- NormalizeData(x)
#   x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
#   return(x)
# })
# 
# # Select integration features
# features <- SelectIntegrationFeatures(object.list = seurat_list)
# # Scale and PCA
# seurat_list <- lapply(seurat_list, function(x) {
#   x <- ScaleData(x, features = features)
#   x <- RunPCA(x, features = features)
#   return(x)
# })
# # Save other assays before integration
# original_assays <- P565_4_merged@assays
# non_rna_assays <- original_assays[!names(original_assays) %in% c("RNA")]
# # Integration
# anchors <- FindIntegrationAnchors(object.list = seurat_list, anchor.features = features)
# P565_4_merged_integrated <- IntegrateData(
#   anchorset = anchors,
#   new.assay.name = "integrated",
#   features = integration_features, # the variable features used for finding anchors
#   layer = "data"  # Specify which layer to pull from
# )
# # Reattach original non-RNA assays
# for (assay_name in names(non_rna_assays)) {
#   P565_4_merged_integrated[[assay_name]] <- non_rna_assays[[assay_name]]
# }
# # Replace P565_4_merged with the enriched object
# P565_4_merged <- P565_4_merged_integrated
# 


saveRDS(P565_4_merged, file = "Donor_0/D0_data.rds")
saveRDS(D0_data, file = "Donor_0/D0_data.rds")
saveRDS(D0_data_CMV, file = "Donor_0/D0_data_CMV.rds")

D0_data <- readRDS (file = "Donor_0/D0_data.rds")
D0_data_CMV <- readRDS (file = "Donor_0/D0_data_CMV.rds")


D0_data_CMV <- subset(D0_data, subset = specificity %in% c('pp65', 'IE1'))
D0_data_Bulk <- subset(D0_data, subset = specificity == 'Bulk')


# Visualization

#DefaultAssay(D0_data_CMV) <- "RNA"
D0_data_Bulk <- NormalizeData(D0_data_Bulk)
D0_data_Bulk <- FindVariableFeatures(D0_data_Bulk, selection.method = "vst")
D0_data_Bulk <- ScaleData(D0_data_Bulk, verbose = FALSE)
D0_data_Bulk <- RunPCA(D0_data_Bulk, npcs = 30, verbose = FALSE)
D0_data_Bulk <- RunUMAP(D0_data_Bulk, reduction = "pca", dims = 1:30, reduction.name = "umap")
D0_data_Bulk <- FindNeighbors(D0_data_Bulk, reduction = "pca", dims = 1:30)
D0_data_Bulk <- FindClusters(D0_data_Bulk, resolution = 0.5)

DefaultAssay(D0_data_CMV) <- "FB"
D0_data_CMV <- NormalizeData(D0_data_CMV)

DimPlot(D0_data, reduction = "umap", group.by = "specificity")
DimPlot(D0_data_CMV, reduction = "umap", group.by = "timepoint")
DimPlot(D0_data_Bulk, reduction = "umap", group.by = "seurat_clusters")
FeaturePlot(D0_data_CMV, features = "anti-human CD57", reduction = "umap",
            cols = c("lightgrey", "darkred"))

FeaturePlot(D0_data, features = "KLRC2", reduction = "umap",
            cols = c("lightgrey", "darkred"))


# Define named vector of colors
cluster_colors <- c(
  "0" = "#C49A00",  
  "1" = "#53B400",  
  "2" = "#8B0000",  
  "3" = "#1976D2",  
  "4" = "#F8766D",  
  "5" = "#A58AFF",  
  "6" = "#FFDB58",  
  "7" = "#00C094"   
)


# Set the timepoint factor levels in desired order
D0_data_CMV$timepoint <- factor(D0_data_CMV$timepoint, levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))

# Generate individual plots for each timepoint in order
plots <- lapply(levels(D0_data_CMV$timepoint), function(tp) {
  DimPlot(D0_data_CMV, group.by = "seurat_clusters", cells = WhichCells(D0_data_CMV, expression = timepoint == tp), cols = cluster_colors) +
    ggtitle(tp) +
    theme(plot.title = element_text(hjust = 0.5))
})

# Combine plots side-by-side
combined_plot <- wrap_plots(plots, nrow = 1)

# Display the plot
print(combined_plot)

# Save to file (optional)
ggsave("Donor_0/exports/cluster_umaps_by_timepoint.png", combined_plot, width = 16, height = 4)


#
# Only pp65- or IE-specific cells
#
# Generate individual plots for each timepoint in order
plots <- lapply(levels(D0_data_CMV$timepoint), function(tp) {
  DimPlot(
    D0_data_CMV,
    group.by = "specificity",
    cells = WhichCells(
      D0_data_CMV,
      expression = timepoint == tp & specificity == "IE1"
    )
  ) +
    scale_color_manual(values = c(
      "pp65" = "#0047FF",
      "IE1"  = "#FF1A1A"
    )) +
    ggtitle(paste0(tp, " - IE1")) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none"
    )
})


# Combine plots side-by-side
combined_plot <- wrap_plots(plots, nrow = 1)

# Save to file (optional)
ggsave("Donor_0/exports/specificity_umaps_by_timepoint_IE1.png", combined_plot, width = 16, height = 4)




# Tetramer-negative data
# Set the timepoint factor levels in desired order
D0_data_Bulk$timepoint <- factor(D0_data_Bulk$timepoint, levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))

# Generate individual plots for each timepoint in order
plots <- lapply(levels(D0_data_Bulk$timepoint), function(tp) {
  DimPlot(D0_data_Bulk, group.by = "seurat_clusters", cells = WhichCells(D0_data_Bulk, expression = timepoint == tp) ) + #, cols = cluster_colors) +
    ggtitle(tp) +
    theme(plot.title = element_text(hjust = 0.5))
})

# Combine plots side-by-side
combined_plot <- wrap_plots(plots, nrow = 1)

# Display the plot
print(combined_plot)

# Save to file (optional)
ggsave("Donor_0/exports/Bulk_data_only_cluster_umaps_by_timepoint.png", combined_plot, width = 16, height = 4)



# CLONOTYPE ANALYSIS
# Extract numbers for each TCRgraph clonotype
#
# Edited to extract Global_Clone_IDs 
#
# Step 1: Extract metadata
meta <- D0_data_CMV@meta.data

# Step 2: Define specificity values to loop over
specificities <- c("pp65", "IE1")

# Step 3: Create output directory
output_dir <- "Donor_0/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Step 4: Loop over each specificity
for (spec in specificities) {
  
  # Filter metadata
  meta_spec <- meta %>%
    filter(specificity == spec, !is.na(cloneId), !grepl("notCalled", cloneId))
  
  # Total cells per timepoint for this specificity
  timepoint_totals <- meta_spec %>%
    group_by(timepoint) %>%
    summarize(n_total = n(), .groups = "drop")
  
  # Count cells per cloneId per timepoint
  clone_counts <- meta_spec %>%
    group_by(cloneId, timepoint) %>%
    summarize(count = n(), .groups = "drop") %>%
    pivot_wider(names_from = timepoint, values_from = count, values_fill = 0)
  
  # Ensure Visit columns are in the correct order
  visit_cols <- c("Visit 1", "Visit 2", "Visit 3", "Visit 4")
  for (col in visit_cols) {
    if (!col %in% colnames(clone_counts)) {
      clone_counts[[col]] <- 0
    }
  }
  clone_counts <- clone_counts[, c("cloneId", visit_cols)]
  
  # Add total for sorting only
  clone_counts$Total <- rowSums(clone_counts[, visit_cols])
  clone_counts <- clone_counts %>% arrange(desc(Total))
  
  # Calculate percentage columns
  for (tp in visit_cols) {
    total_cells <- timepoint_totals$n_total[timepoint_totals$timepoint == tp]
    pct_col <- paste0(tp, "_pct")
    if (length(total_cells) == 0 || total_cells == 0) {
      clone_counts[[pct_col]] <- NA
    } else {
      clone_counts[[pct_col]] <- 100 * clone_counts[[tp]] / total_cells
    }
  }
  
  # === Add Global_Clone_ID mapping, including NAs ===
  clone_to_global <- meta_spec %>%
    group_by(cloneId) %>%
    summarize(Global_Clone_ID = first(Global_Clone_ID), .groups = "drop")  # includes NAs
  
  # Merge into the results
  clone_counts <- left_join(clone_counts, clone_to_global, by = "cloneId") %>%
    relocate(Global_Clone_ID, .after = cloneId)
  
  # Drop Total before writing
  clone_counts <- clone_counts %>% select(-Total)
  
  # Write to CSV
  output_file <- file.path(output_dir, paste0(spec, "_global_clone_counts.csv"))
  write_csv(clone_counts, output_file, na = "NA")  # explicitly write "NA" in output
}


# CLONOTYPE ANALYSIS
# Extract TCR sequences for each clonotype

# Extract metadata
metadata <- D0_data_CMV@meta.data

# Initialize an empty data frame to store results
results <- data.frame(Specificity = character(), Clone_ID = character(), stringsAsFactors = FALSE)

# Function to update or add a row
update_results <- function(results, clone_id, cdr3_aa, specificity) {
  
  # Check if Clone_ID already exists
  row_index <- which(results$Clone_ID == clone_id)
  
  if (length(row_index) == 0) {
    # Clone_ID not found, add a new row with specificity
    new_row <- data.frame(Specificity = specificity, Clone_ID = clone_id, stringsAsFactors = FALSE)
    new_row["CDR3_aa-sequence_1"] <- cdr3_aa
    new_row["Count_1"] <- 1
    results <- bind_rows(results, new_row)
    
  } else {
    # Clone_ID exists, check for CDR3 sequence
    row_index <- row_index[1]  # Ensure single index
    
    # Get the columns relevant to this row
    existing_sequences <- grep("CDR3_aa-sequence_", names(results), value = TRUE)
    
    # Check if the cdr3_aa is already in the row
    match_found <- FALSE
    for (seq_col in existing_sequences) {
      count_col <- gsub("CDR3_aa-sequence_", "Count_", seq_col)
      
      if (!is.na(results[row_index, seq_col]) && results[row_index, seq_col] == cdr3_aa) {
        # Increment count
        results[row_index, count_col] <- results[row_index, count_col] + 1
        match_found <- TRUE
        break
      }
    }
    
    if (!match_found) {
      # Find the next available column for this specific row
      used_cols <- sum(!is.na(results[row_index, existing_sequences]))
      next_index <- used_cols + 1  # The next available index for this row
      
      new_seq_col <- paste0("CDR3_aa-sequence_", next_index)
      new_count_col <- paste0("Count_", next_index)
      
      # Ensure columns exist before assigning values
      if (!(new_seq_col %in% names(results))) {
        results[[new_seq_col]] <- NA
        results[[new_count_col]] <- NA
      }
      
      # Assign new sequence and count in the correct spot
      results[row_index, new_seq_col] <- cdr3_aa
      results[row_index, new_count_col] <- 1
    }
  }
  
  return(results)
}

# Iterate through Seurat metadata
for (i in 1:nrow(metadata)) {
  clone_id <- metadata$cloneId[i]
  cdr3_aa <- metadata$t_cdr3s_aa[i]
  specificity <- metadata$specificity[i]  # Extract specificity
  
  if (!is.na(clone_id) && !is.na(cdr3_aa)) {
    results <- update_results(results, clone_id, cdr3_aa, specificity)
  }
}

# Step 1: Add "Total_count" column (sum of all Count_x)
results <- results %>%
  mutate(Total_count = rowSums(select(., starts_with("Count_")), na.rm = TRUE)) %>%
  relocate(Total_count, .after = Clone_ID)

# Step 2: Filter out rows where Total_count < 5
filtered_results <- results %>%
  filter(Total_count >= 5)

# Step 3: Remove rows where Clone_ID contains "notCalled"
filtered_results <- filtered_results %>%
  filter(!grepl("notCalled", Clone_ID))

# Step 3.5: Only keep pp65 clones
filtered_results <- filtered_results %>%
  filter(Specificity == "IE1")

# Step 4: Remove completely empty columns
filtered_results <- filtered_results %>%
  select(where(~ !all(is.na(.))))

# Step 5: Export to CSV
write.csv(filtered_results, "Donor_0/exports/filtered_TCR_sequences/IE1_clonotypes_TCR_sequences_filtered.csv", row.names = FALSE)


#
# Matching aa sequence values??
#
find_matching_values <- function(file_path) {
  
  # Read the CSV file
  df <- read.csv(file_path, stringsAsFactors = FALSE)
  
  # Normalize column names in case R has converted "_" to "." (e.g., "CDR3_aa-sequence_1" → "CDR3.aa.sequence.1")
  colnames(df) <- gsub("\\.", "-", colnames(df))  # Replace all dots with underscores
  
  # Identify columns that match the pattern "CDR3_aa-sequence_x" (x = 1 to 100)
  scan_cols <- grep("^CDR3_aa-sequence_[1-9][0-9]?$|^CDR3_aa-sequence_100$", colnames(df), value = TRUE)
  
  if (length(scan_cols) == 0) {
    stop("Error: No 'CDR3_aa-sequence_x' columns found in the dataset.")
  }
  
  # Initialize results dataframe
  results <- data.frame(Match_1_gclone = character(),
                        Match_1_cloneID = character(),
                        Match_2_gclone = character(),
                        Match_2_cloneID = character(),
                        Matching_aa_sequence = character(),
                        stringsAsFactors = FALSE)
  
  # Create a list of all (row, column) pairs in the scanned range
  scanned_cells <- expand.grid(Row = 1:nrow(df), Column = scan_cols, stringsAsFactors = FALSE)
  
  # Iterate over all unique pairs of scanned cells
  for (idx1 in 1:(nrow(scanned_cells) - 1)) {
    for (idx2 in (idx1 + 1):nrow(scanned_cells)) {
      row1 <- scanned_cells$Row[idx1]
      col1 <- scanned_cells$Column[idx1]
      row2 <- scanned_cells$Row[idx2]
      col2 <- scanned_cells$Column[idx2]
      
      val1 <- df[row1, col1]
      val2 <- df[row2, col2]
      
      # Check if both values are non-NA before proceeding
      if (!is.na(val1) && !is.na(val2)) {
        # Check if both values are non-empty and match
        if (val1 != "" && val2 != "" && val1 == val2) {
          results <- rbind(results, data.frame(
            Match_1_gclone = df$GClone_name[row1],
            Match_1_cloneID = df$Clone_ID[row1],
            Match_2_gclone = df$GClone_name[row2],
            Match_2_cloneID = df$Clone_ID[row2],
            Matching_aa_sequence = val1,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }
  
  # Create output filename
  file_name <- tools::file_path_sans_ext(basename(file_path))
  output_file <- paste0(dirname(file_path), "/", file_name, "_AA-SEQUENCE-MATCHES.csv")
  
  # Export results if there are matches
  if (nrow(results) > 0) {
    write.csv(results, output_file, row.names = FALSE)
    message("Results saved to: ", output_file)
  } else {
    message("No matches found. No output file created.")
  }
  
  return(results)
}


find_matching_values("Donor_0/exports/filtered_TCR_sequences/pp65_GClone_TCR_sequences.csv")



#
# R script to ADD GLOBAL CLONES and generate clonotype plots
#

# Load the CSV file for adding global clone names to the object
csv_path <- "Donor_0/D0_Global_Clone_IDs.csv"
global_clones <- read.csv(csv_path, stringsAsFactors = FALSE)

# Extract metadata from Seurat object
metadata <- D0_data@meta.data

# Initialize new metadata column with NA
metadata$Global_Clone_ID <- NA

# Create a named lookup vector for fast matching
clone_lookup <- setNames(global_clones$GClone_name, global_clones$Clone_ID)

# Update 'Global_Clone_ID' based on matches in 'cloneId'
metadata$Global_Clone_ID <- clone_lookup[metadata$cloneId]

# Assign updated metadata back to Seurat object
D0_data@meta.data <- metadata


# Open PDF file
pdf(file = "Donor_0/exports/Clonotype_plots.pdf", width = 8, height = 6)

# Loop through IE1 clones
for(i in 1:10) {
  clonotype_id <- paste0("D0_IE1_clone_", i)
  
  # Check if there are any cells with this clonotype
  cells_to_highlight <- WhichCells(D0_data_CMV, expression = Global_Clone_ID == clonotype_id)
  
  # Only plot if there are cells to highlight
  if (length(cells_to_highlight) > 0) {
    # Generate the DimPlot for each clonotype
    plot <- DimPlot(D0_data_CMV, reduction = 'umap', 
                    cells.highlight = cells_to_highlight, 
                    sizes.highlight = 1.25) + ggtitle(clonotype_id)
    
    # Print each plot to the PDF file
    print(plot)
  }
}

#Loop through pp65 clones
for(i in 1:50) {
  clonotype_id <- paste0("D0_pp65_clone_", i)
  
  # Check if there are any cells with this clonotype
  cells_to_highlight <- WhichCells(D0_data_CMV, expression = Global_Clone_ID == clonotype_id)
  
  # Only plot if there are cells to highlight
  if (length(cells_to_highlight) > 0) {
    # Generate the DimPlot for each clonotype
    plot <- DimPlot(D0_data_CMV, reduction = 'umap', 
                    cells.highlight = cells_to_highlight, 
                    sizes.highlight = 1.25) + ggtitle(clonotype_id)
    
    # Print each plot to the PDF file
    print(plot)
  }
}

# Close the PDF file
dev.off()



#
# Find markers for all clusters
#
#install.packages('devtools')
devtools::install_github('immunogenomics/presto')

Idents(D0_data_CMV) <- "seurat_clusters"
Idents(D0_data_CMV) <- "specificity"

cluster_markers <- FindAllMarkers(D0_data_CMV, only.pos = TRUE, logfc.threshold = 0.25)
filtered_markers <- cluster_markers %>% filter(pct.1 > 0.1)

# Plot a heatmap of the top 10 markers per cluster
top10 <- cluster_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- filtered_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top15 <- filtered_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top20 <- filtered_markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

DoHeatmap(D0_data_CMV, features = top15$gene) + scale_fill_gradientn(colors = c("blue", "white", "red"))



#
# List of protein-coding genes to filter DEGs
#

library(biomaRt)

# Your gene list
genes <- rownames(filtered_markers)  # or any character vector of gene symbols

# Connect to Ensembl
mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

# Retrieve gene biotypes
gene_info <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters = "hgnc_symbol",
  values = genes,
  mart = mart
)

# Keep protein-coding genes only
protein_coding_genes <- gene_info$hgnc_symbol[
  gene_info$gene_biotype == "protein_coding"
]

# Remove empty symbols
protein_coding_genes <- protein_coding_genes[protein_coding_genes != ""]

# Subset Seurat object
filtered_markers <- filtered_markers[
  rownames(filtered_markers) %in% protein_coding_genes,
]



#
#
# EXPORTS
#
#

ggsave("Donor_0/exports/D0_CMV_UMAP_ISG+_score.png", plot = 
         FeaturePlot(D0_data_CMV, features = "ISG+_score", reduction = "umap",
                     cols = c("lightgrey", "darkred")), 
       width = 6, height = 5, dpi = 300)

ggsave("Donor_0/exports/D0_CMV_UMAP_by_metacluster.png", plot = 
         DimPlot(D0_data_CMV, group.by = 'metacluster', reduction = "umap"), 
       width = 6, height = 5, dpi = 300)

ggsave(
  "Donor_0/exports/D0_CMV_DEGs_per_cluster.png",
  plot = DoHeatmap(
    D0_data_CMV,
    features = top15$gene,
    group.colors = cluster_colors
  ) +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +  # expression values
    theme(
      axis.text.y = element_text(size = 18)  # gene names
      # legend.text and legend.title can be adjusted here if needed
    ),
  width = 10, height = 15, dpi = 300
)



#
# Stacked bar graphs for each donor split by timepoint
#
# Define export directory
export_dir <- "Donor_0/exports/clone_cluster_stacked_bars"

# Extract metadata and parse fields
meta <- D0_data_CMV@meta.data %>%
  mutate(
    donor = str_extract(Global_Clone_ID, "D\\d") %>% str_remove("D") %>% as.integer(),
    specificity = str_extract(Global_Clone_ID, "(pp65|IE1)"),
    clone_number = str_extract(Global_Clone_ID, "clone_\\d+") %>% str_remove("clone_") %>% as.integer(),
    cluster = seurat_clusters
  )

# Define consistent cluster colors
# all_clusters <- sort(unique(meta$cluster))
# cluster_colors <- setNames(
#  RColorBrewer::brewer.pal(n = max(3, min(12, length(all_clusters))), "Set3")[1:length(all_clusters)],
#  all_clusters
# )


## Extract default Seurat cluster colors from DimPlot
# default_cluster_colors <- scales::hue_pal()(
#   length(unique(D0_data_CMV$seurat_clusters))
# )
# names(default_cluster_colors) <- sort(unique(D0_data_CMV$seurat_clusters))

# cluster_colors should be already defined (at line ~3270 in this code)

# Replace cluster_colors in the plotting function
generate_timepoint_plot <- function(df, d, time_label) {
  # Define all possible clone labels
  ie1_labels <- paste0("IE1 Clone ", 1:7)
  pp65_labels <- paste0("pp65 Clone ", 1:40)
  
  dummy_label <- " "  
  all_labels <- c(ie1_labels, dummy_label, pp65_labels)
  
  # Summarize clone sizes
  clone_summary <- df %>%
    count(specificity, clone_number, name = "n_cells_clone") %>%
    group_by(specificity) %>%
    mutate(clone_size_percent = n_cells_clone / sum(n_cells_clone) * 100) %>%
    ungroup() %>%
    mutate(clone_label = ifelse(specificity == "IE1",
                                paste0("IE1 Clone ", clone_number),
                                paste0("pp65 Clone ", clone_number)))
  
  comp_data <- df %>%
    count(specificity, clone_number, cluster, name = "n_cells") %>%
    group_by(specificity, clone_number) %>%
    mutate(freq = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    mutate(clone_label = ifelse(specificity == "IE1",
                                paste0("IE1 Clone ", clone_number),
                                paste0("pp65 Clone ", clone_number))) %>%
    complete(clone_label = c(ie1_labels, pp65_labels), cluster = sort(unique(D0_data_CMV$seurat_clusters)), fill = list(n_cells = 0, freq = 0)) %>%
    left_join(clone_summary %>% select(clone_label, clone_size_percent), by = "clone_label") %>%
    mutate(clone_size_percent = replace_na(clone_size_percent, 0))
  
  dummy_row <- expand.grid(
    clone_label = dummy_label,
    cluster = sort(unique(D0_data_CMV$seurat_clusters)),
    n_cells = 0,
    freq = 0,
    clone_size_percent = 0
  )
  comp_data <- bind_rows(comp_data, dummy_row)
  
  comp_data$label_color <- ifelse(grepl("IE1", comp_data$clone_label), "red",
                                  ifelse(grepl("pp65", comp_data$clone_label), "blue", "black"))
  
  comp_data$clone_label <- factor(comp_data$clone_label, levels = rev(all_labels))
  
  dummy_position <- which(levels(comp_data$clone_label) == dummy_label)
  ie1_break_position <- dummy_position
  
  # ---- Plot ----
  ggplot(comp_data, aes(x = freq, y = clone_label, fill = factor(cluster))) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(
      data = comp_data %>% filter(clone_label != dummy_label) %>% group_by(clone_label) %>% slice(1),
      aes(x = 1.02, y = clone_label, label = sprintf("%.1f%%", clone_size_percent), color = label_color),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3
    ) +
    scale_fill_manual(values = cluster_colors, drop = FALSE) +  # << change here
    scale_color_identity() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::percent_format(accuracy = 1)) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 8, color = rev(c(rep("red", 7), "black", rep("blue", 40)))),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.margin = margin(10, 30, 10, 10)
    ) +
    guides(fill = guide_legend(title = "Metacluster")) +
    geom_hline(yintercept = ie1_break_position, linetype = "dashed") +
    ggtitle(paste("Donor", d, "-", time_label, "Timepoint"))
}

# Loop through donors and generate plots
for (d in sort(unique(meta$donor))) {
  donor_meta <- meta %>% filter(donor == d)
  
  V1_data <- donor_meta %>% filter(timepoint == "Visit 1")
  V2_data <- donor_meta %>% filter(timepoint == "Visit 2")
  V3_data <- donor_meta %>% filter(timepoint == "Visit 3")
  V4_data <- donor_meta %>% filter(timepoint == "Visit 4")
  
  p_1 <- generate_timepoint_plot(V1_data, d, "Visit 1")
  p_2 <- generate_timepoint_plot(V2_data, d, "Visit 2")
  p_3 <- generate_timepoint_plot(V3_data, d, "Visit 3")
  p_4 <- generate_timepoint_plot(V4_data, d, "Visit 4")
  
  # Combine and save
  if (!is.null(p_1) && !is.null(p_2) && !is.null(p_3) && !is.null(p_4)) {
    combined_plot <- p_1 + p_2 + p_3 + p_4 + plot_layout(ncol = 4, widths = c(1, 1))
    ggsave(
      filename = file.path(export_dir, paste0("Donor", d, "_seurat_clusters_stacked_bars.png")),
      plot = combined_plot,
      width = 21,
      height = 0.3 * length(unique(donor_meta$Global_Clone_ID)) + 2
    )
    message("Saved side-by-side plot for Donor ", d)
  }
}


# Function to generate a plot for one donor and one timepoint (pp65 cutoff at 20)
generate_timepoint_plot_cutoff <- function(df, d, time_label) {
  # Define clone labels with pp65 cutoff
  ie1_labels <- paste0("IE1 Clone ", 1:7)
  pp65_labels <- paste0("pp65 Clone ", 1:20)   # cutoff at 20
  
  # Add a dummy label between IE1 and pp65
  dummy_label <- " "
  all_labels <- c(ie1_labels, dummy_label, pp65_labels)
  
  # Summarize clone sizes
  clone_summary <- df %>%
    count(specificity, clone_number, name = "n_cells_clone") %>%
    group_by(specificity) %>%
    mutate(clone_size_percent = n_cells_clone / sum(n_cells_clone) * 100) %>%
    ungroup() %>%
    mutate(clone_label = ifelse(specificity == "IE1",
                                paste0("IE1 Clone ", clone_number),
                                paste0("pp65 Clone ", clone_number)))
  
  # Cluster composition per clone
  comp_data <- df %>%
    count(specificity, clone_number, cluster, name = "n_cells") %>%
    group_by(specificity, clone_number) %>%
    mutate(freq = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    mutate(clone_label = ifelse(specificity == "IE1",
                                paste0("IE1 Clone ", clone_number),
                                paste0("pp65 Clone ", clone_number))) %>%
    # Add missing clones with 0 cells (cutoff applied here too)
    complete(clone_label = c(ie1_labels, pp65_labels), cluster = sort(unique(D0_data_CMV$seurat_clusters)),
             fill = list(n_cells = 0, freq = 0)) %>%
    left_join(clone_summary %>% select(clone_label, clone_size_percent), by = "clone_label") %>%
    mutate(clone_size_percent = replace_na(clone_size_percent, 0))
  
  # Add dummy row
  dummy_row <- expand.grid(
    clone_label = dummy_label,
    cluster = sort(unique(D0_data_CMV$seurat_clusters)),
    n_cells = 0,
    freq = 0,
    clone_size_percent = 0
  )
  comp_data <- bind_rows(comp_data, dummy_row)
  
  # Label colors
  comp_data$label_color <- ifelse(grepl("IE1", comp_data$clone_label), "red",
                                  ifelse(grepl("pp65", comp_data$clone_label), "blue", "black"))
  
  # Factor ordering with dummy row in middle
  comp_data$clone_label <- factor(comp_data$clone_label, levels = rev(all_labels))
  
  # Dashed line position
  dummy_position <- which(levels(comp_data$clone_label) == dummy_label)
  ie1_break_position <- dummy_position
  
  # ---- Plot ----
  ggplot(comp_data, aes(x = freq, y = clone_label, fill = factor(cluster))) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(
      data = comp_data %>% filter(clone_label != dummy_label) %>% group_by(clone_label) %>% slice(1),
      aes(x = 1.02, y = clone_label, label = sprintf("%.1f%%", clone_size_percent), color = label_color),
      inherit.aes = FALSE,
      hjust = 0,
      size = 3
    ) +
    scale_fill_manual(values = cluster_colors, drop = FALSE) +
    scale_color_identity() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15)), labels = scales::percent_format(accuracy = 1)) +
    theme_classic() +
    theme(
      axis.title = element_blank(),
      axis.text.y = element_text(size = 8, color = rev(c(rep("red", 7), "black", rep("blue", 20)))),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position = "right",
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 8),
      plot.margin = margin(10, 30, 10, 10)
    ) +
    guides(fill = guide_legend(title = "Metacluster")) +
    geom_hline(yintercept = ie1_break_position, linetype = "dashed") +
    ggtitle(paste("Donor", d, "-", time_label, "Timepoint"))
}


# ---- Loop through donors and generate plots with pp65 cutoff ----
for (d in sort(unique(meta$donor))) {
  donor_meta <- meta %>% filter(donor == d)
  
  # Apply cutoff: keep all IE1, pp65 <= 20
  donor_meta_cut <- donor_meta %>%
    filter(specificity == "IE1" | (specificity == "pp65" & clone_number <= 20))
  
  V1_data <- donor_meta_cut %>% filter(timepoint == "Visit 1")
  V2_data <- donor_meta_cut %>% filter(timepoint == "Visit 2")
  V3_data <- donor_meta_cut %>% filter(timepoint == "Visit 3")
  V4_data <- donor_meta_cut %>% filter(timepoint == "Visit 4")
  
  p_1 <- generate_timepoint_plot_cutoff(V1_data, d, "Visit 1")
  p_2 <- generate_timepoint_plot_cutoff(V2_data, d, "Visit 2")
  p_3 <- generate_timepoint_plot_cutoff(V3_data, d, "Visit 3")
  p_4 <- generate_timepoint_plot_cutoff(V4_data, d, "Visit 4")
  
  # Combine and save
  if (!is.null(p_1) && !is.null(p_2) && !is.null(p_3) && !is.null(p_4)) {
    combined_plot <- p_1 + p_2 + p_3 + p_4 + plot_layout(ncol = 4, widths = c(1, 1))
    ggsave(
      filename = file.path(export_dir, paste0("Donor", d, "_seurat_clusters_stacked_bars_pp65_cutoff.png")),
      plot = combined_plot,
      width = 21,
      height = 0.3 * length(unique(donor_meta_cut$Global_Clone_ID)) + 2
    )
    message("Saved side-by-side plot (pp65 cutoff) for Donor ", d)
  }
}





#
# Script to plot donor ages and tetramer % data
#

# Load data
df <- read_csv("Donor_age_tet_percents.csv")

# Rename messy columns to convenient names (adjust if your headers differ)
colnames(df)[2:4] <- c("Age", "pp65", "IE1")

# --- Manually populate Donor and Timepoint based on row order ---
# You stated the 10 rows are ordered:
# Donor 1 - ET, Donor 1 - LT, Donor 2 - ET, Donor 2 - LT, ... Donor 5 - LT
df <- df %>%
  mutate(
    Donor = rep(paste0("Donor ", 1:5), each = 2),
    Timepoint = rep(c("ET", "LT"), times = 5)
  ) %>%
  mutate(
    Donor = factor(Donor, levels = paste0("Donor ", 1:5)),
    Timepoint = factor(Timepoint, levels = c("ET", "LT")),
    Donor_num = as.numeric(Donor)  # y positions 1..5 for plotting
  )

# --- Color palette for donors ---
donors <- levels(df$Donor)
pal <- brewer.pal(max(3, length(donors)), "Set1")   # Set1 is a good qualitative palette
pal <- rep(pal, length.out = length(donors))       # ensure length
donor_cols <- setNames(pal, donors)

# --- Plot 1: Age (x = Age, y = donor row 1..5) ---
p_age <- ggplot() +
  # horizontal line connecting the two timepoints per donor
  geom_line(data = df, aes(x = Age, y = Donor_num, group = Donor, color = Donor), size = 0.9) +
  # ET = empty circle (white fill), outline = donor color
  geom_point(data = filter(df, Timepoint == "ET"),
             aes(x = Age, y = Donor_num, color = Donor),
             shape = 21, fill = "white", size = 4, stroke = 1) +
  # LT = filled circle (fill = donor color), outline = donor color
  geom_point(data = filter(df, Timepoint == "LT"),
             aes(x = Age, y = Donor_num, color = Donor, fill = Donor),
             shape = 21, size = 4, stroke = 1) +
  scale_y_continuous(breaks = 1:5, labels = donors, expand = expansion(add = c(0.5, 0.5)), trans = "reverse") +
  scale_color_manual(values = donor_cols) +
  scale_fill_manual(values = donor_cols, guide = FALSE) + # hide duplicate fill legend
  labs(x = "Age", y = "", title = "Donor ages at ET and LT") +
  theme_bw()

# Save age plot
ggsave("exports/integrated_v5/Donor_age_graphs/Donor_Age_ET_LT.png", p_age, width = 5.5, height = 4, dpi = 300)

# --- Plot 2: Tetramer comparison (%pp65 on x, %IE1 on y) ---

# Main plot with donor points only
main_plot <- ggplot() +
  geom_line(data = df, aes(x = log10(pp65), y = log10(IE1), group = Donor, color = Donor), size = 0.9) +
  geom_point(data = filter(df, Timepoint == "ET"),
             aes(x = log10(pp65), y = log10(IE1), color = Donor),
             shape = 21, fill = "white", size = 4, stroke = 1) +
  geom_point(data = filter(df, Timepoint == "LT"),
             aes(x = log10(pp65), y = log10(IE1), color = Donor, fill = Donor),
             shape = 21, size = 4, stroke = 1) +
  scale_color_manual(values = donor_cols) +
  scale_fill_manual(values = donor_cols) +
  scale_x_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1)) +
  scale_y_continuous(limits = c(-1,1), breaks = c(-1,-0.5,0,0.5,1)) +
  labs(x = "% pp65-tetramer+ of CD8 (log10)",
       y = "% IE1-tetramer+ of CD8 (log10)",
       title = "Tetramer+ (ET vs LT) per donor") +
  theme_bw() +
  theme(legend.title = element_blank()) +
  guides(color = guide_legend(override.aes = list(shape = 21, fill = donor_cols, color = donor_cols, size = 4)))

# Dummy plot for ET vs LT legend
legend_df <- tibble(
  Timepoint = factor(c("Early Timepoint", "Late Timepoint")),
  x = 1:2,
  y = 1:2,
  fillcol = c("white", "black")
)

legend_plot <- ggplot(legend_df, aes(x = x, y = y, shape = Timepoint, fill = fillcol)) +
  geom_point(color = "black", size = 4) +
  scale_shape_manual(values = c("Early Timepoint" = 21, "Late Timepoint" = 21)) +
  scale_fill_identity() +
  theme_void() +
  theme(legend.position = "bottom")

# Extract the legend
etlt_legend <- cowplot::get_legend(legend_plot)

# Combine main plot and ET/LT legend
final_plot <- cowplot::plot_grid(main_plot, etlt_legend, ncol = 1, rel_heights = c(1, 0.1))

# Save tetramer plot
ggsave("exports/integrated_v5/Donor_age_graphs/Tetramer_Percent_ET_LT.png", final_plot, width = 6.5, height = 5, dpi = 300)




#
# Graph proportion in each metacluster for each donor/timepoint combo, for full dataset
#

# 1. Extract relevant metadata
meta <- mdata@meta.data %>%
  dplyr::select(donor_identity, timepoint, metacluster)

# 2. Count cells per cluster for each donor-timepoint combination
cluster_counts <- meta %>%
  group_by(donor_identity, timepoint, metacluster) %>%
  tally(name = "cell_count") %>%
  ungroup()

# 3. Compute percentages
cluster_counts <- cluster_counts %>%
  group_by(donor_identity, timepoint) %>%
  mutate(percent = cell_count / sum(cell_count) * 100) %>%
  ungroup()

# 4. Create donor-timepoint labels
cluster_counts$donor_timepoint <- paste(cluster_counts$donor_identity, cluster_counts$timepoint, sep = " - ")

# 5. Insert blank rows for gaps between donors
donor_order <- c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5")
timepoints <- c("Early", "Late")

# Build desired x-axis order including gaps
x_levels <- c()
for(d in donor_order){
  for(t in timepoints){
    x_levels <- c(x_levels, paste(d, t, sep = " - "))
  }
  x_levels <- c(x_levels, paste0(d, " - Gap"))  # add gap after each donor
}

# Add a dummy row for each gap to ensure ggplot leaves space
gap_rows <- data.frame(
  donor_identity = NA,
  timepoint = NA,
  metacluster = NA,
  cell_count = 0,
  percent = 0,
  donor_timepoint = paste0(donor_order, " - Gap")
)
cluster_counts <- bind_rows(cluster_counts, gap_rows)

# Set factor levels for x-axis
cluster_counts$donor_timepoint <- factor(cluster_counts$donor_timepoint, levels = x_levels)

# 6. Optional: order clusters consistently
cluster_counts$metacluster <- factor(cluster_counts$metacluster, 
                                         levels = sort(unique(meta$metacluster)))

# 7. Plot stacked bar graph
p <- ggplot(cluster_counts, aes(x = donor_timepoint, y = percent, fill = metacluster)) +
  geom_bar(stat = "identity") +
  labs(x = "Donor - Timepoint", y = "Percent of Cells", fill = "Seurat Cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 8. Export plot
ggsave("exports/integrated_v5/total_UMAP_stacked_bar_metaclusters_per_donor_timepoint.png", p, width = 18, height = 6)


# Proportion in each metacluster for each donor/timepoint/specificity combo (full dataset)

# 1. Extract relevant metadata (include specificity)
meta <- mdata@meta.data %>%
  select(donor_identity, timepoint, specificity, metacluster)

# 2. Count cells per cluster for each donor-timepoint-specificity combination
cluster_counts <- meta %>%
  group_by(donor_identity, timepoint, specificity, metacluster) %>%
  tally(name = "cell_count") %>%
  ungroup()

# 3. Compute percentages
cluster_counts <- cluster_counts %>%
  group_by(donor_identity, timepoint, specificity) %>%
  mutate(percent = cell_count / sum(cell_count) * 100) %>%
  ungroup()

# 4. Create donor-timepoint-specificity labels
cluster_counts$donor_condition <- paste(cluster_counts$donor_identity,
                                        cluster_counts$timepoint,
                                        cluster_counts$specificity,
                                        sep = " - ")

# 5. Build desired x-axis order
donor_order <- c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5")
timepoints <- c("Early", "Late")
specificities <- c("Bulk", "pp65", "IE1")

x_levels <- c()
donor_gap_positions <- c()   # track only donor-wide gaps
for(d in donor_order){
  # Early specificities
  for(s in specificities){
    x_levels <- c(x_levels, paste(d, "Early", s, sep = " - "))
  }
  # small gap after Early (no dashed line)
  early_gap_label <- paste0(d, " - EarlyGap")
  x_levels <- c(x_levels, early_gap_label)
  
  # Late specificities
  for(s in specificities){
    x_levels <- c(x_levels, paste(d, "Late", s, sep = " - "))
  }
  
  # large donor-wide gap (dashed line)
  donor_gap_label <- paste0(d, " - DonorGap")
  x_levels <- c(x_levels, donor_gap_label)
  donor_gap_positions <- c(donor_gap_positions, donor_gap_label)
}

# 6. Add dummy rows for gaps
gap_rows <- data.frame(
  donor_identity = NA,
  timepoint = NA,
  specificity = NA,
  metacluster = NA,
  cell_count = 0,
  percent = 0,
  donor_condition = setdiff(x_levels, unique(cluster_counts$donor_condition))
)

cluster_counts <- bind_rows(cluster_counts, gap_rows)

# 7. Set factor levels for x-axis
cluster_counts$donor_condition <- factor(cluster_counts$donor_condition, 
                                         levels = x_levels)

# 8. Explicit metacluster order
metacluster_order <- c("Naive", "Central Memory", "Effector Memory", 
                       "NK-like", "ISG+", "MAIT")
cluster_counts$metacluster <- factor(cluster_counts$metacluster,
                                     levels = metacluster_order)

# 9. Plot stacked bar graph
p <- ggplot(cluster_counts, aes(x = donor_condition, y = percent, fill = metacluster)) +
  geom_bar(stat = "identity") +
  labs(x = "Donor - Timepoint - Specificity", 
       y = "Percent of Cells", 
       fill = "Metacluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 10. Add dashed lines at donor-wide gaps only
gap_indices <- which(levels(cluster_counts$donor_condition) %in% donor_gap_positions)
p <- p + geom_vline(xintercept = gap_indices, linetype = "dashed", color = "black")

# 11. Strip labels for all gap columns
labels_clean <- levels(cluster_counts$donor_condition)
labels_clean[grepl("Gap", labels_clean)] <- ""
p <- p + scale_x_discrete(labels = labels_clean)

# 12. Export plot
ggsave("exports/integrated_v5/total_UMAP_stacked_bar_metaclusters_per_donor_timepoint_specificity.png", 
       p, width = 22, height = 6)




#
# Graph proportion in each metacluster for each donor/specificity combo, for CMV dataset
#

# 1. Extract relevant metadata
meta <- mdata_CMV@meta.data %>%
  select(donor_identity, specificity, metacluster)

# 2. Count cells per cluster for each donor-specificity combination
cluster_counts <- meta %>%
  group_by(donor_identity, specificity, metacluster) %>%
  tally(name = "cell_count") %>%
  ungroup()

# 3. Compute percentages
cluster_counts <- cluster_counts %>%
  group_by(donor_identity, specificity) %>%
  mutate(percent = cell_count / sum(cell_count) * 100) %>%
  ungroup()

# 4. Create donor-specificity labels
cluster_counts$donor_specificity <- paste(cluster_counts$donor_identity, cluster_counts$specificity, sep = " - ")

# 5. Insert blank rows for gaps between donors
donor_order <- c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5")
specificities <- c("pp65", "IE1")

# Build desired x-axis order including gaps
x_levels <- c()
for(d in donor_order){
  for(t in specificities){
    x_levels <- c(x_levels, paste(d, t, sep = " - "))
  }
  x_levels <- c(x_levels, paste0(d, " - Gap"))  # add gap after each donor
}

# Add a dummy row for each gap to ensure ggplot leaves space
gap_rows <- data.frame(
  donor_identity = NA,
  specificity = NA,
  metacluster = NA,
  cell_count = 0,
  percent = 0,
  donor_specificity = paste0(donor_order, " - Gap")
)
cluster_counts <- bind_rows(cluster_counts, gap_rows)

# Set factor levels for x-axis
cluster_counts$donor_specificity <- factor(cluster_counts$donor_specificity, levels = x_levels)

# 6. Explicit metacluster order
metacluster_order <- c("Tcm", "Tem", "T-CMC1", 
                       "NK-like", "ISG+")
cluster_counts$metacluster <- factor(cluster_counts$metacluster,
                                     levels = metacluster_order)

# 7. Plot stacked bar graph
p <- ggplot(cluster_counts, aes(x = donor_specificity, y = percent, fill = metacluster)) +
  geom_bar(stat = "identity") +
  labs(x = "Donor - Timepoint", y = "Percent of Cells", fill = "Seurat Cluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 8. Export plot
ggsave("exports/integrated_v5/CMV_stacked_bar_metaclusters_per_donor_timepoint.png", p, width = 12, height = 6)


# Proportion in each metacluster for each donor/timepoint/specificity combo (CMV dataset)

# 1. Extract relevant metadata (include specificity)
meta <- mdata_CMV@meta.data %>%
  select(donor_identity, timepoint, specificity, metacluster)

# 2. Count cells per cluster for each donor-timepoint-specificity combination
cluster_counts <- meta %>%
  group_by(donor_identity, timepoint, specificity, metacluster) %>%
  tally(name = "cell_count") %>%
  ungroup()

# 3. Compute percentages
cluster_counts <- cluster_counts %>%
  group_by(donor_identity, timepoint, specificity) %>%
  mutate(percent = cell_count / sum(cell_count) * 100) %>%
  ungroup()

# 4. Create donor-timepoint-specificity labels
cluster_counts$donor_condition <- paste(cluster_counts$donor_identity,
                                        cluster_counts$timepoint,
                                        cluster_counts$specificity,
                                        sep = " - ")

# 5. Build desired x-axis order
donor_order <- c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5")
timepoints <- c("Early", "Late")
specificities <- c("pp65", "IE1")

x_levels <- c()
donor_gap_positions <- c()   # track only donor-wide gaps
for(d in donor_order){
  # Early specificities
  for(s in specificities){
    x_levels <- c(x_levels, paste(d, "Early", s, sep = " - "))
  }
  # small gap after Early (no dashed line)
  early_gap_label <- paste0(d, " - EarlyGap")
  x_levels <- c(x_levels, early_gap_label)
  
  # Late specificities
  for(s in specificities){
    x_levels <- c(x_levels, paste(d, "Late", s, sep = " - "))
  }
  
  # large donor-wide gap (dashed line)
  donor_gap_label <- paste0(d, " - DonorGap")
  x_levels <- c(x_levels, donor_gap_label)
  donor_gap_positions <- c(donor_gap_positions, donor_gap_label)
}

# 6. Add dummy rows for gaps
gap_rows <- data.frame(
  donor_identity = NA,
  timepoint = NA,
  specificity = NA,
  metacluster = NA,
  cell_count = 0,
  percent = 0,
  donor_condition = setdiff(x_levels, unique(cluster_counts$donor_condition))
)

cluster_counts <- bind_rows(cluster_counts, gap_rows)

# 7. Set factor levels for x-axis
cluster_counts$donor_condition <- factor(cluster_counts$donor_condition, 
                                         levels = x_levels)

# 8. Explicit metacluster order
metacluster_order <- c("Tcm", "Tem", "T-CMC1", 
                       "NK-like", "ISG+")
cluster_counts$metacluster <- factor(cluster_counts$metacluster,
                                     levels = metacluster_order)

# 9. Plot stacked bar graph
p <- ggplot(cluster_counts, aes(x = donor_condition, y = percent, fill = metacluster)) +
  geom_bar(stat = "identity") +
  labs(x = "Donor - Timepoint - Specificity", 
       y = "Percent of Cells", 
       fill = "Metacluster") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 10. Add dashed lines at donor-wide gaps only
gap_indices <- which(levels(cluster_counts$donor_condition) %in% donor_gap_positions)
p <- p + geom_vline(xintercept = gap_indices, linetype = "dashed", color = "black")

# 11. Strip labels for all gap columns
labels_clean <- levels(cluster_counts$donor_condition)
labels_clean[grepl("Gap", labels_clean)] <- ""
p <- p + scale_x_discrete(labels = labels_clean)

# 12. Export plot
ggsave("exports/integrated_v5/CMV_stacked_bar_metaclusters_per_donor_timepoint_specificity.png", 
       p, width = 18, height = 6)




#
# Stacked bar graph of percent metacluster for pp65 and IE1
#
# Extract metadata
meta <- mdata_CMV@meta.data

# Get metacluster levels in the order used by DimPlot
metacluster_levels <- levels(factor(meta$metacluster))

# Example DimPlot to extract colors
dp <- DimPlot(mdata_CMV, group.by = "metacluster", reduction = "umap")

# Extract the mapping of colors from the ggplot object
# In Seurat v5, the color scale is stored in dp$layers[[1]]$mapping$colour
# But safer is to use the scales via the ggplot_build object

ggb <- ggplot_build(dp)  # build the plot

# The first layer contains the fill colors actually used
# Extract unique colors assigned to each metacluster
layer_data <- ggb$data[[1]]
color_mapping <- unique(layer_data$colour)  # colors actually used in the plot

# Make sure colors are in order of metacluster levels
metacluster_levels <- levels(factor(mdata_CMV$metacluster))
names(color_mapping) <- metacluster_levels  # assign names

# Reverse colors
reversed_colors <- rev(color_mapping)

# Prepare plotting data
plot_data <- meta %>%
  group_by(specificity, metacluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(specificity) %>%
  mutate(percent = n / sum(n) * 100)

# Set desired order
plot_data$metacluster <- factor(
  plot_data$metacluster,
  levels = c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")
)

# Combined stacked bar plot using exact DimPlot colors reversed
combined_plot <- ggplot(plot_data, aes(x = specificity, y = percent, fill = metacluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.6) +
  labs(title = "Cell Phenotype Composition - IE1 vs pp65",
       x = "Specificity", y = "Percent of cells", fill = "Cell Phenotypes") +
  theme_classic(base_size = 14)

# Save to PNG
ggsave("exports/integrated_v5/CMV_metacluster_by_specificity.png", combined_plot,
       width = 6, height = 5, dpi = 300)



#
# Calculate phenotype (metacluster) heterogeneity within each clone (both timepoints combined)
#

# --- Extract and parse metadata ---
meta_parsed <- mdata_CMV@meta.data %>%
  as.data.frame() %>%
  filter(!is.na(Global_Clone_ID)) %>%
  rename(clone_id = Global_Clone_ID) %>%
  mutate(
    donor_identity = str_extract(clone_id, "^D\\d+"),
    specificity    = str_extract(clone_id, "(?<=_)pp65|IE1"),
    clone_number   = as.integer(str_extract(clone_id, "(?<=clone-)\\d+"))
  )

# --- Sanity check: clone_id must map to exactly one donor & specificity ---
stopifnot(
  meta_parsed %>%
    distinct(clone_id, donor_identity, specificity) %>%
    count(clone_id) %>%
    filter(n > 1) %>%
    nrow() == 0
)

# --- Count cells per clone × metacluster ---
clone_counts <- meta_parsed %>%
  count(clone_id, metacluster, name = "n_cells")

# --- Define Gini–Simpson heterogeneity metric ---
gini_simpson <- function(counts) {
  p <- counts / sum(counts)
  1 - sum(p^2)
}

# --- Compute clone-level heterogeneity ---
clone_metrics <- clone_counts %>%
  group_by(clone_id) %>%
  summarise(
    total_cells    = sum(n_cells),
    n_metaclusters = n(),
    gini_simpson   = gini_simpson(n_cells),
    .groups = "drop"
  ) %>%
  left_join(
    meta_parsed %>%
      distinct(clone_id, donor_identity, specificity, clone_number),
    by = "clone_id"
  )

# --- Compute percent of donor × specificity stratum ---
clone_metrics <- clone_metrics %>%
  group_by(donor_identity, specificity) %>%
  mutate(
    pct_of_stratum = total_cells / sum(total_cells) * 100
  ) %>%
  ungroup()

# --- Set specificity order for plotting ---
clone_metrics$specificity <- factor(
  clone_metrics$specificity,
  levels = c("pp65", "IE1")
)

# --- Final sanity check: one row per clone ---
stopifnot(
  nrow(clone_metrics) == n_distinct(clone_metrics$clone_id)
)

# --- Write output ---
write.csv(
  clone_metrics,
  file = "exports/integrated_v5/clone_heterogeneity_plots/clone_metrics.csv",
  row.names = FALSE
)




# --- Scatter plot: normalized Shannon entropy vs clone size ---
Pielous_J_plot <- ggplot() +
  # triangles, slightly smaller
  geom_point(data = subset(clone_metrics, donor_identity %in% c("Donor 3", "Donor 4")),
             aes(x = pct_of_stratum, y = entropy,
                 color = specificity, fill = specificity,
                 shape = donor_identity),
             size = 2.8, alpha = 0.9, stroke = 0.6) +
  # other donors, normal size
  geom_point(data = subset(clone_metrics, !donor_identity %in% c("Donor 3", "Donor 4")),
             aes(x = pct_of_stratum, y = entropy,
                 color = specificity, fill = specificity,
                 shape = donor_identity),
             size = 3.5, alpha = 0.9, stroke = 0.6) +
  scale_color_manual(name = "Specificity",
                     values = c("IE1" = "red", "pp65" = "blue")) +
  scale_fill_manual(name = "Specificity",
                    values = c("IE1" = "red", "pp65" = "blue")) +
  scale_shape_manual(name = "Donor",
                     values = c("Donor 1" = 16,
                                "Donor 2" = 15,
                                "Donor 3" = 17,
                                "Donor 4" = 25,
                                "Donor 5" = 18)) +
  guides(
    shape = guide_legend(override.aes = list(
      fill = c("black","black","black","black","black"),
      color = c("black","black","black","black","black"),
      size = c(3.5, 3.5, 2.8, 2.8, 3.5)
    ))
  ) +
  labs(x = "Clone size (% of pp65- or IE1-specific major clones)",
       y = "Normalized Shannon entropy (Pielou's J)") +
  theme_classic(base_size = 14)

# Export as PNG
ggsave("exports/integrated_v5/clone_heterogeneity_plots/Normalized_Shannon_entropy_Pielou's_J.png", Pielous_J_plot,
       width = 8, height = 6, dpi = 300)

# --- Scatter plot: Gini-Simpson index vs clone size ---
Gini_Simpson_plot <- ggplot() +
  # triangles, slightly smaller
  geom_point(data = subset(clone_metrics, donor_identity %in% c("Donor 3", "Donor 4")),
             aes(x = pct_of_stratum, y = gini_simpson,
                 color = specificity, fill = specificity,
                 shape = donor_identity),
             size = 2.8, alpha = 0.9, stroke = 0.6) +
  # other donors, normal size
  geom_point(data = subset(clone_metrics, !donor_identity %in% c("Donor 3", "Donor 4")),
             aes(x = pct_of_stratum, y = gini_simpson,
                 color = specificity, fill = specificity,
                 shape = donor_identity),
             size = 3.5, alpha = 0.9, stroke = 0.6) +
  scale_color_manual(name = "Specificity",
                     values = c("IE1" = "red", "pp65" = "blue")) +
  scale_fill_manual(name = "Specificity",
                    values = c("IE1" = "red", "pp65" = "blue")) +
  scale_shape_manual(name = "Donor",
                     values = c("Donor 1" = 16,
                                "Donor 2" = 15,
                                "Donor 3" = 17,
                                "Donor 4" = 25,
                                "Donor 5" = 18)) +
  guides(
    shape = guide_legend(override.aes = list(
      fill = rep("black", 5),
      color = rep("black", 5),
      size = c(3.5, 3.5, 2.8, 2.8, 3.5)
    ))
  ) +
  labs(x = "Clone size (% of pp65- or IE1-specific major clones)",
       y = "Gini-Simpson index") +
  theme_classic(base_size = 14)

# Export as PNG
ggsave("exports/integrated_v5/clone_heterogeneity_plots/Gini_Simpson_index.png",
       Gini_Simpson_plot, width = 8, height = 6, dpi = 300)


# --- Scatter plot: raw Shannon entropy vs clone size ---
Shannon_entropy_plot <- ggplot() +
  # triangles, slightly smaller
  geom_point(data = subset(clone_metrics, donor_identity %in% c("Donor 3", "Donor 4")),
             aes(x = pct_of_stratum, y = entropy * log(length(unique(meta_filtered$metacluster))),
                 color = specificity, fill = specificity,
                 shape = donor_identity),
             size = 2.8, alpha = 0.9, stroke = 0.6) +
  # other donors, normal size
  geom_point(data = subset(clone_metrics, !donor_identity %in% c("Donor 3", "Donor 4")),
             aes(x = pct_of_stratum, y = entropy * log(length(unique(meta_filtered$metacluster))),
                 color = specificity, fill = specificity,
                 shape = donor_identity),
             size = 3.5, alpha = 0.9, stroke = 0.6) +
  scale_color_manual(name = "Specificity",
                     values = c("IE1" = "red", "pp65" = "blue")) +
  scale_fill_manual(name = "Specificity",
                    values = c("IE1" = "red", "pp65" = "blue")) +
  scale_shape_manual(name = "Donor",
                     values = c("Donor 1" = 16,
                                "Donor 2" = 15,
                                "Donor 3" = 17,
                                "Donor 4" = 25,
                                "Donor 5" = 18)) +
  guides(
    shape = guide_legend(override.aes = list(
      fill = rep("black", 5),
      color = rep("black", 5),
      size = c(3.5, 3.5, 2.8, 2.8, 3.5)
    ))
  ) +
  labs(x = "Clone size (% of pp65- or IE1-specific major clones)",
       y = "Shannon entropy (raw)") +
  theme_classic(base_size = 14)

# Export as PNG
ggsave("exports/integrated_v5/clone_heterogeneity_plots/Shannon_entropy_raw.png",
       Shannon_entropy_plot, width = 8, height = 6, dpi = 300)


# --- Scatter plot: Gini–Simpson index grouped by specificity ---

Gini_Simpson_grouped_plot <- ggplot(clone_metrics, 
                                    aes(x = specificity, 
                                        y = gini_simpson,
                                        color = specificity, 
                                        fill = specificity,
                                        shape = donor_identity,
                                        size = pct_of_stratum)) +
  geom_jitter(width = 0.2, alpha = 0.9, stroke = 0.6) +
  scale_color_manual(name = "specificity",
                     values = c("IE1" = "red", "pp65" = "blue")) +
  scale_fill_manual(name = "specificity",
                    values = c("IE1" = "red", "pp65" = "blue")) +
  scale_shape_manual(name = "donor_identity",
                     values = c("D1" = 16,
                                "D2" = 15,
                                "D3" = 17,
                                "D4" = 25,
                                "D5" = 18)) +
  scale_size_continuous(name = "pct_of_stratum",
                        range = c(2, 8)) +
  guides(
    shape = guide_legend(override.aes = list(
      fill = rep("black", 5),
      color = rep("black", 5),
      size = 4
    )),
    size = guide_legend(order = 3)
  ) +
  labs(x = "Specificity",
       y = "Gini–Simpson index") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )

# --- Export as PNG ---
ggsave("exports/integrated_v5/clone_heterogeneity_plots/Gini_Simpson_grouped_by_specificity.png",
       Gini_Simpson_grouped_plot, width = 6, height = 7, dpi = 300)


# --- Scatter plot: Gini–Simpson index grouped by specificity (no size scaling) ---

Gini_Simpson_grouped_no_sizing <- ggplot(clone_metrics, 
                                         aes(x = specificity, 
                                             y = gini_simpson,
                                             color = specificity, 
                                             fill = specificity,
                                             shape = donor_identity)) +
  geom_jitter(width = 0.2, size = 3.5, alpha = 0.9, stroke = 0.6) +
  scale_color_manual(name = "Specificity",
                     values = c("IE1" = "red", "pp65" = "blue")) +
  scale_fill_manual(name = "Specificity",
                    values = c("IE1" = "red", "pp65" = "blue")) +
  scale_shape_manual(name = "Donor",
                     values = c("Donor 1" = 16,
                                "Donor 2" = 15,
                                "Donor 3" = 17,
                                "Donor 4" = 25,
                                "Donor 5" = 18)) +
  guides(
    shape = guide_legend(override.aes = list(
      fill = rep("black", 5),
      color = rep("black", 5),
      size = 3.5
    ))
  ) +
  labs(x = "Specificity",
       y = "Gini–Simpson index") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )

# --- Export as PNG ---
ggsave("exports/integrated_v5/clone_heterogeneity_plots/Gini_Simpson_grouped_no_sizing.png",
       Gini_Simpson_grouped_no_sizing, width = 8, height = 6, dpi = 300)





#
# Jensen-Shannon divergence plot
#

# --- Import CSV ---
metacluster_stats <- read.csv("exports/integrated_v5/clone_cluster_pie_charts/metacluster_statistics/metacluster_statistics.csv")

# Convert numeric donor to factor with proper labels
metacluster_stats$donor <- factor(metacluster_stats$donor,
                                  levels = 1:5,
                                  labels = c("Donor 1", "Donor 2", "Donor 3", "Donor 4", "Donor 5"))

# --- Calculate clone frequency per donor & specificity ---
metacluster_stats <- metacluster_stats %>%
  group_by(donor, specificity) %>%
  mutate(total_cells = sum(n_cells_early + n_cells_late),
         clone_frequency = (n_cells_early + n_cells_late) / total_cells * 100) %>%
  ungroup() %>%
  mutate(one_minus_jsd = 1 - js_divergence)   # <-- new column

# Set order of specificity: pp65 on left, IE1 on right
metacluster_stats$specificity <- factor(metacluster_stats$specificity,
                                        levels = c("pp65", "IE1"))

# Create scatter plot
JS_plot <- ggplot(metacluster_stats,
                  aes(x = clone_frequency, y = one_minus_jsd,
                      color = specificity, fill = specificity,
                      shape = donor)) +
  geom_point(size = ifelse(metacluster_stats$donor %in% c("Donor 3","Donor 4"), 2.8, 3.5),
             alpha = 0.9, stroke = 0.6) +
  scale_shape_manual(name = "Donor",
                     values = c("Donor 1" = 16,
                                "Donor 2" = 15,
                                "Donor 3" = 17,
                                "Donor 4" = 25,
                                "Donor 5" = 18)) +
  scale_color_manual(name = "Specificity",
                     values = c("IE1" = "red", "pp65" = "blue")) +
  scale_fill_manual(name = "Specificity",
                    values = c("IE1" = "red", "pp65" = "blue")) +
  guides(
    shape = guide_legend(override.aes = list(
      fill = c("black","black","black","black","black"),
      color = c("black","black","black","black","black"),
      size = c(3.5, 3.5, 2.8, 2.8, 3.5)
    ))
  ) +
  labs(x = "Clone frequency (% of total pp65- or IE1-specific major clones)",
       y = "1 - Jensen-Shannon divergence (stability)") +
  theme_classic(base_size = 14)

# Export as PNG
ggsave("exports/integrated_v5/timepoint_stability_plots/Clonal_phenotype_Jensen-Shannon_divergence.png",
       JS_plot, width = 8, height = 6, dpi = 300)


# --- Create grouped plot ---
JS_grouped_plot <- ggplot(metacluster_stats,
                          aes(x = specificity,
                              y = one_minus_jsd,
                              color = specificity,
                              fill = specificity,
                              shape = donor,
                              size = clone_frequency)) +
  geom_jitter(width = 0.2, alpha = 0.9, stroke = 0.6) +
  scale_color_manual(name = "Specificity",
                     values = c("IE1" = "red", "pp65" = "blue")) +
  scale_fill_manual(name = "Specificity",
                    values = c("IE1" = "red", "pp65" = "blue")) +
  scale_shape_manual(name = "Donor",
                     values = c("Donor 1" = 16,
                                "Donor 2" = 15,
                                "Donor 3" = 17,
                                "Donor 4" = 25,
                                "Donor 5" = 18)) +
  scale_size_continuous(name = "Clone frequency",
                        range = c(2, 8)) +
  guides(
    shape = guide_legend(override.aes = list(
      fill = rep("black", 5),
      color = rep("black", 5),
      size = 4
    )),
    size = guide_legend(order = 3)
  ) +
  labs(x = "Specificity",
       y = "1 - Jensen–Shannon divergence (stability)") +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 12)
  )

JS_grouped_plot <- JS_grouped_plot +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  )

# --- Export as PNG ---
ggsave("exports/integrated_v5/timepoint_stability_plots/Clonal_phenotype_Jensen-Shannon_grouped_by_specificity.png",
       JS_grouped_plot, width = 6, height = 6, dpi = 300)



#
# Jensen-Shannon divergence for donor/specificity clonal dominance over time
#

# --- Step 0: subset metadata to valid specificity and clones ---
data <- mdata@meta.data %>%
  filter(!is.na(Global_Clone_ID),
         specificity %in% c("pp65", "IE1")) %>%
  # ensure clone names match specificity
  filter(
    (specificity == "pp65" & grepl("_pp65-clone-", Global_Clone_ID)) |
      (specificity == "IE1"  & grepl("_IE1-clone-",  Global_Clone_ID))
  )

# --- Step 1: count total cells per donor/specificity/timepoint ---
total_cells <- data %>%
  group_by(donor_identity, specificity, timepoint) %>%
  summarise(total_cells = n(), .groups = "drop")

# --- Step 2: count cells per donor/specificity/clone/timepoint ---
clone_counts <- data %>%
  group_by(donor_identity, specificity, Global_Clone_ID, timepoint) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  left_join(total_cells, by = c("donor_identity", "specificity", "timepoint")) %>%
  mutate(freq = n_cells / total_cells) %>%
  dplyr::select(-total_cells)

# --- Step 3: pivot to wide format ---
clone_wide <- clone_counts %>%
  pivot_wider(
    id_cols = c(donor_identity, specificity, Global_Clone_ID),
    names_from = timepoint,
    values_from = c(n_cells, freq),
    values_fill = 0
  ) %>%
  rename(
    n_cells_Early = n_cells_Early,
    n_cells_Late  = n_cells_Late,
    freq_Early    = freq_Early,
    freq_Late     = freq_Late
  )

# --- Step 4: compute JSD per donor/specificity ---
calc_jsd <- function(early, late) {
  mat <- rbind(early, late)
  JSD(mat, unit = "log2")
}

jsd_donor_spec <- clone_wide %>%
  group_by(donor_identity, specificity) %>%
  summarise(jsd = calc_jsd(freq_Early, freq_Late), .groups = "drop") %>%
  mutate(one_minus_jsd = 1 - jsd)   # <-- new column

# For new plot
plot_df <- jsd_donor_spec %>%
  mutate(
    donor = donor_identity,
    specificity = factor(specificity, levels = c("pp65", "IE1"))
  )



# Create donor-spec factor with spacing tweaks
plot_data <- jsd_donor_spec %>%
  mutate(
    # replace 0 for plotting
    stability_plot = ifelse(one_minus_jsd == 0, 1e-3, one_minus_jsd),
    # create factor for x-axis with a "dummy" space between pp65 and IE1
    donor_spec = paste(donor_identity, specificity, sep = "-"),
    donor_spec = factor(
      donor_spec,
      levels = c(
        "Donor 1-pp65", "Donor 2-pp65", "Donor 3-pp65", "Donor 4-pp65", "Donor 5-pp65",
        "Donor 1-IE1",  "Donor 2-IE1",  "Donor 3-IE1",  "Donor 4-IE1",  "Donor 5-IE1"
      )
    ),
    spacing = case_when(
      specificity == "pp65" ~ as.numeric(donor_spec),
      specificity == "IE1"  ~ as.numeric(donor_spec) + 0.5
    )
  )

# Define colors and legend order
spec_colors <- c("pp65" = "blue", "IE1" = "red")
plot_data$specificity <- factor(plot_data$specificity, levels = c("pp65", "IE1"))

# Plot
JS_plot <- ggplot(plot_data, aes(x = spacing, y = stability_plot, fill = specificity)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = spec_colors, guide = guide_legend(nrow = 2)) +
  scale_x_continuous(
    breaks = plot_data$spacing,
    labels = gsub("-pp65|-IE1", "", plot_data$donor_spec)
  ) +
  labs(
    x = NULL,
    y = "1 - Jensen-Shannon Divergence (stability)",
    fill = "Specificity"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )

# Export as PNG
ggsave("exports/integrated_v5/timepoint_stability_plots/Donor_clone_spread_Jensen-Shannon_divergence.png",
       JS_plot, width = 8, height = 6, dpi = 300)


# order donor/specificity as D1-pp65, D1-IE1, D2-pp65, D2-IE1, etc.
plot_data <- plot_data %>%
  arrange(donor_identity, specificity) %>%
  mutate(order = row_number())

# add spacing: give each donor 2 slots (pp65, IE1) with a gap between donors
plot_data <- plot_data %>%
  group_by(donor_identity) %>%
  mutate(spacing = (cur_group_id() - 1) * 3 + row_number()) %>%
  ungroup()

# plot
JS_plot <- ggplot(plot_data, aes(x = spacing, y = stability_plot, fill = specificity)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = spec_colors, guide = guide_legend(nrow = 2)) +
  scale_x_continuous(
    breaks = tapply(plot_data$spacing, plot_data$donor_identity, mean),  # center tick per donor
    labels = unique(plot_data$donor_identity)
  ) +
  labs(
    x = NULL,
    y = "1 - Jensen-Shannon Divergence (stability)",
    fill = "Specificity"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    legend.position = "right"
  )

# Export as PNG
ggsave("exports/integrated_v5/timepoint_stability_plots/Donor_clone_spread_Jensen-Shannon_divergence_split_by_donor.png",
       JS_plot, width = 8, height = 6, dpi = 300)


# New plot
JS_grouped_plot <- ggplot(
  plot_df,
  aes(
    x = specificity,
    y = one_minus_jsd,
    color = specificity,
    fill  = specificity,
    shape = donor
  )
) +
  geom_jitter(
    width = 0.15,
    alpha = 0.9,
    stroke = 0.6,
    size = 4
  ) +
  scale_color_manual(
    name = "Specificity",
    values = c("IE1" = "red", "pp65" = "blue")
  ) +
  scale_fill_manual(
    name = "Specificity",
    values = c("IE1" = "red", "pp65" = "blue")
  ) +
  scale_shape_manual(
    name = "Donor",
    values = c(
      "Donor 1" = 16,
      "Donor 2" = 15,
      "Donor 3" = 17,
      "Donor 4" = 25,
      "Donor 5" = 18
    )
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(
        fill  = rep("black", 5),
        color = rep("black", 5),
        size  = 4
      )
    )
  ) +
  labs(
    x = "Specificity",
    y = "1 - Jensen–Shannon divergence (stability)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 12)
  )

JS_grouped_plot <- JS_grouped_plot +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  )

ggsave("exports/integrated_v5/timepoint_stability_plots/Donor_clone_spread_Jensen-Shannon_divergence_dots_plot.png",
       JS_grouped_plot, width = 5.9, height = 5.9, dpi = 300)




#Scatter plot of clone frequency change

# Step 1: compute frequency change as percent
clone_scatter <- clone_wide %>%
  mutate(
    freq_change = abs((freq_Late - freq_Early) / freq_Early) * 100
  ) %>%
  filter(!is.na(Global_Clone_ID))

# Plot
freq_change_plot <- ggplot(clone_scatter,
                           aes(x = freq_Early, y = freq_change,
                               color = specificity, fill = specificity,
                               shape = donor_identity)) +
  geom_point(size = ifelse(clone_scatter$donor_identity %in% c("Donor 3","Donor 4"), 2.8, 3.5),
             alpha = 0.9, stroke = 0.6) +
  scale_shape_manual(name = "Donor",
                     values = c("Donor 1" = 16,
                                "Donor 2" = 15,
                                "Donor 3" = 17,
                                "Donor 4" = 25,
                                "Donor 5" = 18)) +
  scale_color_manual(name = "Specificity",
                     values = c("IE1" = "red", "pp65" = "blue")) +
  scale_fill_manual(name = "Specificity",
                    values = c("IE1" = "red", "pp65" = "blue")) +
  guides(
    shape = guide_legend(override.aes = list(
      fill = c("black","black","black","black","black"),
      color = c("black","black","black","black","black"),
      size = c(3.5, 3.5, 2.8, 2.8, 3.5)
    ))
  ) +
  labs(x = "Clone frequency (Early, % of donor-specificity pool)",
       y = "Clonal frequency change (% of Early frequency)") +
  theme_classic(base_size = 14) +
  scale_y_break(c(350, 700))   # <-- key line: creates the split

# Export as PNG
ggsave("exports/integrated_v5/timepoint_stability_plots/clone_frequency_change_plot.png",
       freq_change_plot, width = 8, height = 6, dpi = 300)




#
# Bubble plots for gene expression
#

# Total UMAP, by metacluster

# Define gene list by assay
genes_integrated <- c("CCR7", "SELL", "KLRG1", "GZMK", "GZMB", "KLRC2", 
                      "TYROBP", "IFIT3", "OASL", "KLRB1", "ZBTB16", "ZNF683", "IKZF2")
genes_FB <- c("CD45RA", "CD57")

# DotPlot can't directly handle two assays at once, so we'll stitch them
p1 <- DotPlot(mdata, features = genes_integrated, group.by = "metacluster", assay = "integrated")
p2 <- DotPlot(mdata, features = genes_FB, group.by = "metacluster", assay = "FB")

# Combine gene expression info into one DotPlot object
plot_data <- rbind(p1$data, p2$data)

# Order metaclusters (reversed for y-axis display)
metacluster_order <- c("Naive", "Central Memory", "Effector Memory", 
                       "NK-like", "ISG+", "MAIT")
plot_data$id <- factor(plot_data$id, levels = rev(metacluster_order))

# Order genes (kept left-to-right as written)
gene_order <- c("CD45RA", "SELL", "CCR7", "GZMK", "GZMB", "CD57",
                "KLRG1", "KLRC2", "TYROBP", "ZNF683", "IKZF2", "IFIT3", "OASL", "ZBTB16", "KLRB1")
plot_data$features.plot <- factor(plot_data$features.plot, levels = gene_order)

# Make the bubble plot (genes on x, metaclusters on y)
p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank()   # removes axis titles
  )

# Save
ggsave("exports/integrated_v5/bubble_plots/Total_bubbleplot_by_metacluster_flipped.png", p, 
       width = 10, height = 2.5, dpi = 300)


# Total UMAP, by seurat_cluster

# Define gene list by assay
genes_integrated <- c("CCR7", "SELL", "KLRG1", "GZMK", "GZMB", "KLRC2", 
                      "TYROBP", "IFIT3", "KLRB1")
genes_FB <- c("CD45RA", "CD57")

# DotPlot can't directly handle two assays at once, so stitch them
p1 <- DotPlot(mdata, features = genes_integrated, group.by = "seurat_clusters", assay = "integrated")
p2 <- DotPlot(mdata, features = genes_FB, group.by = "seurat_clusters", assay = "FB")

# Combine gene expression info into one DotPlot object
plot_data <- rbind(p1$data, p2$data)

# Order clusters
cluster_order <- c("1", "5", "2", "8", "9", "0", "3", "10", "4", "7", "6")
plot_data$id <- factor(plot_data$id, levels = cluster_order)

# Order genes
gene_order <- c("CD45RA", "SELL", "CCR7", "GZMK", "GZMB", "CD57",
                "KLRG1", "KLRC2", "TYROBP", "IFIT3", "KLRB1")
plot_data$features.plot <- factor(plot_data$features.plot, levels = rev(gene_order))

# Make the bubble plot
p <- ggplot(plot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("lightgrey", "blue", "red")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  labs(x = "Seurat cluster", y = "Gene")

# Save
ggsave("exports/integrated_v5/Total_bubbleplot_by_seurat_cluster.png", p, 
       width = 8, height = 6, dpi = 300)



# CMV UMAP, by metacluster

# Define gene list by assay
genes_integrated <- c("CCR7", "SELL", "KLRG1", "GZMK", "GZMB", "KLRC2", 
                      "TYROBP", "IFIT3", "KLRB1", "ZNF683", "IKZF2", 
                      "LGALS1", "CMC1", "CD160", "XCL2", "HLA-DRA", "OASL", "IFIT1")
#genes_FB <- c("CD57")

# DotPlot can't directly handle two assays at once, so we'll stitch them
p1 <- DotPlot(mdata_CMV, features = genes_integrated, group.by = "metacluster", assay = "integrated")
#p2 <- DotPlot(mdata_CMV, features = genes_FB, group.by = "metacluster", assay = "FB")

# Combine gene expression info into one DotPlot object
#plot_data <- rbind(p1$data, p2$data)
plot_data <- p1$data

# Order metaclusters
metacluster_order <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")
plot_data$id <- factor(plot_data$id, levels = rev(metacluster_order))

# Order genes
gene_order <- c("CD45RA", "SELL", "CCR7", "GZMK", 
                "GZMB", #"CD57", 
                "KLRG1", "LGALS1",
                "CMC1", "CD160", "KLRB1", "XCL2", "HLA-DRA",
                "KLRC2", "TYROBP", "ZNF683", "IKZF2", 
                "IFIT1", "IFIT3", "OASL")
plot_data$features.plot <- factor(plot_data$features.plot, levels = gene_order)

# Make the bubble plot (genes on x, metaclusters on y)
p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank()   # removes axis titles
  )

# Save
ggsave("exports/integrated_v5/bubble_plots/CMV_bubbleplot_by_metacluster_flipped_tall_v2.png", p, 
       width = 8, height = 2, dpi = 300)


# CMV UMAP, by seurat_cluster

# Define gene list by assay
genes_integrated <- c("CCR7", "SELL", "KLRG1", "GZMK", "GZMB", "KLRC2", 
                      "TYROBP", "IFIT3", "KLRB1", "ZNF683", "IKZF2", 
                      "LGALS1", "CMC1", "CD160", "IFIT1", "OASL")
#genes_FB <- c("CD45RA", "CD57")

# DotPlot can't directly handle two assays at once, so stitch them
p1 <- DotPlot(mdata_CMV, features = genes_integrated, group.by = "seurat_clusters", assay = "integrated")
#p2 <- DotPlot(mdata_CMV, features = genes_FB, group.by = "seurat_clusters", assay = "FB")

# Combine gene expression info into one DotPlot object
#plot_data <- rbind(p1$data, p2$data)
plot_data <- p1$data

# Order clusters
cluster_order <- c("3", "4", "7", "0", "2", "8", "1", "5", "6")
plot_data$id <- factor(plot_data$id, levels = cluster_order)

# Order genes
gene_order <- c(#"CD45RA", 
                "SELL", "CCR7", "GZMK", 
                "GZMB", #"CD57", 
                "KLRG1", "LGALS1",
                "CMC1", "CD160", "KLRB1",
                "KLRC2", "TYROBP", "ZNF683", "IKZF2", "IFIT1", "IFIT3", "OASL")
plot_data$features.plot <- factor(plot_data$features.plot, levels = rev(gene_order))

# Make the bubble plot
p <- ggplot(plot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_text(size = 12)
  ) +
  labs(x = "Seurat cluster", y = "Gene")

# Save
ggsave("exports/integrated_v5/bubble_plots/CMV_bubbleplot_by_seurat_cluster_tall.png", p, 
       width = 8, height = 7.5, dpi = 300)




#
#
# Recalculating clonal dominance (pie charts etc) based on cloneId
#
#

# Gini-Simpson heterogeneity of clone frequency

# --- Step 1: Extract metadata and exclude "notCalled" cloneIds ---
meta_df <- mdata_CMV@meta.data %>%
  select(donor_identity, specificity, cloneId) %>%
  filter(!grepl("notCalled", cloneId))

# --- Step 2: Count cells per clone for each donor × specificity ---
clone_counts <- meta_df %>%
  group_by(donor_identity, specificity, cloneId) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 3: Compute Gini-Simpson index for each donor × specificity ---
gini_simpson <- clone_counts %>%
  group_by(donor_identity, specificity) %>%
  summarise(
    total_cells = sum(cell_count),
    gini_simpson = 1 - sum((cell_count / sum(cell_count))^2),
    .groups = "drop"
  ) %>%
  # --- Make specificity factor so pp65 is left, IE1 right ---
  mutate(specificity = factor(specificity, levels = c("pp65", "IE1")))

# --- Step 4: Plot with custom shapes/colors and save ---
gini_plot <- ggplot(gini_simpson,
                    aes(x = specificity,
                        y = gini_simpson,
                        color = specificity,
                        fill = specificity,
                        shape = donor_identity,
                        group = donor_identity)) +
  geom_line(color = "black", linewidth = 0.6, alpha = 0.7) +
  geom_point(size = ifelse(gini_simpson$donor_identity %in% c("Donor 3","Donor 4"), 2.8, 3.5),
             alpha = 0.9, stroke = 0.6) +
  scale_shape_manual(name = "Donor",
                     values = c("Donor 1" = 16,
                                "Donor 2" = 15,
                                "Donor 3" = 17,
                                "Donor 4" = 25,
                                "Donor 5" = 18)) +
  scale_color_manual(name = "Specificity",
                     values = c("IE1" = "red", "pp65" = "blue")) +
  scale_fill_manual(name = "Specificity",
                    values = c("IE1" = "red", "pp65" = "blue")) +
  guides(
    shape = guide_legend(override.aes = list(
      fill = c("black","black","black","black","black"),
      color = c("black","black","black","black","black"),
      size = c(3.5, 3.5, 2.8, 2.8, 3.5)
    ))
  ) +
  labs(
    x = NULL,
    y = "Gini-Simpson Index",
    title = "Clonal Diversity"
  ) +
  theme_classic(base_size = 14)

# Save plot
ggsave("exports/integrated_v5/updated_clone_frequency_plots/both_timepoints_gini_simpson.png",
       plot = gini_plot, width = 4, height = 5, dpi = 300, bg = "white")


#
# Gini-Simpson heterogeneity change between timepoints
#

# --- Step 1: Extract metadata and exclude "notCalled" cloneIds ---
meta_df <- mdata_CMV@meta.data %>%
  select(donor_identity, specificity, timepoint, cloneId) %>%
  filter(!grepl("notCalled", cloneId))

# --- Step 2: Count cells per clone for each donor × specificity × timepoint ---
clone_counts <- meta_df %>%
  group_by(donor_identity, specificity, timepoint, cloneId) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 3: Compute Gini-Simpson index for each donor × specificity × timepoint ---
gini_simpson <- clone_counts %>%
  group_by(donor_identity, specificity, timepoint) %>%
  summarise(
    total_cells = sum(cell_count),
    gini_simpson = 1 - sum((cell_count / sum(cell_count))^2),
    .groups = "drop"
  ) %>%
  mutate(
    specificity = factor(specificity, levels = c("pp65", "IE1")),
    timepoint = factor(timepoint, levels = c("Early", "Late"))
  )

# --- Common donor shape mapping ---
donor_shapes <- c(
  "Donor 1" = 16,  # solid circle
  "Donor 2" = 15,  # solid square
  "Donor 3" = 17,  # solid triangle (up)
  "Donor 4" = 25,  # inverted triangle (fillable)
  "Donor 5" = 18   # solid diamond
)

# --- Function to generate one specificity plot ---
make_gini_plot <- function(spec_name, point_color) {
  
  df <- gini_simpson %>%
    filter(specificity == spec_name)
  
  ggplot(df,
         aes(x = timepoint,
             y = gini_simpson,
             group = donor_identity,
             shape = donor_identity)) +
    
    geom_line(color = "black", linewidth = 0.6, alpha = 0.7) +
    
    geom_point(
      color = point_color,   # controls 16,15,17,18
      fill  = point_color,   # ensures 25 is filled
      size = 3.5,
      stroke = 0.5           # keeps inverted triangle looking solid
    ) +
    
    scale_shape_manual(values = donor_shapes, name = "Donor") +
    
    labs(
      x = NULL,
      y = "Gini-Simpson Index",
      title = paste0(spec_name, " Clonal Diversity")
    ) +
    
    theme_classic(base_size = 14)
}

# --- Generate plots ---
gini_pp65 <- make_gini_plot("pp65", "blue")
gini_IE1  <- make_gini_plot("IE1", "red")

# --- Save plots ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/pp65_gini_simpson_timepoint.png",
       plot = gini_pp65, width = 4, height = 5, dpi = 300, bg = "white")

ggsave("exports/integrated_v5/updated_clone_frequency_plots/IE1_gini_simpson_timepoint.png",
       plot = gini_IE1, width = 4, height = 5, dpi = 300, bg = "white")



# Pie charts for both timepoints, based on cloneId

# 1% cutoff

# --- Step 1: Extract metadata and exclude "notCalled" ---
clone_counts <- mdata_CMV@meta.data %>%
  select(donor_identity, specificity, cloneId) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(donor_identity, specificity, cloneId) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute proportions ---
clone_props <- clone_counts %>%
  group_by(donor_identity, specificity) %>%
  mutate(total_cells = sum(cell_count),
         prop = cell_count / total_cells) %>%
  ungroup()

# --- Step 3: Collapse <1% clones into "Other (<1%)" ---
clone_props_collapsed <- clone_props %>%
  mutate(Clone_temp = ifelse(prop < 0.01, "Other (<1%)", cloneId)) %>%
  group_by(donor_identity, specificity, Clone_temp) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(donor_identity, specificity) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

# --- Step 4: Rank real clones and append Other ---
clone_ranked <- clone_props_collapsed %>%
  group_by(donor_identity, specificity) %>%
  mutate(is_other = Clone_temp == "Other (<1%)") %>%
  group_modify(~ {
    df_real <- filter(.x, !is_other) %>%
      arrange(desc(prop)) %>%
      mutate(Clone = paste0("Clone ", row_number()))
    
    df_other <- filter(.x, is_other) %>%
      mutate(Clone = "Other (<1%)")
    
    bind_rows(df_real, df_other)
  }) %>%
  ungroup() %>%
  select(-Clone_temp, -is_other)

# --- Step 0b: Factor specificity for facet order (pp65 top, IE1 bottom) ---
clone_ranked$specificity <- factor(clone_ranked$specificity, levels = c("pp65", "IE1"))

# --- Step 5: Factor Clone for plotting ---
clone_levels_fixed <- c(paste0("Clone ", 1:20), "Other (<1%)")
clone_ranked$Clone <- factor(clone_ranked$Clone, levels = clone_levels_fixed)

# --- Step 6: Reverse factor for clockwise slices ---
clone_ranked <- clone_ranked %>%
  group_by(donor_identity, specificity) %>%
  mutate(Clone = factor(Clone, levels = rev(levels(Clone)))) %>%
  ungroup()

# --- Step 7: Define consistent colors ---
clone_colors_fixed <- c(
  "Clone 1"  = "#1b9e77",
  "Clone 2"  = "#d95f02",
  "Clone 3"  = "#7570b3",
  "Clone 4"  = "#e7298a",
  "Clone 5"  = "#66a61e",
  "Clone 6"  = "#e6ab02",
  "Clone 7"  = "#a6761d",
  "Clone 8"  = "#666666",
  "Clone 9"  = "#1f78b4",
  "Clone 10" = "#b2df8a",
  "Clone 11" = "#33a02c",
  "Clone 12" = "#fb9a99",
  "Clone 13" = "#e31a1c",
  "Clone 14" = "#fdbf6f",
  "Clone 15" = "#ff7f00",
  "Clone 16" = "#cab2d6",
  "Clone 17" = "#6a3d9a",
  "Clone 18" = "#ffff99",
  "Clone 19" = "#b15928",
  "Clone 20" = "#8dd3c7",
  "Other (<1%)" = "grey95"  # <-- updated color
)

# --- Step 8: Plot donut pie charts ---
pie_plot <- ggplot(clone_ranked,
                   aes(x = 2, y = prop, fill = Clone)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(specificity ~ donor_identity) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 9: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/clone_piecharts_1_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")


# 2% cutoff

# --- Step 1: Extract metadata and exclude "notCalled" ---
clone_counts <- mdata_CMV@meta.data %>%
  select(donor_identity, specificity, cloneId) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(donor_identity, specificity, cloneId) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute proportions ---
clone_props <- clone_counts %>%
  group_by(donor_identity, specificity) %>%
  mutate(total_cells = sum(cell_count),
         prop = cell_count / total_cells) %>%
  ungroup()

# --- Step 3: Collapse <2% clones into "Other (<2%)" ---
clone_props_collapsed <- clone_props %>%
  mutate(Clone_temp = ifelse(prop < 0.02, "Other (<2%)", cloneId)) %>%
  group_by(donor_identity, specificity, Clone_temp) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(donor_identity, specificity) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

# --- Step 4: Rank real clones and append Other ---
clone_ranked <- clone_props_collapsed %>%
  group_by(donor_identity, specificity) %>%
  mutate(is_other = Clone_temp == "Other (<2%)") %>%
  group_modify(~ {
    df_real <- filter(.x, !is_other) %>%
      arrange(desc(prop)) %>%
      mutate(Clone = paste0("Clone ", row_number()))
    
    df_other <- filter(.x, is_other) %>%
      mutate(Clone = "Other (<2%)")
    
    bind_rows(df_real, df_other)
  }) %>%
  ungroup() %>%
  select(-Clone_temp, -is_other)

# --- Step 0b: Factor specificity for facet order (pp65 top, IE1 bottom) ---
clone_ranked$specificity <- factor(clone_ranked$specificity, levels = c("pp65", "IE1"))

# --- Step 5: Factor Clone for plotting ---
clone_levels_fixed <- c(paste0("Clone ", 1:20), "Other (<2%)")
clone_ranked$Clone <- factor(clone_ranked$Clone, levels = clone_levels_fixed)

# --- Step 6: Reverse factor for clockwise slices ---
clone_ranked <- clone_ranked %>%
  group_by(donor_identity, specificity) %>%
  mutate(Clone = factor(Clone, levels = rev(levels(Clone)))) %>%
  ungroup()

# --- Step 7: Define consistent colors ---
clone_colors_fixed <- c(
  "Clone 1"  = "#1b9e77",
  "Clone 2"  = "#d95f02",
  "Clone 3"  = "#7570b3",
  "Clone 4"  = "#e7298a",
  "Clone 5"  = "#66a61e",
  "Clone 6"  = "#e6ab02",
  "Clone 7"  = "#a6761d",
  "Clone 8"  = "#666666",
  "Clone 9"  = "#1f78b4",
  "Clone 10" = "#b2df8a",
  "Clone 11" = "#33a02c",
  "Clone 12" = "#fb9a99",
  "Clone 13" = "#e31a1c",
  "Clone 14" = "#fdbf6f",
  "Clone 15" = "#ff7f00",
  "Clone 16" = "#cab2d6",
  "Clone 17" = "#6a3d9a",
  "Clone 18" = "#ffff99",
  "Clone 19" = "#b15928",
  "Clone 20" = "#8dd3c7",
  "Other (<2%)" = "grey95"  # <-- updated color
)

# --- Step 8: Plot donut pie charts ---
pie_plot <- ggplot(clone_ranked,
                   aes(x = 2, y = prop, fill = Clone)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(specificity ~ donor_identity) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 9: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/clone_piecharts_2_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")


# 0.5% cutoff

# --- Step 1: Extract metadata and exclude "notCalled" ---
clone_counts <- mdata_CMV@meta.data %>%
  select(donor_identity, specificity, cloneId) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(donor_identity, specificity, cloneId) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute proportions ---
clone_props <- clone_counts %>%
  group_by(donor_identity, specificity) %>%
  mutate(total_cells = sum(cell_count),
         prop = cell_count / total_cells) %>%
  ungroup()

# --- Step 3: Collapse <2% clones into "Other (<0.5%)" ---
clone_props_collapsed <- clone_props %>%
  mutate(Clone_temp = ifelse(prop < 0.005, "Other (<0.5%)", cloneId)) %>%
  group_by(donor_identity, specificity, Clone_temp) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(donor_identity, specificity) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

# --- Step 4: Rank real clones and append Other ---
clone_ranked <- clone_props_collapsed %>%
  group_by(donor_identity, specificity) %>%
  mutate(is_other = Clone_temp == "Other (<0.5%)") %>%
  group_modify(~ {
    df_real <- filter(.x, !is_other) %>%
      arrange(desc(prop)) %>%
      mutate(Clone = paste0("Clone ", row_number()))
    
    df_other <- filter(.x, is_other) %>%
      mutate(Clone = "Other (<0.5%)")
    
    bind_rows(df_real, df_other)
  }) %>%
  ungroup() %>%
  select(-Clone_temp, -is_other)

# --- Step 0b: Factor specificity for facet order (pp65 top, IE1 bottom) ---
clone_ranked$specificity <- factor(clone_ranked$specificity, levels = c("pp65", "IE1"))

# --- Step 5: Factor Clone for plotting ---
clone_levels_fixed <- c(paste0("Clone ", 1:20), "Other (<0.5%)")
clone_ranked$Clone <- factor(clone_ranked$Clone, levels = clone_levels_fixed)

# --- Step 6: Reverse factor for clockwise slices ---
clone_ranked <- clone_ranked %>%
  group_by(donor_identity, specificity) %>%
  mutate(Clone = factor(Clone, levels = rev(levels(Clone)))) %>%
  ungroup()

# --- Step 7: Define consistent colors ---
clone_colors_fixed <- c(
  "Clone 1"  = "#1b9e77",
  "Clone 2"  = "#d95f02",
  "Clone 3"  = "#7570b3",
  "Clone 4"  = "#e7298a",
  "Clone 5"  = "#66a61e",
  "Clone 6"  = "#e6ab02",
  "Clone 7"  = "#a6761d",
  "Clone 8"  = "#666666",
  "Clone 9"  = "#1f78b4",
  "Clone 10" = "#b2df8a",
  "Clone 11" = "#33a02c",
  "Clone 12" = "#fb9a99",
  "Clone 13" = "#e31a1c",
  "Clone 14" = "#fdbf6f",
  "Clone 15" = "#ff7f00",
  "Clone 16" = "#cab2d6",
  "Clone 17" = "#6a3d9a",
  "Clone 18" = "#ffff99",
  "Clone 19" = "#b15928",
  "Clone 20" = "#8dd3c7",
  "Other (<0.5%)" = "grey95"  # <-- updated color
)

# --- Step 8: Plot donut pie charts ---
pie_plot <- ggplot(clone_ranked,
                   aes(x = 2, y = prop, fill = Clone)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(specificity ~ donor_identity) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 9: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/clone_piecharts_0.5_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")


# 2% cutoff, 1% cutoff for minor clones

# --- Step 1: Extract metadata and exclude "notCalled" ---
clone_counts <- mdata_CMV@meta.data %>%
  select(donor_identity, specificity, cloneId) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(donor_identity, specificity, cloneId) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute proportions relative to clones >=1% ---
clone_props <- clone_counts %>%
  group_by(donor_identity, specificity) %>%
  mutate(total_cells_all = sum(cell_count),
         prop_total = cell_count / total_cells_all) %>%
  ungroup() %>%
  filter(prop_total >= 0.01)  # drop clones <1%

# --- Step 3: Assign Other (1–2%) for clones between 1–2%, real clones >2% ---
clone_props_collapsed <- clone_props %>%
  group_by(donor_identity, specificity) %>%
  mutate(Clone_temp = ifelse(prop_total < 0.02, "Other (1–2%)", cloneId)) %>%
  group_by(donor_identity, specificity, Clone_temp) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(donor_identity, specificity) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%  # proportions sum to clones >=1%
  ungroup()

# --- Step 4: Rank real clones (>2%) and append Other (1–2%) ---
clone_ranked <- clone_props_collapsed %>%
  group_by(donor_identity, specificity) %>%
  mutate(is_other = Clone_temp == "Other (1–2%)") %>%
  group_modify(~ {
    df_real <- filter(.x, !is_other) %>%
      arrange(desc(prop)) %>%
      mutate(Clone = paste0("Clone ", row_number()))
    
    df_other <- filter(.x, is_other) %>%
      mutate(Clone = "Other (1–2%)")
    
    bind_rows(df_real, df_other)
  }) %>%
  ungroup() %>%
  select(-Clone_temp, -is_other)

# --- Step 0b: Factor specificity for facet order (pp65 top, IE1 bottom) ---
clone_ranked$specificity <- factor(clone_ranked$specificity, levels = c("pp65", "IE1"))

# --- Step 5: Factor Clone for plotting ---
clone_levels_fixed <- c(paste0("Clone ", 1:20), "Other (1–2%)")
clone_ranked$Clone <- factor(clone_ranked$Clone, levels = clone_levels_fixed)

# --- Step 6: Reverse factor for clockwise slices ---
clone_ranked <- clone_ranked %>%
  group_by(donor_identity, specificity) %>%
  mutate(Clone = factor(Clone, levels = rev(levels(Clone)))) %>%
  ungroup()

# --- Step 7: Define consistent colors ---
clone_colors_fixed <- c(
  "Clone 1"  = "#1b9e77",
  "Clone 2"  = "#d95f02",
  "Clone 3"  = "#7570b3",
  "Clone 4"  = "#e7298a",
  "Clone 5"  = "#66a61e",
  "Clone 6"  = "#e6ab02",
  "Clone 7"  = "#a6761d",
  "Clone 8"  = "#666666",
  "Clone 9"  = "#1f78b4",
  "Clone 10" = "#b2df8a",
  "Clone 11" = "#33a02c",
  "Clone 12" = "#fb9a99",
  "Clone 13" = "#e31a1c",
  "Clone 14" = "#fdbf6f",
  "Clone 15" = "#ff7f00",
  "Clone 16" = "#cab2d6",
  "Clone 17" = "#6a3d9a",
  "Clone 18" = "#ffff99",
  "Clone 19" = "#b15928",
  "Clone 20" = "#8dd3c7",
  "Other (1–2%)" = "grey95"
)

# --- Step 8: Plot donut pie charts ---
pie_plot <- ggplot(clone_ranked,
                   aes(x = 2, y = prop, fill = Clone)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(specificity ~ donor_identity) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 9: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/clone_piecharts_minor_1_major_2_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")


# 2% cutoff, 0.5% cutoff for minor clones

# --- Step 1: Extract metadata and exclude "notCalled" ---
clone_counts <- mdata_CMV@meta.data %>%
  select(donor_identity, specificity, cloneId) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(donor_identity, specificity, cloneId) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute proportions relative to clones >=0.5% ---
clone_props <- clone_counts %>%
  group_by(donor_identity, specificity) %>%
  mutate(total_cells_all = sum(cell_count),
         prop_total = cell_count / total_cells_all) %>%
  ungroup() %>%
  filter(prop_total >= 0.005)  # drop clones <0.5%

# --- Step 3: Assign Other (0.5–2%) for clones between 0.5–2%, real clones >2% ---
clone_props_collapsed <- clone_props %>%
  group_by(donor_identity, specificity) %>%
  mutate(Clone_temp = ifelse(prop_total < 0.02, "Other (0.5–2%)", cloneId)) %>%
  group_by(donor_identity, specificity, Clone_temp) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(donor_identity, specificity) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%  # proportions sum to clones >=0.5%
  ungroup()

# --- Step 4: Rank real clones (>2%) and append Other (0.5–2%) ---
clone_ranked <- clone_props_collapsed %>%
  group_by(donor_identity, specificity) %>%
  mutate(is_other = Clone_temp == "Other (0.5–2%)") %>%
  group_modify(~ {
    df_real <- filter(.x, !is_other) %>%
      arrange(desc(prop)) %>%
      mutate(Clone = paste0("Clone ", row_number()))
    
    df_other <- filter(.x, is_other) %>%
      mutate(Clone = "Other (0.5–2%)")
    
    bind_rows(df_real, df_other)
  }) %>%
  ungroup() %>%
  select(-Clone_temp, -is_other)

# --- Step 0b: Factor specificity for facet order (pp65 top, IE1 bottom) ---
clone_ranked$specificity <- factor(clone_ranked$specificity, levels = c("pp65", "IE1"))

# --- Step 5: Factor Clone for plotting ---
clone_levels_fixed <- c(paste0("Clone ", 1:20), "Other (0.5–2%)")
clone_ranked$Clone <- factor(clone_ranked$Clone, levels = clone_levels_fixed)

# --- Step 6: Reverse factor for clockwise slices ---
clone_ranked <- clone_ranked %>%
  group_by(donor_identity, specificity) %>%
  mutate(Clone = factor(Clone, levels = rev(levels(Clone)))) %>%
  ungroup()

# --- Step 7: Define consistent colors ---
clone_colors_fixed <- c(
  "Clone 1"  = "#1b9e77",
  "Clone 2"  = "#d95f02",
  "Clone 3"  = "#7570b3",
  "Clone 4"  = "#e7298a",
  "Clone 5"  = "#66a61e",
  "Clone 6"  = "#e6ab02",
  "Clone 7"  = "#a6761d",
  "Clone 8"  = "#666666",
  "Clone 9"  = "#1f78b4",
  "Clone 10" = "#b2df8a",
  "Clone 11" = "#33a02c",
  "Clone 12" = "#fb9a99",
  "Clone 13" = "#e31a1c",
  "Clone 14" = "#fdbf6f",
  "Clone 15" = "#ff7f00",
  "Clone 16" = "#cab2d6",
  "Clone 17" = "#6a3d9a",
  "Clone 18" = "#ffff99",
  "Clone 19" = "#b15928",
  "Clone 20" = "#8dd3c7",
  "Other (0.5–2%)" = "grey95"
)

# --- Step 8: Plot donut pie charts ---
pie_plot <- ggplot(clone_ranked,
                   aes(x = 2, y = prop, fill = Clone)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(specificity ~ donor_identity) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 9: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/clone_piecharts_minor_0.5_major_2_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")




#
# Updated ET vs LT clone pie charts
#

# pp65

# --- Step 1: Extract metadata for pp65 and exclude "notCalled" ---
clone_counts <- mdata_CMV@meta.data %>%
  filter(specificity == "pp65") %>%
  select(donor_identity, timepoint, cloneId) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(donor_identity, timepoint, cloneId) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute global proportions per donor across both timepoints ---
donor_totals <- clone_counts %>%
  group_by(donor_identity, cloneId) %>%
  summarise(total_count = sum(cell_count), .groups = "drop") %>%
  group_by(donor_identity) %>%
  mutate(total_cells = sum(total_count),
         global_prop = total_count / total_cells) %>%
  ungroup()

# --- Step 3: Build clone mapping with ranking by global proportion ---
clone_mapping <- donor_totals %>%
  mutate(Clone_temp = ifelse(global_prop < 0.01, "Other (<1%)", cloneId)) %>%
  group_by(donor_identity, Clone_temp) %>%
  summarise(global_prop = sum(global_prop), .groups = "drop") %>%
  group_by(donor_identity) %>%
  mutate(is_other = Clone_temp == "Other (<1%)") %>%
  group_modify(~ {
    df_real <- filter(.x, !is_other) %>%
      arrange(desc(global_prop)) %>%
      mutate(Clone = paste0("Clone ", row_number()))
    df_other <- filter(.x, is_other) %>%
      mutate(Clone = "Other (<1%)")
    bind_rows(df_real, df_other)
  }) %>%
  ungroup() %>%
  select(donor_identity, Clone_temp, Clone)

# --- Step 4: Apply mapping back to timepoint-level data ---
clone_props <- clone_counts %>%
  left_join(donor_totals %>% select(donor_identity, cloneId, global_prop),
            by = c("donor_identity", "cloneId")) %>%
  mutate(Clone_temp = ifelse(global_prop < 0.01, "Other (<1%)", cloneId)) %>%
  left_join(clone_mapping, by = c("donor_identity", "Clone_temp")) %>%
  group_by(donor_identity, timepoint, Clone) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(donor_identity, timepoint) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

# --- Step 5: Factor Clone for plotting ---
clone_levels_fixed <- c(paste0("Clone ", 1:40), "Other (<1%)")  # allow up to 40 if needed
clone_props$Clone <- factor(clone_props$Clone, levels = clone_levels_fixed)

# --- Step 6: Reverse factor for clockwise slices ---
clone_props <- clone_props %>%
  group_by(donor_identity, timepoint) %>%
  mutate(Clone = factor(Clone, levels = rev(levels(Clone)))) %>%
  ungroup()

# --- Step 7: Factor donor/timepoint for facet order ---
clone_props$donor_identity <- factor(clone_props$donor_identity,
                                     levels = c("Donor 1","Donor 2","Donor 3","Donor 4","Donor 5"))
clone_props$timepoint <- factor(clone_props$timepoint, levels = c("Early","Late"))

# --- Step 8: Define consistent colors ---
clone_colors_fixed <- c(
  "Clone 1"  = "#1b9e77","Clone 2"  = "#d95f02","Clone 3"  = "#7570b3","Clone 4"  = "#e7298a",
  "Clone 5"  = "#66a61e","Clone 6"  = "#e6ab02","Clone 7"  = "#a6761d","Clone 8"  = "#666666",
  "Clone 9"  = "#1f78b4","Clone 10" = "#b2df8a","Clone 11" = "#33a02c","Clone 12" = "#fb9a99",
  "Clone 13" = "#e31a1c","Clone 14" = "#fdbf6f","Clone 15" = "#ff7f00","Clone 16" = "#cab2d6",
  "Clone 17" = "#6a3d9a","Clone 18" = "#ffff99","Clone 19" = "#b15928","Clone 20" = "#8dd3c7",
  "Other (<1%)" = "grey95"
)

# --- Step 9: Plot donut pie charts ---
pie_plot <- ggplot(clone_props, aes(x = 2, y = prop, fill = Clone)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(timepoint ~ donor_identity) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 10: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/pp65_clone_piecharts_1_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")


# IE1

# --- Step 1: Extract metadata for IE1 and exclude "notCalled" ---
clone_counts <- mdata_CMV@meta.data %>%
  filter(specificity == "IE1") %>%
  select(donor_identity, timepoint, cloneId) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(donor_identity, timepoint, cloneId) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute global proportions per donor across both timepoints ---
donor_totals <- clone_counts %>%
  group_by(donor_identity, cloneId) %>%
  summarise(total_count = sum(cell_count), .groups = "drop") %>%
  group_by(donor_identity) %>%
  mutate(total_cells = sum(total_count),
         global_prop = total_count / total_cells) %>%
  ungroup()

# --- Step 3: Build clone mapping with ranking by global proportion ---
clone_mapping <- donor_totals %>%
  mutate(Clone_temp = ifelse(global_prop < 0.01, "Other (<1%)", cloneId)) %>%
  group_by(donor_identity, Clone_temp) %>%
  summarise(global_prop = sum(global_prop), .groups = "drop") %>%
  group_by(donor_identity) %>%
  mutate(is_other = Clone_temp == "Other (<1%)") %>%
  group_modify(~ {
    df_real <- filter(.x, !is_other) %>%
      arrange(desc(global_prop)) %>%
      mutate(Clone = paste0("Clone ", row_number()))
    df_other <- filter(.x, is_other) %>%
      mutate(Clone = "Other (<1%)")
    bind_rows(df_real, df_other)
  }) %>%
  ungroup() %>%
  select(donor_identity, Clone_temp, Clone)

# --- Step 4: Apply mapping back to timepoint-level data ---
clone_props <- clone_counts %>%
  left_join(donor_totals %>% select(donor_identity, cloneId, global_prop),
            by = c("donor_identity", "cloneId")) %>%
  mutate(Clone_temp = ifelse(global_prop < 0.01, "Other (<1%)", cloneId)) %>%
  left_join(clone_mapping, by = c("donor_identity", "Clone_temp")) %>%
  group_by(donor_identity, timepoint, Clone) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(donor_identity, timepoint) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

# --- Step 5: Factor Clone for plotting ---
clone_levels_fixed <- c(paste0("Clone ", 1:40), "Other (<1%)")  # allow up to 40 if needed
clone_props$Clone <- factor(clone_props$Clone, levels = clone_levels_fixed)

# --- Step 6: Reverse factor for clockwise slices ---
clone_props <- clone_props %>%
  group_by(donor_identity, timepoint) %>%
  mutate(Clone = factor(Clone, levels = rev(levels(Clone)))) %>%
  ungroup()

# --- Step 7: Factor donor/timepoint for facet order ---
clone_props$donor_identity <- factor(clone_props$donor_identity,
                                     levels = c("Donor 1","Donor 2","Donor 3","Donor 4","Donor 5"))
clone_props$timepoint <- factor(clone_props$timepoint, levels = c("Early","Late"))

# --- Step 8: Define consistent colors ---
clone_colors_fixed <- c(
  "Clone 1"  = "#1b9e77","Clone 2"  = "#d95f02","Clone 3"  = "#7570b3","Clone 4"  = "#e7298a",
  "Clone 5"  = "#66a61e","Clone 6"  = "#e6ab02","Clone 7"  = "#a6761d","Clone 8"  = "#666666",
  "Clone 9"  = "#1f78b4","Clone 10" = "#b2df8a","Clone 11" = "#33a02c","Clone 12" = "#fb9a99",
  "Clone 13" = "#e31a1c","Clone 14" = "#fdbf6f","Clone 15" = "#ff7f00","Clone 16" = "#cab2d6",
  "Clone 17" = "#6a3d9a","Clone 18" = "#ffff99","Clone 19" = "#b15928","Clone 20" = "#8dd3c7",
  "Other (<1%)" = "grey95"
)

# --- Step 9: Plot donut pie charts ---
pie_plot <- ggplot(clone_props, aes(x = 2, y = prop, fill = Clone)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(timepoint ~ donor_identity) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 10: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/IE1_clone_piecharts_1_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")



#
# Updated pie charts for Donor 0
#

# 1% cutoff

# --- Step 1: Extract metadata and exclude "notCalled" ---
clone_counts <- D0_data_CMV@meta.data %>%
  select(timepoint, specificity, cloneId, Global_Clone_ID) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(timepoint, specificity, Global_Clone_ID) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute global proportions per specificity ---
global_totals <- clone_counts %>%
  group_by(specificity, Global_Clone_ID) %>%
  summarise(global_count = sum(cell_count), .groups = "drop") %>%
  group_by(specificity) %>%
  mutate(global_prop = global_count / sum(global_count)) %>%
  ungroup()

# --- Step 3: Identify clones below 1% globally ---
low_clones <- global_totals %>%
  filter(global_prop < 0.01) %>%
  pull(Global_Clone_ID) %>%
  as.character()

# --- Step 4: Collapse low clones and NAs into "Other (<1%)" ---
clone_counts <- clone_counts %>%
  mutate(Global_Clone_ID = as.character(Global_Clone_ID),
         Global_Clone_ID = ifelse(is.na(Global_Clone_ID) | Global_Clone_ID %in% low_clones,
                                  "Other (<1%)",
                                  Global_Clone_ID))

# --- Step 5: Compute proportions per timepoint for plotting ---
clone_props <- clone_counts %>%
  group_by(timepoint, specificity) %>%
  mutate(total_cells_timepoint = sum(cell_count),
         prop = cell_count / total_cells_timepoint) %>%
  ungroup()

# --- Step 5b: Collapse all "Other (<1%)" rows into one per specificity/timepoint ---
clone_props <- clone_props %>%
  group_by(timepoint, specificity, Global_Clone_ID) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(timepoint, specificity) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

# --- Step 6: Factor Global_Clone_ID for plotting (Other last, numeric order) ---
clone_levels <- clone_props %>%
  filter(Global_Clone_ID != "Other (<1%)") %>%
  distinct(Global_Clone_ID) %>%
  mutate(clone_num = as.numeric(sub(".*_(\\d+)$", "\\1", Global_Clone_ID))) %>%
  arrange(clone_num) %>%
  pull(Global_Clone_ID)

clone_levels <- c(clone_levels, "Other (<1%)")
clone_props$Global_Clone_ID <- factor(clone_props$Global_Clone_ID, levels = clone_levels)

# --- Step 7: Reverse factor for clockwise slices ---
clone_props <- clone_props %>%
  group_by(timepoint, specificity) %>%
  mutate(Global_Clone_ID = factor(Global_Clone_ID, levels = rev(levels(Global_Clone_ID)))) %>%
  ungroup()

# --- Step 8: Define consistent colors manually ---
clone_colors_fixed <- c(
  "D0_pp65_clone_1"  = "#1b9e77","D0_pp65_clone_2"  = "#d95f02","D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a","D0_pp65_clone_5"  = "#66a61e","D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d","D0_pp65_clone_8"  = "#666666","D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a","D0_pp65_clone_11" = "#33a02c","D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c","D0_pp65_clone_14" = "#fdbf6f","D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6","D0_pp65_clone_17" = "#6a3d9a","D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928","D0_pp65_clone_20" = "#8dd3c7","D0_pp65_clone_21" = "#a6cee3",
  "D0_pp65_clone_22" = "#1f78b4","D0_pp65_clone_23" = "#b2df8a","D0_pp65_clone_24" = "#fb9a99",
  "D0_pp65_clone_25" = "#fdbf6f","D0_pp65_clone_26" = "#ff7f00","D0_pp65_clone_27" = "#cab2d6",
  "D0_pp65_clone_28" = "#6a3d9a","D0_pp65_clone_29" = "#ffff99","D0_pp65_clone_30" = "#b15928",
  "D0_IE1_clone_1"  = "#1b9e77","D0_IE1_clone_2"  = "#d95f02","D0_IE1_clone_3"  = "#7570b3",
  "D0_IE1_clone_4"  = "#e7298a","D0_IE1_clone_5"  = "#66a61e","D0_IE1_clone_6"  = "#e6ab02",
  "D0_IE1_clone_7"  = "#a6761d","D0_IE1_clone_8"  = "#666666","D0_IE1_clone_9"  = "#1f78b4",
  "D0_IE1_clone_10" = "#b2df8a",
  "Other (<1%)" = "grey95"
)


# --- Step 0: Factor specificity for facet order (pp65 top, IE1 bottom) ---
clone_props$specificity <- factor(clone_props$specificity, levels = c("pp65", "IE1"))

# --- Step 9: Plot donut pie charts ---
pie_plot <- ggplot(clone_props,
                   aes(x = 2, y = prop, fill = Global_Clone_ID)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(specificity ~ timepoint) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 10: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/donor_0/clone_piecharts_1_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")



# 2% cutoff

# --- Step 1: Extract metadata and exclude "notCalled" ---
clone_counts <- D0_data_CMV@meta.data %>%
  select(timepoint, specificity, cloneId, Global_Clone_ID) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(timepoint, specificity, Global_Clone_ID) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute global proportions per specificity ---
global_totals <- clone_counts %>%
  group_by(specificity, Global_Clone_ID) %>%
  summarise(global_count = sum(cell_count), .groups = "drop") %>%
  group_by(specificity) %>%
  mutate(global_prop = global_count / sum(global_count)) %>%
  ungroup()

# --- Step 3: Identify clones below 2% globally ---
low_clones <- global_totals %>%
  filter(global_prop < 0.02) %>%
  pull(Global_Clone_ID) %>%
  as.character()

# --- Step 4: Collapse low clones and NAs into "Other (<2%)" ---
clone_counts <- clone_counts %>%
  mutate(Global_Clone_ID = as.character(Global_Clone_ID),
         Global_Clone_ID = ifelse(is.na(Global_Clone_ID) | Global_Clone_ID %in% low_clones,
                                  "Other (<2%)",
                                  Global_Clone_ID))

# --- Step 5: Compute proportions per timepoint for plotting ---
clone_props <- clone_counts %>%
  group_by(timepoint, specificity) %>%
  mutate(total_cells_timepoint = sum(cell_count),
         prop = cell_count / total_cells_timepoint) %>%
  ungroup()

# --- Step 5b: Collapse all "Other (<2%)" rows into one per specificity/timepoint ---
clone_props <- clone_props %>%
  group_by(timepoint, specificity, Global_Clone_ID) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(timepoint, specificity) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

# --- Step 6: Factor Global_Clone_ID for plotting (Other last, numeric order) ---
clone_levels <- clone_props %>%
  filter(Global_Clone_ID != "Other (<2%)") %>%
  distinct(Global_Clone_ID) %>%
  mutate(clone_num = as.numeric(sub(".*_(\\d+)$", "\\1", Global_Clone_ID))) %>%
  arrange(clone_num) %>%
  pull(Global_Clone_ID)

clone_levels <- c(clone_levels, "Other (<2%)")
clone_props$Global_Clone_ID <- factor(clone_props$Global_Clone_ID, levels = clone_levels)

# --- Step 7: Reverse factor for clockwise slices ---
clone_props <- clone_props %>%
  group_by(timepoint, specificity) %>%
  mutate(Global_Clone_ID = factor(Global_Clone_ID, levels = rev(levels(Global_Clone_ID)))) %>%
  ungroup()

# --- Step 8: Define consistent colors manually ---
clone_colors_fixed <- c(
  "D0_pp65_clone_1"  = "#1b9e77","D0_pp65_clone_2"  = "#d95f02","D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a","D0_pp65_clone_5"  = "#66a61e","D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d","D0_pp65_clone_8"  = "#666666","D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a","D0_pp65_clone_11" = "#33a02c","D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c","D0_pp65_clone_14" = "#fdbf6f","D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6","D0_pp65_clone_17" = "#6a3d9a","D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928","D0_pp65_clone_20" = "#8dd3c7","D0_pp65_clone_21" = "#a6cee3",
  "D0_pp65_clone_22" = "#1f78b4","D0_pp65_clone_23" = "#b2df8a","D0_pp65_clone_24" = "#fb9a99",
  "D0_pp65_clone_25" = "#fdbf6f","D0_pp65_clone_26" = "#ff7f00","D0_pp65_clone_27" = "#cab2d6",
  "D0_pp65_clone_28" = "#6a3d9a","D0_pp65_clone_29" = "#ffff99","D0_pp65_clone_30" = "#b15928",
  "D0_IE1_clone_1"  = "#1b9e77","D0_IE1_clone_2"  = "#d95f02","D0_IE1_clone_3"  = "#7570b3",
  "D0_IE1_clone_4"  = "#e7298a","D0_IE1_clone_5"  = "#66a61e","D0_IE1_clone_6"  = "#e6ab02",
  "D0_IE1_clone_7"  = "#a6761d","D0_IE1_clone_8"  = "#666666","D0_IE1_clone_9"  = "#1f78b4",
  "D0_IE1_clone_10" = "#b2df8a",
  "Other (<2%)" = "grey95"
)

# --- Step 0: Factor specificity for facet order (pp65 top, IE1 bottom) ---
clone_props$specificity <- factor(clone_props$specificity, levels = c("pp65", "IE1"))

# --- Step 9: Plot donut pie charts ---
pie_plot <- ggplot(clone_props,
                   aes(x = 2, y = prop, fill = Global_Clone_ID)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(specificity ~ timepoint) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 10: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/donor_0/clone_piecharts_2_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")



# 0.5% cutoff

# --- Step 1: Extract metadata and exclude "notCalled" ---
clone_counts <- D0_data_CMV@meta.data %>%
  select(timepoint, specificity, cloneId, Global_Clone_ID) %>%
  filter(!grepl("notCalled", cloneId)) %>%
  group_by(timepoint, specificity, Global_Clone_ID) %>%
  summarise(cell_count = n(), .groups = "drop")

# --- Step 2: Compute global proportions per specificity ---
global_totals <- clone_counts %>%
  group_by(specificity, Global_Clone_ID) %>%
  summarise(global_count = sum(cell_count), .groups = "drop") %>%
  group_by(specificity) %>%
  mutate(global_prop = global_count / sum(global_count)) %>%
  ungroup()

# --- Step 3: Identify clones below 0.5% globally ---
low_clones <- global_totals %>%
  filter(global_prop < 0.005) %>%
  pull(Global_Clone_ID) %>%
  as.character()

# --- Step 4: Collapse low clones and NAs into "Other (<0.5%)" ---
clone_counts <- clone_counts %>%
  mutate(Global_Clone_ID = as.character(Global_Clone_ID),
         Global_Clone_ID = ifelse(is.na(Global_Clone_ID) | Global_Clone_ID %in% low_clones,
                                  "Other (<0.5%)",
                                  Global_Clone_ID))

# --- Step 5: Compute proportions per timepoint for plotting ---
clone_props <- clone_counts %>%
  group_by(timepoint, specificity) %>%
  mutate(total_cells_timepoint = sum(cell_count),
         prop = cell_count / total_cells_timepoint) %>%
  ungroup()

# --- Step 5b: Collapse all "Other (<0.5%)" rows into one per specificity/timepoint ---
clone_props <- clone_props %>%
  group_by(timepoint, specificity, Global_Clone_ID) %>%
  summarise(cell_count = sum(cell_count), .groups = "drop") %>%
  group_by(timepoint, specificity) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

# --- Step 6: Factor Global_Clone_ID for plotting (Other last, numeric order) ---
clone_levels <- clone_props %>%
  filter(Global_Clone_ID != "Other (<0.5%)") %>%
  distinct(Global_Clone_ID) %>%
  mutate(clone_num = as.numeric(sub(".*_(\\d+)$", "\\1", Global_Clone_ID))) %>%
  arrange(clone_num) %>%
  pull(Global_Clone_ID)

clone_levels <- c(clone_levels, "Other (<0.5%)")
clone_props$Global_Clone_ID <- factor(clone_props$Global_Clone_ID, levels = clone_levels)

# --- Step 7: Reverse factor for clockwise slices ---
clone_props <- clone_props %>%
  group_by(timepoint, specificity) %>%
  mutate(Global_Clone_ID = factor(Global_Clone_ID, levels = rev(levels(Global_Clone_ID)))) %>%
  ungroup()

# --- Step 8: Define consistent colors manually ---
clone_colors_fixed <- c(
  # pp65 clones (1–40)
  "D0_pp65_clone_1"  = "#1b9e77","D0_pp65_clone_2"  = "#d95f02","D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a","D0_pp65_clone_5"  = "#66a61e","D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d","D0_pp65_clone_8"  = "#666666","D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a","D0_pp65_clone_11" = "#33a02c","D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c","D0_pp65_clone_14" = "#fdbf6f","D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6","D0_pp65_clone_17" = "#6a3d9a","D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928","D0_pp65_clone_20" = "#8dd3c7","D0_pp65_clone_21" = "#a6cee3",
  "D0_pp65_clone_22" = "#1f78b4","D0_pp65_clone_23" = "#b2df8a","D0_pp65_clone_24" = "#fb9a99",
  "D0_pp65_clone_25" = "#fdbf6f","D0_pp65_clone_26" = "#ff7f00","D0_pp65_clone_27" = "#cab2d6",
  "D0_pp65_clone_28" = "#6a3d9a","D0_pp65_clone_29" = "#ffff99","D0_pp65_clone_30" = "#b15928",
  "D0_pp65_clone_31" = "#8dd3c7","D0_pp65_clone_32" = "#a6cee3","D0_pp65_clone_33" = "#1f78b4",
  "D0_pp65_clone_34" = "#b2df8a","D0_pp65_clone_35" = "#fb9a99","D0_pp65_clone_36" = "#fdbf6f",
  "D0_pp65_clone_37" = "#ff7f00","D0_pp65_clone_38" = "#cab2d6","D0_pp65_clone_39" = "#6a3d9a",
  "D0_pp65_clone_40" = "#ffff99",
  
  # IE1 clones (1–10)
  "D0_IE1_clone_1"  = "#1b9e77","D0_IE1_clone_2"  = "#d95f02","D0_IE1_clone_3"  = "#7570b3",
  "D0_IE1_clone_4"  = "#e7298a","D0_IE1_clone_5"  = "#66a61e","D0_IE1_clone_6"  = "#e6ab02",
  "D0_IE1_clone_7"  = "#a6761d","D0_IE1_clone_8"  = "#666666","D0_IE1_clone_9"  = "#1f78b4",
  "D0_IE1_clone_10" = "#b2df8a",
  
  # Other
  "Other (<0.5%)" = "grey95"
)

# --- Step 0: Factor specificity for facet order (pp65 top, IE1 bottom) ---
clone_props$specificity <- factor(clone_props$specificity, levels = c("pp65", "IE1"))

# --- Step 9: Plot donut pie charts ---
pie_plot <- ggplot(clone_props,
                   aes(x = 2, y = prop, fill = Global_Clone_ID)) +
  geom_col(color = "black", width = 1) +
  coord_polar(theta = "y") +
  xlim(0.5, 2.5) +
  facet_grid(specificity ~ timepoint) +
  scale_fill_manual(values = clone_colors_fixed) +
  theme_void(base_size = 14) +
  theme(
    legend.position = "none",
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA)
  )

# --- Step 10: Save plot ---
ggsave("exports/integrated_v5/updated_clone_frequency_plots/donor_0/clone_piecharts_0.5_percent.png",
       plot = pie_plot, width = 12, height = 6, dpi = 300, bg = "white")



#
# Donor 0 pp65 clone frequency over time
#

# --- Define the clone colors ---
clone_colors <- c(
  "D0_pp65_clone_1"  = "#1b9e77","D0_pp65_clone_2"  = "#d95f02","D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a","D0_pp65_clone_5"  = "#66a61e","D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d","D0_pp65_clone_8"  = "#666666","D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a","D0_pp65_clone_11" = "#33a02c","D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c","D0_pp65_clone_14" = "#fdbf6f","D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6","D0_pp65_clone_17" = "#6a3d9a","D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928","D0_pp65_clone_20" = "#8dd3c7"
)

# --- Extract metadata ---
meta_df <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID) %>%
  dplyr::filter(specificity == "pp65", 
                Global_Clone_ID %in% names(clone_colors))

# --- Calculate percentages per timepoint ---
percent_df <- meta_df %>%
  group_by(timepoint, Global_Clone_ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(percent = count / sum(count) * 100)

# --- Ensure timepoints are in the desired order ---
percent_df$timepoint <- factor(percent_df$timepoint, 
                               levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))

# --- Create the plot ---
pp65_plot <- ggplot(percent_df, aes(x = timepoint, y = percent, group = Global_Clone_ID, color = Global_Clone_ID)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = clone_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "Clone Frequency (log scale)") +      # updated y-axis label
  ggtitle("pp65 Major Clones") +               # add plot title
  theme_classic() +
  theme(
    legend.position = "none",                      
    axis.text.x = element_text(size = 14),         
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # center title
  )

# --- Export the plot ---
ggsave("Donor_0/exports/clone_frequency_plots/pp65_clone_freq_over_time_log.png", 
       plot = pp65_plot, width = 8, height = 6, dpi = 300)



#
# Donor 0 IE1 clone frequency over time
#

# --- Define the clone colors ---
clone_colors <- c(
  "D0_IE1_clone_1"  = "#1b9e77","D0_IE1_clone_2"  = "#d95f02","D0_IE1_clone_3"  = "#7570b3",
  "D0_IE1_clone_4"  = "#e7298a","D0_IE1_clone_5"  = "#66a61e","D0_IE1_clone_6"  = "#e6ab02",
  "D0_IE1_clone_7"  = "#a6761d","D0_IE1_clone_8"  = "#666666","D0_IE1_clone_9"  = "#1f78b4",
  "D0_IE1_clone_10" = "#b2df8a"
)

# --- Extract metadata ---
meta_df <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID) %>%
  dplyr::filter(specificity == "IE1", 
                Global_Clone_ID %in% names(clone_colors))

# --- Calculate percentages per timepoint ---
percent_df <- meta_df %>%
  group_by(timepoint, Global_Clone_ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(percent = count / sum(count) * 100)

# --- Ensure timepoints are in the desired order ---
percent_df$timepoint <- factor(percent_df$timepoint, 
                               levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))

# --- Create the plot ---
IE1_plot <- ggplot(percent_df, aes(x = timepoint, y = percent, group = Global_Clone_ID, color = Global_Clone_ID)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = clone_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "Clone Frequency (log scale)") +      # updated y-axis label
  ggtitle("IE1 Major Clones") +               # add plot title
  theme_classic() +
  theme(
    legend.position = "none",                      
    axis.text.x = element_text(size = 14),         
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)  # center title
  )

# --- Export the plot ---
ggsave("Donor_0/exports/clone_frequency_plots/IE1_clone_freq_over_time_log.png", 
       plot = IE1_plot, width = 8, height = 6, dpi = 300)



#
# Donor 0 pp65 clone size change over time (log 2 fold change relative to Visit 1)
#

# --- Define the clone colors ---
clone_colors <- c(
  "D0_pp65_clone_1"  = "#1b9e77","D0_pp65_clone_2"  = "#d95f02","D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a","D0_pp65_clone_5"  = "#66a61e","D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d","D0_pp65_clone_8"  = "#666666","D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a","D0_pp65_clone_11" = "#33a02c","D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c","D0_pp65_clone_14" = "#fdbf6f","D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6","D0_pp65_clone_17" = "#6a3d9a","D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928","D0_pp65_clone_20" = "#8dd3c7"
)

# --- Extract metadata ---
meta_df <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID) %>%
  dplyr::filter(specificity == "pp65", 
                Global_Clone_ID %in% names(clone_colors))

# --- Calculate percentages per timepoint ---
percent_df <- meta_df %>%
  group_by(timepoint, Global_Clone_ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(percent = count / sum(count) * 100)

# --- Ensure timepoints are ordered ---
percent_df$timepoint <- factor(percent_df$timepoint, 
                               levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))

# --- Compute natural log of fold change relative to Visit 1 ---
# Add a small pseudocount (1e-6) to avoid log(0)
percent_df <- percent_df %>%
  group_by(Global_Clone_ID) %>%
  mutate(
    baseline = percent[timepoint == "Visit 1"],
    log2FC = log2((percent + 1e-6) / (baseline + 1e-6))
  ) %>%
  ungroup()

# --- Create the plot ---
pp65_log2FC_plot <- ggplot(percent_df, aes(x = timepoint, y = log2FC, group = Global_Clone_ID, color = Global_Clone_ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = clone_colors) +
  labs(x = NULL, y = "log2(Fold Change vs Visit 1)") +
  ggtitle("pp65 Major Clones — Relative Change (log2 Fold Change)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# --- Export the plot ---
ggsave("Donor_0/exports/clone_frequency_plots/pp65_clone_log2FC_over_time.png", 
       plot = pp65_log2FC_plot, width = 8, height = 6, dpi = 300)


#
# Donor 0 — pp65 and IE1 clone frequency over time
# (% of total CD8s, log scale)
#

# --- CD8 fraction lookup (your supplied numbers, as decimals) ---
cd8_frac <- tibble::tribble(
  ~specificity, ~timepoint,  ~spec_frac_of_cd8,
  "pp65",       "Visit 1",   0.27/100,
  "pp65",       "Visit 2",   0.25/100,
  "pp65",       "Visit 3",   0.20/100,
  "pp65",       "Visit 4",   0.10/100,
  "IE1",        "Visit 1",   3.58/100,
  "IE1",        "Visit 2",   3.00/100,
  "IE1",        "Visit 3",   2.94/100,
  "IE1",        "Visit 4",   0.77/100
)


#
# Donor 0 pp65 clone frequency over time (% of total CD8s)
#

# --- Define the clone colors ---
clone_colors <- c(
  "D0_pp65_clone_1"  = "#1b9e77","D0_pp65_clone_2"  = "#d95f02","D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a","D0_pp65_clone_5"  = "#66a61e","D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d","D0_pp65_clone_8"  = "#666666","D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a","D0_pp65_clone_11" = "#33a02c","D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c","D0_pp65_clone_14" = "#fdbf6f","D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6","D0_pp65_clone_17" = "#6a3d9a","D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928","D0_pp65_clone_20" = "#8dd3c7"
)

# --- Extract metadata ---
meta_df <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID) %>%
  dplyr::filter(specificity == "pp65",
                Global_Clone_ID %in% names(clone_colors))

# --- Calculate percent of total CD8s per timepoint ---
# pct_of_CD8 = (% of specificity for this clone) × (% of CD8s for this specificity at this visit)
percent_df <- meta_df %>%
  group_by(timepoint, Global_Clone_ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(pct_of_spec = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(specificity = "pp65") %>%
  left_join(cd8_frac, by = c("specificity", "timepoint")) %>%
  mutate(percent = pct_of_spec * spec_frac_of_cd8)   # pct_of_spec (%) × frac → %

# --- Order timepoints ---
percent_df$timepoint <- factor(percent_df$timepoint,
                               levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))

# --- Plot ---
pp65_plot <- ggplot(percent_df, aes(x = timepoint, y = percent,
                                    group = Global_Clone_ID, color = Global_Clone_ID)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = clone_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "Clone Frequency (% of total CD8s, log scale)") +
  ggtitle("pp65 Major Clones") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 14),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    plot.title   = element_text(size = 18, face = "bold", hjust = 0.5)
  )

ggsave("Donor_0/exports/clone_frequency_plots/pp65_clone_freq_over_time_pctCD8_log.png",
       plot = pp65_plot, width = 8, height = 6, dpi = 300)


#
# Donor 0 IE1 clone frequency over time (% of total CD8s)
#

clone_colors <- c(
  "D0_IE1_clone_1"  = "#1b9e77","D0_IE1_clone_2"  = "#d95f02","D0_IE1_clone_3"  = "#7570b3",
  "D0_IE1_clone_4"  = "#e7298a","D0_IE1_clone_5"  = "#66a61e","D0_IE1_clone_6"  = "#e6ab02",
  "D0_IE1_clone_7"  = "#a6761d","D0_IE1_clone_8"  = "#666666","D0_IE1_clone_9"  = "#1f78b4",
  "D0_IE1_clone_10" = "#b2df8a"
)

meta_df <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID) %>%
  dplyr::filter(specificity == "IE1",
                Global_Clone_ID %in% names(clone_colors))

percent_df <- meta_df %>%
  group_by(timepoint, Global_Clone_ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(pct_of_spec = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(specificity = "IE1") %>%
  left_join(cd8_frac, by = c("specificity", "timepoint")) %>%
  mutate(percent = pct_of_spec * spec_frac_of_cd8)

percent_df$timepoint <- factor(percent_df$timepoint,
                               levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))

IE1_plot <- ggplot(percent_df, aes(x = timepoint, y = percent,
                                   group = Global_Clone_ID, color = Global_Clone_ID)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = clone_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "Clone Frequency (% of total CD8s, log scale)") +
  ggtitle("IE1 Major Clones") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 14),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    plot.title   = element_text(size = 18, face = "bold", hjust = 0.5)
  )

ggsave("Donor_0/exports/clone_frequency_plots/IE1_clone_freq_over_time_pctCD8_log.png",
       plot = IE1_plot, width = 8, height = 6, dpi = 300)



#
# Donor 0 pp65 clone size change over time
# (log2 fold change of % of total CD8s, relative to Visit 1)
#

# --- Define the clone colors ---
clone_colors <- c(
  "D0_pp65_clone_1"  = "#1b9e77","D0_pp65_clone_2"  = "#d95f02","D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a","D0_pp65_clone_5"  = "#66a61e","D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d","D0_pp65_clone_8"  = "#666666","D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a","D0_pp65_clone_11" = "#33a02c","D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c","D0_pp65_clone_14" = "#fdbf6f","D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6","D0_pp65_clone_17" = "#6a3d9a","D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928","D0_pp65_clone_20" = "#8dd3c7"
)

# (Assumes cd8_frac tibble from previous block is already in your environment)

# --- Extract metadata ---
meta_df <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID) %>%
  dplyr::filter(specificity == "pp65",
                Global_Clone_ID %in% names(clone_colors))

# --- Calculate percent of total CD8s per timepoint ---
percent_df <- meta_df %>%
  group_by(timepoint, Global_Clone_ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(pct_of_spec = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(specificity = "pp65") %>%
  left_join(cd8_frac, by = c("specificity", "timepoint")) %>%
  mutate(percent = pct_of_spec * spec_frac_of_cd8)

# --- Order timepoints ---
percent_df$timepoint <- factor(percent_df$timepoint,
                               levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))

# --- Compute log2 fold change of pct-of-CD8 relative to Visit 1 ---
percent_df <- percent_df %>%
  group_by(Global_Clone_ID) %>%
  mutate(
    baseline = percent[timepoint == "Visit 1"],
    log2FC   = log2((percent + 1e-6) / (baseline + 1e-6))
  ) %>%
  ungroup()

# --- Plot ---
pp65_log2FC_plot <- ggplot(percent_df,
                           aes(x = timepoint, y = log2FC,
                               group = Global_Clone_ID, color = Global_Clone_ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = clone_colors) +
  labs(x = NULL, y = "log2(Fold Change of % of CD8s vs Visit 1)") +
  ggtitle("pp65 Major Clones — Relative Change in % of CD8s (log2 Fold Change)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 14),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5)
  )

ggsave("Donor_0/exports/clone_frequency_plots/pp65_clone_log2FC_pctCD8_over_time.png",
       plot = pp65_log2FC_plot, width = 8, height = 6, dpi = 300)


#
# Donor 0 IE1 clone size change over time
# (log2 fold change of % of total CD8s, relative to Visit 1)
#

clone_colors <- c(
  "D0_IE1_clone_1"  = "#1b9e77","D0_IE1_clone_2"  = "#d95f02","D0_IE1_clone_3"  = "#7570b3",
  "D0_IE1_clone_4"  = "#e7298a","D0_IE1_clone_5"  = "#66a61e","D0_IE1_clone_6"  = "#e6ab02",
  "D0_IE1_clone_7"  = "#a6761d","D0_IE1_clone_8"  = "#666666","D0_IE1_clone_9"  = "#1f78b4",
  "D0_IE1_clone_10" = "#b2df8a"
)

meta_df <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID) %>%
  dplyr::filter(specificity == "IE1",
                Global_Clone_ID %in% names(clone_colors))

percent_df <- meta_df %>%
  group_by(timepoint, Global_Clone_ID) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(timepoint) %>%
  mutate(pct_of_spec = count / sum(count) * 100) %>%
  ungroup() %>%
  mutate(specificity = "IE1") %>%
  left_join(cd8_frac, by = c("specificity", "timepoint")) %>%
  mutate(percent = pct_of_spec * spec_frac_of_cd8)

percent_df$timepoint <- factor(percent_df$timepoint,
                               levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))

percent_df <- percent_df %>%
  group_by(Global_Clone_ID) %>%
  mutate(
    baseline = percent[timepoint == "Visit 1"],
    log2FC   = log2((percent + 1e-6) / (baseline + 1e-6))
  ) %>%
  ungroup()

IE1_log2FC_plot <- ggplot(percent_df,
                          aes(x = timepoint, y = log2FC,
                              group = Global_Clone_ID, color = Global_Clone_ID)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = clone_colors) +
  labs(x = NULL, y = "log2(Fold Change of % of CD8s vs Visit 1)") +
  ggtitle("IE1 Major Clones — Relative Change in % of CD8s (log2 Fold Change)") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 14),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    plot.title   = element_text(size = 16, face = "bold", hjust = 0.5)
  )

ggsave("Donor_0/exports/clone_frequency_plots/IE1_clone_log2FC_pctCD8_over_time.png",
       plot = IE1_log2FC_plot, width = 8, height = 6, dpi = 300)




#
# Remaking clonotype highlight plots
#

# Define output directory
output_dir <- "exports/integrated_v5/Clonotype_plots/"

# Loop through donors
for (x in 1:5) {
  
  # --- IE1 clones (red highlights) ---
  for (i in 1:20) {
    clonotype_id <- paste0("D", x, "_IE1-clone-", i)
    
    # Find cells with this clonotype
    cells_to_highlight <- WhichCells(mdata_CMV, expression = Global_Clone_ID == clonotype_id)
    
    if (length(cells_to_highlight) > 0) {
      # Create the plot
      plot_title <- paste0("Donor ", x, " - IE1 clone ", i)
      p <- DimPlot(
        mdata_CMV, 
        reduction = "umap", 
        cells.highlight = cells_to_highlight, 
        sizes.highlight = 1.25,
        cols.highlight = "red"
      ) + 
        ggtitle(plot_title) +
        theme(legend.position = "none")
      
      # Export to PNG
      file_name <- paste0(output_dir, clonotype_id, ".png")
      ggsave(filename = file_name, plot = p, width = 6, height = 6, dpi = 300)
    }
  }
  
  # --- pp65 clones (blue highlights) ---
  for (i in 1:20) {
    clonotype_id <- paste0("D", x, "_pp65-clone-", i)
    
    # Find cells with this clonotype
    cells_to_highlight <- WhichCells(mdata_CMV, expression = Global_Clone_ID == clonotype_id)
    
    if (length(cells_to_highlight) > 0) {
      # Create the plot
      plot_title <- paste0("Donor ", x, " - pp65 clone ", i)
      p <- DimPlot(
        mdata_CMV, 
        reduction = "umap", 
        cells.highlight = cells_to_highlight, 
        sizes.highlight = 1.25,
        cols.highlight = "blue"
      ) + 
        ggtitle(plot_title) +
        theme(legend.position = "none")
      
      # Export to PNG
      file_name <- paste0(output_dir, clonotype_id, ".png")
      ggsave(filename = file_name, plot = p, width = 6, height = 6, dpi = 300)
    }
  }
}




#
#
# Gene set analysis - correlating 5 phenotypes back to Donor 0
#
#

Idents(mdata_CMV) <- "metacluster"

metacluster_markers <- FindAllMarkers(
  mdata_CMV,
  only.pos = TRUE,
  logfc.threshold = 0.25,
  min.pct = 0.25,
  test.use = "wilcox"
)

gene_sets <- metacluster_markers %>%
  filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  summarise(genes = list(gene)) %>%
  deframe()

gene_sets_filtered <- lapply(gene_sets, function(genes) {
  intersect(genes, rownames(D0_data_CMV))
})

D0_data_CMV <- AddModuleScore(
  D0_data_CMV,
  features = gene_sets_filtered,
  name = names(gene_sets_filtered)
)

# Get gene set names in correct order
set_names <- names(gene_sets_filtered)

# Construct the exact column names Seurat created
old_names <- paste0(set_names, seq_along(set_names))

# Construct desired new names
new_names <- paste0(set_names, "_score")

# Rename metadata columns
colnames(D0_data_CMV@meta.data)[
  match(old_names, colnames(D0_data_CMV@meta.data))
] <- new_names


DimPlot(D0_data_CMV, reduction = "umap", group.by = "seurat_clusters")
FeaturePlot(D0_data_CMV, features = "ISG+_score")



#
#
# Jensen-Shannon Divergence of clone phenotypes for Donor 0
#
#

# ----------------------------
# Jensen–Shannon divergence
# ----------------------------

js_divergence <- function(p, q) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  
  epsilon <- 1e-10
  p <- p + epsilon
  q <- q + epsilon
  
  p <- p / sum(p)
  q <- q / sum(q)
  
  m <- 0.5 * (p + q)
  0.5 * sum(p * log2(p / m)) + 0.5 * sum(q * log2(q / m))
}

# ----------------------------
# Helper function to compute JSD timecourse
# ----------------------------

compute_jsd_timecourse <- function(meta_df) {
  
  meta_df$timepoint <- factor(meta_df$timepoint,
                              levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))
  
  # Cluster proportions
  cluster_dist <- meta_df %>%
    group_by(Global_Clone_ID, timepoint, seurat_clusters) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(Global_Clone_ID, timepoint) %>%
    mutate(prop = n_cells / sum(n_cells)) %>%
    ungroup()
  
  all_clusters <- sort(unique(meta_df$seurat_clusters))
  
  cluster_dist <- cluster_dist %>%
    complete(
      Global_Clone_ID,
      timepoint,
      seurat_clusters = all_clusters,
      fill = list(n_cells = 0, prop = 0)
    )
  
  # Safe JSD computation via wide format
  jsd_df <- cluster_dist %>%
    group_by(Global_Clone_ID) %>%
    group_modify(~{
      
      df_clone <- .x
      
      wide_mat <- df_clone %>%
        select(timepoint, seurat_clusters, prop) %>%
        pivot_wider(
          names_from = seurat_clusters,
          values_from = prop,
          values_fill = 0
        )
      
      baseline_vec <- wide_mat %>%
        filter(timepoint == "Visit 1") %>%
        select(-timepoint) %>%
        as.numeric()
      
      wide_mat %>%
        rowwise() %>%
        mutate(
          JSD = js_divergence(
            c_across(-timepoint),
            baseline_vec
          )
        ) %>%
        ungroup() %>%
        select(timepoint, JSD)
      
    }) %>%
    ungroup() %>%
    mutate(one_minus_jsd = 1 - JSD)
  
  return(jsd_df)
}

# ----------------------------
# Extract metadata separately
# ----------------------------

meta_pp65 <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID, seurat_clusters) %>%
  dplyr::filter(specificity == "pp65",
                !is.na(Global_Clone_ID))

meta_IE1 <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID, seurat_clusters) %>%
  dplyr::filter(specificity == "IE1",
                !is.na(Global_Clone_ID))

# ----------------------------
# Compute JSD separately
# ----------------------------

jsd_pp65 <- compute_jsd_timecourse(meta_pp65)
jsd_IE1  <- compute_jsd_timecourse(meta_IE1)

# ----------------------------
# Plotting function
# ----------------------------

plot_jsd_timecourse <- function(df, spec_label, output_file) {
  
  p <- ggplot(df,
              aes(x = timepoint,
                  y = one_minus_jsd,
                  group = Global_Clone_ID,
                  color = Global_Clone_ID)) +
    geom_line(size = 1) +
    geom_point(size = 2.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    labs(x = NULL,
         y = "1 - Jensen-Shannon divergence (vs Visit 1)",
         color = "Clone ID",
         title = paste0(spec_label, " Clones — Phenotypic Stability Over Time")) +
    theme_classic() +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 9),
      legend.text = element_text(size = 7),
      legend.key.height = unit(0.4, "cm"),
      legend.key.width = unit(0.4, "cm"),
      axis.text.x = element_text(size = 14),
      axis.title.y = element_text(size = 16, margin = margin(r = 15)),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
    )
  
  ggsave(output_file,
         plot = p,
         width = 9, height = 6, dpi = 300)
  
  return(p)
}

# ----------------------------
# Generate plots
# ----------------------------

clone_colors_pp65 <- c(
  "D0_pp65_clone_1"  = "#1b9e77","D0_pp65_clone_2"  = "#d95f02","D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a","D0_pp65_clone_5"  = "#66a61e","D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d","D0_pp65_clone_8"  = "#666666","D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a","D0_pp65_clone_11" = "#33a02c","D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c","D0_pp65_clone_14" = "#fdbf6f","D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6","D0_pp65_clone_17" = "#6a3d9a","D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928","D0_pp65_clone_20" = "#8dd3c7"
)

jsd_pp65_top20 <- jsd_pp65 %>%
  filter(Global_Clone_ID %in% names(clone_colors_pp65))

pp65_plot <- ggplot(jsd_pp65_top20,
                    aes(x = timepoint,
                        y = one_minus_jsd,
                        group = Global_Clone_ID,
                        color = Global_Clone_ID)) +
  geom_line(size = 1) +
  geom_point(size = 2.5) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = clone_colors_pp65) +
  labs(x = NULL,
       y = "1 - Jensen-Shannon divergence (vs Visit 1)",
       color = "Clone ID",
       title = "pp65 Major Clones — Phenotypic Stability Over Time") +
  theme_classic() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9),
    legend.text = element_text(size = 7),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(0.4, "cm"),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5)
  )

ggsave("Donor_0/exports/pp65_clone_phenotype_stability_over_time.png",
       plot = pp65_plot,
       width = 9, height = 6, dpi = 300)


IE1_plot <- plot_jsd_timecourse(
  jsd_IE1,
  "IE1",
  "Donor_0/exports/IE1_clone_phenotype_stability_over_time.png"
)


#
# 1-JSD over time plots for Donor 0, both specificities combined
#

# ----------------------------
# Jensen–Shannon divergence
# ----------------------------

js_divergence <- function(p, q) {
  p <- as.numeric(p)
  q <- as.numeric(q)
  
  epsilon <- 1e-10
  p <- p + epsilon
  q <- q + epsilon
  
  p <- p / sum(p)
  q <- q / sum(q)
  
  m <- 0.5 * (p + q)
  0.5 * sum(p * log2(p / m)) + 0.5 * sum(q * log2(q / m))
}

# ----------------------------
# Helper function to compute JSD timecourse
# ----------------------------

compute_jsd_timecourse <- function(meta_df, min_cells = 10) {
  
  meta_df$timepoint <- factor(meta_df$timepoint,
                              levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4"))
  
  # Retain specificity label (one value per clone)
  specificity_map <- meta_df %>%
    distinct(Global_Clone_ID, specificity)
  
  # Count total cells per clone per timepoint BEFORE completing the grid
  cell_counts <- meta_df %>%
    group_by(Global_Clone_ID, timepoint) %>%
    summarise(total_cells = n(), .groups = "drop")
  
  # Drop clones whose Visit 1 is below the cell threshold — baseline is unreliable
  valid_clones <- cell_counts %>%
    filter(timepoint == "Visit 1", total_cells >= min_cells) %>%
    pull(Global_Clone_ID)
  
  meta_df <- meta_df %>% filter(Global_Clone_ID %in% valid_clones)
  cell_counts <- cell_counts %>% filter(Global_Clone_ID %in% valid_clones)
  
  # Cluster proportions
  all_clusters <- sort(unique(meta_df$seurat_clusters))
  
  cluster_dist <- meta_df %>%
    group_by(Global_Clone_ID, timepoint, seurat_clusters) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    group_by(Global_Clone_ID, timepoint) %>%
    mutate(prop = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    complete(
      Global_Clone_ID,
      timepoint,
      seurat_clusters = all_clusters,
      fill = list(n_cells = 0, prop = 0)
    )
  
  # Compute JSD vs Visit 1 for each clone
  jsd_df <- cluster_dist %>%
    group_by(Global_Clone_ID) %>%
    group_modify(~{
      
      df_clone <- .x
      
      wide_mat <- df_clone %>%
        select(timepoint, seurat_clusters, prop) %>%
        pivot_wider(
          names_from  = seurat_clusters,
          values_from = prop,
          values_fill = 0
        )
      
      baseline_vec <- wide_mat %>%
        filter(timepoint == "Visit 1") %>%
        select(-timepoint) %>%
        as.numeric()
      
      wide_mat %>%
        rowwise() %>%
        mutate(
          JSD = js_divergence(c_across(-timepoint), baseline_vec)
        ) %>%
        ungroup() %>%
        select(timepoint, JSD)
      
    }) %>%
    ungroup() %>%
    mutate(one_minus_jsd = 1 - JSD)
  
  # Join cell counts and mask timepoints below the threshold
  jsd_df <- jsd_df %>%
    left_join(cell_counts, by = c("Global_Clone_ID", "timepoint")) %>%
    mutate(one_minus_jsd = ifelse(total_cells >= min_cells, one_minus_jsd, NA)) %>%
    left_join(specificity_map, by = "Global_Clone_ID")
  
  return(jsd_df)
}

# ----------------------------
# Extract metadata separately
# ----------------------------

meta_pp65 <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID, seurat_clusters) %>%
  dplyr::filter(specificity == "pp65", !is.na(Global_Clone_ID))

meta_IE1 <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID, seurat_clusters) %>%
  dplyr::filter(specificity == "IE1",  !is.na(Global_Clone_ID))

# ----------------------------
# Compute JSD separately, then combine
# ----------------------------

jsd_pp65 <- compute_jsd_timecourse(meta_pp65, min_cells = 10)
jsd_IE1  <- compute_jsd_timecourse(meta_IE1,  min_cells = 10)

jsd_combined <- bind_rows(jsd_pp65, jsd_IE1)

# ----------------------------
# Combined plot — pp65 blue, IE1 red
# ----------------------------

combined_plot <- ggplot(jsd_combined,
                        aes(x     = timepoint,
                            y     = one_minus_jsd,
                            group = Global_Clone_ID,
                            color = specificity)) +
  geom_line(size = 0.8, na.rm = TRUE) +
  geom_point(size = 2.5, na.rm = TRUE) +           # NA points are simply omitted
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("pp65" = "#2166ac", "IE1" = "#d6604d"),
    labels = c("pp65" = "pp65", "IE1" = "IE1")
  ) +
  labs(
    x     = NULL,
    y     = "1 - Jensen-Shannon divergence (vs Visit 1)",
    color = "Specificity",
    title = "pp65 & IE1 Clones — Phenotypic Stability Over Time\n(clones with < 10 cells at a timepoint excluded)"
  ) +
  theme_classic() +
  theme(
    legend.position  = "right",
    legend.title     = element_text(size = 11),
    legend.text      = element_text(size = 10),
    axis.text.x      = element_text(size = 14),
    axis.title.y     = element_text(size = 14, margin = margin(r = 15)),
    plot.title       = element_text(size = 16, face = "bold", hjust = 0.5)
  )

ggsave("Donor_0/exports/combined_clone_phenotype_stability_over_time.png",
       plot   = combined_plot,
       width  = 9, height = 6, dpi = 300)


# Paneled plot

# ----------------------------
# Clone color definitions
# ----------------------------

clone_colors_pp65 <- c(
  "D0_pp65_clone_1"  = "#1b9e77", "D0_pp65_clone_2"  = "#d95f02", "D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a", "D0_pp65_clone_5"  = "#66a61e", "D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d", "D0_pp65_clone_8"  = "#666666", "D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a", "D0_pp65_clone_11" = "#33a02c", "D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c", "D0_pp65_clone_14" = "#fdbf6f", "D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6", "D0_pp65_clone_17" = "#6a3d9a", "D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928", "D0_pp65_clone_20" = "#8dd3c7"
)

clone_colors_IE1 <- c(
  "D0_IE1_clone_1"  = "#1b9e77","D0_IE1_clone_2"  = "#d95f02","D0_IE1_clone_3"  = "#7570b3",
  "D0_IE1_clone_4"  = "#e7298a","D0_IE1_clone_5"  = "#66a61e","D0_IE1_clone_6"  = "#e6ab02",
  "D0_IE1_clone_7"  = "#a6761d","D0_IE1_clone_8"  = "#666666","D0_IE1_clone_9"  = "#1f78b4",
  "D0_IE1_clone_10" = "#b2df8a"
)

clone_colors_all <- c(clone_colors_pp65, clone_colors_IE1)

# ----------------------------
# Filter to only clones that have a defined color
# ----------------------------

jsd_combined_filtered <- jsd_combined %>%
  filter(Global_Clone_ID %in% names(clone_colors_all))

# ----------------------------
# Paneled plot with per-clone colors and shared y-axis
# ----------------------------

combined_plot_paneled <- ggplot(jsd_combined_filtered,
                                aes(x     = timepoint,
                                    y     = one_minus_jsd,
                                    group = Global_Clone_ID,
                                    color = Global_Clone_ID)) +
  geom_line(size = 0.8, na.rm = TRUE, show.legend = FALSE) +
  geom_point(size = 2.5, na.rm = TRUE, show.legend = FALSE) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
  scale_color_manual(values = clone_colors_all, guide = "none") +
  facet_grid(
    rows = vars(factor(specificity, levels = c("IE1", "pp65")))
  ) +
  labs(
    x     = NULL,
    y     = "1 - Jensen-Shannon divergence (vs Visit 1)",
    title = "All Clones — Phenotypic Stability Over Time"
  ) +
  theme_classic() +
  theme(
    axis.text.x      = element_text(size = 14),
    axis.title.y     = element_text(size = 14, margin = margin(r = 15)),
    plot.title       = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.background = element_blank(),
    strip.text       = element_text(size = 13, face = "bold"),
    panel.spacing    = unit(0, "cm"),
    panel.border     = element_rect(color = "grey40", fill = NA, linewidth = 0.5)
  )

ggsave("Donor_0/exports/combined_clone_phenotype_stability_paneled.png",
       plot   = combined_plot_paneled,
       width  = 9, height = 8, dpi = 300)


#
#
# Stacked cluster bar graphs for Donor 0 over timepoints
#
#

############################################################
# Define metacluster color scheme
############################################################

cluster_colors <- c(
  "0" = "#C49A00",
  "1" = "#53B400",
  "2" = "#8B0000",
  "3" = "#1976D2",
  "4" = "#F8766D",
  "5" = "#A58AFF",
  "6" = "#FFDB58",
  "7" = "#00C094"
)

meta <- D0_data_CMV@meta.data

# Ensure clusters are ordered consistently
meta$seurat_clusters <- factor(
  as.character(meta$seurat_clusters),
  levels = names(cluster_colors)
)

############################################################
# PLOT 1: All cells — stacked % per timepoint
############################################################

cluster_counts <- meta %>%
  group_by(timepoint, seurat_clusters) %>%
  summarise(n_cells = n(), .groups = "drop")

timepoint_totals <- meta %>%
  group_by(timepoint) %>%
  summarise(total_cells = n(), .groups = "drop")

plot_df_all <- cluster_counts %>%
  left_join(timepoint_totals, by = "timepoint") %>%
  mutate(percent = 100 * n_cells / total_cells)

p_all <- ggplot(plot_df_all, aes(x = timepoint, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cluster_colors, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ylab("Percent of cells") +
  xlab("Timepoint") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Export
ggsave(
  filename = "Donor_0/exports/D0_all_cells_cluster_composition_over_time.png",
  plot = p_all,
  width = 6,
  height = 5,
  dpi = 300,
  bg = "white"
)


############################################################
# PLOT 2: Faceted specificity cluster composition
############################################################

meta_sub <- meta %>%
  filter(specificity %in% c("pp65", "IE1"))

# Count cells per specificity/timepoint/cluster
cluster_counts <- meta_sub %>%
  group_by(specificity, timepoint, seurat_clusters) %>%
  summarise(n_cells = n(), .groups = "drop")

# Totals per specificity/timepoint
totals <- meta_sub %>%
  group_by(specificity, timepoint) %>%
  summarise(total_cells = n(), .groups = "drop")

# Compute percentages
plot_df <- cluster_counts %>%
  left_join(totals, by = c("specificity", "timepoint")) %>%
  mutate(percent = 100 * n_cells / total_cells)

# Generate faceted plot
p_spec_facet <- ggplot(plot_df,
                       aes(x = timepoint,
                           y = percent,
                           fill = seurat_clusters)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cluster_colors, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  facet_wrap(~ specificity, nrow = 1) +
  ylab("Percent of cells") +
  xlab("Visit") +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold")
  )

# Export
ggsave(
  filename = "Donor_0/exports/cluster_composition_over_time_by_specificity.png",
  plot = p_spec_facet,
  width = 8,
  height = 5,
  dpi = 300,
  bg = "white"
)



#
#
# Gene set scores for standard gene sets
#
#

# ── 1. Define gene sets ──────────────────────────────────────────────────────

# NOTE: "LEFT" in the image has been interpreted as "LEF1" - correct if needed
Tcm_genes <- list(c("CCR7", "SELL", "IL7R", "TCF7", "EOMES", "CD27", "CD28",
                    "BCL2", "BACH2", "ZEB1", "LEF1", "KLF2", "S1PR1", "BCL6", "ID3"))

Tem_genes <- list(c("GZMB", "PRF1", "GZMH", "IFNG", "KLRG1", "ZEB2", "GNLY",
                    "NKG7", "TBX21", "CX3CR1", "CD244", "CXCR3", "HLA-DRA",
                    "FCGR3A", "IL2RB"))

# ── 2. Add module scores ─────────────────────────────────────────────────────

D0_data_CMV <- AddModuleScore(D0_data_CMV,
                              features = Tcm_genes,
                              name     = "Standard_Tcm_score")
# AddModuleScore appends "1" to the name; rename to clean version
D0_data_CMV$Standard_Tcm_score <- D0_data_CMV$Standard_Tcm_score1
D0_data_CMV$Standard_Tcm_score1 <- NULL

D0_data_CMV <- AddModuleScore(D0_data_CMV,
                              features = Tem_genes,
                              name     = "Standard_Tem_score")
D0_data_CMV$Standard_Tem_score <- D0_data_CMV$Standard_Tem_score1
D0_data_CMV$Standard_Tem_score1 <- NULL

# ── 3. Helper: FeaturePlot theme matching your example ───────────────────────

make_feature_plot <- function(obj, feature) {
  FeaturePlot(obj,
              features  = feature,
              reduction = "umap",
              cols      = c("blue", "red"))
}

# ── 4. UMAP plots – all cells ────────────────────────────────────────────────

ggsave("Donor_0/exports/D0_CMV_UMAP_ISG+_score.png",
       plot   = make_feature_plot(D0_data_CMV, "ISG+_score"),
       width  = 6, height = 5, dpi = 300)

ggsave("Donor_0/exports/D0_CMV_UMAP_Standard_Tem_score.png",
       plot   = make_feature_plot(D0_data_CMV, "Standard_Tem_score"),
       width  = 6, height = 5, dpi = 300)

# ── 5. UMAP plots – Visit 4 cells only ──────────────────────────────────────

D0_CMV_V4 <- subset(D0_data_CMV, subset = timepoint == "Visit 4")

ggsave("Donor_0/exports/D0_CMV_UMAP_V4_ISG+_score.png",
       plot   = make_feature_plot(D0_CMV_V4, "ISG+_score"),
       width  = 6, height = 5, dpi = 300)

ggsave("Donor_0/exports/D0_CMV_UMAP_V4_Standard_Tcm_score.png",
       plot   = make_feature_plot(D0_CMV_V4, "Standard_Tcm_score"),
       width  = 6, height = 5, dpi = 300)

# ── 6. Heatmap helper ────────────────────────────────────────────────────────

score_cols <- c("Tem_score", "Tcm_score", "T-CMC1_score", "NK-like_score",
                "ISG+_score", "Standard_Tcm_score", "Standard_Tem_score")

make_specificity_heatmap <- function(obj, title, filename) {
  
  meta <- obj@meta.data
  
  # Mean score per specificity group (pp65 vs IE1)
  pp65_means <- meta %>%
    filter(specificity == "pp65") %>%
    summarise(across(all_of(score_cols), mean, na.rm = TRUE)) %>%
    unlist()
  
  IE1_means <- meta %>%
    filter(specificity == "IE1") %>%
    summarise(across(all_of(score_cols), mean, na.rm = TRUE)) %>%
    unlist()
  
  heatmap_mat <- cbind(pp65 = pp65_means, IE1 = IE1_means)
  rownames(heatmap_mat) <- score_cols
  
  # Scale rows so colour reflects relative enrichment across the two groups
  # heatmap_mat_scaled <- t(scale(t(heatmap_mat)))
  
  pheatmap(heatmap_mat,
           main          = title,
           cluster_rows  = FALSE,
           cluster_cols  = FALSE,
           color         = colorRampPalette(c("blue", "red"))(50),
           border_color  = "black",
           fontsize_row  = 11,
           fontsize_col  = 12,
           filename      = filename,
           width         = 5,
           height        = 6)
}

# ── 7. Heatmaps ──────────────────────────────────────────────────────────────

# All cells
make_specificity_heatmap(D0_data_CMV,
                         title    = "Gene Set Scores: pp65 vs IE1",
                         filename = "Donor_0/exports/D0_CMV_Heatmap_pp65_IE1.png")

# Visit 4 only
make_specificity_heatmap(D0_CMV_V4,
                         title    = "Gene Set Scores: pp65 vs IE1 (Visit 4)",
                         filename = "Donor_0/exports/D0_CMV_Heatmap_pp65_IE1_V4.png")


#
#
# Version for heatmap of gene sets by cluster
#
#


# ── 6. Heatmap helper ────────────────────────────────────────────────────────
score_cols <- c("Tem_score", "Tcm_score", "T-CMC1_score", "NK-like_score",
                "ISG+_score", "Standard_Tcm_score", "Standard_Tem_score")

make_cluster_heatmap <- function(obj, title, filename) {
  
  meta <- obj@meta.data
  
  # Get all unique clusters (sorted numerically)
  clusters <- sort(unique(meta$seurat_clusters))
  
  # Compute mean score for each cluster
  cluster_means <- lapply(clusters, function(cl) {
    meta %>%
      filter(seurat_clusters == cl) %>%
      summarise(across(all_of(score_cols), mean, na.rm = TRUE)) %>%
      unlist()
  })
  
  # Build matrix: rows = scores, columns = clusters
  heatmap_mat <- do.call(cbind, cluster_means)
  colnames(heatmap_mat) <- paste0("Cluster ", clusters)
  rownames(heatmap_mat) <- score_cols
  
  pheatmap(heatmap_mat,
           main         = title,
           cluster_rows = FALSE,
           cluster_cols = FALSE,
           color        = colorRampPalette(c("blue", "red"))(50),
           border_color = "black",
           fontsize_row = 11,
           fontsize_col = 11,
           filename     = filename,
           width        = 8,   # wider to accommodate more columns
           height       = 6)
}

# ── 7. Heatmaps ──────────────────────────────────────────────────────────────

# All cells
make_cluster_heatmap(D0_data_CMV,
                     title    = "Gene Set Scores by Cluster",
                     filename = "Donor_0/exports/D0_CMV_Heatmap_clusters.png")

# Visit 4 only
make_cluster_heatmap(D0_CMV_V4,
                     title    = "Gene Set Scores by Cluster (Visit 4)",
                     filename = "Donor_0/exports/D0_CMV_Heatmap_clusters_V4.png")


#
#
# Stacked bar graphs for donor 0, split by timepoint, with cell # cutoffs
#
#

export_dir <- "Donor_0/exports/clone_cluster_stacked_bars"

# --- Minimum cells per clone per timepoint required to display a bar ---
# Clones below this threshold will still show their label and percentage,
# but the stacked bar area will be left blank.
MIN_CELLS_THRESHOLD <- 10  # <-- adjust this value to change the cutoff

# Extract metadata and parse fields
meta <- D0_data_CMV@meta.data %>%
  mutate(
    specificity  = str_extract(Global_Clone_ID, "(pp65|IE1)"),
    clone_number = str_extract(Global_Clone_ID, "clone_\\d+") %>% str_remove("clone_") %>% as.integer(),
    cluster      = seurat_clusters
  )

# Named color vector for clusters
cluster_colors <- c(
  "0" = "#C49A00",
  "1" = "#53B400",
  "2" = "#8B0000",
  "3" = "#1976D2",
  "4" = "#F8766D",
  "5" = "#A58AFF",
  "6" = "#FFDB58",
  "7" = "#00C094"
)

# Clone label sets (pp65 capped at clone 20)
ie1_labels  <- paste0("IE1 Clone ",  1:7)
pp65_labels <- paste0("pp65 Clone ", 1:20)
dummy_label <- " "
all_labels  <- c(ie1_labels, dummy_label, pp65_labels)


generate_timepoint_plot <- function(df, time_label, show_legend = FALSE) {
  
  clone_summary <- df %>%
    count(specificity, clone_number, name = "n_cells_clone") %>%
    group_by(specificity) %>%
    mutate(clone_size_percent = n_cells_clone / sum(n_cells_clone) * 100) %>%
    ungroup() %>%
    mutate(clone_label = ifelse(specificity == "IE1",
                                paste0("IE1 Clone ",  clone_number),
                                paste0("pp65 Clone ", clone_number)))
  
  comp_data <- df %>%
    count(specificity, clone_number, cluster, name = "n_cells") %>%
    group_by(specificity, clone_number) %>%
    mutate(freq = n_cells / sum(n_cells)) %>%
    ungroup() %>%
    mutate(clone_label = ifelse(specificity == "IE1",
                                paste0("IE1 Clone ",  clone_number),
                                paste0("pp65 Clone ", clone_number))) %>%
    complete(clone_label = c(ie1_labels, pp65_labels),
             cluster     = sort(unique(D0_data_CMV$seurat_clusters)),
             fill        = list(n_cells = 0, freq = 0)) %>%
    left_join(clone_summary %>% select(clone_label, clone_size_percent, n_cells_clone),
              by = "clone_label") %>%
    mutate(
      clone_size_percent = replace_na(clone_size_percent, 0),
      n_cells_clone      = replace_na(n_cells_clone, 0)
    )
  
  # --- Apply cell threshold: zero out bars for clones below MIN_CELLS_THRESHOLD ---
  # Labels/percentages will still display; bars will simply be invisible (freq = 0).
  comp_data <- comp_data %>%
    mutate(freq = ifelse(n_cells_clone >= MIN_CELLS_THRESHOLD, freq, 0))  # <-- threshold used here
  
  dummy_row <- expand.grid(
    clone_label        = dummy_label,
    cluster            = sort(unique(D0_data_CMV$seurat_clusters)),
    n_cells            = 0,
    freq               = 0,
    clone_size_percent = 0,
    n_cells_clone      = 0,
    stringsAsFactors   = FALSE
  )
  comp_data <- bind_rows(comp_data, dummy_row)
  
  comp_data <- comp_data %>%
    mutate(label_color = case_when(
      grepl("IE1",  clone_label) ~ "red",
      grepl("pp65", clone_label) ~ "blue",
      TRUE                       ~ "black"
    ))
  
  comp_data$clone_label <- factor(comp_data$clone_label, levels = rev(all_labels))
  
  dummy_position     <- which(levels(comp_data$clone_label) == dummy_label)
  ie1_break_position <- dummy_position
  
  label_data <- comp_data %>%
    filter(clone_label != dummy_label) %>%
    group_by(clone_label) %>%
    slice(1) %>%
    ungroup()
  
  ggplot(comp_data, aes(x = freq, y = clone_label, fill = factor(cluster))) +
    geom_bar(stat = "identity", width = 0.8) +
    geom_text(
      data        = label_data,
      aes(x = 1.02, y = clone_label,
          label = sprintf("%.1f%%", clone_size_percent),
          color = label_color),
      inherit.aes = FALSE,
      hjust       = 0,
      size        = 3
    ) +
    scale_fill_manual(values = cluster_colors, drop = FALSE) +
    scale_color_identity() +
    scale_x_continuous(
      expand = expansion(mult = c(0, 0.15)),
      labels = scales::percent_format(accuracy = 1)
    ) +
    theme_classic() +
    theme(
      axis.title         = element_blank(),
      axis.text.y        = element_text(
        size  = 8,
        color = rev(c(rep("red", 7), "black", rep("blue", 20)))
      ),
      axis.ticks.y       = element_blank(),
      panel.grid.major.y = element_blank(),
      legend.position    = if (show_legend) "right" else "none",  # <-- controlled here
      legend.title       = element_text(size = 10),
      legend.text        = element_text(size = 8),
      plot.margin        = margin(10, 30, 10, 10)
    ) +
    guides(fill = guide_legend(title = "Metacluster")) +
    geom_hline(yintercept = ie1_break_position, linetype = "dashed") +
    ggtitle(paste("Donor 0 -", time_label, "Timepoint"))
}


# --- Filter to pp65 clones 1–20 only ---
meta_cut <- meta %>%
  filter(specificity == "IE1" | (specificity == "pp65" & clone_number <= 20))

# --- Generate one plot per timepoint ---
p_V1 <- generate_timepoint_plot(meta_cut %>% filter(timepoint == "Visit 1"), "Visit 1")
p_V2 <- generate_timepoint_plot(meta_cut %>% filter(timepoint == "Visit 2"), "Visit 2")
p_V3 <- generate_timepoint_plot(meta_cut %>% filter(timepoint == "Visit 3"), "Visit 3")
p_V4 <- generate_timepoint_plot(meta_cut %>% filter(timepoint == "Visit 4"), "Visit 4", show_legend = TRUE)

# --- Combine and save ---
combined_plot <- p_V1 + p_V2 + p_V3 + p_V4 + plot_layout(ncol = 4, widths = c(1, 1))

ggsave(
  filename = file.path(export_dir, "Donor0_seurat_clusters_stacked_bars_pp65_cutoff_10_cell_threshold.png"),
  plot     = combined_plot,
  width    = 21,
  height   = 0.3 * length(unique(meta_cut$Global_Clone_ID)) + 2
)


#
#
# Gene set analysis from public TIL dataset
#
#

# Download the metadata (contains cluster/cell type labels for all cells)
download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE156728&format=file&file=GSE156728_metadata.txt.gz",
  destfile = "GSE156728_data/metadata.txt.gz"
)

cancer_types <- c("BC", "BCL", "ESCA", "MM", "PACA", "RC", "THCA", "UCEC")

for (ct in cancer_types) {
  url <- paste0(
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE156728&format=file&file=GSE156728_",
    ct, "_10X.CD8.counts.txt.gz"
  )
  download.file(url, destfile = paste0("GSE156728_data/", ct, "_CD8_counts.txt.gz"))
}

# First, make sure metadata rownames are set correctly after fread
meta <- fread("GSE156728_data/metadata.txt.gz", data.table = FALSE)
rownames(meta) <- meta[, 1]  # set first column (cell barcodes) as rownames

# Then use this corrected lapply
seurat_list <- lapply(cancer_types, function(ct) {
  
  # Load counts
  counts <- fread(paste0("GSE156728_data/", ct, "_CD8_counts.txt.gz"), 
                  data.table = FALSE)
  rownames(counts) <- counts[, 1]
  counts <- counts[, -1]
  
  # Convert to sparse matrix BEFORE creating Seurat object
  counts_mat <- as(as.matrix(counts), "dgCMatrix")
  
  # Find shared cells between counts columns and metadata rows
  shared_cells <- intersect(colnames(counts_mat), rownames(meta))
  
  if (length(shared_cells) == 0) {
    message("No shared cells found for cancer type: ", ct)
    return(NULL)
  }
  
  # Subset both to shared cells
  counts_sub <- counts_mat[, shared_cells]
  meta_sub   <- meta[shared_cells, , drop = FALSE]
  
  # Confirm rownames are intact before creating object
  stopifnot(all(colnames(counts_sub) == rownames(meta_sub)))
  
  CreateSeuratObject(
    counts    = counts_sub,
    meta.data = meta_sub,
    project   = ct
  )
})

# Remove any NULLs (cancer types with no shared cells)
seurat_list <- Filter(Negate(is.null), seurat_list)


seurat_merged <- merge(seurat_list[[1]], y = seurat_list[-1])


seurat_merged <- NormalizeData(seurat_merged)
seurat_merged <- FindVariableFeatures(seurat_merged, nfeatures = 2000)
seurat_merged <- ScaleData(seurat_merged)
seurat_merged <- RunPCA(seurat_merged)
seurat_merged <- RunUMAP(seurat_merged, dims = 1:20)

# Visualize using the published cluster labels from the metadata
DimPlot(seurat_merged, group.by = "meta.cluster", label = TRUE)



# Set cluster identity
Idents(seurat_merged) <- "meta.cluster"

seurat_merged <- JoinLayers(seurat_merged)

# Find marker genes for each cluster
til_markers <- FindAllMarkers(
  seurat_merged,
  only.pos    = TRUE,
  min.pct     = 0.25,
  logfc.threshold = 0.5,
  test.use    = "wilcox"
)

# Review the top markers per cluster
top_markers <- til_markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 50)

# Export as a named list of gene sets
gene_sets <- split(top_markers$gene, top_markers$cluster)

# Load into my object
mdata_CMV <- AddModuleScore_UCell(
  mdata_CMV,
  features = gene_sets
)

Idents(mdata_CMV) <- "metacluster"

jpeg("exports/integrated_v5/TIL_UCell_scores.jpeg", 
     width   = 5600, 
     height  = 4200, 
     res     = 300,
     quality = 95)

VlnPlot(
  mdata_CMV,
  features = paste0(names(gene_sets), "_UCell"),
  ncol     = 4,    # 4 cols x 5 rows = 20 slots, fits 17 panels
  pt.size  = 0
)

dev.off()

# Heatmap
# Step 1: Pull UCell scores + metacluster labels into a data frame
ucell_cols <- paste0(names(gene_sets), "_UCell")

score_df <- mdata_CMV@meta.data[, c("metacluster", ucell_cols)]

# Step 2: Average each UCell score per metacluster
avg_scores <- score_df %>%
  group_by(metacluster) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  tibble::column_to_rownames("metacluster")

# Step 3: Scale across gene sets (z-score per column) so 
# differences in magnitude don't dominate the color scale
avg_scaled <- scale(avg_scores)

# Step 4: Clean up column names for readability (remove "_UCell" suffix)
colnames(avg_scaled) <- gsub("_UCell$", "", colnames(avg_scaled))

p <- pheatmap(
  t(avg_scaled),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale        = "none",
  color        = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "TIL Gene Set UCell Scores by CMV Metacluster"
)

ggsave(
  filename  = "exports/integrated_v5/TIL_UCell_heatmap.jpeg",
  plot      = p,
  width     = 10,
  height    = 7,
  dpi       = 300,
  device    = "jpeg"
)

# ── Heatmap: raw average UCell scores (unscaled) ─────────────────────────────

avg_scores_raw <- avg_scores  # already computed above; rows = metaclusters, cols = gene sets

# Clean up column names
colnames(avg_scores_raw) <- gsub("_UCell$", "", colnames(avg_scores_raw))

p_raw <- pheatmap(
  t(avg_scores_raw),
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale        = "none",
  color        = colorRampPalette(c("navy", "white", "firebrick3"))(50),
  border_color = NA,
  fontsize_row = 10,
  fontsize_col = 10,
  main         = "TIL Gene Set UCell Scores by CMV Metacluster (Raw)"
)

ggsave(
  filename = "exports/integrated_v5/TIL_UCell_heatmap_raw.jpeg",
  plot     = p_raw,
  width    = 10,
  height   = 7,
  dpi      = 300,
  device   = "jpeg"
)

# UMAP plots with gene set scores

ucell_cols <- paste0(names(gene_sets), "_UCell")

# Generate one FeaturePlot per gene set
plot_list <- lapply(ucell_cols, function(feature) {
  FeaturePlot(
    mdata_CMV,
    features   = feature,
    reduction  = "umap",
    pt.size    = 0.3,
    order      = TRUE      # plots high-scoring cells on top
  ) +
    scale_colour_gradientn(
      colours = c("lightgrey", "navy", "firebrick3")
    ) +
    labs(title = gsub("_UCell$", "", feature)) +
    theme(
      plot.title   = element_text(size = 9, face = "bold"),
      legend.key.size = unit(0.3, "cm"),
      legend.text  = element_text(size = 7),
      axis.title   = element_blank(),
      axis.text    = element_blank(),
      axis.ticks   = element_blank()
    )
})

# Combine with patchwork into a 4-column grid
combined <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(
    title   = "TIL Gene Set UCell Scores on CMV CD8 T Cell UMAP",
    theme   = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  )

# Export
ggsave(
  filename  = "exports/integrated_v5/TIL_UCell_UMAP_panels.jpeg",
  plot      = combined,
  width     = 15,       # 4 cols * 5 inches
  height    = 30,       # 5 rows * 5 inches (fits 17 in a 4x5 grid)
  dpi       = 300,
  device    = "jpeg",
  limitsize = FALSE
)



#
#
# Gene sets from CMV metaclusters applied to public TIL dataset
#
#

# ── 1. Derive gene sets from CMV metaclusters ─────────────────────────────────

Idents(mdata_CMV) <- "metacluster"

metacluster_markers <- FindAllMarkers(
  mdata_CMV,
  only.pos        = TRUE,
  logfc.threshold = 0.25,
  min.pct         = 0.25,
  test.use        = "wilcox"
)

gene_sets_cmv <- metacluster_markers %>%
  filter(p_val_adj < 0.01) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  slice_head(n = 50) %>%
  summarise(genes = list(gene)) %>%
  deframe()

# ── 2. Apply to seurat_merged with UCell ──────────────────────────────────────

# Filter each gene set to genes present in the TIL object
gene_sets_cmv_filtered <- lapply(gene_sets_cmv, function(genes) {
  intersect(genes, rownames(seurat_merged))
})

seurat_merged <- AddModuleScore_UCell(
  seurat_merged,
  features = gene_sets_cmv_filtered
)

ucell_cols_cmv <- paste0(names(gene_sets_cmv_filtered), "_UCell")

# ── 3. Helper: build average score matrix and plot heatmap ───────────────────

make_til_heatmap <- function(obj, group_var, title, filename,
                             scale = c("none", "row", "column")) {
  
  scale <- match.arg(scale)
  
  score_df <- obj@meta.data[, c(group_var, ucell_cols_cmv)]
  
  avg_mat <- score_df %>%
    group_by(across(all_of(group_var))) %>%
    summarise(across(everything(), mean, na.rm = TRUE)) %>%
    tibble::column_to_rownames(group_var)
  
  colnames(avg_mat) <- gsub("_UCell$", "", colnames(avg_mat))
  
  # Choose color scale: diverging for z-scored, sequential-ish for raw
  color_pal <- if (scale == "row") {
    colorRampPalette(c("navy", "white", "firebrick3"))(50)
  } else {
    colorRampPalette(c("navy", "white", "firebrick3"))(50)
  }
  
  p <- pheatmap(
    t(avg_mat),
    cluster_rows = TRUE,
    cluster_cols = TRUE,
    scale        = scale,        # "row" → z-score each gene set across clusters
    color        = color_pal,
    border_color = NA,
    fontsize_row = 10,
    fontsize_col = 10,
    main         = title
  )
  
  ggsave(
    filename = filename,
    plot     = p,
    width    = 10,
    height   = 6,
    dpi      = 300,
    device   = "jpeg"
  )
}

# ── 4. Heatmap by TIL cluster ─────────────────────────────────────────────────

make_til_heatmap(
  seurat_merged,
  group_var = "meta.cluster",
  title     = "CMV Metacluster Gene Sets in TIL Dataset (by Cluster)",
  filename  = "exports/integrated_v5/CMV_genesets_TIL_by_cluster.jpeg"
)

# Z-scored version
make_til_heatmap(
  seurat_merged,
  group_var = "meta.cluster",
  title     = "CMV Metacluster Gene Sets in TIL Dataset (by Cluster, row Z-score)",
  filename  = "exports/integrated_v5/CMV_genesets_TIL_by_cluster_zscore.jpeg",
  scale     = "row"
)

# ── 5. Heatmap by cancer type ─────────────────────────────────────────────────

make_til_heatmap(
  seurat_merged,
  group_var = "cancerType",
  title     = "CMV Metacluster Gene Sets in TIL Dataset (by Cancer Type)",
  filename  = "exports/integrated_v5/CMV_genesets_TIL_by_cancerType.jpeg"
)

# ── 6. Heatmap by cancer type × tissue location ───────────────────────────────

make_til_heatmap_split <- function(obj, group_vars, title, filename,
                                   width = 16, height = 7,
                                   order_by = c("cancer", "location")) {
  
  order_by <- match.arg(order_by)
  
  score_df <- obj@meta.data[, c(group_vars, ucell_cols_cmv)]
  
  score_df <- score_df %>%
    mutate(group = paste(!!!syms(group_vars), sep = "_"))
  
  avg_mat <- score_df %>%
    group_by(group) %>%
    summarise(across(all_of(ucell_cols_cmv), mean, na.rm = TRUE), .groups = "drop") %>%
    tibble::column_to_rownames("group")
  
  colnames(avg_mat) <- gsub("_UCell$", "", colnames(avg_mat))
  
  col_anno <- data.frame(
    CancerType = sub("_[^_]+$", "", rownames(avg_mat)),
    Location   = sub("^.*_",    "", rownames(avg_mat)),
    row.names  = rownames(avg_mat)
  )
  
  loc_order <- c("T", "N", "P")
  
  col_order <- col_anno %>%
    tibble::rownames_to_column("group") %>%
    mutate(Location = factor(Location, levels = loc_order))
  
  # Order columns: by Location first (all T, then N, then P) or by CancerType first
  col_order <- if (order_by == "location") {
    col_order %>% arrange(Location, CancerType) %>% pull(group)
  } else {
    col_order %>% arrange(CancerType, Location) %>% pull(group)
  }
  
  avg_mat  <- avg_mat[col_order, ]
  col_anno <- col_anno[col_order, , drop = FALSE]
  
  anno_colors <- list(
    Location = c("T" = "firebrick3", "N" = "steelblue3", "P" = "forestgreen"),
    CancerType = setNames(
      colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(
        length(unique(col_anno$CancerType))),
      unique(col_anno$CancerType)
    )
  )
  
  p <- pheatmap(
    t(avg_mat),
    annotation_col    = col_anno,
    annotation_colors = anno_colors,
    cluster_rows      = TRUE,
    cluster_cols      = FALSE,
    scale             = "row",
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(50),
    border_color      = NA,
    fontsize_row      = 10,
    fontsize_col      = 8,
    angle_col         = 45,
    main              = title
  )
  
  ggsave(
    filename  = filename,
    plot      = p,
    width     = width,
    height    = height,
    dpi       = 300,
    device    = "jpeg",
    limitsize = FALSE
  )
}

# Location-first ordering: all T, then all N, then all P
make_til_heatmap_split(
  seurat_merged,
  group_vars = c("cancerType", "loc"),
  title      = "CMV Metacluster Gene Sets in TIL Dataset (by Location × Cancer Type)",
  filename   = "exports/integrated_v5/CMV_genesets_TIL_by_loc_cancerType.jpeg",
  order_by   = "location"
)

make_til_heatmap_split(
  seurat_merged,
  group_vars = c("cancerType", "loc"),
  title      = "CMV Metacluster Gene Sets in TIL Dataset (by Cancer Type × Location)",
  filename   = "exports/integrated_v5/CMV_genesets_TIL_by_cancerType_loc.jpeg"
)



# =============================================================================
# TCR CDR3 Motif Analysis from Seurat Object
# -----------------------------------------------------------------------------
# Input:  Seurat object `mdata_CMV` with metadata columns:
#           - t_cdr3s_aa  : CDR3 AA sequences (e.g. "TRA:CAVXX;TRB:CASSXX")
#           - donor       : donor ID (5 donors)
#           - specificity : "pp65" or "IE1"
# Output: Sequence logos, summary tables, and (optionally) motif clusters
# =============================================================================

# ── 1. PACKAGES ───────────────────────────────────────────────────────────────
# Install any missing packages before loading
required_pkgs <- c(
  "Seurat",       # Seurat object handling
  "dplyr",        # Data wrangling
  "tidyr",        # Reshaping
  "stringr",      # String parsing
  "ggplot2",      # Base plotting
  "ggseqlogo",    # Sequence logo plots  ← primary motif viz
  "patchwork",    # Combine plots
  "viridis"       # Colour scales
)

# Optional but recommended for deeper analysis
optional_pkgs <- c(
  "immunarch",    # Immune repertoire statistics (diversity, clonality)
  "Biostrings",   # AA sequence handling / PWM construction
  "msa"           # Multiple sequence alignment (Bioconductor)
)

install_if_missing <- function(pkgs, bioc = FALSE) {
  missing <- pkgs[!pkgs %in% installed.packages()[, "Package"]]
  if (length(missing)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(missing, ask = FALSE)
    } else {
      install.packages(missing)
    }
  }
}

install_if_missing(required_pkgs)
install_if_missing(c("Biostrings", "msa"), bioc = TRUE)   # Bioconductor

suppressPackageStartupMessages({
  lapply(required_pkgs, library, character.only = TRUE)
})

# Load optional packages quietly (won't error if absent)
for (pkg in optional_pkgs) {
  suppressPackageStartupMessages(
    tryCatch(library(pkg, character.only = TRUE),
             error = function(e) message("Optional package not loaded: ", pkg))
  )
}


# ── 2. EXTRACT METADATA FROM SEURAT ──────────────────────────────────────────

#' Pull the relevant metadata columns out of the Seurat object into a tidy
#' data frame, then parse the CDR3 string into separate TRA / TRB columns.
#'
#' The `t_cdr3s_aa` field (produced by scRepertoire / Seurat integration)
#' has the format:  "TRA:CAVXXXXF;TRB:CASSXXXXF"
#' Some cells carry only one chain; some carry multiple TRA alleles
#' (separated by ";TRA:").  We take the FIRST TRA and FIRST TRB per cell.

extract_cdr3_df <- function(seurat_obj,
                            cdr3_col    = "t_cdr3s_aa",
                            extra_cols  = c("donor", "specificity")) {
  
  # Pull metadata
  meta <- seurat_obj@meta.data %>%
    tibble::rownames_to_column("cell_barcode") %>%
    dplyr::select(cell_barcode,
                  all_of(c(cdr3_col, extra_cols))) %>%
    dplyr::rename(cdr3_raw = all_of(cdr3_col))
  
  # ── Parse TRA and TRB from the combined string ──────────────────────────
  # Handles formats:
  #   "TRA:CAVXX;TRB:CASSXX"
  #   "TRB:CASSXX"            (TRA missing)
  #   "TRA:CAVXX"             (TRB missing)
  #   NA / "None"             (no VDJ recovered)
  
  parse_chain <- function(raw_string, chain = c("TRA", "TRB")) {
    chain <- match.arg(chain)
    ifelse(
      is.na(raw_string) | raw_string %in% c("None", "none", ""),
      NA_character_,
      {
        # Extract first occurrence of the requested chain
        pattern <- paste0("(?<=", chain, ":)[^;]+")
        m <- regmatches(raw_string,
                        regexpr(pattern, raw_string, perl = TRUE))
        ifelse(length(m) == 0 | m == "", NA_character_, m)
      }
    )
  }
  
  meta <- meta %>%
    dplyr::mutate(
      TRA = mapply(parse_chain, cdr3_raw, MoreArgs = list(chain = "TRA")),
      TRB = mapply(parse_chain, cdr3_raw, MoreArgs = list(chain = "TRB"))
    )
  
  return(meta)
}

# Run extraction (assumes `mdata_CMV` is already loaded in your environment)
cdr3_df <- extract_cdr3_df(
  seurat_obj  = mdata_CMV,
  cdr3_col    = "t_cdr3s_aa",
  extra_cols  = c("donor", "specificity")
)

message(sprintf(
  "Extracted %d cells | TRA recovered: %d | TRB recovered: %d",
  nrow(cdr3_df),
  sum(!is.na(cdr3_df$TRA)),
  sum(!is.na(cdr3_df$TRB))
))


# ── 3. HELPER: STANDARDISE SEQUENCES FOR MOTIF ANALYSIS ──────────────────────

#' Prepare a character vector of CDR3 sequences for sequence-logo / motif work.
#'
#' @param seqs   Character vector of AA sequences (NAs allowed)
#' @param min_n  Minimum number of valid sequences required (returns NULL if fewer)
#' @param trim_c Remove the conserved C at position 1? (FALSE by default)
#' @param trim_f Remove the conserved F/W at the last position? (FALSE by default)
#' @return Cleaned character vector, or NULL

prepare_sequences <- function(seqs,
                              min_n   = 5,
                              trim_c  = FALSE,
                              trim_f  = FALSE) {
  seqs <- seqs[!is.na(seqs) & nchar(seqs) > 0]
  seqs <- seqs[grepl("^[ACDEFGHIKLMNPQRSTVWY]+$", seqs)]  # standard AA only
  
  if (length(seqs) < min_n) return(NULL)
  
  if (trim_c) seqs <- substr(seqs, 2, nchar(seqs))         # drop leading C
  if (trim_f) seqs <- substr(seqs, 1, nchar(seqs) - 1)     # drop trailing F/W
  
  return(seqs)
}


# ── 4. CORE MOTIF ANALYSIS FUNCTION ──────────────────────────────────────────

#' Build positional weight matrices and sequence logos for a set of CDR3s.
#'
#' ggseqlogo handles variable-length sequences by aligning from the N-terminus;
#' for CDR3s it is often more informative to align from both ends.  We therefore
#' produce three views: full sequence, N-terminal 4 AA, C-terminal 4 AA.
#'
#' @param seqs       Character vector of cleaned CDR3 AA sequences
#' @param title      Title string for plots
#' @param n_term_len Length of N-terminal window (default 4)
#' @param c_term_len Length of C-terminal window (default 4)
#' @return A named list: $full, $nterm, $cterm  (ggplot objects)

make_logo_plots <- function(seqs,
                            title      = "",
                            n_term_len = 4,
                            c_term_len = 4) {
  
  # ── Full-length logo (variable length — ggseqlogo stacks by position)
  p_full <- ggseqlogo(seqs, method = "prob") +
    ggtitle(paste0(title, "\nFull CDR3  (n=", length(seqs), ")")) +
    theme(plot.title = element_text(size = 9, face = "bold"))
  
  # ── N-terminal window
  nterm_seqs <- substr(seqs, 1, n_term_len)
  nterm_seqs <- nterm_seqs[nchar(nterm_seqs) == n_term_len]
  p_nterm <- ggseqlogo(nterm_seqs, method = "prob") +
    ggtitle(paste0("N-term [1:", n_term_len, "]")) +
    theme(plot.title = element_text(size = 9))
  
  # ── C-terminal window
  cterm_seqs <- substr(seqs, nchar(seqs) - c_term_len + 1, nchar(seqs))
  cterm_seqs <- cterm_seqs[nchar(cterm_seqs) == c_term_len]
  p_cterm <- ggseqlogo(cterm_seqs, method = "prob") +
    ggtitle(paste0("C-term [-", c_term_len, ":]")) +
    theme(plot.title = element_text(size = 9))
  
  return(list(full = p_full, nterm = p_nterm, cterm = p_cterm))
}


# ── 5. LENGTH DISTRIBUTION HELPER ─────────────────────────────────────────────

#' CDR3 length distribution bar chart

plot_length_dist <- function(seqs, title = "") {
  data.frame(len = nchar(seqs)) %>%
    dplyr::count(len) %>%
    ggplot(aes(x = len, y = n, fill = n)) +
    geom_col(colour = "white", width = 0.8) +
    scale_fill_viridis_c(guide = "none") +
    labs(title   = title,
         x       = "CDR3 length (AA)",
         y       = "Cell count") +
    theme_minimal(base_size = 9) +
    theme(plot.title = element_text(size = 9, face = "bold"))
}


# ── 6. MAIN ANALYSIS LOOP: donor × specificity ────────────────────────────────

#' Iterate over every unique donor × specificity combination, generate motif
#' plots for TRA and TRB, and save one PDF per combination.

run_motif_analysis <- function(df,
                               output_dir = "exports/integrated_v5/tcr_motif_output",
                               chain      = c("TRA", "TRB"),
                               min_cells  = 5) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  chain <- match.arg(chain, several.ok = TRUE)
  
  # All unique combinations present in the data
  combos <- df %>%
    dplyr::filter(!is.na(donor), !is.na(specificity)) %>%
    dplyr::distinct(donor, specificity) %>%
    dplyr::arrange(donor, specificity)
  
  message(sprintf("Found %d donor × specificity combinations", nrow(combos)))
  
  # Storage for summary statistics
  summary_list <- list()
  
  for (i in seq_len(nrow(combos))) {
    
    don  <- combos$donor[i]
    spec <- combos$specificity[i]
    tag  <- paste0(don, "_", spec)
    
    message(sprintf("  Processing: %s", tag))
    
    subset_df <- df %>%
      dplyr::filter(donor == don, specificity == spec)
    
    plot_list <- list()  # collect all panels for this combo
    
    for (ch in chain) {
      
      seqs <- prepare_sequences(subset_df[[ch]], min_n = min_cells)
      
      if (is.null(seqs)) {
        message(sprintf("    %s: fewer than %d sequences — skipping", ch, min_cells))
        next
      }
      
      # ── Sequence logos
      logos <- make_logo_plots(
        seqs  = seqs,
        title = sprintf("%s | %s | %s", ch, don, spec)
      )
      
      # ── Length distribution
      p_len <- plot_length_dist(
        seqs  = seqs,
        title = paste0(ch, " length distribution")
      )
      
      # Assemble: length bar | full logo | nterm | cterm
      combined <- (p_len | logos$full) / (logos$nterm | logos$cterm) +
        plot_annotation(
          title    = sprintf("Donor: %s  |  Specificity: %s  |  Chain: %s", don, spec, ch),
          subtitle = sprintf("n = %d cells", length(seqs)),
          theme    = theme(plot.title    = element_text(size = 12, face = "bold"),
                           plot.subtitle = element_text(size = 10))
        )
      
      plot_list[[ch]] <- combined
      
      # ── Per-chain summary stats
      summary_list[[paste(tag, ch, sep = "_")]] <- data.frame(
        donor         = don,
        specificity   = spec,
        chain         = ch,
        n_cells       = length(seqs),
        n_unique      = length(unique(seqs)),
        clonality     = 1 - (length(unique(seqs)) / length(seqs)),   # simple clonality index
        median_length = median(nchar(seqs)),
        modal_length  = as.integer(names(sort(table(nchar(seqs)),
                                              decreasing = TRUE)[1]))
      )
    }
    
    # ── Save plots to PDF
    if (length(plot_list) > 0) {
      pdf_path <- file.path(output_dir, paste0(tag, "_motif_analysis.pdf"))
      pdf(pdf_path, width = 14, height = 8)
      for (p in plot_list) print(p)
      dev.off()
      message(sprintf("    Saved → %s", pdf_path))
    }
  }
  
  # Combine summary table
  summary_df <- dplyr::bind_rows(summary_list)
  write.csv(summary_df,
            file = file.path(output_dir, "motif_summary_stats.csv"),
            row.names = FALSE)
  message("\nSummary statistics written to motif_summary_stats.csv")
  
  return(invisible(summary_df))
}


# ── 7. BONUS: CROSS-DONOR OVERLAY LOGOS ──────────────────────────────────────
#' For each specificity group, pool all donors and plot a combined logo.
#' Useful for identifying conserved motifs across donors.

run_pooled_logos <- function(df,
                             output_dir = "exports/integrated_v5/tcr_motif_output",
                             chain      = c("TRA", "TRB"),
                             min_cells  = 10) {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  chain <- match.arg(chain, several.ok = TRUE)
  
  specificities <- unique(na.omit(df$specificity))
  
  for (spec in specificities) {
    plot_list <- list()
    for (ch in chain) {
      seqs <- df %>%
        dplyr::filter(specificity == spec) %>%
        dplyr::pull(ch) %>%
        prepare_sequences(min_n = min_cells)
      
      if (is.null(seqs)) next
      
      logos <- make_logo_plots(
        seqs  = seqs,
        title = sprintf("ALL DONORS | %s | %s  (n=%d)", ch, spec, length(seqs))
      )
      p_len <- plot_length_dist(seqs, title = "Length dist (all donors)")
      
      combined <- (p_len | logos$full) / (logos$nterm | logos$cterm) +
        plot_annotation(
          title = sprintf("Pooled — Specificity: %s  |  Chain: %s", spec, ch),
          theme = theme(plot.title = element_text(size = 12, face = "bold"))
        )
      plot_list[[ch]] <- combined
    }
    
    if (length(plot_list) > 0) {
      pdf_path <- file.path(output_dir,
                            paste0("POOLED_", spec, "_motif_analysis.pdf"))
      pdf(pdf_path, width = 14, height = 8)
      for (p in plot_list) print(p)
      dev.off()
      message(sprintf("Pooled logo saved → %s", pdf_path))
    }
  }
}


# ── 8. BONUS: IMMUNARCH CLONOTYPE STATISTICS (if package is available) ────────
#' Computes diversity indices (Shannon, Gini-Simpson), top clonotype tables,
#' and V-gene usage per donor × specificity group.

run_immunarch_stats <- function(df,
                                output_dir = "tcr_motif_output",
                                chain      = "TRB") {
  
  if (!requireNamespace("immunarch", quietly = TRUE)) {
    message("immunarch not installed — skipping repertoire statistics.")
    return(invisible(NULL))
  }
  library(immunarch)
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # immunarch expects a named list of data frames, each with a "CDR3.aa" column
  groups <- df %>%
    dplyr::filter(!is.na(donor), !is.na(specificity)) %>%
    dplyr::mutate(group = paste0(donor, "_", specificity)) %>%
    dplyr::filter(!is.na(.data[[chain]]))
  
  immdata_list <- split(groups, groups$group)
  immdata_list <- lapply(immdata_list, function(d) {
    data.frame(
      CDR3.aa = d[[chain]],
      Clones  = 1L   # each row = one cell
    )
  })
  
  # Build immunarch object
  imm <- list(data = immdata_list,
              meta = data.frame(Sample = names(immdata_list)))
  
  # Clonality / diversity
  div <- repDiversity(imm$data, .method = "gini.simp")
  print(div)
  
  # Top-10 clonotypes per group
  top_clones <- repClonality(imm$data, .method = "top", .head = 10)
  
  # Save
  saveRDS(imm, file.path(output_dir, "immunarch_object.rds"))
  write.csv(div,        file.path(output_dir, "diversity_gini_simpson.csv"))
  write.csv(top_clones, file.path(output_dir, "top10_clonotypes.csv"))
  
  message("immunarch stats saved.")
  return(invisible(imm))
}


# ── 9. RUN EVERYTHING ─────────────────────────────────────────────────────────

OUT_DIR <- "exports/integrated_v5/tcr_motif_output"   # ← change to your preferred output path

# Main per-donor × specificity analysis
summary_stats <- run_motif_analysis(
  df         = cdr3_df,
  output_dir = OUT_DIR,
  chain      = c("TRA", "TRB"),
  min_cells  = 5
)

# Cross-donor pooled logos per specificity
run_pooled_logos(
  df         = cdr3_df,
  output_dir = OUT_DIR,
  chain      = c("TRA", "TRB"),
  min_cells  = 10
)

# immunarch repertoire statistics (TRB only — richer diversity signal)
run_immunarch_stats(
  df         = cdr3_df,
  output_dir = OUT_DIR,
  chain      = "TRB"
)

# Quick peek at the summary table
print(summary_stats)




# Adding TCR sequences straight from the contig files

# Extract everything BEFORE the 16-char 10x barcode (format: XXXXXXXXXXXXXXXX-1)
barcodes <- colnames(mdata_CMV)

# The raw 10x barcode is always 16 uppercase letters + "-1" (+ optional "_2" suffix)
# The prefix is everything before it
prefixes <- gsub("_?[ACGT]{16}-1.*$", "", barcodes)

# Show all unique prefixes and how many cells each has
sort(table(prefixes), decreasing = TRUE)


# Check suffixes for project3 barcodes
barcodes <- colnames(mdata_CMV)
p2_barcodes <- barcodes[grepl("^project2", barcodes)]
# Extract what comes after the 16-char barcode
suffixes <- gsub(".*[ACGT]{16}-1", "", p2_barcodes)
table(suffixes)


###############################################################################
# add_cdr3_from_contigs.R
#
# Reads all filtered_contig_annotations CSV files and patches t_cdr3s_aa
# (filling NAs for donors 4 & 5) and adds t_cdr3s_nt (all donors) into
# mdata_CMV by matching on cell barcodes.
#
# Matching strategy
# -----------------
# project1, project2_1, project2_2:
#   Each lane contains one donor. Match key = raw 10x barcode (e.g. ACGT...-1)
#
# project3_donor4 / project3_donor5:
#   Both donors are mixed across both lanes (_1 and _2). Donor assignment was
#   already determined by hashtag demultiplexing when the Seurat object was
#   built. We therefore match on raw_barcode + lane_suffix (e.g. ACGT...-1_1)
#   which uniquely identifies a cell within the pool, and let the Seurat
#   barcode's donor prefix carry the donor identity.
#
# Confirmed Seurat barcode formats:
#   project1_AAACCTGCACCACGTG-1          <- pool565-1_1
#   project2_1_AAACCTGAGTACGACG-1        <- pool565-2_1
#   project2_2_AAACCTGAGAGGTAGA-1        <- pool565-2_2
#   project3_donor4_AAACCAAAGAAAGTCA-1_1 <- pool565-3_1 (lane _1)
#   project3_donor4_AAACCAAAGAAAGTCA-1_2 <- pool565-3_2 (lane _2)
#   project3_donor5_GTTCTATGTAGGATAG-1_1 <- pool565-3_1 (lane _1)
#   project3_donor5_GTTCTATGTAGGATAG-1_2 <- pool565-3_2 (lane _2)
###############################################################################


# ── 0. PACKAGES ───────────────────────────────────────────────────────────────
library(dplyr)
library(tidyr)
library(readr)


# ── 1. FILE CONFIGURATION ─────────────────────────────────────────────────────
# Edit file paths to match your system.
#
# lane_suffix: appended to raw barcodes to form the lookup key.
#   ""   for project1 and project2 (no suffix in Seurat barcodes)
#   "_1" for pool565-3_1 (lane 1 cells have _1 in their Seurat barcode)
#   "_2" for pool565-3_2 (lane 2 cells have _2 in their Seurat barcode)

contig_config <- list(
  
  list(
    file        = "scTCRseq_contig_files/pool565-1_1_filtered_contig_annotations.csv",
    lane_suffix = ""
  ),
  
  list(
    file        = "scTCRseq_contig_files/pool565-2_1_filtered_contig_annotations.csv",
    lane_suffix = ""
  ),
  
  list(
    file        = "scTCRseq_contig_files/pool565-2_2_filtered_contig_annotations.csv",
    lane_suffix = ""
  ),
  
  list(
    file        = "scTCRseq_contig_files/pool565-3_1_filtered_contig_annotations.csv",
    lane_suffix = "_1"
  ),
  
  list(
    file        = "scTCRseq_contig_files/pool565-3_2_filtered_contig_annotations.csv",
    lane_suffix = "_2"
  )
)


# ── 2. READ AND AGGREGATE ONE CONTIG FILE ─────────────────────────────────────

read_contigs <- function(file, lane_suffix) {
  
  df <- read_csv(file, show_col_types = FALSE)
  
  # pool565-3 files have an extra leading `sample` column — drop it
  if ("sample" %in% colnames(df)) {
    df <- select(df, -sample)
  }
  
  df <- df %>%
    select(barcode, chain, productive, high_confidence, cdr3, cdr3_nt, umis) %>%
    filter(
      productive      == TRUE,
      high_confidence == TRUE,
      chain           %in% c("TRA", "TRB"),
      !is.na(cdr3),    cdr3    != "", cdr3    != "None",
      !is.na(cdr3_nt), cdr3_nt != ""
    )
  
  if (nrow(df) == 0) {
    warning(sprintf("No productive TRA/TRB contigs found in %s", basename(file)))
    return(NULL)
  }
  
  # For cells with multiple contigs of the same chain, keep the highest-UMI one
  df <- df %>%
    group_by(barcode, chain) %>%
    slice_max(order_by = umis, n = 1, with_ties = FALSE) %>%
    ungroup()
  
  # Pivot to one row per cell
  df_wide <- df %>%
    select(barcode, chain, cdr3, cdr3_nt) %>%
    pivot_wider(
      names_from  = chain,
      values_from = c(cdr3, cdr3_nt),
      names_glue  = "{chain}_{.value}"
      # produces: TRA_cdr3, TRB_cdr3, TRA_cdr3_nt, TRB_cdr3_nt
    )
  
  # Build combined "TRA:XXX;TRB:YYY" strings matching the existing format
  df_wide <- df_wide %>%
    mutate(
      t_cdr3s_aa = case_when(
        !is.na(TRA_cdr3) & !is.na(TRB_cdr3) ~ paste0("TRA:", TRA_cdr3,    ";TRB:", TRB_cdr3),
        !is.na(TRA_cdr3)                     ~ paste0("TRA:", TRA_cdr3),
        !is.na(TRB_cdr3)                     ~ paste0("TRB:", TRB_cdr3),
        TRUE                                 ~ NA_character_
      ),
      t_cdr3s_nt = case_when(
        !is.na(TRA_cdr3_nt) & !is.na(TRB_cdr3_nt) ~ paste0("TRA:", TRA_cdr3_nt, ";TRB:", TRB_cdr3_nt),
        !is.na(TRA_cdr3_nt)                        ~ paste0("TRA:", TRA_cdr3_nt),
        !is.na(TRB_cdr3_nt)                        ~ paste0("TRB:", TRB_cdr3_nt),
        TRUE                                        ~ NA_character_
      ),
      # Lookup key = raw barcode + lane suffix
      # For project1/2: "AAACCTGCACCACGTG-1"
      # For project3:   "AAACCAAAGAAAGTCA-1_1"  or  "GTTCTATGTAGGATAG-1_2"
      lookup_key = paste0(barcode, lane_suffix)
    ) %>%
    select(lookup_key, t_cdr3s_aa, t_cdr3s_nt)
  
  message(sprintf("  %-50s -> %d cells (lane suffix: '%s')",
                  basename(file), nrow(df_wide), lane_suffix))
  return(df_wide)
}


# ── 3. BUILD MASTER LOOKUP TABLE ──────────────────────────────────────────────

lookup_list <- lapply(contig_config, function(cfg) {
  read_contigs(cfg$file, cfg$lane_suffix)
})

lookup_list <- Filter(Negate(is.null), lookup_list)

# Combine; if the same key appears twice (shouldn't happen) keep the first
cdr3_lookup <- bind_rows(lookup_list) %>%
  distinct(lookup_key, .keep_all = TRUE)

message(sprintf("\nTotal unique keys in CDR3 lookup: %d\n", nrow(cdr3_lookup)))


# ── 4. BUILD LOOKUP KEY FROM SEURAT BARCODES ──────────────────────────────────
# Strip the donor/project prefix from each Seurat barcode to get the lookup key.
# The raw 10x barcode is always 16 uppercase ACGT chars + "-1" + optional "_1"/"_2".
# Everything before that is the prefix we want to remove.
#
# e.g. "project3_donor4_AAACCAAAGAAAGTCA-1_1"
#       -> strip "project3_donor4_"
#       -> lookup key = "AAACCAAAGAAAGTCA-1_1"

extract_lookup_key <- function(seurat_barcodes) {
  # Match: 16-char barcode + "-1" + optional "_1" or "_2"
  # Return just that matched portion (the lookup key)
  regmatches(
    seurat_barcodes,
    regexpr("[ACGT]{16}-1(_[12])?$", seurat_barcodes)
  )
}


# ── 5. DIAGNOSTIC: CHECK MATCH RATES BEFORE PATCHING ─────────────────────────

seurat_barcodes <- colnames(mdata_CMV)
lookup_keys     <- extract_lookup_key(seurat_barcodes)

n_match <- sum(lookup_keys %in% cdr3_lookup$lookup_key)
message(sprintf(
  "Overall barcode match: %d / %d Seurat cells (%.1f%%)\n",
  n_match, length(seurat_barcodes),
  100 * n_match / length(seurat_barcodes)
))

# Per-donor breakdown
overlap_check <- mdata_CMV@meta.data %>%
  tibble::rownames_to_column("cell") %>%
  mutate(
    lookup_key = extract_lookup_key(cell),
    matched    = lookup_key %in% cdr3_lookup$lookup_key
  ) %>%
  group_by(donor_identity) %>%
  summarise(
    total_cells   = n(),
    matched_cells = sum(matched),
    pct_matched   = round(100 * matched_cells / total_cells, 1),
    .groups = "drop"
  )

message("── Match rates by donor ──")
print(overlap_check)

# Uncomment to hard-stop if any donor matches fewer than 5% of cells:
# if (any(overlap_check$pct_matched < 5)) {
#   stop("Low match rate for one or more donors — check contig_config above.")
# }


# ── 6. PATCH SEURAT METADATA ──────────────────────────────────────────────────

aa_lut <- setNames(cdr3_lookup$t_cdr3s_aa, cdr3_lookup$lookup_key)
nt_lut <- setNames(cdr3_lookup$t_cdr3s_nt, cdr3_lookup$lookup_key)

meta        <- mdata_CMV@meta.data
cells       <- rownames(meta)
cell_keys   <- extract_lookup_key(cells)

# Ensure columns exist
if (!"t_cdr3s_aa" %in% colnames(meta)) meta$t_cdr3s_aa <- NA_character_
if (!"t_cdr3s_nt" %in% colnames(meta)) meta$t_cdr3s_nt <- NA_character_

matched_aa <- aa_lut[cell_keys]
matched_nt <- nt_lut[cell_keys]

# t_cdr3s_aa: only fill currently-NA cells (preserves existing donors 1-3 data)
meta$t_cdr3s_aa <- ifelse(
  is.na(meta$t_cdr3s_aa) & !is.na(matched_aa),
  matched_aa,
  meta$t_cdr3s_aa
)

# t_cdr3s_nt: fill where currently NA (new column for all donors)
meta$t_cdr3s_nt <- ifelse(
  is.na(meta$t_cdr3s_nt) & !is.na(matched_nt),
  matched_nt,
  meta$t_cdr3s_nt
)

mdata_CMV@meta.data <- meta


# ── 7. VERIFY ─────────────────────────────────────────────────────────────────

message("\n── Post-patch CDR3 coverage by donor ──")

mdata_CMV@meta.data %>%
  group_by(donor_identity) %>%
  summarise(
    total   = n(),
    has_aa  = sum(!is.na(t_cdr3s_aa)),
    has_nt  = sum(!is.na(t_cdr3s_nt)),
    pct_aa  = round(100 * has_aa / total, 1),
    pct_nt  = round(100 * has_nt / total, 1),
    .groups = "drop"
  ) %>%
  print()

message("\nSample rows for donors 4 & 5:")
mdata_CMV@meta.data %>%
  filter(grepl("Donor 4|Donor 5", donor_identity), !is.na(t_cdr3s_aa)) %>%
  select(donor_identity, cloneId, t_cdr3s_aa, t_cdr3s_nt) %>%
  head(6) %>%
  print()


# ── 8. APPLY TO mdata TOO (optional) ─────────────────────────────────────────
# Uncomment if you want the parent object updated as well.

# meta_full     <- mdata@meta.data
# cells_full    <- rownames(meta_full)
# cell_keys_full <- extract_lookup_key(cells_full)
# if (!"t_cdr3s_aa" %in% colnames(meta_full)) meta_full$t_cdr3s_aa <- NA_character_
# if (!"t_cdr3s_nt" %in% colnames(meta_full)) meta_full$t_cdr3s_nt <- NA_character_
# meta_full$t_cdr3s_aa <- ifelse(is.na(meta_full$t_cdr3s_aa) & !is.na(aa_lut[cell_keys_full]), aa_lut[cell_keys_full], meta_full$t_cdr3s_aa)
# meta_full$t_cdr3s_nt <- ifelse(is.na(meta_full$t_cdr3s_nt) & !is.na(nt_lut[cell_keys_full]), nt_lut[cell_keys_full], meta_full$t_cdr3s_nt)
# mdata@meta.data <- meta_full


# ── 9. SAVE ───────────────────────────────────────────────────────────────────

saveRDS(mdata_CMV, file = "final_mdata_CMV.rds")
# saveRDS(mdata, file = "final_mdata.rds")



# Export metacluster composition per Global_Clone_ID
# Add this to your main Seurat notebook and run with mdata_CMV loaded

library(dplyr)
library(tidyr)

# Define the expected metacluster levels (sets column order in output)
metacluster_levels <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")

metacluster_composition <- mdata_CMV@meta.data %>%
  filter(!is.na(Global_Clone_ID), !is.na(metacluster)) %>%
  group_by(Global_Clone_ID, metacluster) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  group_by(Global_Clone_ID) %>%
  mutate(
    total_cells = sum(n_cells),
    pct         = n_cells / total_cells * 100
  ) %>%
  ungroup() %>%
  select(Global_Clone_ID, metacluster, pct) %>%
  pivot_wider(
    names_from   = metacluster,
    values_from  = pct,
    values_fill  = 0,
    names_prefix = "pct_"
  ) %>%
  # Ensure all 5 metacluster columns are present even if a cluster has 0 cells
  # for all clones (unlikely but safe)
  { missing <- setdiff(paste0("pct_", metacluster_levels), names(.))
  for (col in missing) .[[col]] <- 0
  . } %>%
  select(Global_Clone_ID, paste0("pct_", metacluster_levels)) %>%
  # Sort by donor then clone number, matching the ordering used in TCRseq analysis
  mutate(
    donor      = str_extract(Global_Clone_ID, "^D\\d+"),
    clone_num  = as.integer(str_extract(Global_Clone_ID, "\\d+$"))
  ) %>%
  arrange(donor, clone_num) %>%
  select(-donor, -clone_num)

write.csv(
  metacluster_composition,
  file      = "exports/Global_Clone_metacluster_composition.csv",
  row.names = FALSE
)

message(sprintf(
  "Exported metacluster composition for %d global clones",
  nrow(metacluster_composition)
))

# Quick sanity check — each row should sum to 100
row_sums <- metacluster_composition %>%
  select(starts_with("pct_")) %>%
  rowSums()
stopifnot(all(abs(row_sums - 100) < 0.01))




#
#
# Adding Emily Ford's TCR clusters to the Donor 0 data
#
#

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(readr)

# ── 1. Load cluster CSV files ────────────────────────────────────────────────

df_et48 <- read_csv("/box/Campbell/Sam K/CMV_project/####_RNAseq/# Donor 0/cmv_seurat_dataset_paired_clusters_et48_032626.csv")
df_et36 <- read_csv("/box/Campbell/Sam K/CMV_project/####_RNAseq/# Donor 0/cmv_seurat_dataset_paired_clusters_et36_032626.csv")

# ── 2. TCR parsing and cluster assignment ────────────────────────────────────

assign_tcr_cluster <- function(cdr3_vec, cluster_df) {
  
  cluster_df <- cluster_df[!is.na(cluster_df$cluster), ]
  
  sapply(cdr3_vec, function(x) {
    
    if (is.na(x)) return(NA_character_)
    
    parts    <- strsplit(x, ";")[[1]]
    trb_seqs <- sub("^TRB:", "", parts[grepl("^TRB:", parts)])
    tra_seqs <- sub("^TRA:", "", parts[grepl("^TRA:", parts)])
    
    if (length(trb_seqs) > 0) {
      hits <- cluster_df[cluster_df$cdr3_b_aa %in% trb_seqs, ]
      if (nrow(hits) == 0) return(NA_character_)
      
      if (length(tra_seqs) > 0) {
        hits_narrow <- hits[hits$cdr3_a_aa %in% tra_seqs, ]
        if (nrow(hits_narrow) > 0) hits <- hits_narrow
      }
      
    } else if (length(tra_seqs) > 0) {
      hits <- cluster_df[cluster_df$cdr3_a_aa %in% tra_seqs, ]
      if (nrow(hits) == 0) return(NA_character_)
      
    } else {
      return(NA_character_)
    }
    
    best <- hits[which.max(hits$count), ]
    return(as.character(best$cluster))
    
  }, USE.NAMES = FALSE)
}

# ── 3. Add metadata to Seurat object ─────────────────────────────────────────

D0_data_CMV$TCR_cluster_et48 <- assign_tcr_cluster(
  D0_data_CMV$t_cdr3s_aa, df_et48
)

D0_data_CMV$TCR_cluster_et36 <- assign_tcr_cluster(
  D0_data_CMV$t_cdr3s_aa, df_et36
)

for (tag in c("et48", "et36")) {
  col        <- paste0("TCR_cluster_", tag)
  n_assigned <- sum(!is.na(D0_data_CMV@meta.data[[col]]))
  n_total    <- ncol(D0_data_CMV)
  message(sprintf("[%s] Assigned: %d / %d cells (%.1f%%)",
                  tag, n_assigned, n_total, 100 * n_assigned / n_total))
}

# ── 4. Color maps ─────────────────────────────────────────────────────────────

clone_colors_pp65 <- c(
  "Clone 1"     = "#1b9e77",
  "Clone 2"     = "#d95f02",
  "Clone 3"     = "#7570b3",
  "Clone 4"     = "#e7298a",
  "Clone 5"     = "#66a61e",
  "Clone 6"     = "#e6ab02",
  "Clone 7"     = "#a6761d",
  "Clone 8"     = "#666666",
  "Clone 9"     = "#1f78b4",
  "Clone 10"    = "#b2df8a",
  "Clone 11"    = "#33a02c",
  "Clone 12"    = "#fb9a99",
  "Clone 13"    = "#e31a1c",
  "Clone 14"    = "#fdbf6f",
  "Clone 15"    = "#ff7f00",
  "Clone 16"    = "#cab2d6",
  "Clone 17"    = "#6a3d9a",
  "Clone 18"    = "#ffff99",
  "Clone 19"    = "#b15928",
  "Clone 20"    = "#8dd3c7",
  "Other (<1%)" = "grey70"
)

clone_colors_ie1 <- c(
  "Clone 1"     = "#1b9e77",
  "Clone 2"     = "#d95f02",
  "Clone 3"     = "#7570b3",
  "Clone 4"     = "#e7298a",
  "Clone 5"     = "#66a61e",
  "Clone 6"     = "#e6ab02",
  "Clone 7"     = "#a6761d",
  "Other (<1%)" = "grey70"
)

spec_color_maps <- list(
  pp65 = clone_colors_pp65,
  IE1  = clone_colors_ie1
)

# ── 5. Visualisation helpers ──────────────────────────────────────────────────

get_clone_label <- function(global_clone_ids, spec) {
  pattern <- paste0("^D0_", spec, "_clone_(\\d+)$")
  ifelse(
    grepl(pattern, global_clone_ids, perl = TRUE),
    paste0("Clone ", sub(pattern, "\\1", global_clone_ids, perl = TRUE)),
    "Other (<1%)"
  )
}

make_legend_plot <- function(color_map) {
  dummy_df <- data.frame(
    x     = seq_along(color_map),
    y     = 1,
    label = factor(names(color_map), levels = names(color_map))
  )
  ggplot(dummy_df, aes(x = x, y = y, color = label)) +
    geom_point(size = 3) +
    scale_color_manual(values = color_map, name = "Global Clone ID",
                       drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
    theme_void() +
    theme(legend.position = "right")
}

get_top_clusters <- function(seurat_obj, meta_col, target_spec, n_top) {
  seurat_obj@meta.data %>%
    filter(specificity == target_spec,
           !is.na(.data[[meta_col]])) %>%
    count(.data[[meta_col]], name = "n_cells") %>%
    arrange(desc(n_cells)) %>%
    slice_head(n = n_top) %>%
    pull(.data[[meta_col]])
}

# timepoint_filter = NULL means use all cells
make_cluster_umap_colored <- function(seurat_obj, meta_col, cluster_id,
                                      spec, color_map,
                                      timepoint_filter = NULL) {
  
  md <- seurat_obj@meta.data
  
  if (!is.null(timepoint_filter)) {
    cells_keep <- rownames(md)[md$timepoint %in% timepoint_filter]
    seurat_obj <- subset(seurat_obj, cells = cells_keep)
    md         <- seurat_obj@meta.data
  }
  
  umap_df      <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  umap_cols    <- colnames(umap_df)
  umap_df$cell <- rownames(umap_df)
  
  in_cluster <- !is.na(md[[meta_col]]) & md[[meta_col]] == as.character(cluster_id)
  n          <- sum(in_cluster)
  
  umap_df$clone_label <- NA_character_
  umap_df$clone_label[in_cluster] <- get_clone_label(
    md$Global_Clone_ID[in_cluster], spec
  )
  umap_df$clone_label <- factor(umap_df$clone_label, levels = names(color_map))
  
  bg_df <- umap_df[!in_cluster, ]
  fg_df <- umap_df[ in_cluster, ]
  
  ggplot() +
    geom_point(
      data  = bg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]]),
      color = "lightgrey", size = 0.3
    ) +
    geom_point(
      data = fg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]],
          color = clone_label),
      size = 1.4
    ) +
    scale_color_manual(values = color_map, drop = FALSE) +
    ggtitle(paste0("Cluster ", cluster_id, "\n(n = ", n, ")")) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 7, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin     = margin(2, 2, 2, 2)
    )
}

plot_cluster_grid_colored <- function(seurat_obj, meta_col, clusters, spec,
                                      color_map, timepoint_filter = NULL,
                                      ncol = 5, main_title = NULL) {
  
  plots <- lapply(clusters, function(cl) {
    make_cluster_umap_colored(seurat_obj, meta_col, cl, spec,
                              color_map, timepoint_filter)
  })
  
  shared_legend <- cowplot::get_legend(make_legend_plot(color_map))
  panel_grid    <- wrap_plots(plots, ncol = ncol)
  
  if (!is.null(main_title)) {
    panel_grid <- panel_grid + plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(size = 12, face = "bold",
                                              hjust = 0.5))
    )
  }
  
  cowplot::plot_grid(panel_grid, shared_legend,
                     rel_widths = c(1, 0.12), nrow = 1)
}

# ── 6. Build and export UMAP plots ───────────────────────────────────────────

export_dir   <- "Donor_0/exports/TCR_clusters"
combined_dir <- file.path(export_dir, "all_timepoints_combined")
dir.create(export_dir,   showWarnings = FALSE, recursive = TRUE)
dir.create(combined_dir, showWarnings = FALSE, recursive = TRUE)

timepoint_groups <- list(
  "V1_V2_V3" = c("Visit 1", "Visit 2", "Visit 3"),
  "V4"        = c("Visit 4")
)

specs <- list(
  pp65 = list(n_top = 30, color_map = spec_color_maps$pp65),
  IE1  = list(n_top = 10, color_map = spec_color_maps$IE1)
)

for (tag in c("et48", "et36")) {
  
  meta_col <- paste0("TCR_cluster_", tag)
  
  for (spec_name in names(specs)) {
    
    n_top     <- specs[[spec_name]]$n_top
    color_map <- specs[[spec_name]]$color_map
    
    top_clusters <- get_top_clusters(D0_data_CMV, meta_col, spec_name, n_top)
    message(sprintf("[%s | %s] %d clusters identified",
                    tag, spec_name, length(top_clusters)))
    
    n_rows <- ceiling(length(top_clusters) / 5)
    grid_h <- n_rows * (22 / 5) + 1.5
    
    # ── All timepoints combined ────────────────────────────────────────────
    p_all <- plot_cluster_grid_colored(
      seurat_obj       = D0_data_CMV,
      meta_col         = meta_col,
      clusters         = top_clusters,
      spec             = spec_name,
      color_map        = color_map,
      timepoint_filter = NULL,
      ncol             = 5,
      main_title       = sprintf("Top %d %s-specific TCR clusters — %s (all visits)",
                                 n_top, spec_name, tag)
    )
    
    ggsave(
      filename = file.path(combined_dir,
                           sprintf("%s_top%d_clusters_%s_all_visits.png",
                                   spec_name, n_top, tag)),
      plot = p_all, width = 22, height = grid_h, dpi = 150, bg = "white"
    )
    message(sprintf("[%s | %s | all] Saved to %s", tag, spec_name, combined_dir))
    
    # ── Per-timepoint split ────────────────────────────────────────────────
    for (tp_label in names(timepoint_groups)) {
      
      tp_filter  <- timepoint_groups[[tp_label]]
      main_title <- sprintf("Top %d %s-specific TCR clusters — %s (%s)",
                            n_top, spec_name, tag, gsub("_", "/", tp_label))
      
      p <- plot_cluster_grid_colored(
        seurat_obj       = D0_data_CMV,
        meta_col         = meta_col,
        clusters         = top_clusters,
        spec             = spec_name,
        color_map        = color_map,
        timepoint_filter = tp_filter,
        ncol             = 5,
        main_title       = main_title
      )
      
      ggsave(
        filename = file.path(export_dir,
                             sprintf("%s_top%d_clusters_%s_%s.png",
                                     spec_name, n_top, tag, tp_label)),
        plot = p, width = 22, height = grid_h, dpi = 150, bg = "white"
      )
      message(sprintf("[%s | %s | %s] Saved to %s", tag, spec_name, tp_label, export_dir))
    }
  }
}



# ── 6. Export cell-level CSV files ──────────────────────────────────────────
#
# Four files: {pp65, IE1} × {et48, et36}
# One row per unique Global_Clone_ID, columns for each TCR cluster it appears in,
# values = cell counts. Rows ordered by descending total cells per clone.

for (tag in c("et48", "et36")) {
  
  meta_col <- paste0("TCR_cluster_", tag)
  
  for (spec in c("pp65", "IE1")) {
    
    df_out <- D0_data_CMV@meta.data %>%
      filter(specificity == spec,
             !is.na(.data[[meta_col]])) %>%
      select(Global_Clone_ID, TCR_cluster = all_of(meta_col)) %>%
      # Count cells per clone × cluster combination
      count(Global_Clone_ID, TCR_cluster, name = "n_cells") %>%
      # Pivot so each cluster becomes its own column
      tidyr::pivot_wider(
        names_from  = TCR_cluster,
        values_from = n_cells,
        values_fill = 0,
        names_prefix = "cluster_"
      ) %>%
      # Add a total column and sort by it
      mutate(total_cells = rowSums(across(starts_with("cluster_")))) %>%
      arrange(desc(total_cells))
    
    out_path <- file.path(
      export_dir,
      paste0(spec, "_clones_TCR_clusters_", tag, ".csv")
    )
    
    write_csv(df_out, out_path)
    message(sprintf("[%s | %s] %d clones written to %s",
                    tag, spec, nrow(df_out), out_path))
  }
}

# ── 7. Export cluster-level summary CSV files ────────────────────────────────
#
# Four files: {pp65, IE1} × {et48, et36}
# One row per unique TCR cluster, columns for each Global_Clone_ID,
# values = cell counts. Rows ordered by descending total cells per cluster.

for (tag in c("et48", "et36")) {
  
  meta_col <- paste0("TCR_cluster_", tag)
  
  for (spec in c("pp65", "IE1")) {
    
    df_out <- D0_data_CMV@meta.data %>%
      filter(specificity == spec,
             !is.na(.data[[meta_col]])) %>%
      select(Global_Clone_ID, TCR_cluster = all_of(meta_col)) %>%
      count(TCR_cluster, Global_Clone_ID, name = "n_cells") %>%
      tidyr::pivot_wider(
        names_from  = Global_Clone_ID,
        values_from = n_cells,
        values_fill = 0
      ) %>%
      mutate(total_cells = rowSums(across(-TCR_cluster))) %>%
      arrange(desc(total_cells))
    
    out_path <- file.path(
      export_dir,
      paste0(spec, "_TCR_clusters_by_clone_", tag, ".csv")
    )
    
    write_csv(df_out, out_path)
    message(sprintf("[%s | %s] %d clusters written to %s",
                    tag, spec, nrow(df_out), out_path))
  }
}




#
# TCR cluster analysis for donors 1-5
#

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(cowplot)
library(readr)

# ── 1. Load CSV ──────────────────────────────────────────────────────────────

df_et64 <- read_csv("/box/Campbell/Sam K/CMV_project/####_RNAseq/# Cross-experiment analyses/### Seurat analyses/cmv_d1-5_dataset_paired_clusters_et64_042826.csv")

# ── 2. TCR cluster assignment ────────────────────────────────────────────────

assign_tcr_cluster <- function(cdr3_vec, cluster_df) {
  
  count_cols <- c("IE1_Early", "Bulk_Early", "pp65_Early",
                  "Bulk_Late",  "IE1_Late",  "pp65_Late")
  
  cluster_df <- cluster_df[!is.na(cluster_df$cluster), ]
  cluster_df$total_count <- rowSums(cluster_df[, count_cols], na.rm = TRUE)
  
  sapply(cdr3_vec, function(x) {
    if (is.na(x)) return(NA_character_)
    
    parts    <- strsplit(x, ";")[[1]]
    trb_seqs <- sub("^TRB:", "", parts[grepl("^TRB:", parts)])
    tra_seqs <- sub("^TRA:", "", parts[grepl("^TRA:", parts)])
    
    if (length(trb_seqs) > 0) {
      hits <- cluster_df[cluster_df$cdr3_b_aa %in% trb_seqs, ]
      if (nrow(hits) == 0) return(NA_character_)
      if (length(tra_seqs) > 0) {
        hits_narrow <- hits[hits$cdr3_a_aa %in% tra_seqs, ]
        if (nrow(hits_narrow) > 0) hits <- hits_narrow
      }
    } else if (length(tra_seqs) > 0) {
      hits <- cluster_df[cluster_df$cdr3_a_aa %in% tra_seqs, ]
      if (nrow(hits) == 0) return(NA_character_)
    } else {
      return(NA_character_)
    }
    
    best <- hits[which.max(hits$total_count), ]
    return(as.character(best$cluster))
    
  }, USE.NAMES = FALSE)
}

mdata_CMV$TCR_cluster_et64 <- assign_tcr_cluster(mdata_CMV$t_cdr3s_aa, df_et64)

n_assigned <- sum(!is.na(mdata_CMV$TCR_cluster_et64))
n_total    <- ncol(mdata_CMV)
message(sprintf("[et64] Assigned: %d / %d cells (%.1f%%)",
                n_assigned, n_total, 100 * n_assigned / n_total))

# ── 3. Export directory ──────────────────────────────────────────────────────

export_dir <- "exports/TCR_clusters"
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

# ── 4. Export summary CSVs ───────────────────────────────────────────────────

for (spec in c("pp65", "IE1")) {
  
  base_df <- mdata_CMV@meta.data %>%
    filter(specificity == spec, !is.na(TCR_cluster_et64)) %>%
    select(Global_Clone_ID, TCR_cluster = TCR_cluster_et64)
  
  # Clone-level: one row per Global_Clone_ID, columns = TCR clusters
  clone_df <- base_df %>%
    count(Global_Clone_ID, TCR_cluster, name = "n_cells") %>%
    pivot_wider(names_from = TCR_cluster, values_from = n_cells,
                values_fill = 0, names_prefix = "cluster_") %>%
    mutate(total_cells = rowSums(across(starts_with("cluster_")))) %>%
    arrange(desc(total_cells))
  
  write_csv(clone_df,
            file.path(export_dir, paste0(spec, "_clones_TCR_clusters_et64.csv")))
  message(sprintf("[et64 | %s] Clone-level: %d clones written", spec, nrow(clone_df)))
  
  # Cluster-level: one row per TCR cluster, columns = Global_Clone_IDs
  cluster_df_out <- base_df %>%
    count(TCR_cluster, Global_Clone_ID, name = "n_cells") %>%
    pivot_wider(names_from = Global_Clone_ID, values_from = n_cells,
                values_fill = 0) %>%
    mutate(total_cells = rowSums(across(-TCR_cluster))) %>%
    arrange(desc(total_cells))
  
  write_csv(cluster_df_out,
            file.path(export_dir, paste0(spec, "_TCR_clusters_by_clone_et64.csv")))
  message(sprintf("[et64 | %s] Cluster-level: %d clusters written",
                  spec, nrow(cluster_df_out)))
}

# ── All-donors combined UMAP plots ───────────────────────────────────────────

# Donor-level clone labeling — all clones from same donor get same label
get_clone_label_by_donor <- function(global_clone_ids, spec) {
  pattern <- paste0("^(D\\d)_", spec, "-clone-\\d+$")
  ifelse(
    grepl(pattern, as.character(global_clone_ids), perl = TRUE),
    sub(pattern, "\\1", as.character(global_clone_ids), perl = TRUE),
    "Other"
  )
}

# Donor-level color map
donor_color_map <- c(
  "D1"    = "#1b9e77",
  "D2"    = "#d95f02",
  "D3"    = "#7570b3",
  "D4"    = "#e7298a",
  "D5"    = "#66a61e",
  "Other" = "grey70"
)

make_legend_plot <- function(color_map) {
  dummy_df <- data.frame(
    x     = seq_along(color_map),
    y     = 1,
    label = factor(names(color_map), levels = names(color_map))
  )
  ggplot(dummy_df, aes(x = x, y = y, color = label)) +
    geom_point(size = 3) +
    scale_color_manual(values = color_map, name = "Global Clone ID",
                       drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
    theme_void() +
    theme(legend.position = "right")
}

make_cluster_umap_colored <- function(seurat_obj, meta_col, cluster_id,
                                      spec, color_map,
                                      timepoint_filter = NULL) {
  
  md <- seurat_obj@meta.data
  
  if (!is.null(timepoint_filter)) {
    cells_keep <- rownames(md)[md$timepoint %in% timepoint_filter]
    seurat_obj <- subset(seurat_obj, cells = cells_keep)
    md         <- seurat_obj@meta.data
  }
  
  umap_df      <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  umap_cols    <- colnames(umap_df)
  umap_df$cell <- rownames(umap_df)
  
  in_cluster <- !is.na(md[[meta_col]]) & md[[meta_col]] == as.character(cluster_id)
  n          <- sum(in_cluster)
  
  umap_df$clone_label <- NA_character_
  umap_df$clone_label[in_cluster] <- get_clone_label_by_donor(
    md$Global_Clone_ID[in_cluster], spec
  )
  umap_df$clone_label <- factor(umap_df$clone_label, levels = names(color_map))
  
  bg_df <- umap_df[!in_cluster, ]
  fg_df <- umap_df[ in_cluster, ]
  
  ggplot() +
    geom_point(
      data  = bg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]]),
      color = "lightgrey", size = 0.3
    ) +
    geom_point(
      data = fg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]],
          color = clone_label),
      size = 1.4
    ) +
    scale_color_manual(values = color_map, drop = FALSE) +
    ggtitle(paste0("Cluster ", cluster_id, "\n(n = ", n, ")")) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 7, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin     = margin(2, 2, 2, 2)
    )
}

plot_cluster_grid_colored <- function(seurat_obj, meta_col, clusters, spec,
                                      color_map, timepoint_filter = NULL,
                                      ncol = 5, main_title = NULL) {
  
  plots <- lapply(clusters, function(cl) {
    make_cluster_umap_colored(seurat_obj, meta_col, cl, spec,
                              color_map, timepoint_filter)
  })
  
  shared_legend <- cowplot::get_legend(make_legend_plot(color_map))
  panel_grid    <- wrap_plots(plots, ncol = ncol)
  
  if (!is.null(main_title)) {
    panel_grid <- panel_grid + plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(size = 12, face = "bold",
                                              hjust = 0.5))
    )
  }
  
  cowplot::plot_grid(panel_grid, shared_legend,
                     rel_widths = c(1, 0.12), nrow = 1)
}

get_top_clusters <- function(seurat_obj, meta_col, target_spec, n_top) {
  seurat_obj@meta.data %>%
    filter(specificity == target_spec, !is.na(.data[[meta_col]])) %>%
    count(.data[[meta_col]], name = "n_cells") %>%
    arrange(desc(n_cells)) %>%
    slice_head(n = n_top) %>%
    pull(.data[[meta_col]])
}

# ── Generate plots ────────────────────────────────────────────────────────────

timepoint_groups <- list(
  "all"   = NULL,
  "Early" = "Early",
  "Late"  = "Late"
)

specs <- list(
  pp65 = list(n_top = 30, grid_h = 24),
  IE1  = list(n_top = 30, grid_h = 24)
)

for (spec_name in names(specs)) {
  
  n_top        <- specs[[spec_name]]$n_top
  top_clusters <- get_top_clusters(mdata_CMV, "TCR_cluster_et64", spec_name, n_top)
  
  message(sprintf("[et64 | %s] %d clusters identified", spec_name, length(top_clusters)))
  
  n_rows <- ceiling(length(top_clusters) / 5)
  grid_h <- n_rows * (22 / 5) + 1.5
  
  for (tp_label in names(timepoint_groups)) {
    
    tp_filter  <- timepoint_groups[[tp_label]]
    tp_display <- if (tp_label == "all") "All timepoints" else tp_label
    main_title <- sprintf("Top %d %s-specific TCR clusters — et64 (%s)",
                          n_top, spec_name, tp_display)
    
    p <- plot_cluster_grid_colored(
      seurat_obj       = mdata_CMV,
      meta_col         = "TCR_cluster_et64",
      clusters         = top_clusters,
      spec             = spec_name,
      color_map        = donor_color_map,
      timepoint_filter = tp_filter,
      ncol             = 5,
      main_title       = main_title
    )
    
    out_file <- file.path(
      export_dir,
      sprintf("%s_top%d_clusters_et64_%s.png", spec_name, n_top, tp_label)
    )
    
    ggsave(filename = out_file, plot = p,
           width = 22, height = grid_h, dpi = 150, bg = "white")
    
    message(sprintf("[et64 | %s | %s] Saved to %s", spec_name, tp_label, out_file))
  }
}


# ── By-donor UMAP plots ──────────────────────────────────────────────────────

donor_export_dir <- file.path(export_dir, "by_donor")
dir.create(donor_export_dir, showWarnings = FALSE, recursive = TRUE)

# Fixed qualitative palette — same across all donors, up to 10 clones + Other
clone_palette <- c(
  "Clone 1"  = "#1b9e77",
  "Clone 2"  = "#d95f02",
  "Clone 3"  = "#7570b3",
  "Clone 4"  = "#e7298a",
  "Clone 5"  = "#66a61e",
  "Clone 6"  = "#e6ab02",
  "Clone 7"  = "#a6761d",
  "Clone 8"  = "#386cb0",
  "Clone 9"  = "#f0027f",
  "Clone 10" = "#bf5b17",
  "Other"    = "grey70"
)

# Label clone IDs relative to a single donor
get_clone_label_donor <- function(global_clone_ids, donor, spec) {
  own_pattern <- paste0("^", donor, "_", spec, "-clone-(\\d+)$")
  any_clone_pattern <- "^D\\d_"  # matches any well-formed donor clone ID
  
  dplyr::case_when(
    grepl(own_pattern, as.character(global_clone_ids), perl = TRUE) ~
      paste0("Clone ", sub(own_pattern, "\\1",
                           as.character(global_clone_ids), perl = TRUE)),
    grepl(any_clone_pattern, as.character(global_clone_ids), perl = TRUE) ~
      "Other donor",
    TRUE ~ "No clone"
  )
}

# Legend plot using the fixed palette
make_legend_plot_donor <- function(color_map) {
  dummy_df <- data.frame(
    x     = seq_along(color_map),
    y     = 1,
    label = factor(names(color_map), levels = names(color_map))
  )
  ggplot(dummy_df, aes(x = x, y = y, color = label)) +
    geom_point(size = 3) +
    scale_color_manual(values = color_map, name = "Global Clone ID",
                       drop = FALSE) +
    guides(color = guide_legend(override.aes = list(size = 3), ncol = 1)) +
    theme_void() +
    theme(legend.position = "right")
}

# Single UMAP panel for one cluster, colored by clone within a single donor.
# Background = all cells across all donors (grey). Foreground = cells in this
# cluster; donor-matching clones get palette colors, all others get "Other".
make_cluster_umap_donor <- function(seurat_obj, meta_col, cluster_id,
                                    donor, spec, color_map) {
  
  md <- seurat_obj@meta.data
  
  umap_df      <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  umap_cols    <- colnames(umap_df)
  umap_df$cell <- rownames(umap_df)
  
  in_cluster <- !is.na(md[[meta_col]]) & md[[meta_col]] == as.character(cluster_id)
  n_donor    <- sum(in_cluster & grepl(paste0("^", donor, "_"), md$Global_Clone_ID))
  
  umap_df$clone_label <- NA_character_
  umap_df$clone_label[in_cluster] <- get_clone_label_donor(
    md$Global_Clone_ID[in_cluster], donor, spec
  )
  umap_df$clone_label <- factor(umap_df$clone_label, levels = names(color_map))
  
  bg_df <- umap_df[!in_cluster, ]
  fg_df <- umap_df[ in_cluster, ]
  
  ggplot() +
    geom_point(
      data  = bg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]]),
      color = "lightgrey", size = 0.3
    ) +
    geom_point(
      data = fg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]],
          color = clone_label),
      size = 1.4
    ) +
    scale_color_manual(values = color_map, drop = FALSE) +
    ggtitle(paste0("Cluster ", cluster_id, "\n(n = ", n_donor, ")")) +
    theme_void() +
    theme(
      plot.title      = element_text(size = 7, face = "bold", hjust = 0.5),
      legend.position = "none",
      plot.margin     = margin(2, 2, 2, 2)
    )
}

plot_cluster_grid_donor <- function(seurat_obj, meta_col, clusters,
                                    donor, spec, color_map,
                                    ncol = 5, main_title = NULL) {
  
  plots <- lapply(clusters, function(cl) {
    make_cluster_umap_donor(seurat_obj, meta_col, cl, donor, spec, color_map)
  })
  
  shared_legend <- cowplot::get_legend(make_legend_plot_donor(color_map))
  panel_grid    <- wrap_plots(plots, ncol = ncol)
  
  if (!is.null(main_title)) {
    panel_grid <- panel_grid + plot_annotation(
      title = main_title,
      theme = theme(plot.title = element_text(size = 12, face = "bold",
                                              hjust = 0.5))
    )
  }
  
  cowplot::plot_grid(panel_grid, shared_legend,
                     rel_widths = c(1, 0.12), nrow = 1)
}

# Rank top clusters using only cells from a specific donor
get_top_clusters_donor <- function(seurat_obj, meta_col, target_spec,
                                   donor, n_top) {
  seurat_obj@meta.data %>%
    filter(
      specificity == target_spec,
      !is.na(.data[[meta_col]]),
      grepl(paste0("^", .env$donor, "_"), as.character(Global_Clone_ID))
    ) %>%
    count(.data[[meta_col]], name = "n_cells") %>%
    arrange(desc(n_cells)) %>%
    slice_head(n = n_top) %>%
    pull(.data[[meta_col]])
}

# ── Generate plots ────────────────────────────────────────────────────────────

donors <- paste0("D", 1:5)

for (spec_name in c("pp65", "IE1")) {
  for (donor in donors) {
    
    top_clusters <- get_top_clusters_donor(
      mdata_CMV, "TCR_cluster_et64", spec_name, donor, n_top = 30
    )
    
    message(sprintf("[et64 | %s | %s] %d clusters identified",
                    spec_name, donor, length(top_clusters)))
    
    if (length(top_clusters) == 0) {
      message(sprintf("[et64 | %s | %s] No clusters found, skipping", spec_name, donor))
      next
    }
    
    # Trim palette to only clones that exist for this donor/spec,
    # but always keep "Other"
    pattern      <- paste0("^", donor, "_", spec_name, "-clone-(\\d+)$")
    present_nums <- mdata_CMV@meta.data %>%
      filter(grepl(paste0("^", .env$donor, "_", spec_name, "-clone-(\\d+)$"),
                   as.character(Global_Clone_ID), perl = TRUE)) %>%
      pull(Global_Clone_ID) %>%
      as.character() %>%
      sub(paste0("^", donor, "_", spec_name, "-clone-(\\d+)$"), "\\1", ., perl = TRUE) %>%
      as.integer() %>%
      unique() %>%
      sort()
    
    present_labels <- paste0("Clone ", present_nums)
    active_palette <- c(
      clone_palette[present_labels],   # donor's own clone colors
      "Other donor" = "#FFE900",       # bright yellow
      "No clone"    = "grey70"
    )
    
    main_title <- sprintf(
      "Top %d %s-specific TCR clusters — et64 (%s, all timepoints)",
      length(top_clusters), spec_name, donor
    )
    
    # Grid height scales with number of rows needed
    n_rows      <- ceiling(length(top_clusters) / 5)
    panel_width <- 22 / 5          # each panel's share of the total width
    grid_h      <- n_rows * panel_width + 1.5  # +1.5 for title headroom
    
    p <- plot_cluster_grid_donor(
      seurat_obj = mdata_CMV,
      meta_col   = "TCR_cluster_et64",
      clusters   = top_clusters,
      donor      = donor,
      spec       = spec_name,
      color_map  = active_palette,
      ncol       = 5,
      main_title = main_title
    )
    
    out_file <- file.path(
      donor_export_dir,
      sprintf("%s_%s_top%d_clusters_et64.png",
              spec_name, donor, length(top_clusters))
    )
    
    ggsave(filename = out_file, plot = p,
           width = 22, height = grid_h, dpi = 150, bg = "white")
    
    message(sprintf("[et64 | %s | %s] Saved to %s", spec_name, donor, out_file))
  }
}



# ── TCR cluster composition heatmap ─────────────────────────────────────────
#
# One plot per specificity (pp65, IE1)
# Rows    = top 30 clusters by total cell count (descending)
# Columns = Global_Clone_ID grouped by donor (D1→D5), plus D_minor-clones
#           per donor; clones ordered by clone number within donor
# Fill    = % of cluster total
# Label   = raw cell count (non-zero only)

library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)

heatmap_export_dir <- file.path(export_dir, "cluster_heatmaps")
dir.create(heatmap_export_dir, showWarnings = FALSE, recursive = TRUE)

# ── Helper: build the long-form data frame for one specificity ───────────────

build_heatmap_df <- function(seurat_obj, meta_col, target_spec, n_top = 30) {
  
  md <- seurat_obj@meta.data %>%
    filter(specificity == target_spec, !is.na(.data[[meta_col]])) %>%
    select(cluster   = all_of(meta_col),
           clone_id  = Global_Clone_ID,
           donor_raw = donor_identity)
  
  # Standardise donor labels to D1–D5
  md <- md %>%
    mutate(donor = str_extract(donor_raw, "\\d") %>% paste0("D", .))
  
  # Top N clusters by total cell count
  top_clusters <- md %>%
    count(cluster, name = "total_cells") %>%
    arrange(desc(total_cells)) %>%
    slice_head(n = n_top) %>%
    pull(cluster)
  
  md <- md %>% filter(cluster %in% top_clusters)
  
  # Classify each cell as a named clone or donor-specific minor-clones
  md <- md %>%
    mutate(
      col_id = case_when(
        !is.na(clone_id) ~ clone_id,
        TRUE             ~ paste0(donor, "_minor-clones")
      )
    )
  
  # Count cells per cluster × col_id
  counts <- md %>%
    count(cluster, col_id, name = "n_cells")
  
  # Cluster totals for % calculation
  cluster_totals <- md %>%
    count(cluster, name = "total_cells")
  
  counts <- counts %>%
    left_join(cluster_totals, by = "cluster") %>%
    mutate(pct = 100 * n_cells / total_cells)
  
  # Attach cluster total for row ordering
  counts <- counts %>%
    left_join(
      cluster_totals %>% rename(row_order_n = total_cells),
      by = "cluster"
    )
  
  list(long_df = counts, top_clusters = top_clusters,
       cluster_totals = cluster_totals)
}

# ── Helper: build ordered column levels ──────────────────────────────────────

build_col_order <- function(long_df, seurat_obj, target_spec) {
  
  donors <- paste0("D", 1:5)
  
  col_levels <- lapply(donors, function(d) {
    
    # Named clones: extract clone number, sort numerically
    pattern     <- paste0("^", d, "_", target_spec, "-clone-(\\d+)$")
    named_cols  <- long_df %>%
      filter(grepl(pattern, col_id, perl = TRUE)) %>%
      pull(col_id) %>%
      unique()
    
    clone_nums  <- as.integer(
      str_match(named_cols, paste0(d, "_", target_spec, "-clone-(\\d+)"))[, 2]
    )
    named_cols  <- named_cols[order(clone_nums)]
    
    # Minor-clones column for this donor (include even if absent in top clusters)
    minor_col   <- paste0(d, "_minor-clones")
    
    c(named_cols, minor_col)
  })
  
  unlist(col_levels)
}

# ── Helper: donor group lines & labels for x-axis annotation ─────────────────

build_donor_bands <- function(col_levels) {
  
  donors <- paste0("D", 1:5)
  
  band_df <- lapply(donors, function(d) {
    idxs <- which(str_starts(col_levels, d))
    if (length(idxs) == 0) return(NULL)
    data.frame(donor = d, xmin = min(idxs) - 0.5, xmax = max(idxs) + 0.5)
  }) %>%
    bind_rows()
  
  band_df$xmid <- (band_df$xmin + band_df$xmax) / 2
  band_df
}

# ── Helper: short x-axis labels (strip donor prefix) ─────────────────────────

make_short_labels <- function(col_levels, target_spec) {
  col_levels %>%
    str_remove(paste0("^D\\d_", target_spec, "-")) %>%  # remove "D1_pp65-"
    str_remove("^D\\d_")                                # remove "D1_" from minor
}

# ── Main plot function ────────────────────────────────────────────────────────

make_cluster_heatmap <- function(seurat_obj, meta_col, target_spec,
                                 n_top = 30) {
  
  res         <- build_heatmap_df(seurat_obj, meta_col, target_spec, n_top)
  long_df     <- res$long_df
  top_clusters <- res$top_clusters
  
  col_levels  <- build_col_order(long_df, seurat_obj, target_spec)
  band_df     <- build_donor_bands(col_levels)
  short_labs  <- make_short_labels(col_levels, target_spec)
  
  # Row order: cluster sorted by total cells descending
  row_levels  <- res$cluster_totals %>%
    filter(cluster %in% top_clusters) %>%
    arrange(desc(total_cells)) %>%
    pull(cluster)
  
  # Complete grid (fill missing combos with 0)
  plot_df <- long_df %>%
    complete(cluster = top_clusters, col_id = col_levels,
             fill = list(n_cells = 0, pct = 0)) %>%
    mutate(
      cluster = factor(cluster, levels = rev(row_levels)),  # top = highest row
      col_id  = factor(col_id,  levels = col_levels)
    )
  
  # Label: show count only when > 0
  plot_df <- plot_df %>%
    mutate(label = if_else(n_cells > 0, as.character(n_cells), ""))
  
  # Donor band background rectangles
  donor_colors <- c(D1 = "#e8f4f8", D2 = "#fef3e2", D3 = "#f0f4e8",
                    D4 = "#fce4ec", D5 = "#ede7f6")
  
  n_rows <- length(row_levels)
  
  p <- ggplot(plot_df,
              aes(x = col_id, y = cluster, fill = pct)) +
    
    # Donor band backgrounds
    geom_rect(
      data = band_df,
      aes(xmin = xmin, xmax = xmax, ymin = 0.5, ymax = n_rows + 0.5,
          fill = NULL),
      fill  = donor_colors[band_df$donor],
      alpha = 0.25,
      inherit.aes = FALSE
    ) +
    
    # Heatmap tiles
    geom_tile(color = "white", linewidth = 0.4) +
    
    # Cell count labels
    geom_text(aes(label = label), size = 2.2, color = "black") +
    
    # Fill scale
    scale_fill_gradientn(
      colours = c("white", "#ffffcc", "#fd8d3c", "#800026"),
      values  = rescale(c(0, 1, 10, max(plot_df$pct, na.rm = TRUE))),
      limits  = c(0, max(plot_df$pct, na.rm = TRUE)),
      name    = "% of cluster"
    ) +
    
    # Axes
    scale_x_discrete(labels = short_labs, expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    
    # Donor labels above the plot via band midpoints
    annotate(
      "text",
      x     = band_df$xmid,
      y     = n_rows + 1.1,
      label = band_df$donor,
      size  = 3.5, fontface = "bold"
    ) +
    
    labs(
      title = sprintf("TCR cluster composition — %s (et64, top %d clusters)",
                      target_spec, n_top),
      x = NULL, y = "TCR cluster"
    ) +
    
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y      = element_text(size = 8),
      panel.grid       = element_blank(),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
      legend.position  = "right",
      plot.margin      = margin(t = 30, r = 10, b = 10, l = 10)
    )
  
  p
}

# ── Generate and save ─────────────────────────────────────────────────────────

for (spec_name in c("pp65", "IE1")) {
  
  p <- make_cluster_heatmap(
    seurat_obj  = mdata_CMV,
    meta_col    = "TCR_cluster_et64",
    target_spec = spec_name,
    n_top       = 30
  )
  
  # Width scales with number of clone columns (~0.35 in per column + legend)
  n_cols   <- length(build_col_order(
    build_heatmap_df(mdata_CMV, "TCR_cluster_et64", spec_name)$long_df,
    mdata_CMV, spec_name
  ))
  
  plot_w <- max(16, n_cols * 0.38 + 3)
  plot_h <- 12
  
  out_file <- file.path(
    heatmap_export_dir,
    sprintf("%s_cluster_composition_heatmap_et64.png", spec_name)
  )
  
  ggsave(filename = out_file, plot = p,
         width = plot_w, height = plot_h, dpi = 150, bg = "white")
  
  message(sprintf("[et64 | %s] Heatmap saved to %s", spec_name, out_file))
}


# ── Cluster card plots ────────────────────────────────────────────────────────
#
# One PNG per cluster with >10 cells (specificity-agnostic)
# Layout (2×2):
#   top-left    : text summary — cluster ID, total cells, per-donor breakdown
#   top-right   : UMAP — background (grey) + foreground colored by
#                 Global_Clone_ID or donor_minor-clones group
#   bottom-left : TCR CDR3 logo (TRA over TRB)
#   bottom-right: metacluster stacked bar chart, one bar per clone/minor group
#
# Fresh qualitative palette assigned per plot for the UMAP groups.
# Fixed metacluster palette shared across all cards.

library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggseqlogo)

# Export directory for the donors 1-5 (et64) cards. Defined explicitly so the
# script does not depend on an `export_dir` left in the session by other code.
export_dir      <- "exports/TCR_clusters"
card_export_dir <- file.path(export_dir, "cluster_cards")
dir.create(card_export_dir, showWarnings = FALSE, recursive = TRUE)

# ── Palette helpers ───────────────────────────────────────────────────────────

# Qualitative colors to draw from for the per-plot UMAP group palette
QUAL_PALETTE <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#386cb0", "#f0027f", "#bf5b17",
  "#4daf4a", "#e41a1c", "#984ea3", "#ff7f00", "#a65628",
  "#f781bf", "#999999", "#8dd3c7", "#fb8072", "#80b1d3"
)

assign_palette <- function(groups) {
  n <- length(groups)
  colors <- QUAL_PALETTE[seq_len(min(n, length(QUAL_PALETTE)))]
  if (n > length(QUAL_PALETTE)) {
    colors <- colorRampPalette(QUAL_PALETTE)(n)
  }
  setNames(colors, groups)
}

# Fixed metacluster levels & palette (shared across all cards)
metacluster_levels <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")

metacluster_palette <- c(
  "Tcm"     = "#F8766D",
  "Tem"     = "#A3A500",
  "T-CMC1"  = "#00BF7D",
  "NK-like" = "#00B0F6",
  "ISG+"    = "#E76BF3"
)

# ── ggseqlogo "chemistry" color scheme (default for amino acids) ─────────────
# Used to build a standalone legend; logos themselves have legends suppressed.

AA_CHEMISTRY <- data.frame(
  letter = c("G","S","T","Y","C","N","Q","K","R","H","D","E","P","A","W","F","L","I","M","V"),
  group  = c(rep("Polar", 5), rep("Neutral", 2), rep("Basic", 3),
             rep("Acidic", 2), rep("Hydrophobic", 8)),
  col    = c(rep("#109648", 5), rep("#5E239D", 2), rep("#255C99", 3),
             rep("#D62839", 2), rep("#221E22", 8)),
  stringsAsFactors = FALSE
)

# Build & export a single legend PNG covering all amino acids present across
# every cluster's logos. `present_letters` is the set of AAs that actually
# appear in any TRA/TRB modal-length logo across all cards.
export_aa_legend <- function(present_letters, out_file) {

  leg_df <- AA_CHEMISTRY %>%
    filter(letter %in% present_letters)

  if (nrow(leg_df) == 0) return(invisible(NULL))

  # Group -> color mapping (one swatch per chemistry group), with the
  # member amino acids listed in each label
  group_levels <- c("Polar", "Neutral", "Basic", "Acidic", "Hydrophobic")

  group_summary <- leg_df %>%
    group_by(group) %>%
    summarise(
      col   = dplyr::first(col),
      aas   = paste(sort(letter), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    arrange(group) %>%
    mutate(label = sprintf("%s (%s)", group, aas))

  color_map <- setNames(group_summary$col, group_summary$label)

  dummy <- data.frame(
    x     = seq_len(nrow(group_summary)),
    y     = 1,
    label = factor(group_summary$label, levels = group_summary$label)
  )

  leg_plot <- ggplot(dummy, aes(x = x, y = y, fill = label)) +
    geom_tile() +
    scale_fill_manual(values = color_map, name = "Amino acid chemistry") +
    guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10),
      legend.key.size = unit(0.6, "cm")
    )

  legend_grob <- cowplot::get_legend(leg_plot)

  ggsave(filename = out_file,
         plot     = cowplot::ggdraw(legend_grob),
         width    = 3.5, height = 3, dpi = 150, bg = "white")

  message(sprintf("[et64] AA chemistry legend saved to %s", out_file))
}

# ── Build per-cluster summary data ───────────────────────────────────────────

build_cluster_summary <- function(seurat_obj, meta_col, min_cells = 10) {

  md <- seurat_obj@meta.data %>%
    filter(!is.na(.data[[meta_col]])) %>%
    transmute(
      cluster     = .data[[meta_col]],
      clone_id    = Global_Clone_ID,
      donor       = paste0("D", str_extract(donor_identity, "\\d")),
      spec        = specificity,
      cdr3        = t_cdr3s_aa,
      metacluster = metacluster
    )

  # Clusters with >10 cells
  eligible <- md %>%
    count(cluster, name = "total_cells") %>%
    filter(total_cells > min_cells) %>%
    arrange(desc(total_cells))

  list(md = md, eligible = eligible)
}

# ── Top-left: text panel ──────────────────────────────────────────────────────

make_text_panel <- function(cluster_id, cluster_md) {

  total <- nrow(cluster_md)

  lines <- c(
    sprintf("Cluster %s", cluster_id),
    sprintf("n = %d cells", total),
    ""
  )

  donors_present <- sort(unique(cluster_md$donor))

  for (d in donors_present) {
    d_md  <- cluster_md %>% filter(donor == d)
    d_n   <- nrow(d_md)

    lines <- c(lines, sprintf("%s  (n = %d)", d, d_n))

    # Named clones: extract spec + clone number, sort by spec then clone number
    named <- d_md %>%
      filter(!is.na(clone_id)) %>%
      count(clone_id, spec, name = "n_cells") %>%
      mutate(clone_num = as.integer(str_extract(clone_id, "\\d+$"))) %>%
      arrange(spec, clone_num)

    for (i in seq_len(nrow(named))) {
      # Strip donor prefix for brevity: D1_pp65-clone-3 -> pp65-clone-3
      short <- str_remove(named$clone_id[i], paste0("^", d, "_"))
      lines <- c(lines, sprintf("  %s:  %d", short, named$n_cells[i]))
    }

    # Minor clones (NA Global_Clone_ID)
    n_minor <- d_md %>% filter(is.na(clone_id)) %>% nrow()
    if (n_minor > 0) {
      lines <- c(lines, sprintf("  minor:  %d", n_minor))
    }

    lines <- c(lines, "")  # blank line between donors
  }

  text_block <- paste(lines, collapse = "\n")

  ggplot() +
    annotate(
      "text",
      x = 0, y = 1,
      label      = text_block,
      hjust      = 0, vjust = 1,
      size       = 3.2,
      family     = "mono",
      lineheight = 1.4
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(10, 0, 10, 10))
}

# ── Top-right: UMAP panel ─────────────────────────────────────────────────────

make_umap_panel <- function(seurat_obj, cluster_id, meta_col,
                             cluster_md, color_map) {

  umap_df      <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  umap_cols    <- colnames(umap_df)
  umap_df$cell <- rownames(umap_df)

  md_full <- seurat_obj@meta.data %>%
    mutate(
      donor    = paste0("D", str_extract(donor_identity, "\\d")),
      group_id = case_when(
        !is.na(Global_Clone_ID) ~ as.character(Global_Clone_ID),
        !is.na(.data[[meta_col]]) & .data[[meta_col]] == as.character(cluster_id) ~
          paste0(donor, "_minor-clones"),
        TRUE ~ NA_character_
      )
    )

  in_cluster <- !is.na(seurat_obj@meta.data[[meta_col]]) &
    seurat_obj@meta.data[[meta_col]] == as.character(cluster_id)

  umap_df$group_id <- NA_character_
  umap_df$group_id[in_cluster] <- md_full$group_id[in_cluster]
  umap_df$group_id <- factor(umap_df$group_id, levels = names(color_map))

  bg_df <- umap_df[!in_cluster, ]
  fg_df <- umap_df[ in_cluster, ]

  ggplot() +
    geom_point(
      data  = bg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]]),
      color = "lightgrey", size = 0.3
    ) +
    geom_point(
      data = fg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]],
          color = group_id),
      size = 1.4
    ) +
    scale_color_manual(
      values = color_map,
      name   = "Group",
      drop   = FALSE,
      guide  = guide_legend(
        override.aes = list(size = 3),
        ncol = 1
      )
    ) +
    theme_void() +
    theme(
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.key.size = unit(0.4, "cm"),
      plot.margin     = margin(5, 5, 5, 5)
    )
}

# ── Bottom-left: TCR CDR3 logo panel (TRA over TRB) ──────────────────────────
#
# Parses t_cdr3s_aa for the cluster's cells, splits into TRA / TRB sequence
# lists, and builds a two-row logo (TRA on top, TRB below). Each chain's logo
# uses only the modal (most common) CDR3 length so the logo is rectangular and
# interpretable; the length and sequence count are noted in the panel title.

extract_cdr3 <- function(cdr3_vec, chain_prefix) {
  # chain_prefix is "TRA:" or "TRB:"
  seqs <- unlist(lapply(cdr3_vec, function(x) {
    if (is.na(x)) return(character(0))
    parts <- strsplit(x, ";")[[1]]
    sub(paste0("^", chain_prefix), "",
        parts[grepl(paste0("^", chain_prefix), parts)])
  }))
  seqs[nchar(seqs) > 0]
}

# Keep only sequences of the modal (most common) length so the logo is
# rectangular and interpretable.
filter_modal_length <- function(seqs) {
  if (length(seqs) == 0) return(character(0))
  lens      <- nchar(seqs)
  modal_len <- as.integer(names(sort(table(lens), decreasing = TRUE))[1])
  seqs[lens == modal_len]
}

make_logo_panel <- function(cluster_md_cells_cdr3, weight_mode = "unique") {
  # cluster_md_cells_cdr3 = vector of t_cdr3s_aa strings for cells in cluster
  # weight_mode:
  #   "unique" - each distinct CDR3 contributes once (a dominant clone does not
  #              dominate the position weights)
  #   "cell"   - each cell contributes (logo weighted by clonal abundance)

  pull_seqs <- function(prefix) {
    seqs <- extract_cdr3(cluster_md_cells_cdr3, prefix)
    if (weight_mode == "unique") seqs <- unique(seqs)
    filter_modal_length(seqs)
  }

  tra_seqs <- pull_seqs("TRA:")
  trb_seqs <- pull_seqs("TRB:")

  # Title count label reflects the weighting mode
  count_label <- if (weight_mode == "unique") "unique" else "cells"

  build_one <- function(seqs, title) {
    if (length(seqs) < 1) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = sprintf("%s: no\nsequences", title),
                   size = 3) +
          theme_void()
      )
    }
    # Single sequence: use probability method so it renders as clean
    # full-height letters (bits would also be full height but implies diversity)
    logo_method <- if (length(seqs) == 1) "prob" else "bits"
    ggseqlogo(seqs, method = logo_method) +
      ggtitle(sprintf("%s  (%d %s, len = %d)",
                      title, length(seqs), count_label, nchar(seqs[1]))) +
      theme(
        plot.title      = element_text(size = 8, face = "bold"),
        axis.text       = element_text(size = 6),
        axis.title      = element_text(size = 7),
        legend.position = "none"
      )
  }

  tra_logo <- build_one(tra_seqs, "TRA")
  trb_logo <- build_one(trb_seqs, "TRB")

  # TRA on top, TRB below
  cowplot::plot_grid(tra_logo, trb_logo, ncol = 1)
}

# ── Bottom-right: metacluster stacked bar panel ──────────────────────────────
#
# One bar per group (named clone or donor minor-clones), y = % of that group's
# cells in each metacluster. Bar order matches the text panel: donor (D1→D5),
# then clone number within donor, minor-clones last per donor.

make_metacluster_bar_panel <- function(cluster_md) {

  # Assign each cell a group_id: named clone or "<donor>_minor-clones"
  md_grouped <- cluster_md %>%
    mutate(
      group_id = if_else(!is.na(clone_id),
                         as.character(clone_id),
                         paste0(donor, "_minor-clones"))
    )

  # Keep only groups with >5 cells in the cluster (counting all cells,
  # including those without a metacluster assignment)
  keep_groups <- md_grouped %>%
    count(group_id, name = "total_in_cluster") %>%
    filter(total_in_cluster > 5) %>%
    pull(group_id)

  bar_df <- md_grouped %>%
    filter(group_id %in% keep_groups, !is.na(metacluster)) %>%
    count(group_id, donor, clone_id, metacluster, name = "n") %>%
    group_by(group_id) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(metacluster = factor(metacluster, levels = metacluster_levels))

  # Build ordered group levels: donor first, then clone number, minor last
  group_order <- bar_df %>%
    distinct(group_id, donor, clone_id) %>%
    mutate(
      donor_num = as.integer(str_extract(donor, "\\d")),
      is_minor  = is.na(clone_id),
      clone_num = if_else(is_minor, Inf,
                          as.integer(str_extract(clone_id, "\\d+$"))) %>%
                    as.numeric()
    ) %>%
    arrange(donor_num, is_minor, clone_num) %>%
    pull(group_id)

  bar_df <- bar_df %>%
    mutate(group_id = factor(group_id, levels = group_order))

  # x labels: prepend donor to named clones (e.g. "D1 pp65-clone-3"),
  # and render minor-clone groups as "Dx minor"
  short_label <- function(g) {
    ifelse(
      grepl("_minor-clones$", g),
      str_replace(g, "_minor-clones$", " minor"),
      str_replace(g, "^(D\\d)_", "\\1 ")
    )
  }

  ggplot(bar_df, aes(x = group_id, y = percent, fill = metacluster)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = metacluster_palette, name = "Metacluster",
                      drop = FALSE) +
    scale_x_discrete(labels = short_label) +
    labs(x = NULL, y = "% of cells") +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y     = element_text(size = 7),
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.key.size = unit(0.35, "cm"),
      plot.margin     = margin(5, 5, 5, 5)
    )
}

# ── Assemble and save one card (2×2 layout) ──────────────────────────────────

save_cluster_card <- function(seurat_obj, cluster_id, meta_col,
                               cluster_md, export_dir,
                               weight_mode = "unique") {

  # Determine foreground groups in this cluster for the UMAP palette
  named_groups <- cluster_md %>%
    filter(!is.na(clone_id)) %>%
    pull(clone_id) %>% unique() %>% sort()

  minor_groups <- cluster_md %>%
    filter(is.na(clone_id)) %>%
    pull(donor) %>% unique() %>% sort() %>%
    paste0("_minor-clones")

  all_groups <- c(named_groups, minor_groups)
  color_map  <- assign_palette(all_groups)

  # Build the four panels
  text_panel <- make_text_panel(cluster_id, cluster_md)
  umap_panel <- make_umap_panel(seurat_obj, cluster_id, meta_col,
                                cluster_md, color_map)
  logo_panel <- make_logo_panel(cluster_md$cdr3, weight_mode = weight_mode)
  bar_panel  <- make_metacluster_bar_panel(cluster_md)

  # Top row: text + UMAP
  top_row <- plot_grid(
    text_panel, umap_panel,
    nrow = 1, rel_widths = c(0.3, 0.7)
  )

  # Bottom row: logo + bar chart
  bottom_row <- plot_grid(
    logo_panel, bar_panel,
    nrow = 1, rel_widths = c(0.45, 0.55)
  )

  # Stack the two rows (double height)
  card <- plot_grid(
    top_row, bottom_row,
    ncol = 1, rel_heights = c(0.55, 0.45)
  )

  out_file <- file.path(
    export_dir,
    sprintf("cluster_%s_card.png", cluster_id)
  )

  ggsave(filename = out_file, plot = card,
         width = 10, height = 10, dpi = 150, bg = "white")
}

# ── Run across all eligible clusters ─────────────────────────────────────────

orig_dir <- file.path(card_export_dir, "original_logo_weights")
dir.create(orig_dir, showWarnings = FALSE, recursive = TRUE)

res      <- build_cluster_summary(mdata_CMV, "TCR_cluster_et64", min_cells = 10)
md_all   <- res$md
eligible <- res$eligible

message(sprintf("[et64] %d clusters with >10 cells", nrow(eligible)))

# Accumulate the set of amino acids appearing in any logo (modal-length
# TRA/TRB sequences) so the standalone legend covers exactly what's shown.
present_aa <- character(0)

for (i in seq_len(nrow(eligible))) {

  cl_id <- eligible$cluster[i]
  cl_md <- md_all %>% filter(cluster == cl_id)

  # Track amino acids that will appear in this card's logos. The modal-length
  # AA set is identical for unique- and cell-weighted modes, so tracking on
  # the deduplicated modal-length sequences covers both card versions.
  tra_seqs <- filter_modal_length(unique(extract_cdr3(cl_md$cdr3, "TRA:")))
  trb_seqs <- filter_modal_length(unique(extract_cdr3(cl_md$cdr3, "TRB:")))
  if (length(tra_seqs) >= 1) {
    present_aa <- union(present_aa, unlist(strsplit(tra_seqs, "")))
  }
  if (length(trb_seqs) >= 1) {
    present_aa <- union(present_aa, unlist(strsplit(trb_seqs, "")))
  }

  # Unique-weighted card (each distinct CDR3 weighted equally) -> main folder
  save_cluster_card(
    seurat_obj  = mdata_CMV,
    cluster_id  = cl_id,
    meta_col    = "TCR_cluster_et64",
    cluster_md  = cl_md,
    export_dir  = card_export_dir,
    weight_mode = "unique"
  )

  # Original cell-weighted card (logo weighted by clonal abundance)
  # -> original_logo_weights/ subfolder
  save_cluster_card(
    seurat_obj  = mdata_CMV,
    cluster_id  = cl_id,
    meta_col    = "TCR_cluster_et64",
    cluster_md  = cl_md,
    export_dir  = orig_dir,
    weight_mode = "cell"
  )

  if (i %% 10 == 0) {
    message(sprintf("[et64] %d / %d clusters done (both versions)",
                    i, nrow(eligible)))
  }
}

# Export the standalone AA chemistry legend into both folders (identical key)
export_aa_legend(
  present_letters = present_aa,
  out_file        = file.path(card_export_dir, "aa_chemistry_legend.png")
)
export_aa_legend(
  present_letters = present_aa,
  out_file        = file.path(orig_dir, "aa_chemistry_legend.png")
)

message(sprintf("[et64] All cluster cards saved to %s (+ original_logo_weights/)",
                card_export_dir))



#
#
# Donor 0 - Violin plots for gene set scores – Visit 4 cells
#
#

library(ggplot2)
library(dplyr)
library(tidyr)

D0_CMV_V4 <- subset(D0_data_CMV, subset = timepoint == "Visit 4")

# ── 1. Define the 5 custom score columns (exclude Standard_* sets) ───────────

score_cols_custom <- c("Tem_score", "Tcm_score", "T-CMC1_score",
                       "NK-like_score", "ISG+_score")

# ── 2. Helper: pivot meta.data to long format ─────────────────────────────────

pivot_scores <- function(obj, group_var) {
  obj@meta.data %>%
    select(all_of(c(group_var, score_cols_custom))) %>%
    pivot_longer(
      cols      = all_of(score_cols_custom),
      names_to  = "Score",
      values_to = "Value"
    ) %>%
    mutate(Score = factor(Score, levels = score_cols_custom),
           Group = as.factor(.data[[group_var]]))
}

# ── 3. Shared violin theme ────────────────────────────────────────────────────

violin_theme <- theme_bw() +
  theme(
    strip.text         = element_text(size = 11, face = "bold"),
    axis.text          = element_text(size = 10),
    axis.title         = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    legend.position    = "none"
  )

# ── 4. Plot 1: scores by cluster (Visit 4) ───────────────────────────────────

long_clusters <- pivot_scores(D0_CMV_V4, "seurat_clusters") %>%
  mutate(Group = factor(paste0("Cluster ", Group),
                        levels = paste0("Cluster ", sort(unique(D0_CMV_V4$seurat_clusters)))))

p_clusters <- ggplot(long_clusters, aes(x = Score, y = Value, fill = Score)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white", alpha = 0.6) +
  facet_wrap(~ Group, ncol = 1) +
  labs(
    title = "Gene Set Scores by Cluster (Visit 4)",
    x     = NULL,
    y     = "Module Score"
  ) +
  violin_theme

ggsave("Donor_0/exports/D0_CMV_Violin_clusters_V4.png",
       plot   = p_clusters,
       width  = 7, height = 3 * n_distinct(long_clusters$Group), dpi = 300)

# ── 5. Plot 2: scores by specificity (Visit 4) ───────────────────────────────

long_specificity <- pivot_scores(D0_CMV_V4, "specificity")

p_specificity <- ggplot(long_specificity, aes(x = Score, y = Value, fill = Score)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white", alpha = 0.6) +
  facet_wrap(~ Group, ncol = 1) +
  labs(
    title = "Gene Set Scores by Specificity (Visit 4)",
    x     = NULL,
    y     = "Module Score"
  ) +
  violin_theme

ggsave("Donor_0/exports/D0_CMV_Violin_specificity_V4.png",
       plot   = p_specificity,
       width  = 7, height = 3 * n_distinct(long_specificity$Group), dpi = 300)


# ── 4. Plot 1: scores by cluster (Visit 4) – transposed ──────────────────────

p_clusters_t <- ggplot(long_clusters, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white", alpha = 0.6) +
  facet_wrap(~ Score, ncol = 1) +
  labs(
    title = "Gene Set Scores by Cluster (Visit 4)",
    x     = NULL,
    y     = "Module Score"
  ) +
  violin_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Donor_0/exports/D0_CMV_Violin_clusters_V4_transposed.png",
       plot   = p_clusters_t,
       width  = 3 * n_distinct(long_clusters$Group), height = 16, dpi = 300)

# ── 5. Plot 2: scores by specificity (Visit 4) – transposed ──────────────────

p_specificity_t <- ggplot(long_specificity, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white", alpha = 0.6) +
  facet_wrap(~ Score, ncol = 1) +
  labs(
    title = "Gene Set Scores by Specificity (Visit 4)",
    x     = NULL,
    y     = "Module Score"
  ) +
  violin_theme

ggsave("Donor_0/exports/D0_CMV_Violin_specificity_V4_transposed.png",
       plot   = p_specificity_t,
       width  = 3 * n_distinct(long_specificity$Group), height = 16, dpi = 300)


#
#
# Donor 0 - Violin plots for gene set scores – Visits 1–3 cells (transposed only)
#
#

D0_CMV_V1to3 <- subset(D0_data_CMV, subset = timepoint %in% c("Visit 1", "Visit 2", "Visit 3"))

# ── 1. Define the 5 custom score columns (exclude Standard_* sets) ───────────

score_cols_custom <- c("Tem_score", "Tcm_score", "T-CMC1_score",
                       "NK-like_score", "ISG+_score")

# ── 2. Helper: pivot meta.data to long format ─────────────────────────────────

pivot_scores <- function(obj, group_var) {
  obj@meta.data %>%
    select(all_of(c(group_var, score_cols_custom))) %>%
    pivot_longer(
      cols      = all_of(score_cols_custom),
      names_to  = "Score",
      values_to = "Value"
    ) %>%
    mutate(Score = factor(Score, levels = score_cols_custom),
           Group = as.factor(.data[[group_var]]))
}

# ── 3. Shared violin theme ────────────────────────────────────────────────────

violin_theme <- theme_bw() +
  theme(
    strip.text         = element_text(size = 11, face = "bold"),
    axis.text          = element_text(size = 10),
    axis.title         = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    legend.position    = "none"
  )

# ── 4. Build long-format tables ──────────────────────────────────────────────

long_clusters <- pivot_scores(D0_CMV_V1to3, "seurat_clusters") %>%
  mutate(Group = factor(paste0("Cluster ", Group),
                        levels = paste0("Cluster ", sort(unique(D0_CMV_V1to3$seurat_clusters)))))

long_specificity <- pivot_scores(D0_CMV_V1to3, "specificity")

# ── 5. Plot 1: scores by cluster (Visits 1–3) – transposed ───────────────────

p_clusters_t <- ggplot(long_clusters, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white", alpha = 0.6) +
  facet_wrap(~ Score, ncol = 1) +
  labs(
    title = "Gene Set Scores by Cluster (Visits 1–3)",
    x     = NULL,
    y     = "Module Score"
  ) +
  violin_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Donor_0/exports/D0_CMV_Violin_clusters_V1to3_transposed.png",
       plot   = p_clusters_t,
       width  = 3 * n_distinct(long_clusters$Group), height = 16, dpi = 300)

# ── 6. Plot 2: scores by specificity (Visits 1–3) – transposed ───────────────

p_specificity_t <- ggplot(long_specificity, aes(x = Group, y = Value, fill = Group)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white", alpha = 0.6) +
  facet_wrap(~ Score, ncol = 1) +
  labs(
    title = "Gene Set Scores by Specificity (Visits 1–3)",
    x     = NULL,
    y     = "Module Score"
  ) +
  violin_theme

ggsave("Donor_0/exports/D0_CMV_Violin_specificity_V1to3_transposed.png",
       plot   = p_specificity_t,
       width  = 3 * n_distinct(long_specificity$Group), height = 16, dpi = 300)


#
#
# Donor 0 - Violin plots: pp65 vs IE1, split by visit group, colored by specificity
#
#

# Subset to CMV-specific cells (pp65 or IE1)
D0_CMV_all <- subset(D0_data_CMV, subset = specificity %in% c("pp65", "IE1"))

# ── 1. Define the 5 custom score columns ─────────────────────────────────────

score_cols_custom <- c("Tem_score", "Tcm_score", "T-CMC1_score",
                       "NK-like_score", "ISG+_score")

# ── 2. Build long-format table, keeping specificity for fill ─────────────────

long_df <- D0_CMV_all@meta.data %>%
  mutate(
    visit_group = ifelse(timepoint == "Visit 4", "V4", "V1-3"),
    specificity = factor(specificity, levels = c("pp65", "IE1")),
    Group       = factor(paste(specificity, visit_group),
                         levels = c("pp65 V1-3", "pp65 V4",
                                    "IE1 V1-3",  "IE1 V4"))
  ) %>%
  select(Group, specificity, all_of(score_cols_custom)) %>%
  pivot_longer(
    cols      = all_of(score_cols_custom),
    names_to  = "Score",
    values_to = "Value"
  ) %>%
  mutate(Score = factor(Score, levels = score_cols_custom))

# ── 3. Specificity color palette ─────────────────────────────────────────────

specificity_colors <- c("pp65" = "#1f77b4", "IE1" = "#d62728")

# ── 4. Shared violin theme ───────────────────────────────────────────────────

violin_theme <- theme_bw() +
  theme(
    strip.text         = element_text(size = 11, face = "bold"),
    axis.text          = element_text(size = 10),
    axis.title         = element_text(size = 11),
    panel.grid.major.x = element_blank(),
    legend.position    = "none"
  )

# ── 5. Plot: 4 violins per panel, colored by specificity, dashed separator ──

p <- ggplot(long_df, aes(x = Group, y = Value, fill = specificity)) +
  geom_violin(scale = "width", trim = FALSE, alpha = 0.75) +
  geom_boxplot(width = 0.12, outlier.size = 0.5, fill = "white", alpha = 0.6) +
  geom_vline(xintercept = 2.5, linetype = "dashed", color = "grey40") +
  facet_wrap(~ Score, ncol = 1) +
  scale_fill_manual(values = specificity_colors) +
  labs(
    title = "Gene Set Scores by Specificity and Visit Group",
    x     = NULL,
    y     = "Module Score"
  ) +
  violin_theme +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("Donor_0/exports/D0_CMV_Violin_specificity_visitgroup_transposed.png",
       plot   = p,
       width  = 8, height = 16, dpi = 300)



#
# Export D0 UMAPs for cluster ID comparison
#

cluster_colors <- c(
  "0" = "#C49A00",
  "1" = "#53B400",
  "2" = "#8B0000",
  "3" = "#1976D2",
  "4" = "#F8766D",
  "5" = "#A58AFF",
  "6" = "#FFDB58",
  "7" = "#00C094"
)

p1 <- DimPlot(
  D0_data_CMV,
  group.by  = "seurat_clusters",
  reduction = "umap",
  label     = TRUE,
  pt.size   = 0.5
) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("D0 CMV – All Visits")

p2 <- DimPlot(
  subset(D0_data_CMV, subset = timepoint == "Visit 4"),
  group.by  = "seurat_clusters",
  reduction = "umap",
  label     = TRUE,
  pt.size   = 0.5
) +
  scale_color_manual(values = cluster_colors) +
  ggtitle("D0 CMV – Visit 4 Only")

ggsave(
  filename = "Donor_0/exports/D0_CMV_seurat_clusters_allVisits_UMAP.jpeg",
  plot     = p1,
  width    = 8,
  height   = 7,
  dpi      = 300,
  device   = "jpeg"
)

ggsave(
  filename = "Donor_0/exports/D0_CMV_seurat_clusters_Visit4_UMAP.jpeg",
  plot     = p2,
  width    = 8,
  height   = 7,
  dpi      = 300,
  device   = "jpeg"
)



#
#
# Gini-Simpson index diversity of D0 clones
#
#

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(readr)

# ── 0. Settings ─────────────────────────────────────────────────────────────
min_cells   <- 5  # exclude clone × timepoint combos with fewer cells (GS unreliable)
export_dir  <- "Donor_0/exports/heterogeneity"
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

# Assumes clone_colors_fixed is already defined in your environment
# (20 "Clone N" entries + "Other (<1%)" → grey95)

# ── 1. Gini-Simpson index per clone × timepoint ─────────────────────────────

gini_simpson <- function(clusters) {
  n <- length(clusters)
  if (n < 1) return(NA_real_)
  p <- as.numeric(table(clusters)) / n
  1 - sum(p^2)
}

md <- D0_data_CMV@meta.data %>%
  filter(!is.na(Global_Clone_ID),
         !is.na(seurat_clusters),
         !is.na(timepoint),
         specificity %in% c("pp65", "IE1")) %>%
  mutate(timepoint = factor(timepoint,
                            levels = c("Visit 1", "Visit 2", "Visit 3", "Visit 4")))

gs_df <- md %>%
  group_by(specificity, Global_Clone_ID, timepoint) %>%
  summarise(
    n_cells       = n(),
    gini_simpson  = gini_simpson(seurat_clusters),
    .groups       = "drop"
  ) %>%
  filter(n_cells >= min_cells)

# ── 2. Map Global_Clone_ID → "Clone N" / "Other (<1%)" ──────────────────────

get_clone_label <- function(global_clone_ids, spec) {
  pattern <- paste0("^D0_", spec, "_clone_(\\d+)$")
  ifelse(
    grepl(pattern, global_clone_ids, perl = TRUE),
    paste0("Clone ", sub(pattern, "\\1", global_clone_ids, perl = TRUE)),
    "Other (<1%)"
  )
}

gs_df <- gs_df %>%
  mutate(
    clone_label = case_when(
      specificity == "pp65" & grepl("^D0_pp65_clone_\\d+$", Global_Clone_ID) ~
        paste0("Clone ", sub("^D0_pp65_clone_(\\d+)$", "\\1", Global_Clone_ID)),
      specificity == "IE1"  & grepl("^D0_IE1_clone_\\d+$",  Global_Clone_ID) ~
        paste0("Clone ", sub("^D0_IE1_clone_(\\d+)$",  "\\1", Global_Clone_ID)),
      TRUE ~ "Other (<1%)"
    ),
    clone_label = factor(clone_label, levels = names(clone_colors_fixed))
  )

# Export the underlying data
write_csv(gs_df, file.path(export_dir,
                           "D0_CMV_gini_simpson_by_clone_timepoint.csv"))

# ── 3. Plot ─────────────────────────────────────────────────────────────────

make_gs_plot <- function(df, spec, color_map) {
  df_sub <- df %>% filter(specificity == spec)
  
  ggplot(df_sub,
         aes(x      = timepoint,
             y      = gini_simpson,
             group  = Global_Clone_ID,
             color  = clone_label)) +
    geom_line(linewidth = 0.7, alpha = 0.85) +
    geom_point(size = 2) +
    scale_color_manual(values = color_map, drop = TRUE) +
    scale_y_continuous(
      limits = c(0, 1 - 1/8 + 0.02),
      breaks = seq(0, 0.8, by = 0.2)
    ) +
    labs(
      title    = paste0(spec, "-specific clones"),
      subtitle = "Gini-Simpson index over 8 Seurat clusters",
      x        = NULL,
      y        = "Gini-Simpson index",
      color    = NULL
    ) +
    theme_classic(base_size = 12) +
    theme(
      plot.title      = element_text(face = "bold"),
      plot.subtitle   = element_text(size = 10, color = "grey30"),
      legend.position = "right"
    )
}

p_pp65 <- make_gs_plot(gs_df, "pp65", clone_colors_fixed)
p_ie1  <- make_gs_plot(gs_df, "IE1",  clone_colors_fixed)

# ── 4. Export ───────────────────────────────────────────────────────────────

ggsave(file.path(export_dir, "D0_CMV_gini_simpson_pp65.jpeg"),
       plot = p_pp65, width = 9,  height = 6, dpi = 300, device = "jpeg")

ggsave(file.path(export_dir, "D0_CMV_gini_simpson_IE1.jpeg"),
       plot = p_ie1,  width = 9,  height = 6, dpi = 300, device = "jpeg")

ggsave(file.path(export_dir, "D0_CMV_gini_simpson_combined.jpeg"),
       plot = p_pp65 + p_ie1,
       width = 16, height = 6, dpi = 300, device = "jpeg")


# ── 5. Clonal repertoire heterogeneity per specificity × timepoint ──────────
# Named clones (Clone 1…N) share a bin; each "Other (<1%)" cell is its own bin

cell_df <- md %>%
  mutate(
    clone_label = case_when(
      specificity == "pp65" & grepl("^D0_pp65_clone_\\d+$", Global_Clone_ID) ~
        paste0("Clone ", sub("^D0_pp65_clone_(\\d+)$", "\\1", Global_Clone_ID)),
      specificity == "IE1"  & grepl("^D0_IE1_clone_\\d+$",  Global_Clone_ID) ~
        paste0("Clone ", sub("^D0_IE1_clone_(\\d+)$",  "\\1", Global_Clone_ID)),
      TRUE ~ "Other (<1%)"
    ),
    bin_id = if_else(clone_label == "Other (<1%)",
                     paste0("Other_", row_number()),
                     clone_label)
  )

pool_gs <- cell_df %>%
  group_by(specificity, timepoint) %>%
  summarise(
    n_cells      = n(),
    n_named      = sum(clone_label != "Other (<1%)"),
    n_other      = sum(clone_label == "Other (<1%)"),
    gini_simpson = gini_simpson(bin_id),
    .groups      = "drop"
  )

write_csv(pool_gs, file.path(export_dir,
                             "D0_CMV_gini_simpson_pool_by_timepoint.csv"))

# ── 6. Plot ─────────────────────────────────────────────────────────────────

pool_colors <- c("pp65" = "#1f78b4", "IE1" = "#e31a1c")

p_pool <- ggplot(pool_gs,
                 aes(x = timepoint, y = gini_simpson,
                     group = specificity, color = specificity)) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  scale_color_manual(values = pool_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title    = "Clonal repertoire heterogeneity",
    subtitle = "Gini-Simpson index over clone bins (each 'Other' cell = own bin)",
    x        = NULL,
    y        = "Gini-Simpson index",
    color    = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(size = 10, color = "grey30"),
    legend.position = "right"
  )

ggsave(file.path(export_dir, "D0_CMV_gini_simpson_pool.jpeg"),
       plot = p_pool, width = 7, height = 5, dpi = 300, device = "jpeg")


# ── 7. Bootstrap 95% CIs for pool-level GS ──────────────────────────────────

library(boot)

set.seed(42)
R_boot <- 2000

boot_gs <- function(data, indices) {
  gini_simpson(data$bin_id[indices])
}

pool_gs_ci <- cell_df %>%
  group_by(specificity, timepoint) %>%
  group_modify(~ {
    b  <- boot(data = .x, statistic = boot_gs, R = R_boot)
    ci <- boot.ci(b, type = "perc")$percent[4:5]
    tibble(
      n_cells      = nrow(.x),
      gini_simpson = b$t0,
      ci_lo        = ci[1],
      ci_hi        = ci[2]
    )
  }) %>%
  ungroup()

write_csv(pool_gs_ci, file.path(export_dir,
                                "D0_CMV_gini_simpson_pool_with_CIs.csv"))

# ── 8. Plot with bootstrap CIs ──────────────────────────────────────────────

dodge <- position_dodge(width = 0.15)

p_pool_ci <- ggplot(pool_gs_ci,
                    aes(x     = timepoint,
                        y     = gini_simpson,
                        group = specificity,
                        color = specificity)) +
  geom_line(linewidth = 1, position = dodge) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.15, linewidth = 0.7, position = dodge) +
  geom_point(size = 3, position = dodge) +
  scale_color_manual(values = pool_colors) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title    = "Clonal repertoire heterogeneity",
    subtitle = paste0("Gini-Simpson index with bootstrap 95% CIs (R = ",
                      R_boot, ")"),
    x        = NULL,
    y        = "Gini-Simpson index",
    color    = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(size = 10, color = "grey30"),
    legend.position = "right"
  )

ggsave(file.path(export_dir, "D0_CMV_gini_simpson_pool_CI.jpeg"),
       plot = p_pool_ci, width = 7, height = 5, dpi = 300, device = "jpeg")



#
#
# Pool-level Gini-Simpson diversity: bootstrap CIs + pairwise permutation tests
#
#

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(tibble)
library(boot)

# ── 0. Settings ─────────────────────────────────────────────────────────────

export_dir <- "Donor_0/exports/heterogeneity"
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
R_boot          <- 2000      # bootstrap replicates for CIs
R_perm          <- 10000     # permutation replicates per pairwise test
alpha           <- 0.05
p_adjust_method <- "holm"    # within each specificity (6 comparisons)
adjacent_only   <- TRUE      # only bracket consecutive visits (V1-V2, V2-V3, V3-V4)

tp_levels   <- c("Visit 1", "Visit 2", "Visit 3", "Visit 4")
pool_colors <- c("pp65" = "#1f78b4", "IE1" = "#e31a1c")
dodge_shift <- 0.0375        # horizontal offset per specificity
tip_len     <- 0.018         # bracket tip height

# ── 1. Helpers ──────────────────────────────────────────────────────────────

gini_simpson <- function(clusters) {
  n <- length(clusters)
  if (n < 1) return(NA_real_)
  p <- as.numeric(table(clusters)) / n
  1 - sum(p^2)
}

# fast version on pre-coded integers (used inside the permutation loop)
gs_codes <- function(codes, nbins) {
  n <- length(codes)
  if (n < 1) return(NA_real_)
  1 - sum((tabulate(codes, nbins = nbins) / n)^2)
}

perm_test_gs <- function(bins_a, bins_b, R = R_perm) {
  f     <- factor(c(bins_a, bins_b))
  codes <- as.integer(f)
  nbins <- nlevels(f)
  na    <- length(bins_a)
  n     <- length(codes)

  gs_a <- gs_codes(codes[seq_len(na)],  nbins)
  gs_b <- gs_codes(codes[-seq_len(na)], nbins)
  obs  <- gs_a - gs_b

  null <- vapply(seq_len(R), function(i) {
    idx <- sample.int(n, na)
    gs_codes(codes[idx], nbins) - gs_codes(codes[-idx], nbins)
  }, numeric(1))

  tibble(
    gs_1     = gs_a,
    gs_2     = gs_b,
    delta_gs = obs,
    p_value  = (1 + sum(abs(null) >= abs(obs) - 1e-12)) / (R + 1)
  )
}

# ── 2. Cell-level bins (named clones share a bin; each "Other" cell own bin) ─

md <- D0_data_CMV@meta.data %>%
  filter(!is.na(Global_Clone_ID),
         !is.na(seurat_clusters),
         !is.na(timepoint),
         specificity %in% c("pp65", "IE1")) %>%
  mutate(timepoint = factor(timepoint, levels = tp_levels))

cell_df <- md %>%
  mutate(
    clone_label = case_when(
      specificity == "pp65" & grepl("^D0_pp65_clone_\\d+$", Global_Clone_ID) ~
        paste0("Clone ", sub("^D0_pp65_clone_(\\d+)$", "\\1", Global_Clone_ID)),
      specificity == "IE1"  & grepl("^D0_IE1_clone_\\d+$",  Global_Clone_ID) ~
        paste0("Clone ", sub("^D0_IE1_clone_(\\d+)$",  "\\1", Global_Clone_ID)),
      TRUE ~ "Other (<1%)"
    ),
    bin_id = if_else(clone_label == "Other (<1%)",
                     paste0("Other_", row_number()),
                     paste0(specificity, "_", clone_label))
  )

# ── 3. Point estimates + bootstrap 95% CIs ──────────────────────────────────

boot_gs <- function(data, indices) gini_simpson(data$bin_id[indices])

pool_gs_ci <- cell_df %>%
  group_by(specificity, timepoint) %>%
  group_modify(~ {
    b  <- boot(data = .x, statistic = boot_gs, R = R_boot)
    ci <- boot.ci(b, type = "perc")$percent[4:5]
    tibble(
      n_cells      = nrow(.x),
      n_named      = sum(.x$clone_label != "Other (<1%)"),
      n_other      = sum(.x$clone_label == "Other (<1%)"),
      gini_simpson = b$t0,
      ci_lo        = ci[1],
      ci_hi        = ci[2]
    )
  }) %>%
  ungroup() %>%
  mutate(
    tp_idx = match(as.character(timepoint), tp_levels),
    shift  = if_else(specificity == "pp65", dodge_shift, -dodge_shift),
    x_pos  = tp_idx + shift
  )

write_csv(pool_gs_ci %>% select(-tp_idx, -shift, -x_pos),
          file.path(export_dir, "D0_CMV_gini_simpson_pool_with_CIs.csv"))

# ── 4. Pairwise permutation tests within each specificity ───────────────────

pairwise_gs <- lapply(sort(unique(as.character(cell_df$specificity))), function(spec) {
  d   <- cell_df %>% filter(specificity == spec)
  tps <- tp_levels[tp_levels %in% unique(as.character(d$timepoint))]
  prs <- combn(tps, 2, simplify = FALSE)

  lapply(prs, function(pr) {
    a <- d$bin_id[as.character(d$timepoint) == pr[1]]
    b <- d$bin_id[as.character(d$timepoint) == pr[2]]
    perm_test_gs(a, b) %>%
      mutate(specificity = spec,
             group1 = pr[1], group2 = pr[2],
             n1 = length(a), n2 = length(b),
             .before = 1)
  }) %>% bind_rows()
}) %>% bind_rows()

pairwise_gs <- pairwise_gs %>%
  group_by(specificity) %>%
  mutate(p_adj = p.adjust(p_value, method = p_adjust_method)) %>%
  ungroup() %>%
  mutate(
    signif = case_when(
      p_adj <= 0.0001 ~ "****",
      p_adj <= 0.001  ~ "***",
      p_adj <= 0.01   ~ "**",
      p_adj <= alpha  ~ "*",
      TRUE            ~ "ns"
    )
  ) %>%
  arrange(specificity, group1, group2)

print(as.data.frame(pairwise_gs), digits = 3)

write_csv(pairwise_gs,
          file.path(export_dir, "D0_CMV_gini_simpson_pool_pairwise_tests.csv"))

# ── 5. Build significance brackets ──────────────────────────────────────────

sig_pairs <- pairwise_gs %>%
  filter(signif != "ns") %>%
  mutate(i1 = match(group1, tp_levels),
         i2 = match(group2, tp_levels))

if (adjacent_only) {
  sig_pairs <- sig_pairs %>% filter(abs(i2 - i1) == 1)
}

spec_top <- pool_gs_ci %>%
  group_by(specificity) %>%
  summarise(top = max(ci_hi, na.rm = TRUE), .groups = "drop")

bracket_df <- sig_pairs %>%
  left_join(spec_top, by = "specificity") %>%
  group_by(specificity) %>%
  arrange(i1, i2, .by_group = TRUE) %>%
  mutate(level = row_number(),
         y     = top + 0.035 + (level - 1) * 0.065) %>%
  ungroup() %>%
  mutate(shift = if_else(specificity == "pp65", dodge_shift, -dodge_shift),
         x1    = i1 + shift,
         x2    = i2 + shift,
         id    = paste(specificity, group1, group2, sep = "_"))

bracket_path <- if (nrow(bracket_df) > 0) {
  do.call(rbind, lapply(seq_len(nrow(bracket_df)), function(i) {
    b <- bracket_df[i, ]
    data.frame(
      id          = b$id,
      specificity = b$specificity,
      x           = c(b$x1, b$x1, b$x2, b$x2),
      y           = c(b$y - tip_len, b$y, b$y, b$y - tip_len),
      stringsAsFactors = FALSE
    )
  }))
} else NULL

y_max <- if (nrow(bracket_df) > 0) max(1.02, max(bracket_df$y) + 0.06) else 1.02

# ── 6. Plot ─────────────────────────────────────────────────────────────────

p_pool_ci <- ggplot(pool_gs_ci,
                    aes(x     = x_pos,
                        y     = gini_simpson,
                        group = specificity,
                        color = specificity)) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.08, linewidth = 0.7) +
  geom_point(size = 3) +
  scale_color_manual(values = pool_colors) +
  scale_x_continuous(breaks = seq_along(tp_levels),
                     labels = tp_levels,
                     expand = expansion(add = 0.5)) +
  scale_y_continuous(limits = c(0, y_max), breaks = seq(0, 1, by = 0.2)) +
  labs(
    title    = "Clonal repertoire heterogeneity",
    subtitle = paste0("Gini-Simpson index with bootstrap 95% CIs (R = ", R_boot, ")"),
    caption  = paste0("Pairwise permutation tests on |\u0394 Gini-Simpson| (R = ", R_perm,
                      "), ", p_adjust_method, "-adjusted within specificity; ",
                      "* p < 0.05, ** p < 0.01, *** p < 0.001, **** p < 0.0001"),
    x        = NULL,
    y        = "Gini-Simpson index",
    color    = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title      = element_text(face = "bold"),
    plot.subtitle   = element_text(size = 10, color = "grey30"),
    plot.caption    = element_text(size = 8, color = "grey40", hjust = 0),
    legend.position = "right"
  )

if (!is.null(bracket_path)) {
  p_pool_ci <- p_pool_ci +
    geom_path(data = bracket_path,
              aes(x = x, y = y, group = id, color = specificity),
              inherit.aes = FALSE, linewidth = 0.6, show.legend = FALSE) +
    geom_text(data = bracket_df,
              aes(x = (x1 + x2) / 2, y = y + 0.010,
                  label = signif, color = specificity),
              inherit.aes = FALSE, vjust = 0, size = 5, show.legend = FALSE)
}

# ── 7. Export ───────────────────────────────────────────────────────────────

ggsave(file.path(export_dir, "D0_CMV_gini_simpson_pool_CI.jpeg"),
       plot = p_pool_ci, width = 7.5, height = 5.5, dpi = 300, device = "jpeg")

ggsave(file.path(export_dir, "D0_CMV_gini_simpson_pool_CI.png"),
       plot = p_pool_ci, width = 7.5, height = 5.5, dpi = 300, device = "png")




#
# Donor 0 — Visit 3 / Visit 4 phenotype & dynamics correlation analysis
#

library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)

# ── 0. Settings ─────────────────────────────────────────────────────────────
min_cells   <- 5
fdr_alpha   <- 0.05
top_n_hits  <- 6   # how many top scatter plots to export per specificity
out_dir     <- "Donor_0/exports/V3_V4_correlations"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# (Assumes clone_colors_fixed for the scatter color scheme.)


# ── 1. Build clone × visit feature table ────────────────────────────────────
#   - 8 cluster fractions per (clone, visit)
#   - clone log2FC in % of specificity from V3 → V4

md <- D0_data_CMV@meta.data %>%
  filter(!is.na(Global_Clone_ID),
         !is.na(seurat_clusters),
         !is.na(timepoint),
         specificity %in% c("pp65", "IE1"),
         timepoint   %in% c("Visit 3", "Visit 4")) %>%
  mutate(timepoint       = as.character(timepoint),
         seurat_clusters = as.character(seurat_clusters))

# Per (specificity, clone, visit): cell counts and cluster fractions
n_cv <- md %>%
  group_by(specificity, Global_Clone_ID, timepoint) %>%
  summarise(n_cells = n(), .groups = "drop")

cluster_frac <- md %>%
  group_by(specificity, Global_Clone_ID, timepoint, seurat_clusters) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(specificity, Global_Clone_ID, timepoint) %>%
  mutate(frac = n / sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  complete(
    nesting(specificity, Global_Clone_ID, timepoint),
    seurat_clusters = as.character(0:7),
    fill = list(frac = 0)
  ) %>%
  left_join(n_cv, by = c("specificity", "Global_Clone_ID", "timepoint")) %>%
  filter(n_cells >= min_cells)

# Pivot to wide: one row per (specificity, clone), columns V3_C0…V4_C7
wide_clust <- cluster_frac %>%
  mutate(
    visit_short = if_else(timepoint == "Visit 3", "V3", "V4"),
    col_name    = paste0(visit_short, "_C", seurat_clusters)
  ) %>%
  select(specificity, Global_Clone_ID, col_name, frac) %>%
  pivot_wider(names_from = col_name, values_from = frac)

# Clone-level log2FC in % of specificity (V3 → V4)
clone_pct_spec <- md %>%
  count(specificity, Global_Clone_ID, timepoint, name = "count") %>%
  group_by(specificity, timepoint) %>%
  mutate(pct_of_spec = count / sum(count) * 100) %>%
  ungroup() %>%
  select(specificity, Global_Clone_ID, timepoint, pct_of_spec) %>%
  pivot_wider(names_from = timepoint, values_from = pct_of_spec,
              names_prefix = "pctSpec_") %>%
  mutate(
    log2FC_V3_V4 = log2((`pctSpec_Visit 4` + 1e-6) / (`pctSpec_Visit 3` + 1e-6))
  ) %>%
  select(specificity, Global_Clone_ID, log2FC_V3_V4)

# Final feature table: must have both V3 and V4 cluster data + log2FC
feat <- wide_clust %>%
  inner_join(clone_pct_spec, by = c("specificity", "Global_Clone_ID")) %>%
  filter(if_all(starts_with("V3_C"), ~ !is.na(.)),
         if_all(starts_with("V4_C"), ~ !is.na(.)),
         !is.na(log2FC_V3_V4))

message(sprintf("Clones meeting threshold: pp65 = %d, IE1 = %d",
                sum(feat$specificity == "pp65"),
                sum(feat$specificity == "IE1")))


# ── 2. Compute Spearman correlations ────────────────────────────────────────
# Three families per specificity:
#   (a) V3 cluster fracs   × V4 cluster fracs       (8 × 8 = 64 pairs)
#   (b) V3 cluster fracs   × log2FC_V3_V4           (8 pairs)
#   (c) V4 cluster fracs   × log2FC_V3_V4           (8 pairs)

compute_cor_pairs <- function(df, x_cols, y_cols) {
  expand.grid(x = x_cols, y = y_cols, stringsAsFactors = FALSE) %>%
    rowwise() %>%
    mutate(
      n   = sum(!is.na(df[[x]]) & !is.na(df[[y]])),
      res = list(
        if (n >= 4 && sd(df[[x]], na.rm = TRUE) > 0 &&
            sd(df[[y]], na.rm = TRUE) > 0) {
          suppressWarnings(cor.test(df[[x]], df[[y]],
                                    method = "spearman",
                                    exact = FALSE))
        } else NULL
      ),
      rho     = if (is.null(res)) NA_real_ else unname(res$estimate),
      p_value = if (is.null(res)) NA_real_ else res$p.value
    ) %>%
    ungroup() %>%
    select(-res)
}

cor_results <- list()

for (spec in c("pp65", "IE1")) {
  df_s   <- feat %>% filter(specificity == spec)
  v3_cols <- paste0("V3_C", 0:7)
  v4_cols <- paste0("V4_C", 0:7)
  
  r_v3v4   <- compute_cor_pairs(df_s, v3_cols, v4_cols) %>%
    mutate(family = "V3_vs_V4")
  r_v3_fc  <- compute_cor_pairs(df_s, v3_cols, "log2FC_V3_V4") %>%
    mutate(family = "V3_vs_log2FC")
  r_v4_fc  <- compute_cor_pairs(df_s, v4_cols, "log2FC_V3_V4") %>%
    mutate(family = "V4_vs_log2FC")
  
  all_r <- bind_rows(r_v3v4, r_v3_fc, r_v4_fc) %>%
    mutate(specificity = spec,
           p_adj_BH    = p.adjust(p_value, method = "BH"),
           sig         = !is.na(p_adj_BH) & p_adj_BH < fdr_alpha)
  
  cor_results[[spec]] <- all_r
}

cor_all <- bind_rows(cor_results)
write_csv(cor_all, file.path(out_dir, "V3_V4_correlations_all_pairs.csv"))


# ── 3. Heatmaps of rho per specificity ──────────────────────────────────────

make_rho_heatmap <- function(df, spec, family,
                             x_levels, y_levels,
                             title) {
  d <- df %>%
    filter(specificity == spec, family == !!family) %>%
    mutate(
      x      = factor(x, levels = x_levels),
      y      = factor(y, levels = y_levels),
      stars  = case_when(
        is.na(p_adj_BH)      ~ "",
        p_adj_BH < 0.001     ~ "***",
        p_adj_BH < 0.01      ~ "**",
        p_adj_BH < fdr_alpha ~ "*",
        TRUE                 ~ ""
      )
    )
  
  ggplot(d, aes(x = x, y = y, fill = rho)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(aes(label = stars), size = 4, vjust = 0.78) +
    scale_fill_gradient2(low = "#2166ac", mid = "white", high = "#b2182b",
                         midpoint = 0, limits = c(-1, 1),
                         name = "Spearman ρ") +
    coord_fixed() +
    labs(title = title, x = NULL, y = NULL) +
    theme_classic(base_size = 11) +
    theme(
      axis.text.x   = element_text(angle = 45, hjust = 1),
      plot.title    = element_text(face = "bold")
    )
}

for (spec in c("pp65", "IE1")) {
  n_spec <- sum(feat$specificity == spec)
  
  p_v3v4 <- make_rho_heatmap(
    cor_all, spec, "V3_vs_V4",
    x_levels = paste0("V3_C", 0:7),
    y_levels = paste0("V4_C", 0:7),
    title    = sprintf("%s clones — V3 vs V4 cluster fractions (n = %d)",
                       spec, n_spec)
  )
  
  p_v3fc <- make_rho_heatmap(
    cor_all, spec, "V3_vs_log2FC",
    x_levels = paste0("V3_C", 0:7),
    y_levels = "log2FC_V3_V4",
    title    = sprintf("%s clones — V3 cluster fractions vs log2FC(V3→V4)",
                       spec)
  )
  
  p_v4fc <- make_rho_heatmap(
    cor_all, spec, "V4_vs_log2FC",
    x_levels = paste0("V4_C", 0:7),
    y_levels = "log2FC_V3_V4",
    title    = sprintf("%s clones — V4 cluster fractions vs log2FC(V3→V4)",
                       spec)
  )
  
  ggsave(file.path(out_dir, sprintf("%s_heatmap_V3_vs_V4.png", spec)),
         p_v3v4, width = 7, height = 7, dpi = 300)
  ggsave(file.path(out_dir, sprintf("%s_heatmap_V3_vs_log2FC.png", spec)),
         p_v3fc, width = 7, height = 3, dpi = 300)
  ggsave(file.path(out_dir, sprintf("%s_heatmap_V4_vs_log2FC.png", spec)),
         p_v4fc, width = 7, height = 3, dpi = 300)
}


# ── 4. Scatter plots for top hits ───────────────────────────────────────────
#   Top N (by |rho|) among FDR-significant pairs, per specificity.
#   Falls back to "top N by raw p" if nothing survives FDR.

get_clone_label <- function(global_clone_ids, spec) {
  pattern <- paste0("^D0_", spec, "_clone_(\\d+)$")
  ifelse(
    grepl(pattern, global_clone_ids, perl = TRUE),
    paste0("Clone ", sub(pattern, "\\1", global_clone_ids, perl = TRUE)),
    "Other"
  )
}

make_scatter <- function(df_spec, x_var, y_var, spec, rho, p_adj) {
  d <- df_spec %>%
    select(Global_Clone_ID, all_of(c(x_var, y_var))) %>%
    mutate(clone_label = get_clone_label(Global_Clone_ID, spec))
  
  ggplot(d, aes(x = .data[[x_var]], y = .data[[y_var]],
                color = clone_label)) +
    geom_point(size = 3, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = clone_label),
                             size = 3, max.overlaps = 30,
                             segment.size = 0.3,
                             segment.color = "grey60",
                             show.legend = FALSE) +
    scale_color_manual(values = clone_colors_fixed, guide = "none") +
    labs(
      title    = sprintf("%s: %s vs %s", spec, x_var, y_var),
      subtitle = sprintf("Spearman ρ = %.2f, FDR-adj p = %.3g",
                         rho, p_adj),
      x = x_var, y = y_var
    ) +
    theme_classic(base_size = 12) +
    theme(plot.title    = element_text(face = "bold"),
          plot.subtitle = element_text(size = 10, color = "grey30"))
}

for (spec in c("pp65", "IE1")) {
  df_s <- feat %>% filter(specificity == spec)
  res  <- cor_all %>% filter(specificity == spec, !is.na(rho))
  
  sig_hits <- res %>% filter(sig) %>% arrange(desc(abs(rho)))
  used_fdr <- nrow(sig_hits) > 0
  top      <- if (used_fdr) head(sig_hits, top_n_hits) else
    head(arrange(res, p_value), top_n_hits)
  
  message(sprintf("[%s] %d FDR-significant pairs; exporting top %d scatters%s.",
                  spec, sum(res$sig), nrow(top),
                  if (!used_fdr) " (FDR empty → top by raw p)" else ""))
  
  for (i in seq_len(nrow(top))) {
    h <- top[i, ]
    p <- make_scatter(df_s, h$x, h$y, spec, h$rho, h$p_adj_BH)
    fname <- sprintf("%s_top%02d_%s_vs_%s.png", spec, i, h$x, h$y)
    ggsave(file.path(out_dir, fname), p,
           width = 6, height = 5, dpi = 300)
  }
}

message("Done. Outputs in: ", out_dir)


#
#
# TCR logo plots
#
#

library(Seurat)
library(ggseqlogo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)

# ---- Extract and parse TCR data ----
md <- mdata_CMV@meta.data

tcr_long <- md %>%
  mutate(.cell = rownames(.)) %>%
  select(.cell, metacluster, t_cdr3s_aa) %>%
  filter(!is.na(metacluster), !is.na(t_cdr3s_aa), t_cdr3s_aa != "", t_cdr3s_aa != "NA") %>%
  mutate(part = str_split(t_cdr3s_aa, ";")) %>%
  unnest(part) %>%
  mutate(
    part  = str_trim(part),
    chain = str_match(part, "^(TRA|TRB):")[, 2],
    seq   = str_match(part, "^(?:TRA|TRB):(.+)$")[, 2]
  ) %>%
  filter(!is.na(chain), !is.na(seq), seq != "") %>%
  select(.cell, metacluster, chain, seq)

# ---- Helper: build faceted logos across all CDR3 lengths for one chain ----
make_chain_logo <- function(seqs, chain_label, min_n = 5) {
  if (length(seqs) == 0) return(NULL)

  # Group sequences by length, keep lengths with >= min_n
  len_groups <- split(seqs, nchar(seqs))
  len_groups <- len_groups[sapply(len_groups, length) >= min_n]
  if (length(len_groups) == 0) return(NULL)

  # Order facets by length
  lens <- as.integer(names(len_groups))
  len_groups <- len_groups[order(lens)]
  counts <- sapply(len_groups, length)
  names(len_groups) <- paste0("len ", names(len_groups))

  # Per-facet cell-count annotation
  ann <- data.frame(
    seq_group = factor(names(len_groups), levels = names(len_groups)),
    label = paste0("n = ", counts)
  )

  ggseqlogo(len_groups, method = "bits", ncol = 1) +
    geom_text(
      data = ann, inherit.aes = FALSE,
      aes(x = Inf, y = Inf, label = label),
      hjust = 1.05, vjust = 1.4, size = 3, color = "grey30"
    ) +
    ggtitle(chain_label) +
    theme(
      plot.title = element_text(size = 10, face = "bold"),
      strip.text = element_text(size = 8)
    )
}

# Count length facets (>= min_n) for a set of sequences, for canvas sizing
n_len_facets <- function(seqs, min_n = 5) {
  if (length(seqs) == 0) return(0)
  tab <- table(nchar(seqs))
  sum(tab >= min_n)
}

# ---- Generate one stacked figure (TRA over TRB) per metacluster ----
metaclusters <- sort(unique(tcr_long$metacluster))
dir.create("exports/tcr_logos", recursive = TRUE, showWarnings = FALSE)

min_n <- 5

for (mc in metaclusters) {
  sub <- tcr_long %>% filter(metacluster == mc)
  tra <- sub %>% filter(chain == "TRA") %>% pull(seq)
  trb <- sub %>% filter(chain == "TRB") %>% pull(seq)

  p_tra <- make_chain_logo(tra, "TRA", min_n = min_n)
  p_trb <- make_chain_logo(trb, "TRB", min_n = min_n)

  plots <- Filter(Negate(is.null), list(p_tra, p_trb))
  if (length(plots) == 0) {
    message("Skipping metacluster ", mc, " (insufficient sequences)")
    next
  }

  n_facets <- max(n_len_facets(tra, min_n) + n_len_facets(trb, min_n), 1)

  combined <- wrap_plots(plots, ncol = 1) +
    plot_annotation(title = paste0("Metacluster: ", mc))

  ggsave(
    filename = file.path("exports/tcr_logos", paste0("logo_", mc, ".png")),
    plot = combined, width = 8, height = 1.6 * n_facets + 1, dpi = 300
  )
}

library(Seurat)
library(ggseqlogo)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(patchwork)
library(msa)
library(Biostrings)

# =====================================================================
# Parse TCR data into long form: one row per (cell, chain, sequence)
# =====================================================================
md <- mdata_CMV@meta.data

tcr_long <- md %>%
  mutate(.cell = rownames(.)) %>%
  select(.cell, metacluster, t_cdr3s_aa) %>%
  filter(!is.na(metacluster), !is.na(t_cdr3s_aa), t_cdr3s_aa != "", t_cdr3s_aa != "NA") %>%
  mutate(part = str_split(t_cdr3s_aa, ";")) %>%
  unnest(part) %>%
  mutate(
    part  = str_trim(part),
    chain = str_match(part, "^(TRA|TRB):")[, 2],
    seq   = str_match(part, "^(?:TRA|TRB):(.+)$")[, 2]
  ) %>%
  filter(!is.na(chain), !is.na(seq), seq != "") %>%
  select(.cell, metacluster, chain, seq)

# Standard 20 amino acids (gaps handled separately, never passed as a letter)
AA <- c("A","C","D","E","F","G","H","I","K","L",
        "M","N","P","Q","R","S","T","V","W","Y")

# Build a 20 x L position frequency matrix from a character matrix of
# aligned/gapped sequences. Gaps ('-') contribute nothing; each column is
# normalized over its non-gap residues.
make_pfm <- function(mat) {
  pfm <- apply(mat, 2, function(col) {
    counts <- table(factor(col[col != "-"], levels = AA))
    if (sum(counts) == 0) return(setNames(rep(0, length(AA)), AA))
    counts / sum(counts)
  })
  rownames(pfm) <- AA
  pfm
}

# =====================================================================
# Method 1: flank trimming — keep N-term + C-term residues, drop middle
# =====================================================================
logo_flank <- function(seqs, chain_label, n_left = 4, n_right = 4, min_n = 5) {
  L <- nchar(seqs)
  seqs <- seqs[L >= (n_left + n_right)]
  if (length(seqs) < min_n) return(NULL)

  trimmed <- paste0(
    substr(seqs, 1, n_left),
    substr(seqs, nchar(seqs) - n_right + 1, nchar(seqs))
  )

  ggseqlogo(trimmed, method = "bits") +
    ggtitle(sprintf("%s — %dN+%dC (n=%d)", chain_label, n_left, n_right, length(trimmed))) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

# =====================================================================
# Method 2: IMGT-style center gapping — pad each seq to max length by
# inserting gaps in the middle, anchoring both ends. A support bar strip
# below the logo shows how many cells contribute to each column.
# =====================================================================
center_gap <- function(s, target_len) {
  L <- nchar(s)
  if (L >= target_len) return(s)
  pad  <- target_len - L
  left <- ceiling(L / 2)            # split residues at the midpoint
  paste0(substr(s, 1, left),
         strrep("-", pad),
         substr(s, left + 1, L))
}

logo_centergap <- function(seqs, chain_label, min_n = 5) {
  if (length(seqs) < min_n) return(NULL)
  target <- max(nchar(seqs))
  gapped <- vapply(seqs, center_gap, character(1), target_len = target, USE.NAMES = FALSE)

  mat     <- do.call(rbind, strsplit(gapped, ""))   # rows = seqs, cols = positions
  support <- colSums(mat != "-")
  pfm     <- make_pfm(mat)

  sup_df <- data.frame(pos = seq_len(ncol(pfm)), n = support)

  p_logo <- ggseqlogo(pfm, method = "bits") +
    ggtitle(sprintf("%s — center-gapped (n=%d, len=%d)", chain_label, length(seqs), target)) +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.text.x = element_blank(), axis.title.x = element_blank())

  p_sup <- ggplot(sup_df, aes(pos, n)) +
    geom_col(fill = "grey60", width = 0.9) +
    scale_x_continuous(breaks = sup_df$pos, expand = c(0.01, 0.01)) +
    labs(x = "position", y = "cells") +
    theme_minimal(base_size = 8) +
    theme(panel.grid.minor = element_blank())

  p_logo / p_sup + plot_layout(heights = c(4, 1))
}

# =====================================================================
# Method 3: multiple sequence alignment (Clustal Omega) then logo
# =====================================================================
logo_msa <- function(seqs, chain_label, min_n = 5) {
  if (length(seqs) < min_n) return(NULL)
  if (length(unique(seqs)) < 2) return(NULL)   # MSA degenerate on a single unique seq

  aln <- tryCatch(
    msa(AAStringSet(seqs), method = "ClustalOmega"),
    error = function(e) { message("MSA failed for ", chain_label, ": ", conditionMessage(e)); NULL }
  )
  if (is.null(aln)) return(NULL)
  aln_chr <- as.character(aln)

  mat <- do.call(rbind, strsplit(aln_chr, ""))
  pfm <- make_pfm(mat)

  ggseqlogo(pfm, method = "bits") +
    ggtitle(sprintf("%s — MSA (ClustalOmega, n=%d, len=%d)",
                    chain_label, length(seqs), nchar(aln_chr[1]))) +
    theme(plot.title = element_text(size = 10, face = "bold"))
}

# =====================================================================
# Driver: run all three methods across every metacluster
# =====================================================================
metaclusters <- sort(unique(tcr_long$metacluster))
min_n <- 5

base_dir <- "exports/tcr_logos"
dirs <- list(
  flank     = file.path(base_dir, "flank_trimmed"),
  centergap = file.path(base_dir, "center_gapped"),
  msa       = file.path(base_dir, "msa_aligned")
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

# Build a stacked TRA-over-TRB figure for a given method, then save
build_and_save <- function(builder, sub, mc, out_dir, height_per = 2.6) {
  tra <- sub %>% filter(chain == "TRA") %>% pull(seq)
  trb <- sub %>% filter(chain == "TRB") %>% pull(seq)

  p_tra <- builder(tra, "TRA", min_n = min_n)
  p_trb <- builder(trb, "TRB", min_n = min_n)

  plots <- Filter(Negate(is.null), list(p_tra, p_trb))
  if (length(plots) == 0) {
    message("  skip (insufficient sequences)")
    return(invisible(NULL))
  }

  combined <- wrap_plots(plots, ncol = 1) +
    plot_annotation(title = paste0("Metacluster: ", mc))

  ggsave(
    filename = file.path(out_dir, paste0("logo_", mc, ".png")),
    plot = combined, width = 9, height = height_per * length(plots) + 0.5, dpi = 300
  )
}

for (mc in metaclusters) {
  message("Metacluster: ", mc)
  sub <- tcr_long %>% filter(metacluster == mc)

  message(" flank-trimmed");  build_and_save(logo_flank,     sub, mc, dirs$flank,     height_per = 2.2)
  message(" center-gapped");  build_and_save(logo_centergap, sub, mc, dirs$centergap, height_per = 3.0)
  message(" msa-aligned");    build_and_save(logo_msa,       sub, mc, dirs$msa,       height_per = 2.6)
}

message("Done. Outputs in: ", normalizePath(base_dir, mustWork = FALSE))



# ── Cluster card plots — D0_data_CMV (Donor 0, et48) ─────────────────────────
#
# One PNG per et48 cluster with >10 cells (specificity-agnostic)
# Layout (2×2):
#   top-left    : text summary — cluster ID, total cells, clone breakdown
#   top-right   : UMAP — background (grey) + foreground colored by
#                 Global_Clone_ID or specificity-split minor-clones group
#   bottom-left : TCR CDR3 logo (TRA over TRB), modal-length
#   bottom-right: seurat_clusters stacked bar chart, one bar per clone /
#                 minor-clones group (>5 cells in cluster)
#
# Single donor (Donor 0): no D1–D5 grouping. Clones are D0_<spec>_clone_<n>
# (underscore form); cells with NA Global_Clone_ID are pooled as minor-clones,
# split by per-cell specificity into "minor-clones (pp65)" / "(IE1)".
# Fresh qualitative palette assigned per plot for the UMAP groups.
# Fixed seurat_clusters palette shared across all cards.

library(Seurat)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(ggseqlogo)

# Base export directory for Donor 0 cards (per-tag subfolders created below)
export_dir <- "Donor_0/exports/TCR_clusters"
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

# ── Palette helpers ───────────────────────────────────────────────────────────

# Qualitative colors for the per-plot UMAP group palette
QUAL_PALETTE <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#386cb0", "#f0027f", "#bf5b17",
  "#4daf4a", "#e41a1c", "#984ea3", "#ff7f00", "#a65628",
  "#f781bf", "#999999", "#8dd3c7", "#fb8072", "#80b1d3"
)

assign_palette <- function(groups) {
  n <- length(groups)
  colors <- QUAL_PALETTE[seq_len(min(n, length(QUAL_PALETTE)))]
  if (n > length(QUAL_PALETTE)) {
    colors <- colorRampPalette(QUAL_PALETTE)(n)
  }
  setNames(colors, groups)
}

# Fixed seurat_clusters palette (shared across all cards)
cluster_colors <- c(
  "0" = "#C49A00",
  "1" = "#53B400",
  "2" = "#8B0000",
  "3" = "#1976D2",
  "4" = "#F8766D",
  "5" = "#A58AFF",
  "6" = "#FFDB58",
  "7" = "#00C094"
)
seurat_cluster_levels <- names(cluster_colors)

# ── Build per-cluster summary data ───────────────────────────────────────────

build_cluster_summary <- function(seurat_obj, meta_col, min_cells = 10) {

  md <- seurat_obj@meta.data %>%
    filter(!is.na(.data[[meta_col]])) %>%
    transmute(
      cluster        = .data[[meta_col]],
      clone_id       = Global_Clone_ID,
      spec           = specificity,
      cdr3           = t_cdr3s_aa,
      seurat_cluster = as.character(seurat_clusters)
    )

  # Clusters with >10 cells
  eligible <- md %>%
    dplyr::count(cluster, name = "total_cells") %>%
    filter(total_cells > min_cells) %>%
    arrange(desc(total_cells))

  list(md = md, eligible = eligible)
}

# ── Top-left: text panel ──────────────────────────────────────────────────────

make_text_panel <- function(cluster_id, cluster_md) {

  total <- nrow(cluster_md)

  lines <- c(
    sprintf("Cluster %s", cluster_id),
    sprintf("n = %d cells", total),
    ""
  )

  # Named clones: sort by spec (parsed from clone ID) then clone number
  named <- cluster_md %>%
    filter(!is.na(clone_id)) %>%
    count(clone_id, name = "n_cells") %>%
    mutate(
      spec      = str_extract(clone_id, "pp65|IE1"),
      clone_num = as.integer(str_extract(clone_id, "\\d+$"))
    ) %>%
    arrange(spec, clone_num)

  for (i in seq_len(nrow(named))) {
    # Strip donor prefix for brevity: D0_pp65_clone_3 -> pp65_clone_3
    short <- str_remove(named$clone_id[i], "^D0_")
    lines <- c(lines, sprintf("%s:  %d", short, named$n_cells[i]))
  }

  # Minor clones (NA Global_Clone_ID), split by per-cell specificity
  minor <- cluster_md %>%
    filter(is.na(clone_id)) %>%
    count(spec, name = "n_cells") %>%
    arrange(spec)

  for (i in seq_len(nrow(minor))) {
    lines <- c(lines, sprintf("minor (%s):  %d",
                              minor$spec[i], minor$n_cells[i]))
  }

  text_block <- paste(lines, collapse = "\n")

  ggplot() +
    annotate(
      "text",
      x = 0, y = 1,
      label      = text_block,
      hjust      = 0, vjust = 1,
      size       = 3.2,
      family     = "mono",
      lineheight = 1.4
    ) +
    xlim(0, 1) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin = margin(10, 0, 10, 10))
}

# ── Top-right: UMAP panel ─────────────────────────────────────────────────────

make_umap_panel <- function(seurat_obj, cluster_id, meta_col,
                             cluster_md, color_map) {

  umap_df      <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  umap_cols    <- colnames(umap_df)
  umap_df$cell <- rownames(umap_df)

  md_full <- seurat_obj@meta.data %>%
    mutate(
      group_id = case_when(
        !is.na(Global_Clone_ID) ~ as.character(Global_Clone_ID),
        !is.na(.data[[meta_col]]) & .data[[meta_col]] == as.character(cluster_id) ~
          paste0("minor-clones (", specificity, ")"),
        TRUE ~ NA_character_
      )
    )

  in_cluster <- !is.na(seurat_obj@meta.data[[meta_col]]) &
    seurat_obj@meta.data[[meta_col]] == as.character(cluster_id)

  umap_df$group_id <- NA_character_
  umap_df$group_id[in_cluster] <- md_full$group_id[in_cluster]
  umap_df$group_id <- factor(umap_df$group_id, levels = names(color_map))

  bg_df <- umap_df[!in_cluster, ]
  fg_df <- umap_df[ in_cluster, ]

  ggplot() +
    geom_point(
      data  = bg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]]),
      color = "lightgrey", size = 0.3
    ) +
    geom_point(
      data = fg_df,
      aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]],
          color = group_id),
      size = 1.4
    ) +
    scale_color_manual(
      values = color_map,
      name   = "Group",
      drop   = FALSE,
      guide  = guide_legend(
        override.aes = list(size = 3),
        ncol = 1
      )
    ) +
    theme_void() +
    theme(
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.key.size = unit(0.4, "cm"),
      plot.margin     = margin(5, 5, 5, 5)
    )
}

# ── Bottom-left: TCR CDR3 logo panel (TRA over TRB), modal-length ────────────

extract_cdr3 <- function(cdr3_vec, chain_prefix) {
  # chain_prefix is "TRA:" or "TRB:"
  seqs <- unlist(lapply(cdr3_vec, function(x) {
    if (is.na(x)) return(character(0))
    parts <- strsplit(x, ";")[[1]]
    sub(paste0("^", chain_prefix), "",
        parts[grepl(paste0("^", chain_prefix), parts)])
  }))
  seqs[nchar(seqs) > 0]
}

# Keep only sequences of the modal (most common) length so the logo is
# rectangular and interpretable.
filter_modal_length <- function(seqs) {
  if (length(seqs) == 0) return(character(0))
  lens      <- nchar(seqs)
  modal_len <- as.integer(names(sort(table(lens), decreasing = TRUE))[1])
  seqs[lens == modal_len]
}

make_logo_panel <- function(cluster_md_cells_cdr3, weight_mode = "unique") {
  # cluster_md_cells_cdr3 = vector of t_cdr3s_aa strings for cells in cluster
  # weight_mode:
  #   "unique" - each distinct CDR3 contributes once (a dominant clone does not
  #              dominate the position weights)
  #   "cell"   - each cell contributes (logo weighted by clonal abundance)

  pull_seqs <- function(prefix) {
    seqs <- extract_cdr3(cluster_md_cells_cdr3, prefix)
    if (weight_mode == "unique") seqs <- unique(seqs)
    filter_modal_length(seqs)
  }

  tra_seqs <- pull_seqs("TRA:")
  trb_seqs <- pull_seqs("TRB:")

  # Title count label reflects the weighting mode
  count_label <- if (weight_mode == "unique") "unique" else "cells"

  build_one <- function(seqs, title) {
    if (length(seqs) < 1) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = sprintf("%s: no\nsequences", title),
                   size = 3) +
          theme_void()
      )
    }
    # Single sequence: use probability method so it renders as clean
    # full-height letters (bits would also be full height but implies diversity)
    logo_method <- if (length(seqs) == 1) "prob" else "bits"
    ggseqlogo(seqs, method = logo_method) +
      ggtitle(sprintf("%s  (%d %s, len = %d)",
                      title, length(seqs), count_label, nchar(seqs[1]))) +
      theme(
        plot.title      = element_text(size = 8, face = "bold"),
        axis.text       = element_text(size = 6),
        axis.title      = element_text(size = 7),
        legend.position = "none"
      )
  }

  tra_logo <- build_one(tra_seqs, "TRA")
  trb_logo <- build_one(trb_seqs, "TRB")

  # TRA on top, TRB below
  cowplot::plot_grid(tra_logo, trb_logo, ncol = 1)
}

# ── Bottom-right: seurat_clusters stacked bar panel ──────────────────────────
#
# One bar per group (named clone or "minor-clones"), y = % of that group's
# cells in each seurat_cluster. Groups with >5 cells in the cluster only.
# Bar order: named clones by spec then clone number, minor-clones last.

make_seurat_bar_panel <- function(cluster_md) {

  # Assign each cell a group_id: named clone, or specificity-split minor-clones.
  # Named clones carry their own specificity in the ID; minor clones (NA clone)
  # take specificity from the per-cell `spec` column.
  md_grouped <- cluster_md %>%
    mutate(
      group_id = if_else(
        !is.na(clone_id),
        as.character(clone_id),
        paste0("minor-clones (", spec, ")")
      )
    )

  # Keep only groups with >5 cells in the cluster
  keep_groups <- md_grouped %>%
    count(group_id, name = "total_in_cluster") %>%
    filter(total_in_cluster > 5) %>%
    pull(group_id)

  bar_df <- md_grouped %>%
    filter(group_id %in% keep_groups, !is.na(seurat_cluster)) %>%
    count(group_id, clone_id, seurat_cluster, name = "n") %>%
    group_by(group_id) %>%
    mutate(percent = n / sum(n) * 100) %>%
    ungroup() %>%
    mutate(seurat_cluster = factor(seurat_cluster,
                                   levels = seurat_cluster_levels))

  # Ordered group levels: named clones by spec + clone number, then minor
  # groups (pp65 before IE1) last. For named clones, specificity is parsed from
  # the Global_Clone_ID itself; for minor groups, from the "(spec)" suffix.
  # This guarantees one row per group_id, so factor levels are unique.
  group_order <- bar_df %>%
    distinct(group_id, clone_id) %>%
    mutate(
      is_minor  = is.na(clone_id),
      spec      = if_else(
        is_minor,
        str_extract(group_id, "pp65|IE1"),
        str_extract(clone_id, "pp65|IE1")
      ),
      clone_num = if_else(is_minor, Inf,
                          as.integer(str_extract(clone_id, "\\d+$"))) %>%
                    as.numeric()
    ) %>%
    arrange(is_minor, spec, clone_num) %>%
    pull(group_id)

  bar_df <- bar_df %>%
    mutate(group_id = factor(group_id, levels = group_order))

  # x labels: strip donor prefix from named clones; render minor groups as
  # "minor (pp65)" / "minor (IE1)"
  short_label <- function(g) {
    ifelse(
      grepl("^minor-clones", g),
      str_replace(g, "^minor-clones ", "minor "),
      str_remove(g, "^D0_")
    )
  }

  ggplot(bar_df, aes(x = group_id, y = percent, fill = seurat_cluster)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = cluster_colors, name = "Seurat cluster",
                      drop = FALSE) +
    scale_x_discrete(labels = short_label) +
    labs(x = NULL, y = "% of cells") +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x     = element_text(angle = 45, hjust = 1, size = 7),
      axis.text.y     = element_text(size = 7),
      legend.text     = element_text(size = 7),
      legend.title    = element_text(size = 8, face = "bold"),
      legend.key.size = unit(0.35, "cm"),
      plot.margin     = margin(5, 5, 5, 5)
    )
}

# ── Assemble and save one card (2×2 layout) ──────────────────────────────────

save_cluster_card <- function(seurat_obj, cluster_id, meta_col,
                               cluster_md, export_dir,
                               weight_mode = "unique") {

  # Foreground groups for the UMAP palette
  named_groups <- cluster_md %>%
    filter(!is.na(clone_id)) %>%
    pull(clone_id) %>% unique() %>% sort()

  # Minor clones split by per-cell specificity, e.g. "minor-clones (pp65)"
  minor_groups <- cluster_md %>%
    filter(is.na(clone_id)) %>%
    mutate(group_id = paste0("minor-clones (", spec, ")")) %>%
    pull(group_id) %>% unique() %>% sort()

  all_groups <- c(named_groups, minor_groups)
  color_map  <- assign_palette(all_groups)

  # Build the four panels
  text_panel <- make_text_panel(cluster_id, cluster_md)
  umap_panel <- make_umap_panel(seurat_obj, cluster_id, meta_col,
                                cluster_md, color_map)
  logo_panel <- make_logo_panel(cluster_md$cdr3, weight_mode = weight_mode)
  bar_panel  <- make_seurat_bar_panel(cluster_md)

  # Top row: text + UMAP
  top_row <- plot_grid(
    text_panel, umap_panel,
    nrow = 1, rel_widths = c(0.3, 0.7)
  )

  # Bottom row: logo + bar chart
  bottom_row <- plot_grid(
    logo_panel, bar_panel,
    nrow = 1, rel_widths = c(0.45, 0.55)
  )

  # Stack the two rows
  card <- plot_grid(
    top_row, bottom_row,
    ncol = 1, rel_heights = c(0.55, 0.45)
  )

  out_file <- file.path(
    export_dir,
    sprintf("cluster_%s_card.png", cluster_id)
  )

  ggsave(filename = out_file, plot = card,
         width = 10, height = 10, dpi = 150, bg = "white")
}

# ── AA chemistry legend (ggseqlogo "chemistry" scheme) ───────────────────────

AA_CHEMISTRY <- data.frame(
  letter = c("G","S","T","Y","C","N","Q","K","R","H","D","E","P","A","W","F","L","I","M","V"),
  group  = c(rep("Polar", 5), rep("Neutral", 2), rep("Basic", 3),
             rep("Acidic", 2), rep("Hydrophobic", 8)),
  col    = c(rep("#109648", 5), rep("#5E239D", 2), rep("#255C99", 3),
             rep("#D62839", 2), rep("#221E22", 8)),
  stringsAsFactors = FALSE
)

export_aa_legend <- function(present_letters, out_file) {

  leg_df <- AA_CHEMISTRY %>%
    filter(letter %in% present_letters)

  if (nrow(leg_df) == 0) return(invisible(NULL))

  group_levels <- c("Polar", "Neutral", "Basic", "Acidic", "Hydrophobic")

  group_summary <- leg_df %>%
    group_by(group) %>%
    summarise(
      col   = dplyr::first(col),
      aas   = paste(sort(letter), collapse = ", "),
      .groups = "drop"
    ) %>%
    mutate(group = factor(group, levels = group_levels)) %>%
    arrange(group) %>%
    mutate(label = sprintf("%s (%s)", group, aas))

  color_map <- setNames(group_summary$col, group_summary$label)

  dummy <- data.frame(
    x     = seq_len(nrow(group_summary)),
    y     = 1,
    label = factor(group_summary$label, levels = group_summary$label)
  )

  leg_plot <- ggplot(dummy, aes(x = x, y = y, fill = label)) +
    geom_tile() +
    scale_fill_manual(values = color_map, name = "Amino acid chemistry") +
    guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    theme_void() +
    theme(
      legend.position = "right",
      legend.title    = element_text(size = 11, face = "bold"),
      legend.text     = element_text(size = 10),
      legend.key.size = unit(0.6, "cm")
    )

  legend_grob <- cowplot::get_legend(leg_plot)

  ggsave(filename = out_file,
         plot     = cowplot::ggdraw(legend_grob),
         width    = 3.5, height = 3, dpi = 150, bg = "white")

  message(sprintf("AA chemistry legend saved to %s", out_file))
}

# ── Run across all eligible clusters for one et tag ──────────────────────────

run_cluster_cards <- function(seurat_obj, et_tag, base_export_dir,
                              min_cells = 10) {

  meta_col   <- paste0("TCR_cluster_", et_tag)
  export_dir <- file.path(base_export_dir, paste0("cluster_cards_", et_tag))
  orig_dir   <- file.path(export_dir, "original_logo_weights")
  dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(orig_dir,   showWarnings = FALSE, recursive = TRUE)

  res      <- build_cluster_summary(seurat_obj, meta_col, min_cells)
  md_all   <- res$md
  eligible <- res$eligible

  message(sprintf("[%s] %d clusters with >10 cells", et_tag, nrow(eligible)))

  present_aa <- character(0)

  for (i in seq_len(nrow(eligible))) {

    cl_id <- eligible$cluster[i]
    cl_md <- md_all %>% filter(cluster == cl_id)

    # Track amino acids that will appear in this card's logos. The modal-length
    # AA set is identical for unique- and cell-weighted modes, so tracking on
    # the deduplicated modal-length sequences covers both card versions.
    tra_seqs <- filter_modal_length(unique(extract_cdr3(cl_md$cdr3, "TRA:")))
    trb_seqs <- filter_modal_length(unique(extract_cdr3(cl_md$cdr3, "TRB:")))
    if (length(tra_seqs) >= 1) {
      present_aa <- union(present_aa, unlist(strsplit(tra_seqs, "")))
    }
    if (length(trb_seqs) >= 1) {
      present_aa <- union(present_aa, unlist(strsplit(trb_seqs, "")))
    }

    # Unique-weighted card (each distinct CDR3 weighted equally) -> main folder
    save_cluster_card(
      seurat_obj  = seurat_obj,
      cluster_id  = cl_id,
      meta_col    = meta_col,
      cluster_md  = cl_md,
      export_dir  = export_dir,
      weight_mode = "unique"
    )

    # Original cell-weighted card (logo weighted by clonal abundance)
    # -> original_logo_weights/ subfolder
    save_cluster_card(
      seurat_obj  = seurat_obj,
      cluster_id  = cl_id,
      meta_col    = meta_col,
      cluster_md  = cl_md,
      export_dir  = orig_dir,
      weight_mode = "cell"
    )

    if (i %% 10 == 0) {
      message(sprintf("[%s] %d / %d clusters done (both versions)",
                      et_tag, i, nrow(eligible)))
    }
  }

  # Export the standalone AA chemistry legend into both folders (identical key)
  export_aa_legend(
    present_letters = present_aa,
    out_file        = file.path(export_dir, "aa_chemistry_legend.png")
  )
  export_aa_legend(
    present_letters = present_aa,
    out_file        = file.path(orig_dir, "aa_chemistry_legend.png")
  )

  message(sprintf("[%s] All cluster cards saved to %s (+ original_logo_weights/)",
                  et_tag, export_dir))
}

# ── Generate cards for both et tags ──────────────────────────────────────────

for (et_tag in c("et48", "et36")) {
  run_cluster_cards(
    seurat_obj      = D0_data_CMV,
    et_tag          = et_tag,
    base_export_dir = export_dir,
    min_cells       = 10
  )
}


# ============================================================================
#  D0_data_CMV — Timepoint-split subsetting, re-clustering & analysis
# ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(tidyr)
library(purrr)
library(readr)
library(tibble)

# ── 0. Setup: export folders ────────────────────────────────────────────────
base_export <- "Donor_0/exports/timepoint_split"
subdirs <- c(
  "objects",
  "umaps_clusters",
  "specificity_overlay",
  "stacked_bars",
  "clone_panels/V1_3",
  "clone_panels/V4",
  "v4_feature_plots",
  "v4_heatmap"
)
for (d in subdirs) dir.create(file.path(base_export, d), recursive = TRUE, showWarnings = FALSE)

# Convenience: dynamic ggsave wrapper (JPEG, 300 dpi)
save_jpeg <- function(plot, path, width, height) {
  ggsave(filename = path, plot = plot, width = width, height = height,
         dpi = 300, device = "jpeg", limitsize = FALSE)
}

# ── 1. Subset by timepoint ──────────────────────────────────────────────────
# Sanity check on timepoint values
message("Timepoint values present: ",
        paste(sort(unique(D0_data_CMV$timepoint)), collapse = ", "))

V4_cells   <- colnames(D0_data_CMV)[D0_data_CMV$timepoint == "Visit 4"]
V1_3_cells <- colnames(D0_data_CMV)[D0_data_CMV$timepoint %in%
                                      c("Visit 1", "Visit 2", "Visit 3")]

obj_V4   <- subset(D0_data_CMV, cells = V4_cells)
obj_V1_3 <- subset(D0_data_CMV, cells = V1_3_cells)

message(sprintf("Visit 4 subset:   %d cells", ncol(obj_V4)))
message(sprintf("Visit 1-3 subset: %d cells", ncol(obj_V1_3)))

# ── 2. Fresh re-clustering pipeline ─────────────────────────────────────────
# Full pipeline: LogNormalize -> HVG -> Scale -> PCA -> Neighbors -> Clusters -> UMAP
recluster_obj <- function(obj, npcs = 30, dims = 1:30, resolution = 0.8,
                          n_features = 2000, label = "object") {
  message(sprintf("\n── Re-clustering %s ──", label))

  DefaultAssay(obj) <- "RNA"

  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 1e4, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = n_features, verbose = FALSE)
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- FindClusters(obj, resolution = resolution, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims, verbose = FALSE)

  message(sprintf("%s: %d clusters", label, length(levels(obj$seurat_clusters))))
  obj
}

obj_V4   <- recluster_obj(obj_V4,   resolution = 0.8, label = "Visit 4")
obj_V1_3 <- recluster_obj(obj_V1_3, resolution = 0.8, label = "Visit 1-3")

# Save reclustered objects
saveRDS(obj_V4,   file.path(base_export, "objects", "obj_V4_reclustered.rds"))
saveRDS(obj_V1_3, file.path(base_export, "objects", "obj_V1_3_reclustered.rds"))

# ── Color conventions ───────────────────────────────────────────────────────
spec_colors <- c("pp65" = "#1976D2", "IE1" = "#8B0000")   # blue / red

clone_colors_pp65 <- c(
  "Clone 1"  = "#1b9e77", "Clone 2"  = "#d95f02", "Clone 3"  = "#7570b3",
  "Clone 4"  = "#e7298a", "Clone 5"  = "#66a61e", "Clone 6"  = "#e6ab02",
  "Clone 7"  = "#a6761d", "Clone 8"  = "#666666", "Clone 9"  = "#1f78b4",
  "Clone 10" = "#b2df8a", "Clone 11" = "#33a02c", "Clone 12" = "#fb9a99",
  "Clone 13" = "#e31a1c", "Clone 14" = "#fdbf6f", "Clone 15" = "#ff7f00",
  "Clone 16" = "#cab2d6", "Clone 17" = "#6a3d9a", "Clone 18" = "#ffff99",
  "Clone 19" = "#b15928", "Clone 20" = "#8dd3c7"
)
clone_colors_ie1 <- clone_colors_pp65[1:7]

# ── 3. Cluster UMAPs for each object ────────────────────────────────────────
plot_cluster_umap <- function(obj, title) {
  DimPlot(obj, group.by = "seurat_clusters", reduction = "umap",
          label = TRUE, pt.size = 0.5) +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

p_clu_V4   <- plot_cluster_umap(obj_V4,   "D0 CMV — Visit 4: seurat_clusters")
p_clu_V1_3 <- plot_cluster_umap(obj_V1_3, "D0 CMV — Visit 1-3: seurat_clusters")

save_jpeg(p_clu_V4,
          file.path(base_export, "umaps_clusters", "V4_seurat_clusters_UMAP.jpeg"), 8, 7)
save_jpeg(p_clu_V1_3,
          file.path(base_export, "umaps_clusters", "V1_3_seurat_clusters_UMAP.jpeg"), 8, 7)

# ── 4. Specificity overlay UMAPs ────────────────────────────────────────────
# Background = all cells in subset; colored by specificity (pp65 blue / IE1 red)
plot_specificity_umap <- function(obj, title) {
  DimPlot(obj, group.by = "specificity", reduction = "umap", pt.size = 0.5) +
    scale_color_manual(values = spec_colors, na.value = "grey90") +
    ggtitle(title) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
}

p_spec_V4   <- plot_specificity_umap(obj_V4,   "D0 CMV — Visit 4: pp65 vs IE1")
p_spec_V1_3 <- plot_specificity_umap(obj_V1_3, "D0 CMV — Visit 1-3: pp65 vs IE1")

save_jpeg(p_spec_V4,
          file.path(base_export, "specificity_overlay", "V4_specificity_UMAP.jpeg"), 8, 7)
save_jpeg(p_spec_V1_3,
          file.path(base_export, "specificity_overlay", "V1_3_specificity_UMAP.jpeg"), 8, 7)

# ── 5. Stacked bar graphs: specificity distribution across clusters ─────────
plot_stacked_bar <- function(obj, title) {
  df <- obj@meta.data %>%
    filter(!is.na(specificity), specificity %in% c("pp65", "IE1")) %>%
    count(specificity, seurat_clusters, name = "n") %>%
    group_by(specificity) %>%                       # normalize within each specificity
    mutate(prop = n / sum(n)) %>%
    ungroup()

  n_clusters <- length(levels(obj$seurat_clusters))
  cluster_pal <- setNames(
    scales::hue_pal()(n_clusters),
    levels(obj$seurat_clusters)
  )

  list(
    count = ggplot(df, aes(x = specificity, y = n, fill = seurat_clusters)) +
      geom_col() +
      scale_fill_manual(values = cluster_pal) +
      labs(title = paste0(title, " (counts)"), x = "Specificity",
           y = "n cells", fill = "seurat_cluster") +
      theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5)),
    prop = ggplot(df, aes(x = specificity, y = prop, fill = seurat_clusters)) +
      geom_col() +
      scale_fill_manual(values = cluster_pal) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(title = paste0(title, " (% in each cluster)"), x = "Specificity",
           y = "Percent of specificity", fill = "seurat_cluster") +
      theme_bw() + theme(plot.title = element_text(face = "bold", hjust = 0.5))
  )
}

bars_V4   <- plot_stacked_bar(obj_V4,   "Visit 4")
bars_V1_3 <- plot_stacked_bar(obj_V1_3, "Visit 1-3")

save_jpeg(bars_V4$count,
          file.path(base_export, "stacked_bars", "V4_specificity_counts.jpeg"), 5, 5)
save_jpeg(bars_V4$prop,
          file.path(base_export, "stacked_bars", "V4_specificity_proportion.jpeg"), 5, 5)
save_jpeg(bars_V1_3$count,
          file.path(base_export, "stacked_bars", "V1_3_specificity_counts.jpeg"), 5, 5)
save_jpeg(bars_V1_3$prop,
          file.path(base_export, "stacked_bars", "V1_3_specificity_proportion.jpeg"), 5, 5)

# ── 6. Per-clone highlighted UMAP panels ────────────────────────────────────
# One panel per unique Global_Clone_ID; background = all cells in subset.
# pp65 -> 20 clones, IE1 -> 7 clones.

# Single UMAP panel highlighting one clone
make_clone_panel <- function(obj, clone_id, hl_color) {
  cells_hl <- colnames(obj)[!is.na(obj$Global_Clone_ID) &
                              obj$Global_Clone_ID == clone_id]
  n <- length(cells_hl)

  DimPlot(obj, reduction = "umap",
          cells.highlight = cells_hl,
          cols.highlight  = hl_color,
          cols            = "grey85",
          pt.size         = 0.3,
          sizes.highlight = 1.4) +
    ggtitle(paste0(clone_id, "\n(n = ", n, ")")) +
    theme_bw() +
    theme(
      plot.title      = element_text(size = 8, face = "bold", hjust = 0.5),
      legend.position = "none",
      axis.title      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      plot.margin     = margin(2, 2, 2, 2)
    )
}

# Build & export a clone grid for one specificity in one object
export_clone_grid <- function(obj, specificity, clone_palette, out_subdir, obj_label,
                              max_clone = NULL) {
  prefix <- paste0("D0_", specificity, "_clone_")

  # Clones present in this object for this specificity, ordered numerically
  present <- obj$Global_Clone_ID[!is.na(obj$Global_Clone_ID) &
                                   str_starts(obj$Global_Clone_ID, prefix)]
  clone_ids <- unique(present)
  clone_num <- as.integer(str_remove(clone_ids, prefix))
  clone_ids <- clone_ids[order(clone_num)]
  clone_num <- sort(clone_num)

  # Cap at top `max_clone` by clone number, if specified
  if (!is.null(max_clone)) {
    keep      <- clone_num <= max_clone
    clone_ids <- clone_ids[keep]
  }

  if (length(clone_ids) == 0) {
    message(sprintf("No %s clones found in %s", specificity, obj_label))
    return(invisible(NULL))
  }

  hl_color <- spec_colors[[specificity]]
  panels <- map(clone_ids, ~ make_clone_panel(obj, .x, hl_color))

  ncol_grid <- 5
  nrow_grid <- ceiling(length(panels) / ncol_grid)
  grid <- wrap_plots(panels, ncol = ncol_grid) +
    plot_annotation(
      title = sprintf("%s — %s clones (n = %d)", obj_label, specificity, length(clone_ids)),
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
    )

  out_path <- file.path(base_export, "clone_panels", out_subdir,
                        sprintf("%s_%s_clone_panels.jpeg", obj_label, specificity))
  save_jpeg(grid, out_path, width = ncol_grid * 2.6, height = nrow_grid * 2.8)
  message(sprintf("Exported %d %s clone panels for %s", length(clone_ids), specificity, obj_label))
}

# Visit 4 object
export_clone_grid(obj_V4,   "IE1",  clone_colors_ie1,  "V4",   "V4",   max_clone = 7)
export_clone_grid(obj_V4,   "pp65", clone_colors_pp65, "V4",   "V4",   max_clone = 20)

# Visit 1-3 object
export_clone_grid(obj_V1_3, "IE1",  clone_colors_ie1,  "V1_3", "V1_3", max_clone = 7)
export_clone_grid(obj_V1_3, "pp65", clone_colors_pp65, "V1_3", "V1_3", max_clone = 20)

# ── 7. V4 feature plots: module scores ──────────────────────────────────────
score_cols <- c("Tcm_score", "Tem_score", "T-CMC1_score", "NK-like_score", "ISG+_score")

# Confirm columns exist
missing_scores <- setdiff(score_cols, colnames(obj_V4@meta.data))
if (length(missing_scores) > 0)
  warning("Missing score columns in obj_V4: ", paste(missing_scores, collapse = ", "))
present_scores <- intersect(score_cols, colnames(obj_V4@meta.data))

for (sc in present_scores) {
  p <- FeaturePlot(obj_V4, features = sc, reduction = "umap", pt.size = 0.5) +
    scale_color_gradientn(colors = c("grey90", "#FDBB84", "#B30000")) +
    ggtitle(paste0("Visit 4 — ", sc)) +
    theme_bw() +
    theme(plot.title = element_text(face = "bold", hjust = 0.5))
  safe_name <- str_replace_all(sc, "[^A-Za-z0-9]", "_")
  save_jpeg(p, file.path(base_export, "v4_feature_plots",
                        sprintf("V4_featureplot_%s.jpeg", safe_name)), 7, 6.5)
}

# Combined paneled version
feat_panels <- map(present_scores, function(sc) {
  FeaturePlot(obj_V4, features = sc, reduction = "umap", pt.size = 0.4) +
    scale_color_gradientn(colors = c("grey90", "#FDBB84", "#B30000")) +
    ggtitle(sc) +
    theme_bw() +
    theme(plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.key.size = unit(0.3, "cm"))
})
feat_grid <- wrap_plots(feat_panels, ncol = 3)
save_jpeg(feat_grid, file.path(base_export, "v4_feature_plots",
                              "V4_featureplots_combined.jpeg"),
          width = 15, height = 5 * ceiling(length(present_scores) / 3))

# ── 8. V4 heatmap: mean module score per cluster ────────────────────────────
score_mat <- obj_V4@meta.data %>%
  select(seurat_clusters, all_of(present_scores)) %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(present_scores), mean, na.rm = TRUE), .groups = "drop")

# Export underlying values as CSV
write_csv(score_mat,
          file.path(base_export, "v4_heatmap", "V4_mean_scores_per_cluster.csv"))

# Long format for ggplot heatmap
heat_df <- score_mat %>%
  pivot_longer(-seurat_clusters, names_to = "score", values_to = "mean_value") %>%
  group_by(score) %>%                                    # scale within each score (z-score)
  mutate(z = as.numeric(scale(mean_value))) %>%
  ungroup() %>%
  mutate(score = factor(score, levels = present_scores))

p_heat <- ggplot(heat_df, aes(x = seurat_clusters, y = score, fill = z)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_text(aes(label = sprintf("%.2f", mean_value)), size = 3) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       midpoint = 0, name = "z-score") +
  labs(title = "Visit 4 — mean module score per cluster",
       x = "seurat_cluster", y = NULL) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),
        axis.text.y = element_text(face = "bold"))

save_jpeg(p_heat, file.path(base_export, "v4_heatmap", "V4_score_heatmap.jpeg"),
          width = 1.2 * length(levels(obj_V4$seurat_clusters)) + 3, height = 5)


# ── 9. Cluster-defining gene heatmaps (protein-coding only) ─────────────────
dir.create(file.path(base_export, "marker_heatmaps"), recursive = TRUE, showWarnings = FALSE)

# Protein-coding gene reference (build once)
# BiocManager::install("EnsDb.Hsapiens.v86")   # if not already installed
library(EnsDb.Hsapiens.v86)
coding_genes <- genes(EnsDb.Hsapiens.v86,
                      filter = GeneBiotypeFilter("protein_coding"))$gene_name
coding_genes <- unique(coding_genes[!is.na(coding_genes) & coding_genes != ""])
detach("package:EnsDb.Hsapiens.v86", unload = TRUE)
detach("package:ensembldb", unload = TRUE)
message(sprintf("Protein-coding reference: %d gene symbols", length(coding_genes)))

is_protein_coding <- function(genes) genes %in% coding_genes

export_marker_heatmap <- function(obj, obj_label, n_top = 10, resolution_assay = "RNA") {
  message(sprintf("\n── FindAllMarkers: %s ──", obj_label))

  DefaultAssay(obj) <- resolution_assay
  Idents(obj) <- "seurat_clusters"

  markers <- FindAllMarkers(
    obj,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.25,
    verbose         = FALSE
  )

  # ── Keep only protein-coding genes ──
  n_before <- nrow(markers)
  markers  <- markers %>% filter(is_protein_coding(gene))
  message(sprintf("%s: %d / %d marker rows retained after protein-coding filter",
                  obj_label, nrow(markers), n_before))

  # Export full (coding-only) marker table
  write_csv(markers,
            file.path(base_export, "marker_heatmaps",
                      sprintf("%s_FindAllMarkers_coding.csv", obj_label)))

  # Top N markers per cluster by avg_log2FC
  top_markers <- markers %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = n_top, with_ties = FALSE) %>%
    ungroup()

  write_csv(top_markers,
            file.path(base_export, "marker_heatmaps",
                      sprintf("%s_top%d_markers_coding.csv", obj_label, n_top)))

  # Scale the genes we want to plot (DoHeatmap reads scale.data)
  genes_to_plot <- unique(top_markers$gene)
  obj <- ScaleData(obj, features = genes_to_plot, verbose = FALSE)

  p <- DoHeatmap(
    obj,
    features = genes_to_plot,
    group.by = "seurat_clusters",
    size     = 3.5
  ) +
    scale_fill_gradientn(colors = c("#2166AC", "white", "#B2182B"),
                         na.value = "white") +
    ggtitle(sprintf("%s — top %d cluster markers (protein-coding)", obj_label, n_top)) +
    theme(plot.title  = element_text(face = "bold", hjust = 0.5),
          axis.text.y = element_text(size = 15))

  n_clusters <- length(levels(obj$seurat_clusters))
  out_path <- file.path(base_export, "marker_heatmaps",
                        sprintf("%s_top%d_marker_heatmap_coding.jpeg", obj_label, n_top))
  save_jpeg(p,
            out_path,
            width  = max(8, 1.1 * n_clusters + 4),
            height = max(7, 0.22 * length(genes_to_plot) + 2))

  message(sprintf("%s: %d genes plotted across %d clusters",
                  obj_label, length(genes_to_plot), n_clusters))
  invisible(top_markers)
}

export_marker_heatmap(obj_V4,   "V4",   n_top = 10)
export_marker_heatmap(obj_V1_3, "V1_3", n_top = 10)

# ── 10. Clone panels on the full D0_data_CMV object ─────────────────────────

export_clone_grid(D0_data_CMV, "IE1",  clone_colors_ie1,  "full", "full", max_clone = 7)
export_clone_grid(D0_data_CMV, "pp65", clone_colors_pp65, "full", "full", max_clone = 20)


export_clone_grid_skinny <- function(obj, specificity, clone_palette, out_subdir, obj_label,
                              max_clone = NULL) {
  prefix <- paste0("D0_", specificity, "_clone_")

  # Clones present in this object for this specificity, ordered numerically
  present <- obj$Global_Clone_ID[!is.na(obj$Global_Clone_ID) &
                                   str_starts(obj$Global_Clone_ID, prefix)]
  clone_ids <- unique(present)
  clone_num <- as.integer(str_remove(clone_ids, prefix))
  clone_ids <- clone_ids[order(clone_num)]
  clone_num <- sort(clone_num)

  # Cap at top `max_clone` by clone number, if specified
  if (!is.null(max_clone)) {
    keep      <- clone_num <= max_clone
    clone_ids <- clone_ids[keep]
  }

  if (length(clone_ids) == 0) {
    message(sprintf("No %s clones found in %s", specificity, obj_label))
    return(invisible(NULL))
  }

  hl_color <- spec_colors[[specificity]]
  panels <- map(clone_ids, ~ make_clone_panel(obj, .x, hl_color))

  ncol_grid <- 3
  nrow_grid <- ceiling(length(panels) / ncol_grid)
  grid <- wrap_plots(panels, ncol = ncol_grid) +
    plot_annotation(
      title = sprintf("%s — %s clones (n = %d)", obj_label, specificity, length(clone_ids)),
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
    )

  out_path <- file.path(base_export, "clone_panels", out_subdir,
                        sprintf("%s_%s_clone_panels_skinny.jpeg", obj_label, specificity))
  save_jpeg(grid, out_path, width = ncol_grid * 2.6, height = nrow_grid * 2.8)
  message(sprintf("Exported %d %s clone panels for %s", length(clone_ids), specificity, obj_label))
}

export_clone_grid_skinny(D0_data_CMV, "IE1",  clone_colors_ie1,  "full", "full", max_clone = 7)


# ============================================================================
#  Timepoint-colored clone panels (full D0_data_CMV object)
#  Recreates: full_pp65_clone_panels  and  full_IE1_clone_panels_skinny
#  but colors each dot by `timepoint` (Visit 1-4) instead of specificity.
#
#  Run this AFTER the main D0_data_CMV script, in the same R session
#  (needs: D0_data_CMV, base_export, save_jpeg(), clone_colors_pp65/ie1
#   from the original pipeline for clone ID bookkeeping).
# ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(purrr)

# ── 0. New output subfolder ─────────────────────────────────────────────────
tp_out_dir <- file.path(base_export, "clone_panels", "full", "timepoint_colors")
dir.create(tp_out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Timepoint color palette ──────────────────────────────────────────────
timepoint_colors <- c(
  "Visit 1" = "#4DAF4A",
  "Visit 2" = "#377EB8",
  "Visit 3" = "#FF7F00",
  "Visit 4" = "#E41A1C"
)

# ── 2. Single clone panel, dots colored by timepoint ────────────────────────
make_clone_panel_timepoint <- function(obj, clone_id, timepoint_colors) {
  emb <- as.data.frame(Embeddings(obj, "umap"))
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$clone     <- obj$Global_Clone_ID
  emb$timepoint <- obj$timepoint

  emb$highlight <- !is.na(emb$clone) & emb$clone == clone_id
  n <- sum(emb$highlight)

  emb_bg <- emb[!emb$highlight, ]
  emb_hl <- emb[emb$highlight, ]
  emb_hl$timepoint <- factor(emb_hl$timepoint, levels = names(timepoint_colors))

  ggplot() +
    geom_point(data = emb_bg, aes(UMAP_1, UMAP_2),
               color = "grey85", size = 0.3) +
    geom_point(data = emb_hl, aes(UMAP_1, UMAP_2, color = timepoint),
               size = 1.4) +
    scale_color_manual(values = timepoint_colors, drop = FALSE, name = "Timepoint") +
    ggtitle(paste0(clone_id, "\n(n = ", n, ")")) +
    theme_bw() +
    theme(
      plot.title      = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.title      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      plot.margin     = margin(2, 2, 2, 2),
      legend.position = "right"
    )
}

# ── 3. Grid builder (shared legend via guides = "collect") ─────────────────
export_clone_grid_timepoint <- function(obj, specificity, out_dir, obj_label,
                                        ncol_grid = 5, max_clone = NULL,
                                        skinny = FALSE) {
  prefix <- paste0("D0_", specificity, "_clone_")

  present <- obj$Global_Clone_ID[!is.na(obj$Global_Clone_ID) &
                                   str_starts(obj$Global_Clone_ID, prefix)]
  clone_ids <- unique(present)
  clone_num <- as.integer(str_remove(clone_ids, prefix))
  clone_ids <- clone_ids[order(clone_num)]
  clone_num <- sort(clone_num)

  if (!is.null(max_clone)) {
    keep      <- clone_num <= max_clone
    clone_ids <- clone_ids[keep]
  }

  if (length(clone_ids) == 0) {
    message(sprintf("No %s clones found in %s", specificity, obj_label))
    return(invisible(NULL))
  }

  panels <- map(clone_ids, ~ make_clone_panel_timepoint(obj, .x, timepoint_colors))

  nrow_grid <- ceiling(length(panels) / ncol_grid)
  grid <- wrap_plots(panels, ncol = ncol_grid, guides = "collect") +
    plot_annotation(
      title = sprintf("%s — %s clones by timepoint (n = %d)",
                      obj_label, specificity, length(clone_ids)),
      theme = theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
    )

  suffix   <- if (skinny) "_clone_panels_skinny_timepoint" else "_clone_panels_timepoint"
  out_path <- file.path(out_dir, sprintf("%s_%s%s.jpeg", obj_label, specificity, suffix))

  save_jpeg(grid, out_path,
            width  = ncol_grid * 2.6 + 1.2,   # extra room for shared legend
            height = nrow_grid * 2.8)

  message(sprintf("Exported %d %s clone panels (timepoint colors) for %s",
                  length(clone_ids), specificity, obj_label))
}

# ── 4. Recreate the two requested plots ─────────────────────────────────────
# full_pp65_clone_panels  -> timepoint-colored version, ncol = 5 (wide grid)
export_clone_grid_timepoint(D0_data_CMV, "pp65", tp_out_dir, "full",
                            ncol_grid = 5, max_clone = 20, skinny = FALSE)

# full_IE1_clone_panels_skinny -> timepoint-colored version, ncol = 3 (skinny grid)
export_clone_grid_timepoint(D0_data_CMV, "IE1", tp_out_dir, "full",
                            ncol_grid = 3, max_clone = 7, skinny = TRUE)





#
#
# Percent of cells in each metacluster by donor/specificity
#
#

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

# ── 0. Settings ──────────────────────────────────────────────────────────────
metacluster_order <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")

spec_colors <- c("pp65" = "blue", "IE1" = "red")

# Donor shapes, listed in numeric order for the legend
donor_shapes <- c(
  "Donor 1" = 16,  # solid circle
  "Donor 2" = 15,  # solid square
  "Donor 3" = 17,  # solid triangle (up)
  "Donor 4" = 25,  # inverted triangle (fillable)
  "Donor 5" = 18   # solid diamond
)

export_dir <- "exports/integrated_v5/metacluster_specificity_paired"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Compute per-donor × specificity metacluster percentages ───────────────
pct_df <- mdata_CMV@meta.data %>%
  dplyr::filter(
    !is.na(metacluster),
    specificity %in% c("pp65", "IE1"),
    donor_identity %in% names(donor_shapes)
  ) %>%
  dplyr::count(donor_identity, specificity, metacluster, name = "n_cells") %>%
  dplyr::group_by(donor_identity, specificity) %>%
  dplyr::mutate(percent = n_cells / sum(n_cells) * 100) %>%
  dplyr::ungroup() %>%
  tidyr::complete(
    donor_identity = names(donor_shapes),
    specificity    = c("pp65", "IE1"),
    metacluster    = metacluster_order,
    fill = list(n_cells = 0, percent = 0)
  ) %>%
  dplyr::mutate(
    metacluster    = factor(metacluster, levels = metacluster_order),
    specificity    = factor(specificity, levels = c("pp65", "IE1")),
    donor_identity = factor(donor_identity, levels = names(donor_shapes))
  )

# ── 2. Paired t-tests per metacluster, with BH-FDR correction ────────────────
wide_df <- pct_df %>%
  dplyr::select(donor_identity, specificity, metacluster, percent) %>%
  tidyr::pivot_wider(names_from = specificity, values_from = percent)

ttest_results <- wide_df %>%
  dplyr::group_by(metacluster) %>%
  dplyr::summarise(
    n_pairs = sum(!is.na(pp65) & !is.na(IE1)),
    p_value = tryCatch(
      stats::t.test(pp65, IE1, paired = TRUE)$p.value,
      error = function(e) NA_real_
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    p_label = dplyr::case_when(
      is.na(p_value)              ~ "NA",
      p_value > 0.1               ~ "ns",
      p_value >= 0.05             ~ sprintf("p = %.3g", p_value),  # 0.05–0.1
      p_value >= 0.01             ~ "*",                           # 0.01–0.05
      p_value >= 0.001            ~ "**",                          # 0.001–0.01
      TRUE                        ~ "***"                          # < 0.001
    )
  )

write.csv(ttest_results,
          file.path(export_dir, "metacluster_paired_ttest_results.csv"),
          row.names = FALSE)
print(ttest_results)

# ── 3. Position data for plotting ────────────────────────────────────────────
offset <- 0.18
pct_df <- pct_df %>%
  dplyr::mutate(
    x_base = as.numeric(metacluster),
    x_pos  = x_base + ifelse(specificity == "pp65", -offset, offset)
  )

# Connecting-line data: one row per donor × metacluster (pp65 → IE1)
line_df <- pct_df %>%
  dplyr::select(donor_identity, metacluster, specificity, x_pos, percent) %>%
  tidyr::pivot_wider(
    id_cols     = c(donor_identity, metacluster),
    names_from  = specificity,
    values_from = c(x_pos, percent)
  )

# Bracket + p-value label geometry, per metacluster
y_top <- max(pct_df$percent, na.rm = TRUE)
bracket_df <- pct_df %>%
  dplyr::group_by(metacluster) %>%
  dplyr::summarise(y_max = max(percent, na.rm = TRUE), .groups = "drop") %>%
  dplyr::mutate(
    x_base   = as.numeric(metacluster),
    x_left   = x_base - offset,
    x_right  = x_base + offset,
    y_bar    = y_max + y_top * 0.08,   # horizontal bracket bar height
    y_tick   = y_bar - y_top * 0.02,   # downward tick ends
    y_text   = y_bar + y_top * 0.04    # p-value text height
  ) %>%
  dplyr::left_join(dplyr::select(ttest_results, metacluster, p_label),
                   by = "metacluster")

# Vertical dashed separators between metacluster groups
sep_x <- seq_len(length(metacluster_order) - 1) + 0.5

# ── 4. Plot ──────────────────────────────────────────────────────────────────
paired_plot <- ggplot() +
  # Vertical dashed separators between groups
  geom_vline(xintercept = sep_x, linetype = "dashed",
             color = "grey70", linewidth = 0.4) +
  # Connecting lines per donor (pp65 → IE1)
  geom_segment(
    data = line_df,
    aes(x = x_pos_pp65, xend = x_pos_IE1,
        y = percent_pp65, yend = percent_IE1),
    color = "black", linewidth = 0.4, alpha = 0.7
  ) +
  # Points: triangles (Donors 3 & 4) slightly smaller
  geom_point(
    data = dplyr::filter(pct_df, donor_identity %in% c("Donor 3", "Donor 4")),
    aes(x = x_pos, y = percent, color = specificity,
        fill = specificity, shape = donor_identity),
    size = 2.8, alpha = 0.9, stroke = 0.6
  ) +
  geom_point(
    data = dplyr::filter(pct_df, !donor_identity %in% c("Donor 3", "Donor 4")),
    aes(x = x_pos, y = percent, color = specificity,
        fill = specificity, shape = donor_identity),
    size = 3.5, alpha = 0.9, stroke = 0.6
  ) +
  # Significance brackets (horizontal bar + two downward ticks)
  geom_segment(data = bracket_df,
               aes(x = x_left, xend = x_right, y = y_bar, yend = y_bar),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = bracket_df,
               aes(x = x_left, xend = x_left, y = y_bar, yend = y_tick),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_segment(data = bracket_df,
               aes(x = x_right, xend = x_right, y = y_bar, yend = y_tick),
               inherit.aes = FALSE, linewidth = 0.4) +
  # p-value text above each bracket
  geom_text(data = bracket_df,
            aes(x = x_base, y = y_text, label = p_label),
            inherit.aes = FALSE, size = 4.5) +
  scale_color_manual(name = "Specificity", values = spec_colors) +
  scale_fill_manual(name = "Specificity", values = spec_colors) +
  scale_shape_manual(
    name   = "Donor",
    values = donor_shapes,
    breaks = names(donor_shapes),   # force Donor 1 → Donor 5 order
    limits = names(donor_shapes)
  ) +
  scale_x_continuous(
    breaks = seq_along(metacluster_order),
    labels = metacluster_order
  ) +
  # Headroom above the tallest bracket (Tem) for the p-value text
  expand_limits(y = y_top * 1.18) +
  guides(
    shape = guide_legend(order = 2, override.aes = list(
      fill  = rep("black", 5),
      color = rep("black", 5),
      size  = c(3.5, 3.5, 2.8, 2.8, 3.5)  # matches Donor 1→5
    )),
    color = guide_legend(order = 1),
    fill  = guide_legend(order = 1)
  ) +
  labs(
    x = NULL,
    y = "Percent of specificity-defined cells (%)",
    title = "Phenotype composition: pp65 vs IE1 (paired by donor)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 16),
    legend.position = "right",
    plot.title = element_text(margin = margin(b = 18))  # space title above bracket
  )

ggsave(
  file.path(export_dir, "metacluster_specificity_paired_dotplot.png"),
  plot = paired_plot, width = 9, height = 6, dpi = 300, bg = "white"
)


#
#
# D0 seurat_cluster bubble plots
#
#

# All requested genes are RNA genes, so a single DotPlot call on the RNA assay
# replaces the two-assay stitch used in the metacluster version.
genes_rna <- c("HLA-DRA", "HLA-DRB1", "CX3CR1", "KLRC4", "TNF",
               "CCL3", "CCL4", "XCL1", "IL7R", "CD28",
               "GZMK", "LILRB1", "PLAC8", "PROK2", "ARL15",
               "WWOX", "MAL")

p1 <- DotPlot(D0_data_CMV, features = genes_rna,
              group.by = "seurat_clusters", assay = "RNA")

# Combine into the same plot_data structure as the metacluster version
plot_data <- p1$data

# Order clusters: 0 at top → 7 at bottom.
# ggplot draws the first level at the bottom, so reverse the 0→7 order
# (mirrors the rev(metacluster_order) trick in the metacluster plot).
cluster_order <- as.character(0:7)
plot_data$id <- factor(plot_data$id, levels = rev(cluster_order))

# Order genes (left → right as written)
gene_order <- c("HLA-DRA", "HLA-DRB1", "CX3CR1", "KLRC4", "TNF",
                "CCL3", "CCL4", "XCL1", "IL7R", "CD28",
                "GZMK", "LILRB1", "PLAC8", "PROK2", "ARL15",
                "WWOX", "MAL")
plot_data$features.plot <- factor(plot_data$features.plot, levels = gene_order)

# Make the bubble plot (genes on x, seurat_clusters on y)
p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank()   # removes axis titles
  )

# Save
dir.create("Donor_0/exports/bubble_plots", recursive = TRUE, showWarnings = FALSE)
ggsave("Donor_0/exports/bubble_plots/D0_CMV_bubbleplot_by_cluster_flipped_tall.png", p,
       width = 8, height = 5, dpi = 300)

# D0 CMV, by seurat_cluster — abridged

genes_rna <- c("ARL15", "WWOX", "TNF", "CCL3", "CCL4", "ZFP36",
               "HLA-DRA", "HLA-DRB1", "CX3CR1", "KLRC4", "LILRB1",
               "PLAC8", "PROK2", "IL7R", "CD28", "GZMK")

p1 <- DotPlot(D0_data_CMV, features = genes_rna,
              group.by = "seurat_clusters", assay = "RNA")

plot_data <- p1$data

# Keep only the six requested clusters
clusters_keep <- c("5", "2", "0", "1", "4", "3")
plot_data <- plot_data[plot_data$id %in% clusters_keep, ]

# Order clusters: 5 at top → 4 at bottom.
# ggplot draws the first level at the bottom, so reverse the desired order.
plot_data$id <- factor(plot_data$id, levels = rev(clusters_keep))

# Order genes (left → right as written)
plot_data$features.plot <- factor(plot_data$features.plot, levels = genes_rna)

# Make the bubble plot (genes on x, seurat_clusters on y)
p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank()   # removes axis titles
  )

# Save
dir.create("Donor_0/exports/bubble_plots", recursive = TRUE, showWarnings = FALSE)
ggsave("Donor_0/exports/bubble_plots/D0_CMV_bubbleplot_by_cluster_abridged.png", p,
       width = 8, height = 4, dpi = 300)

# D0 CMV, by seurat_cluster — abridged_2

genes_rna <- c("TNF", "CCL3", "CCL4", "XCL1", "ZFP36",
               "HLA-DRA", "HLA-DRB1", "CX3CR1", "KLRC4", "LILRB1",
               "PLAC8", "PROK2", "IL7R", "CD28", "GZMK")

p1 <- DotPlot(D0_data_CMV, features = genes_rna,
              group.by = "seurat_clusters", assay = "RNA")

plot_data <- p1$data

# Keep only the five requested clusters
clusters_keep <- c("2", "0", "1", "4", "3")
plot_data <- plot_data[plot_data$id %in% clusters_keep, ]

# Order clusters: 2 at top → 3 at bottom.
# ggplot draws the first level at the bottom, so reverse the desired order.
plot_data$id <- factor(plot_data$id, levels = rev(clusters_keep))

# Order genes (left → right as written)
plot_data$features.plot <- factor(plot_data$features.plot, levels = genes_rna)

# Make the bubble plot (genes on x, seurat_clusters on y)
p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title = element_blank()   # removes axis titles
  )

# Save
dir.create("Donor_0/exports/bubble_plots", recursive = TRUE, showWarnings = FALSE)
ggsave("Donor_0/exports/bubble_plots/D0_CMV_bubbleplot_by_cluster_abridged_2.png", p,
       width = 8, height = 4, dpi = 300)



#
#
# Full dataset UMAPs with specificity overlays
#
#

library(Seurat)
library(ggplot2)
library(dplyr)

# ---- Output directory ----
out_dir <- "exports/integrated_v5/total_umap_specificity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Helper to build + export a highlighted UMAP ----
make_spec_umap <- function(obj, spec_value, hl_color, plot_label, file_label,
                           hl_alpha = 0.35) {
  # Pull UMAP embedding + specificity into a plotting data frame
  emb <- as.data.frame(Embeddings(obj, reduction = "umap"))
  colnames(emb)[1:2] <- c("UMAP_1", "UMAP_2")
  emb$specificity <- obj$specificity

  bg <- emb
  fg <- dplyr::filter(emb, specificity == spec_value)

  p <- ggplot() +
    geom_point(data = bg, aes(UMAP_1, UMAP_2),
               color = "grey85", size = 0.3) +
    geom_point(data = fg, aes(UMAP_1, UMAP_2),
               color = hl_color, size = 0.5, alpha = hl_alpha) +
    ggtitle(plot_label) +
    theme_classic() +
    theme(
      plot.title      = element_text(hjust = 0.5, face = "bold"),
      legend.position = "none"
    ) +
    coord_equal()

  ggsave(
    filename = file.path(out_dir, paste0("UMAP_", file_label, ".png")),
    plot     = p,
    width    = 7, height = 6, dpi = 300
  )

  message("Exported: ", file_label, " (", nrow(fg), " cells)")
}

# ---- Generate the three plots ----
make_spec_umap(mdata, "Bulk", "darkgreen", "Tetramer-negative", "Bulk")
make_spec_umap(mdata, "pp65", "blue",      "pp65-specific",     "pp65")
make_spec_umap(mdata, "IE1",  "red",       "IE1-specific",      "IE1")


#
#
# CoNGA setup — export mdata_CMV with cloneId-defined clonotypes
#
#

library(Seurat); library(dplyr); library(stringr); library(readr)
library(DropletUtils)

# ── 0. Settings ──────────────────────────────────────────────────────────────
conga_dir <- "exports/conga/mdata_CMV"
dir.create(file.path(conga_dir, "gex_mtx"), recursive = TRUE, showWarnings = FALSE)

obj <- mdata_CMV   # primary CMV object (pp65 + IE1, donors 1-5)

# ── 1. Map each pool's TCR contig file to its Seurat barcode scheme ───────────
#   pool565-1_1  -> project1_<raw>
#   pool565-2_1  -> project2_1_<raw>
#   pool565-2_2  -> project2_2_<raw>
#   pool565-3_1  -> project3_donor[45]_<raw>_1   (donor from hashtag demux)
#   pool565-3_2  -> project3_donor[45]_<raw>_2
contig_config <- list(
  list(file = "scTCRseq_contig_files/pool565-1_1_filtered_contig_annotations.csv",
       seurat_regex = "^project1_([ACGT]{16}-1)$"),
  list(file = "scTCRseq_contig_files/pool565-2_1_filtered_contig_annotations.csv",
       seurat_regex = "^project2_1_([ACGT]{16}-1)$"),
  list(file = "scTCRseq_contig_files/pool565-2_2_filtered_contig_annotations.csv",
       seurat_regex = "^project2_2_([ACGT]{16}-1)$"),
  list(file = "scTCRseq_contig_files/pool565-3_1_filtered_contig_annotations.csv",
       seurat_regex = "^project3_donor[45]_([ACGT]{16}-1)_1$"),
  list(file = "scTCRseq_contig_files/pool565-3_2_filtered_contig_annotations.csv",
       seurat_regex = "^project3_donor[45]_([ACGT]{16}-1)_2$")
)

seurat_barcodes <- colnames(obj)

remap_one_pool <- function(cfg, seurat_barcodes) {
  hits    <- str_match(seurat_barcodes, cfg$seurat_regex)
  belongs <- !is.na(hits[, 1])
  if (sum(belongs) == 0) {
    message("  No Seurat cells match ", basename(cfg$file), " — skipping"); return(NULL)
  }
  raw2full <- setNames(hits[belongs, 1], hits[belongs, 2])  # raw -> full Seurat name

  df <- read_csv(cfg$file, show_col_types = FALSE)
  if ("sample" %in% colnames(df)) df <- select(df, -sample)  # pool565-3 extra col

  df <- df %>% mutate(barcode = raw2full[barcode]) %>% filter(!is.na(barcode))
  message(sprintf("  %-45s -> %d contigs / %d cells",
                  basename(cfg$file), nrow(df), n_distinct(df$barcode)))
  df
}

combined_contigs <- contig_config %>%
  lapply(remap_one_pool, seurat_barcodes = seurat_barcodes) %>%
  Filter(Negate(is.null), .) %>%
  bind_rows()

## ── 2. Keep productive, high-confidence TRA/TRB and attach clone labels ──────
meta <- obj@meta.data
meta$barcode <- rownames(meta)        # authoritative cell names (a stale 'barcode' col already exists)
meta <- meta %>% select(barcode, cloneId, Global_Clone_ID)

contigs <- combined_contigs %>%
  filter(productive == TRUE, high_confidence == TRUE,
         chain %in% c("TRA", "TRB"),
         !is.na(cdr3),    cdr3    != "", cdr3    != "None",
         !is.na(cdr3_nt), cdr3_nt != "") %>%
  left_join(meta, by = "barcode")

# ── 3. Define canonical clone sequences (canonical beta + dominant alpha) ─────
# Dominant chain of a given type within a clone: ranked by number of cells
# carrying it, ties broken by total UMIs; returns one representative contig row.
pick_dominant_chain <- function(df, chain_type) {
  sub <- df %>% filter(chain == chain_type)
  if (nrow(sub) == 0) return(NULL)
  top <- sub %>%
    group_by(v_gene, j_gene, cdr3, cdr3_nt) %>%
    summarise(n_cells = n_distinct(barcode),
              total_umis = sum(umis, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(n_cells), desc(total_umis)) %>%
    dplyr::slice(1)
  sub %>% filter(v_gene == top$v_gene, j_gene == top$j_gene,
                 cdr3 == top$cdr3, cdr3_nt == top$cdr3_nt) %>% dplyr::slice(1)
}

clustered <- contigs %>% filter(!is.na(cloneId), !grepl("notCalled", cloneId))
uncalled  <- contigs %>% filter(is.na(cloneId) |  grepl("notCalled", cloneId))

clone_ids         <- unique(clustered$cloneId)
clone_out         <- vector("list", length(clone_ids))
multi_beta_clones <- character(0)   # >1 distinct beta — review (shouldn't happen w/ maxB=1)
no_alpha_clones   <- character(0)   # no alpha recovered — CoNGA will drop these
no_beta_clones    <- character(0)   # no beta recovered  — CoNGA will drop these

for (i in seq_along(clone_ids)) {
  cid <- clone_ids[i]
  cdf <- clustered %>% filter(cloneId == cid)

  n_betas <- cdf %>% filter(chain == "TRB") %>%
    distinct(v_gene, j_gene, cdr3, cdr3_nt) %>% nrow()
  if (n_betas > 1) multi_beta_clones <- c(multi_beta_clones, cid)

  canon_tra <- pick_dominant_chain(cdf, "TRA")
  canon_trb <- pick_dominant_chain(cdf, "TRB")
  if (is.null(canon_tra)) no_alpha_clones <- c(no_alpha_clones, cid)
  if (is.null(canon_trb)) no_beta_clones  <- c(no_beta_clones,  cid)

  cells <- unique(cdf$barcode)
  rows  <- list()
  for (cr in list(canon_tra, canon_trb)) {
    if (is.null(cr)) next
    rep_rows <- cr[rep(1, length(cells)), , drop = FALSE]   # broadcast canonical chain
    rep_rows$barcode          <- cells
    rep_rows$raw_clonotype_id <- cid                        # group key = cloneId
    rep_rows$contig_id        <- paste0(cells, "_", cr$chain[1], "_conga")
    rows[[length(rows) + 1]]  <- rep_rows
  }
  clone_out[[i]] <- bind_rows(rows)
}
clustered_canonical <- bind_rows(clone_out)

# notCalled cells: each its own singleton clonotype, original sequence untouched
uncalled_out <- uncalled %>% mutate(raw_clonotype_id = paste0("uncalled-", barcode))

final_contigs <- bind_rows(clustered_canonical, uncalled_out) %>%
  select(-cloneId, -Global_Clone_ID)   # drop helper cols (not 10x format)

contig_out <- file.path(conga_dir, "conga_filtered_contig_annotations.csv")
write_csv(final_contigs, contig_out)

# Diagnostics
message(sprintf("\nClonotypes from cloneId: %d   |   uncalled singletons: %d",
                length(clone_ids), n_distinct(uncalled_out$barcode)))
message(sprintf("  clones with no alpha (dropped by CoNGA): %d", length(no_alpha_clones)))
message(sprintf("  clones with no beta  (dropped by CoNGA): %d", length(no_beta_clones)))
message(sprintf("  clones with >1 distinct beta (REVIEW):   %d", length(multi_beta_clones)))
if (length(multi_beta_clones) > 0) {
  write_csv(data.frame(cloneId = multi_beta_clones),
            file.path(conga_dir, "multi_beta_clones_REVIEW.csv"))
  message("    -> written to multi_beta_clones_REVIEW.csv")
}
message(sprintf("Wrote %d contigs / %d cells to %s",
                nrow(final_contigs), n_distinct(final_contigs$barcode), contig_out))

# ── 4. Export GEX counts to 10x mtx (CoNGA does its own normalization) ────────
DefaultAssay(obj) <- "RNA"
obj <- JoinLayers(obj)
counts <- GetAssayData(obj, assay = "RNA", layer = "counts")
write10xCounts(
  x = counts, path = file.path(conga_dir, "gex_mtx"),
  barcodes = colnames(counts), gene.id = rownames(counts),
  gene.symbol = rownames(counts), version = "3", overwrite = TRUE
)
message("Wrote GEX mtx to ", file.path(conga_dir, "gex_mtx"))

# ── 5. Export per-cell metadata ──────────────────────────────────────────────
# donor/specificity/timepoint/metacluster = CoNGA batch_keys (low cardinality);
# cloneId/Global_Clone_ID = labels for merging CoNGA results back in R later.
md <- obj@meta.data
md$barcode <- rownames(md)
md %>%
  select(barcode, donor_identity, specificity, timepoint, metacluster,
         cloneId, Global_Clone_ID) %>%
  write_csv(file.path(conga_dir, "cell_batch_info.csv"))
message("Wrote metadata to ", file.path(conga_dir, "cell_batch_info.csv"))


#
#
# CoNGA integrated rerun — export per-clone integrated PC centroids
#
#

library(Seurat); library(dplyr); library(stringr); library(readr); library(tibble)

# ── 0. Settings ──────────────────────────────────────────────────────────────
conga_dir <- "exports/conga/mdata_CMV"
int_dir   <- "exports/conga/mdata_CMV_integrated"
dir.create(file.path(int_dir, "out"), recursive = TRUE, showWarnings = FALSE)

obj   <- mdata_CMV
n_pcs <- 30

# ── 1. Per-cell integrated PCs + cloneId ─────────────────────────────────────
pcs <- Embeddings(obj, "pca")[, 1:n_pcs, drop = FALSE]
md  <- obj@meta.data
md$barcode <- rownames(md)
stopifnot(identical(rownames(pcs), md$barcode))

cell_pc <- as.data.frame(pcs) %>%
  setNames(paste0("PC", seq_len(n_pcs))) %>%
  tibble::rownames_to_column("barcode") %>%
  left_join(md %>% select(barcode, cloneId), by = "barcode")

# ── 2. Clone key matching CoNGA's raw_clonotype_id / clone_id ─────────────────
cell_pc <- cell_pc %>%
  mutate(clone_key = if_else(
    is.na(cloneId) | grepl("notCalled", cloneId),
    paste0("uncalled-", barcode),
    cloneId))

# ── 3. Centroid over all cells in a clone (singletons => own coords) ──────────
clone_pc <- cell_pc %>%
  group_by(clone_key) %>%
  summarise(across(all_of(paste0("PC", seq_len(n_pcs))), mean), .groups = "drop")

# ── 4. Verify against CoNGA's clones file ────────────────────────────────────
clones    <- read_tsv(file.path(conga_dir, "conga_clones.tsv"), show_col_types = FALSE)
conga_ids <- unique(clones$clone_id)

missing <- setdiff(conga_ids, clone_pc$clone_key)
extra   <- setdiff(clone_pc$clone_key, conga_ids)
message(sprintf("CoNGA clonotypes: %d | with PCs: %d | missing: %d | extra (dropped): %d",
                length(conga_ids),
                length(intersect(conga_ids, clone_pc$clone_key)),
                length(missing), length(extra)))
if (length(missing) > 0) {
  message("  First few missing: ", paste(head(missing, 5), collapse = ", "))
  write_csv(data.frame(clone_key = missing),
            file.path(int_dir, "clones_missing_pcs_REVIEW.csv"))
}

# ── 5. Write keyed clone × PC matrix (only clonotypes CoNGA kept) ─────────────
clone_pc_out <- clone_pc %>% filter(clone_key %in% conga_ids)
write_csv(clone_pc_out, file.path(int_dir, "clone_integrated_pcs.csv"))
message(sprintf("Wrote %d clone × %d PC rows to %s",
                nrow(clone_pc_out), n_pcs,
                file.path(int_dir, "clone_integrated_pcs.csv")))

# Per-cell integrated PCs (fallback for clones whose key doesn't match)
cell_pc_out <- as.data.frame(pcs) %>%
  setNames(paste0("PC", seq_len(n_pcs))) %>%
  tibble::rownames_to_column("barcode")
write_csv(cell_pc_out, file.path(int_dir, "cell_integrated_pcs.csv"))
message(sprintf("Wrote %d cell × %d PC rows (fallback)", nrow(cell_pc_out), n_pcs))


#
#
# CoNGA results — Stage 1: cross-donor bicluster comparison (default vs integrated)
#
#

library(dplyr); library(tidyr); library(readr); library(stringr); library(tibble)

# ── 0. Settings ──────────────────────────────────────────────────────────────
def_dir <- "exports/conga/mdata_CMV"
int_dir <- "exports/conga/mdata_CMV_integrated"
out_dir <- "exports/conga/conga_exports"          # downstream analyses live here
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj     <- mdata_CMV

# ── 1. Helper: load a run's graph_vs_graph results, attach clone_id + donors ──
# Both runs share ONE clones file (the integrated run reused the default's),
# so clone_index -> clone_id is identical across runs.
clones_file <- file.path(def_dir, "conga_clones.tsv")
clones_idx  <- read_tsv(clones_file, show_col_types = FALSE) %>%
  mutate(clone_index = dplyr::row_number() - 1L) %>%
  select(clone_index, clone_id)

load_conga_gvg <- function(run_dir, prefix) {
  read_tsv(file.path(run_dir, "out", paste0(prefix, "_graph_vs_graph.tsv")),
           show_col_types = FALSE) %>%
    left_join(clones_idx, by = "clone_index")
}

# cloneId -> donor membership
meta <- obj@meta.data
meta$barcode <- rownames(meta)
clone_donors <- meta %>%
  filter(!is.na(cloneId), !grepl("notCalled", cloneId)) %>%
  group_by(cloneId) %>%
  summarise(n_cells       = dplyr::n(),
            donors        = paste(sort(unique(donor_identity)), collapse = ","),
            n_donors      = dplyr::n_distinct(donor_identity),
            specificities = paste(sort(unique(specificity)), collapse = ","),
            .groups = "drop")

# ── 2. Load both runs ────────────────────────────────────────────────────────
gvg_def <- load_conga_gvg(def_dir, "cmv_conga") %>%
  left_join(clone_donors, by = c("clone_id" = "cloneId"))
gvg_int <- load_conga_gvg(int_dir, "cmv_conga_int") %>%
  left_join(clone_donors, by = c("clone_id" = "cloneId"))

# join-health check
cat(sprintf("NA clone_id after join — default: %d / %d | integrated: %d / %d\n",
            sum(is.na(gvg_def$clone_id)), nrow(gvg_def),
            sum(is.na(gvg_int$clone_id)), nrow(gvg_int)))

# ── 3. Per-run summary ───────────────────────────────────────────────────────
summarise_run <- function(gvg, label) {
  gvg %>% summarise(
    run              = label,
    n_hit_clones     = dplyr::n(),
    n_with_donorinfo = sum(!is.na(n_donors)),
    n_cross_donor    = sum(n_donors >= 2, na.rm = TRUE),
    n_donor_private  = sum(n_donors == 1, na.rm = TRUE),
    n_biclusters     = dplyr::n_distinct(paste(gex_cluster, tcr_cluster)))
}
comparison <- bind_rows(summarise_run(gvg_def, "default"),
                        summarise_run(gvg_int, "integrated"))
print(comparison)

# ── 4. Bicluster-level donor composition ─────────────────────────────────────
bicluster_table <- function(gvg, label) {
  gvg %>%
    group_by(gex_cluster, tcr_cluster) %>%
    summarise(run = label,
              n_clones = dplyr::n(),
              donors_in_bicluster = paste(sort(unique(unlist(
                strsplit(paste(donors, collapse=","), ",")))), collapse = ","),
              n_donors = dplyr::n_distinct(unlist(
                strsplit(paste(donors, collapse=","), ","))),
              min_conga_score = min(conga_score, na.rm = TRUE),
              .groups = "drop") %>%
    arrange(desc(n_donors), desc(n_clones))
}
bic_def <- bicluster_table(gvg_def, "default")
bic_int <- bicluster_table(gvg_int, "integrated")
cat("\n=== DEFAULT run biclusters ===\n");   print(bic_def)
cat("\n=== INTEGRATED run biclusters ===\n"); print(bic_int)

# ── 5. Robust vs integrated-only hits ────────────────────────────────────────
hits_def <- unique(gvg_def$clone_id); hits_int <- unique(gvg_int$clone_id)
robust       <- intersect(hits_def, hits_int)
integ_only   <- setdiff(hits_int, hits_def)
default_only <- setdiff(hits_def, hits_int)
cat(sprintf("\nHit clones — robust (both): %d | integrated-only: %d | default-only: %d\n",
            length(robust), length(integ_only), length(default_only)))

integ_only_xdonor <- clone_donors %>%
  filter(cloneId %in% integ_only, n_donors >= 2) %>%
  arrange(desc(n_donors), desc(n_cells))
cat("\n=== Integrated-only CROSS-DONOR hit clones (the headline) ===\n")
print(integ_only_xdonor)

# ── 6. Save downstream outputs ───────────────────────────────────────────────
write_csv(comparison,                   file.path(out_dir, "compare_run_summary.csv"))
write_csv(bind_rows(bic_def, bic_int),  file.path(out_dir, "compare_biclusters.csv"))
write_csv(integ_only_xdonor,            file.path(out_dir, "integrated_only_crossdonor_hits.csv"))
message("Wrote comparison outputs to ", out_dir)


# Direct clone-level comparison via clone_index (stable across runs)
gd <- gvg_def %>% select(clone_index, clone_id, donors, n_donors,
                         conga_score, gex_cluster, tcr_cluster) %>%
                 rename_with(~paste0(.x, "_def"), -clone_index)
gi <- gvg_int %>% select(clone_index, clone_id, donors, n_donors,
                         conga_score, gex_cluster, tcr_cluster) %>%
                 rename_with(~paste0(.x, "_int"), -clone_index)

joined <- full_join(gd, gi, by = "clone_index") %>%
  mutate(in_default    = !is.na(conga_score_def),
         in_integrated = !is.na(conga_score_int),
         status = case_when(in_default & in_integrated ~ "both",
                            in_default                  ~ "default_only",
                            in_integrated               ~ "integrated_only"))

cat("Hit-clone overlap:\n"); print(table(joined$status))

# For clones that are hits in EITHER run, show donor composition + both scores
both_view <- joined %>%
  mutate(clone_id  = coalesce(clone_id_def, clone_id_int),
         donors    = coalesce(donors_def, donors_int),
         n_donors  = coalesce(n_donors_def, n_donors_int)) %>%
  select(clone_index, clone_id, donors, n_donors, status,
         conga_score_def, conga_score_int) %>%
  arrange(status, desc(n_donors))
cat("\nAll hit clones across both runs:\n")
print(both_view, n = 40)

# The key question: of the DEFAULT cross-donor biclusters, are those clones'
# scores preserved or destroyed in the integrated run?
cat("\nDefault-run hit clones — score in default vs integrated:\n")
joined %>%
  filter(in_default) %>%
  mutate(clone_id = coalesce(clone_id_def, clone_id_int)) %>%
  select(clone_id, donors_def, n_donors_def,
         conga_score_def, conga_score_int) %>%
  arrange(desc(n_donors_def)) %>%
  print(n = 40)

write_csv(both_view, "exports/conga/conga_exports/clone_level_run_comparison.csv")


#
#
# CoNGA — bicluster membership: TCR sequences + Global_Clone_ID per bicluster
#
#

library(dplyr); library(tidyr); library(readr); library(stringr)

out_dir <- "exports/conga/conga_exports"
obj     <- mdata_CMV

# ── 1. cloneId -> Global_Clone_ID (+ donor/specificity) from Seurat ──────────
meta <- obj@meta.data
meta$barcode <- rownames(meta)
clone_labels <- meta %>%
  filter(!is.na(cloneId), !grepl("notCalled", cloneId)) %>%
  group_by(cloneId) %>%
  summarise(
    lab_donors        = paste(sort(unique(donor_identity)), collapse = ","),
    lab_specificities = paste(sort(unique(specificity)), collapse = ","),
    global_clone_ids  = {
      g <- unique(Global_Clone_ID[!is.na(Global_Clone_ID)])
      if (length(g) == 0) NA_character_ else paste(g, collapse = ",")
    },
    n_cells = dplyr::n(),
    .groups = "drop")

# ── 2. Builder: strip gvg to CoNGA-native cols, then attach labels fresh ─────
bicluster_members <- function(gvg, label) {
  gvg %>%
    select(clone_id, gex_cluster, tcr_cluster, conga_score,
           va, ja, cdr3a, vb, jb, cdr3b) %>%          # CoNGA-native only
    left_join(clone_labels, by = c("clone_id" = "cloneId")) %>%
    transmute(
      run = label,
      bicluster = paste0("gex", gex_cluster, "_tcr", tcr_cluster),
      gex_cluster, tcr_cluster,
      clone_id,
      global_clone_ids,
      donor = lab_donors,
      specificity = lab_specificities,
      conga_score,
      va, ja, cdr3a, vb, jb, cdr3b,
      n_cells) %>%
    arrange(bicluster, conga_score, donor)
}

memb_def <- bicluster_members(gvg_def, "default")
memb_int <- bicluster_members(gvg_int, "integrated")

# ── 3. Per-bicluster convergence summary ─────────────────────────────────────
convergence_summary <- function(memb, label) {
  memb %>%
    group_by(run, bicluster, gex_cluster, tcr_cluster) %>%
    summarise(
      n_clones         = dplyr::n(),
      donors           = paste(sort(unique(unlist(strsplit(paste(donor, collapse=","), ",")))),
                               collapse = ","),
      n_donors         = dplyr::n_distinct(unlist(strsplit(paste(donor, collapse=","), ","))),
      n_distinct_vb    = dplyr::n_distinct(vb),
      n_distinct_cdr3b = dplyr::n_distinct(cdr3b),
      cdr3b_len_range  = paste0(min(nchar(cdr3b)), "-", max(nchar(cdr3b))),
      shared_vb        = paste(names(which(table(vb) > 1)), collapse = ","),
      n_with_global    = sum(!is.na(global_clone_ids)),
      global_clone_ids = paste(unique(global_clone_ids[!is.na(global_clone_ids)]),
                               collapse = "; "),
      min_score        = min(conga_score),
      .groups = "drop") %>%
    arrange(desc(n_donors), min_score)
}

conv_def <- convergence_summary(memb_def, "default")
conv_int <- convergence_summary(memb_int, "integrated")

cat("\n=== DEFAULT run — bicluster convergence summary ===\n");    print(conv_def, width = Inf)
cat("\n=== INTEGRATED run — bicluster convergence summary ===\n");  print(conv_int, width = Inf)

# ── 4. Full per-clone TCRs for default cross-donor biclusters ────────────────
cat("\n=== DEFAULT cross-donor biclusters — full member TCRs ===\n")
xdonor_bics <- conv_def %>% filter(n_donors >= 2) %>% pull(bicluster)
memb_def %>%
  filter(bicluster %in% xdonor_bics) %>%
  select(bicluster, clone_id, global_clone_ids, donor, specificity,
         conga_score, vb, jb, cdr3b, va, ja, cdr3a) %>%
  print(n = 60, width = Inf)

# ── 5. Save ──────────────────────────────────────────────────────────────────
write_csv(memb_def, file.path(out_dir, "bicluster_members_default.csv"))
write_csv(memb_int, file.path(out_dir, "bicluster_members_integrated.csv"))
write_csv(bind_rows(conv_def, conv_int),
          file.path(out_dir, "bicluster_convergence_summary.csv"))
message("Wrote bicluster membership + convergence tables to ", out_dir)



#
#
# CoNGA default biclusters — Seurat UMAP highlight, colored by cloneId
#
#

library(Seurat); library(dplyr); library(ggplot2); library(readr); library(stringr)

out_dir <- "exports/conga/conga_exports/bicluster_umaps"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj <- mdata_CMV

# ── 0. Bicluster -> member cloneIds (from default gvg, all biclusters) ────────
memb <- gvg_def %>%
  transmute(bicluster = paste0("gex", gex_cluster, "_tcr", tcr_cluster),
            cloneId = clone_id) %>%
  distinct() %>%
  filter(!is.na(cloneId))

all_bics <- memb %>% distinct(bicluster) %>% pull(bicluster)
message("Plotting ", length(all_bics), " biclusters: ", paste(all_bics, collapse = ", "))

# ── 1. Base UMAP coords + metadata, once ─────────────────────────────────────
umap_xy <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap_xy) <- c("UMAP_1", "UMAP_2")
md <- obj@meta.data
umap_xy$cloneId <- md$cloneId

# ── 2. One plot per bicluster ────────────────────────────────────────────────
for (bic in all_bics) {
  member_clones <- memb %>% filter(bicluster == bic) %>% pull(cloneId) %>% unique()

  df <- umap_xy %>%
    mutate(highlight = if_else(cloneId %in% member_clones, cloneId, NA_character_))

  bg <- df %>% filter(is.na(highlight))
  fg <- df %>% filter(!is.na(highlight)) %>%
          mutate(highlight = factor(highlight, levels = member_clones))

  pal <- setNames(scales::hue_pal()(length(member_clones)), member_clones)

  p <- ggplot() +
    geom_point(data = bg, aes(UMAP_1, UMAP_2),
               color = "grey88", size = 0.4, alpha = 0.6) +
    geom_point(data = fg, aes(UMAP_1, UMAP_2, color = highlight),
               size = 1.8, alpha = 0.9) +
    scale_color_manual(values = pal, name = "cloneId") +
    labs(title = paste0("Default CoNGA bicluster: ", bic),
         subtitle = paste0(length(member_clones),
                           if_else(length(member_clones) == 1, " clone", " clones"))) +
    theme_classic() +
    theme(legend.position = "right",
          plot.title = element_text(face = "bold")) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  n_cells_hl <- sum(!is.na(df$highlight))
  ggsave(file.path(out_dir, paste0("bicluster_", bic, "_umap.png")),
         p, width = 8, height = 6, dpi = 300)
  message(sprintf("  %s: %d clones, %d cells highlighted",
                  bic, length(member_clones), n_cells_hl))
}
message("Wrote bicluster UMAPs to ", out_dir)




#
#
# TCRdist-cluster phenotypic coherence — N_PERM=1000, size-filtered, annotated
#
#
library(Seurat); library(dplyr); library(tidyr); library(FNN); library(readr)

out_dir <- "exports/conga/conga_exports/tcrdist_phenotype_test"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj <- mdata_CMV

set.seed(1)
K           <- 100
N_PERM      <- 10000
N_PCS       <- 30
OVERFETCH   <- 1500
MIN_CELLS   <- 20     # cluster cell-count floor
MIN_CLONES  <- 3      # cluster distinct-clone floor (the important one)

# ── 0. Per-cell table ────────────────────────────────────────────────────────
pcs <- Embeddings(obj, "pca")[, 1:N_PCS, drop = FALSE]
md  <- obj@meta.data
cells <- tibble(
  barcode = colnames(obj),
  cluster = as.character(md$TCR_cluster_et64),
  cloneId = as.character(md$cloneId),
  donor   = as.character(md$donor_identity),
  specificity = as.character(md$specificity),
  global  = as.character(md$Global_Clone_ID)
)
keep <- !is.na(cells$cluster) & !is.na(cells$cloneId) & !grepl("notCalled", cells$cloneId)
cells <- cells[keep, ]; pcs <- pcs[keep, , drop = FALSE]

# ── 1. Cluster composition + size filter ─────────────────────────────────────
cluster_info <- cells %>% group_by(cluster) %>%
  dplyr::summarise(
    n_cells  = dplyr::n(),
    n_clones = dplyr::n_distinct(cloneId),
    n_donors = dplyr::n_distinct(donor),
    donors   = paste(sort(unique(donor)), collapse = ","),
    specificities = paste(sort(unique(specificity)), collapse = ","),
    global_clone_ids = {
      g <- unique(global[!is.na(global)])
      if (length(g) == 0) NA_character_ else paste(sort(g), collapse = "; ")
    },
    n_global_clones = dplyr::n_distinct(global[!is.na(global)]),
    .groups = "drop")

keep_clusters <- cluster_info %>%
  filter(n_cells >= MIN_CELLS, n_clones >= MIN_CLONES) %>% pull(cluster)
message(sprintf("Clusters passing filter (>=%d cells, >=%d clones): %d of %d",
                MIN_CELLS, MIN_CLONES, length(keep_clusters), nrow(cluster_info)))

# IMPORTANT: kNN is still computed over ALL cells (neighbors can be any cell);
# we only restrict which clusters get TESTED + corrected.
# ── 2. kNN with same-clone exclusion (over all cells) ────────────────────────
OVERFETCH <- min(OVERFETCH, nrow(cells) - 1)
knn <- get.knn(pcs, k = OVERFETCH)$nn.index
nbr_clone  <- matrix(cells$cloneId[knn], nrow = nrow(knn))
nbr_donor  <- matrix(cells$donor[knn],   nrow = nrow(knn))
same_clone <- nbr_clone == cells$cloneId
nbr_xdonor <- nbr_donor != cells$donor

surv_idx <- vector("list", nrow(cells)); n_dropped <- 0L
for (r in seq_len(nrow(cells))) {
  surv <- which(!same_clone[r, ])
  if (length(surv) < K) n_dropped <- n_dropped + 1L
  surv_idx[[r]] <- if (length(surv) == 0) integer(0) else surv[seq_len(min(K, length(surv)))]
}
message(sprintf("Cells with <K non-same-clone neighbors: %d (%.1f%%)",
                n_dropped, 100*n_dropped/nrow(cells)))

compute_stats <- function(cluster_vec) {
  nbr_clu <- matrix(cluster_vec[knn], nrow = nrow(knn))
  n_same <- numeric(nrow(cells)); n_xdon <- numeric(nrow(cells))
  for (r in seq_len(nrow(cells))) {
    s <- surv_idx[[r]]
    if (length(s) == 0) { n_same[r] <- NA; n_xdon[r] <- NA; next }
    sc <- nbr_clu[r, s] == cluster_vec[r]
    n_same[r] <- sum(sc); n_xdon[r] <- sum(sc & nbr_xdonor[r, s])
  }
  list(n_same = n_same, n_xdon = n_xdon)
}

# ── 3. Aggregators (only over kept clusters) ─────────────────────────────────
agg_cell <- function(df, stat) {
  df %>% filter(cluster %in% keep_clusters) %>% group_by(cluster) %>%
    dplyr::summarise(score = mean(.data[[stat]], na.rm = TRUE), .groups = "drop")
}
agg_clone <- function(df, stat) {
  df %>% filter(cluster %in% keep_clusters) %>% group_by(cluster, cloneId) %>%
    dplyr::summarise(cm = mean(.data[[stat]], na.rm = TRUE), .groups = "drop") %>%
    group_by(cluster) %>%
    dplyr::summarise(score = mean(cm, na.rm = TRUE), .groups = "drop")
}

perm_unstrat <- function() sample(cells$cluster)
perm_strat   <- function() {
  out <- cells$cluster
  for (d in unique(cells$donor)) { ix <- which(cells$donor == d); out[ix] <- sample(cells$cluster[ix]) }
  out
}

run_null <- function(stat, aggregator, perm_fn, obs_cells) {
  obs_tbl <- aggregator(obs_cells, stat) %>% dplyr::rename(obs_score = score)
  clusters <- obs_tbl$cluster
  null_mat <- matrix(NA_real_, nrow = N_PERM, ncol = length(clusters),
                     dimnames = list(NULL, clusters))
  for (p in seq_len(N_PERM)) {
    pl <- perm_fn(); st <- compute_stats(pl)
    tmp <- cells; tmp[[stat]] <- st[[stat]]; tmp$cluster <- pl
    a <- agg_clone(tmp, stat)
    null_mat[p, ] <- a$score[match(clusters, a$cluster)]
    if (p %% 1000 == 0) message(sprintf("    perm %d / %d", p, N_PERM))   # <- add this line
  }
  pvals <- sapply(seq_along(clusters), function(j)
    (1 + sum(null_mat[, j] >= obs_tbl$obs_score[j], na.rm = TRUE)) / (N_PERM + 1))
  obs_tbl %>% mutate(null_mean = colMeans(null_mat, na.rm = TRUE),
                     effect = obs_score - null_mean, p_emp = pvals,
                     p_adj = p.adjust(pvals, "BH"))
}

# ── 4. Observed + run all 8 ──────────────────────────────────────────────────
obs <- compute_stats(cells$cluster)
cells$n_same <- obs$n_same; cells$n_xdon <- obs$n_xdon

combos <- tidyr::crossing(stat = c("n_same","n_xdon"),
                          weighting = c("cell","clone"),
                          null_type = c("unstratified","stratified"))
all_res <- list()
for (i in seq_len(nrow(combos))) {
  cc <- combos[i, ]
  aggr <- if (cc$weighting=="cell") agg_cell else agg_clone
  pfn  <- if (cc$null_type=="unstratified") perm_unstrat else perm_strat
  message(sprintf("[%d/8] %s | %s | %s", i, cc$stat, cc$weighting, cc$null_type))
  all_res[[i]] <- run_null(cc$stat, aggr, pfn, cells) %>%
    mutate(metric = ifelse(cc$stat=="n_same","same_cluster","cross_donor"),
           weighting = cc$weighting, null_type = cc$null_type)
}
results <- bind_rows(all_res) %>%
  left_join(cluster_info, by = "cluster") %>%
  arrange(metric, weighting, null_type, p_adj)

# ── 5. Save ──────────────────────────────────────────────────────────────────
print(results, n = Inf, width = Inf)
write_csv(results, file.path(out_dir, "tcrdist_phenotype_results_n1000.csv"))
write_csv(cluster_info, file.path(out_dir, "cluster_composition_annotated.csv"))
write_csv(cells %>% select(barcode, cluster, cloneId, donor, specificity, global, n_same, n_xdon),
          file.path(out_dir, "per_cell_statistics_n1000.csv"))
message("Wrote n=1000 results to ", out_dir)



#
#
# Cluster cards — kNN bias-test variant (SELF-CONTAINED, et64 / mdata_CMV)
# Layout (2×2): top-left text | top-right UMAP | bottom-left TCR logo |
#               bottom-right = three-row same-cluster kNN bias-test summary.
# All helpers are namespaced (kc_) so this block runs independently of any
# previously-sourced cluster-card block.
#
#
library(Seurat); library(dplyr); library(tidyr); library(stringr)
library(ggplot2); library(cowplot); library(ggseqlogo); library(readr)

# ── 0. Settings ──────────────────────────────────────────────────────────────
knn_export_dir <- "exports/TCR_clusters/cluster_cards/kNN_bias_test"
dir.create(knn_export_dir, recursive = TRUE, showWarnings = FALSE)

# Main readout: same_cluster metric, clone-weighted, unstratified null.
# (Set WEIGHTING <- "cell" for the cell-weighted effect/p instead.)
res_file  <- "exports/conga/conga_exports/tcrdist_phenotype_test/tcrdist_phenotype_results_n1000.csv"
WEIGHTING <- "clone"

knn_lookup <- read_csv(res_file, show_col_types = FALSE) %>%
  dplyr::filter(metric == "same_cluster",
                weighting == WEIGHTING,
                null_type == "unstratified") %>%
  dplyr::transmute(cluster = as.character(cluster), effect, p_adj)
if (nrow(knn_lookup) == 0)
  stop("No rows for metric=same_cluster, weighting=", WEIGHTING,
       ", null_type=unstratified in ", res_file)

# ── 1. Palette + AA-legend helpers ───────────────────────────────────────────
KC_QUAL_PALETTE <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#386cb0", "#f0027f", "#bf5b17",
  "#4daf4a", "#e41a1c", "#984ea3", "#ff7f00", "#a65628",
  "#f781bf", "#999999", "#8dd3c7", "#fb8072", "#80b1d3"
)

kc_assign_palette <- function(groups) {
  n <- length(groups)
  colors <- KC_QUAL_PALETTE[seq_len(min(n, length(KC_QUAL_PALETTE)))]
  if (n > length(KC_QUAL_PALETTE)) colors <- colorRampPalette(KC_QUAL_PALETTE)(n)
  setNames(colors, groups)
}

KC_AA_CHEMISTRY <- data.frame(
  letter = c("G","S","T","Y","C","N","Q","K","R","H","D","E","P","A","W","F","L","I","M","V"),
  group  = c(rep("Polar", 5), rep("Neutral", 2), rep("Basic", 3),
             rep("Acidic", 2), rep("Hydrophobic", 8)),
  col    = c(rep("#109648", 5), rep("#5E239D", 2), rep("#255C99", 3),
             rep("#D62839", 2), rep("#221E22", 8)),
  stringsAsFactors = FALSE
)

kc_export_aa_legend <- function(present_letters, out_file) {
  leg_df <- KC_AA_CHEMISTRY %>% dplyr::filter(letter %in% present_letters)
  if (nrow(leg_df) == 0) return(invisible(NULL))
  group_levels <- c("Polar", "Neutral", "Basic", "Acidic", "Hydrophobic")
  group_summary <- leg_df %>%
    dplyr::group_by(group) %>%
    dplyr::summarise(col = dplyr::first(col),
                     aas = paste(base::sort(letter), collapse = ", "), .groups = "drop") %>%
    dplyr::mutate(group = factor(group, levels = group_levels)) %>%
    dplyr::arrange(group) %>%
    dplyr::mutate(label = sprintf("%s (%s)", group, aas))
  color_map <- setNames(group_summary$col, group_summary$label)
  dummy <- data.frame(x = seq_len(nrow(group_summary)), y = 1,
                      label = factor(group_summary$label, levels = group_summary$label))
  leg_plot <- ggplot(dummy, aes(x = x, y = y, fill = label)) +
    geom_tile() +
    scale_fill_manual(values = color_map, name = "Amino acid chemistry") +
    guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    theme_void() +
    theme(legend.position = "right",
          legend.title = element_text(size = 11, face = "bold"),
          legend.text  = element_text(size = 10),
          legend.key.size = unit(0.6, "cm"))
  legend_grob <- cowplot::get_legend(leg_plot)
  ggsave(filename = out_file, plot = cowplot::ggdraw(legend_grob),
         width = 3.5, height = 3, dpi = 150, bg = "white")
}

# ── 2. Per-cluster summary (et64 schema: cluster, clone_id, donor, spec, cdr3, metacluster) ──
kc_build_cluster_summary <- function(seurat_obj, meta_col, min_cells = 10) {
  md <- seurat_obj@meta.data %>%
    dplyr::filter(!is.na(.data[[meta_col]])) %>%
    dplyr::transmute(
      cluster     = as.character(.data[[meta_col]]),
      clone_id    = Global_Clone_ID,
      donor       = paste0("D", str_extract(donor_identity, "\\d")),
      spec        = specificity,
      cdr3        = t_cdr3s_aa,
      metacluster = metacluster
    )
  eligible <- md %>%
    dplyr::count(cluster, name = "total_cells") %>%
    dplyr::filter(total_cells > min_cells) %>%
    dplyr::arrange(dplyr::desc(total_cells))
  list(md = md, eligible = eligible)
}

# ── 3. Top-left: text panel ──────────────────────────────────────────────────
kc_make_text_panel <- function(cluster_id, cluster_md) {
  total <- nrow(cluster_md)
  lines <- c(sprintf("Cluster %s", cluster_id), sprintf("n = %d cells", total), "")
  for (d in base::sort(base::unique(cluster_md$donor))) {
    d_md <- cluster_md %>% dplyr::filter(donor == d)
    lines <- c(lines, sprintf("%s  (n = %d)", d, nrow(d_md)))
    named <- d_md %>%
      dplyr::filter(!is.na(clone_id)) %>%
      dplyr::count(clone_id, spec, name = "n_cells") %>%
      dplyr::mutate(clone_num = as.integer(str_extract(clone_id, "\\d+$"))) %>%
      dplyr::arrange(spec, clone_num)
    for (i in seq_len(nrow(named)))
      lines <- c(lines, sprintf("  %s:  %d",
                                str_remove(named$clone_id[i], paste0("^", d, "_")),
                                named$n_cells[i]))
    n_minor <- d_md %>% dplyr::filter(is.na(clone_id)) %>% nrow()
    if (n_minor > 0) lines <- c(lines, sprintf("  minor:  %d", n_minor))
    lines <- c(lines, "")
  }
  ggplot() +
    annotate("text", x = 0, y = 1, label = paste(lines, collapse = "\n"),
             hjust = 0, vjust = 1, size = 3.2, family = "mono", lineheight = 1.4) +
    xlim(0, 1) + ylim(0, 1) + theme_void() +
    theme(plot.margin = margin(10, 0, 10, 10))
}

# ── 4. Top-right: UMAP panel ─────────────────────────────────────────────────
kc_make_umap_panel <- function(seurat_obj, cluster_id, meta_col, cluster_md, color_map) {
  umap_df      <- as.data.frame(Embeddings(seurat_obj, reduction = "umap"))
  umap_cols    <- colnames(umap_df)
  umap_df$cell <- rownames(umap_df)
  md_full <- seurat_obj@meta.data %>%
    dplyr::mutate(
      donor    = paste0("D", str_extract(donor_identity, "\\d")),
      group_id = dplyr::case_when(
        !is.na(Global_Clone_ID) ~ as.character(Global_Clone_ID),
        !is.na(.data[[meta_col]]) & .data[[meta_col]] == as.character(cluster_id) ~
          paste0(donor, "_minor-clones"),
        TRUE ~ NA_character_))
  in_cluster <- !is.na(seurat_obj@meta.data[[meta_col]]) &
    seurat_obj@meta.data[[meta_col]] == as.character(cluster_id)
  umap_df$group_id <- NA_character_
  umap_df$group_id[in_cluster] <- md_full$group_id[in_cluster]
  umap_df$group_id <- factor(umap_df$group_id, levels = names(color_map))
  ggplot() +
    geom_point(data = umap_df[!in_cluster, ],
               aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]]),
               color = "lightgrey", size = 0.3) +
    geom_point(data = umap_df[in_cluster, ],
               aes(x = .data[[umap_cols[1]]], y = .data[[umap_cols[2]]], color = group_id),
               size = 1.4) +
    scale_color_manual(values = color_map, name = "Group", drop = FALSE,
                       guide = guide_legend(override.aes = list(size = 3), ncol = 1)) +
    theme_void() +
    theme(legend.text = element_text(size = 7),
          legend.title = element_text(size = 8, face = "bold"),
          legend.key.size = unit(0.4, "cm"),
          plot.margin = margin(5, 5, 5, 5))
}

# ── 5. Bottom-left: TCR CDR3 logo panel (TRA over TRB, modal length) ─────────
kc_extract_cdr3 <- function(cdr3_vec, chain_prefix) {
  seqs <- unlist(lapply(cdr3_vec, function(x) {
    if (is.na(x)) return(character(0))
    parts <- strsplit(x, ";")[[1]]
    sub(paste0("^", chain_prefix), "", parts[grepl(paste0("^", chain_prefix), parts)])
  }))
  seqs[nchar(seqs) > 0]
}
kc_filter_modal_length <- function(seqs) {
  if (length(seqs) == 0) return(character(0))
  lens <- nchar(seqs)
  modal_len <- as.integer(names(base::sort(table(lens), decreasing = TRUE))[1])
  seqs[lens == modal_len]
}
kc_make_logo_panel <- function(cdr3_vec, weight_mode = "unique") {
  pull_seqs <- function(prefix) {
    s <- kc_extract_cdr3(cdr3_vec, prefix)
    if (weight_mode == "unique") s <- base::unique(s)
    kc_filter_modal_length(s)
  }
  count_label <- if (weight_mode == "unique") "unique" else "cells"
  build_one <- function(seqs, title) {
    if (length(seqs) < 1)
      return(ggplot() + annotate("text", x = 0.5, y = 0.5,
                                 label = sprintf("%s: no\nsequences", title), size = 3) +
               theme_void())
    ggseqlogo(seqs, method = if (length(seqs) == 1) "prob" else "bits") +
      ggtitle(sprintf("%s  (%d %s, len = %d)", title, length(seqs), count_label, nchar(seqs[1]))) +
      theme(plot.title = element_text(size = 8, face = "bold"),
            axis.text = element_text(size = 6), axis.title = element_text(size = 7),
            legend.position = "none")
  }
  cowplot::plot_grid(build_one(pull_seqs("TRA:"), "TRA"),
                     build_one(pull_seqs("TRB:"), "TRB"), ncol = 1)
}

# ── 6. Bottom-right: NEW three-row kNN bias-test panel ───────────────────────
kc_make_knn_test_panel <- function(cluster_id, lookup) {
  row <- lookup[lookup$cluster == as.character(cluster_id), ]
  if (nrow(row) == 0) {
    eff_txt <- "Effect size: not tested"; p_txt <- "p (BH-adj): not tested"
  } else {
    eff_txt <- sprintf("Effect size: %.2f", row$effect[1])
    p_txt   <- if (row$p_adj[1] < 1e-3) sprintf("p (BH-adj) = %.1e", row$p_adj[1])
               else                     sprintf("p (BH-adj) = %.3f", row$p_adj[1])
  }
  ggplot() +
    annotate("text", x = 0.5, y = 0.66, label = "Same-cluster kNN bias test",
             fontface = "bold", size = 4.4) +
    annotate("text", x = 0.5, y = 0.46, label = eff_txt, size = 4.0) +
    annotate("text", x = 0.5, y = 0.30, label = p_txt,   size = 4.0) +
    xlim(0, 1) + ylim(0, 1) + theme_void() +
    theme(plot.margin = margin(5, 5, 5, 5))
}

# ── 7. Card assembler ────────────────────────────────────────────────────────
kc_save_card <- function(seurat_obj, cluster_id, meta_col, cluster_md,
                         export_dir, lookup, weight_mode = "unique") {
  named_groups <- cluster_md %>% dplyr::filter(!is.na(clone_id)) %>%
    dplyr::pull(clone_id) %>% base::unique() %>% base::sort()
  minor_groups <- cluster_md %>% dplyr::filter(is.na(clone_id)) %>%
    dplyr::pull(donor) %>% base::unique() %>% base::sort() %>% paste0("_minor-clones")
  color_map <- kc_assign_palette(c(named_groups, minor_groups))

  top_row    <- cowplot::plot_grid(
    kc_make_text_panel(cluster_id, cluster_md),
    kc_make_umap_panel(seurat_obj, cluster_id, meta_col, cluster_md, color_map),
    nrow = 1, rel_widths = c(0.3, 0.7))
  bottom_row <- cowplot::plot_grid(
    kc_make_logo_panel(cluster_md$cdr3, weight_mode = weight_mode),
    kc_make_knn_test_panel(cluster_id, lookup),
    nrow = 1, rel_widths = c(0.45, 0.55))
  card <- cowplot::plot_grid(top_row, bottom_row, ncol = 1, rel_heights = c(0.55, 0.45))

  ggsave(filename = file.path(export_dir, sprintf("cluster_%s_card.png", cluster_id)),
         plot = card, width = 10, height = 10, dpi = 150, bg = "white")
}

# ── 8. Run across all eligible clusters ──────────────────────────────────────
res      <- kc_build_cluster_summary(mdata_CMV, "TCR_cluster_et64", min_cells = 10)
md_all   <- res$md
eligible <- res$eligible
present_aa <- character(0)
message(sprintf("[kNN cards] %d clusters with >10 cells", nrow(eligible)))

for (i in seq_len(nrow(eligible))) {
  cl_id <- eligible$cluster[i]
  cl_md <- md_all %>% dplyr::filter(cluster == cl_id)

  tra <- kc_filter_modal_length(base::unique(kc_extract_cdr3(cl_md$cdr3, "TRA:")))
  trb <- kc_filter_modal_length(base::unique(kc_extract_cdr3(cl_md$cdr3, "TRB:")))
  if (length(tra) >= 1) present_aa <- base::union(present_aa, unlist(strsplit(tra, "")))
  if (length(trb) >= 1) present_aa <- base::union(present_aa, unlist(strsplit(trb, "")))

  kc_save_card(mdata_CMV, cl_id, "TCR_cluster_et64", cl_md,
               knn_export_dir, knn_lookup, weight_mode = "unique")
  if (i %% 10 == 0) message(sprintf("[kNN cards] %d / %d done", i, nrow(eligible)))
}

kc_export_aa_legend(present_aa, file.path(knn_export_dir, "aa_chemistry_legend.png"))
message(sprintf("[kNN cards] All cards saved to %s", knn_export_dir))



#
#
# CoNGA INTEGRATED-run biclusters — Seurat UMAP highlight, colored by cloneId
#   Mirror of the default-run bicluster UMAPs, reading the integrated run instead.
#
#
library(Seurat); library(dplyr); library(ggplot2); library(readr); library(stringr)

int_dir <- "exports/conga/mdata_CMV_integrated"
def_dir <- "exports/conga/mdata_CMV"          # clones file lives here; shared by both runs
out_dir <- "exports/conga/conga_exports/bicluster_umaps_integrated"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj <- mdata_CMV

# ── 0. Load integrated graph_vs_graph, attach clone_id (0-based positional) ───
gvg_int <- read_tsv(file.path(int_dir, "out", "cmv_conga_int_graph_vs_graph.tsv"),
                    show_col_types = FALSE)
clones  <- read_tsv(file.path(def_dir, "conga_clones.tsv"), show_col_types = FALSE) %>%
  dplyr::mutate(clone_index = dplyr::row_number() - 1L)
gvg_int <- gvg_int %>%
  dplyr::left_join(clones %>% dplyr::select(clone_index, clone_id), by = "clone_index")
cat(sprintf("NA clone_id after join: %d / %d\n",
            sum(is.na(gvg_int$clone_id)), nrow(gvg_int)))   # expect 0 or very small

# ── 1. Bicluster -> member cloneIds (all biclusters, incl. single-clone) ──────
memb <- gvg_int %>%
  dplyr::transmute(bicluster = paste0("gex", gex_cluster, "_tcr", tcr_cluster),
                   cloneId = clone_id) %>%
  dplyr::distinct() %>%
  dplyr::filter(!is.na(cloneId))
all_bics <- memb %>% dplyr::distinct(bicluster) %>% dplyr::pull(bicluster)
message("Plotting ", length(all_bics), " integrated biclusters: ",
        paste(all_bics, collapse = ", "))

# ── 2. Base UMAP coords + metadata, once ─────────────────────────────────────
umap_xy <- as.data.frame(Embeddings(obj, "umap"))
colnames(umap_xy) <- c("UMAP_1", "UMAP_2")
umap_xy$cloneId <- obj@meta.data$cloneId

# ── 3. One plot per bicluster ────────────────────────────────────────────────
for (bic in all_bics) {
  member_clones <- memb %>% dplyr::filter(bicluster == bic) %>%
    dplyr::pull(cloneId) %>% base::unique()
  df <- umap_xy %>%
    dplyr::mutate(highlight = dplyr::if_else(cloneId %in% member_clones,
                                             cloneId, NA_character_))
  bg <- df %>% dplyr::filter(is.na(highlight))
  fg <- df %>% dplyr::filter(!is.na(highlight)) %>%
    dplyr::mutate(highlight = factor(highlight, levels = member_clones))
  pal <- setNames(scales::hue_pal()(length(member_clones)), member_clones)

  p <- ggplot() +
    geom_point(data = bg, aes(UMAP_1, UMAP_2),
               color = "grey88", size = 0.4, alpha = 0.6) +
    geom_point(data = fg, aes(UMAP_1, UMAP_2, color = highlight),
               size = 1.8, alpha = 0.9) +
    scale_color_manual(values = pal, name = "cloneId") +
    labs(title = paste0("Integrated CoNGA bicluster: ", bic),
         subtitle = paste0(length(member_clones),
                           dplyr::if_else(length(member_clones) == 1, " clone", " clones"))) +
    theme_classic() +
    theme(legend.position = "right", plot.title = element_text(face = "bold")) +
    guides(color = guide_legend(override.aes = list(size = 3)))

  n_cells_hl <- sum(!is.na(df$highlight))
  ggsave(file.path(out_dir, paste0("bicluster_", bic, "_umap.png")),
         p, width = 8, height = 6, dpi = 300)
  message(sprintf("  %s: %d clones, %d cells highlighted",
                  bic, length(member_clones), n_cells_hl))
}
message("Wrote integrated bicluster UMAPs to ", out_dir)



#
#
# TCR-cluster phenotypic coherence test — specificity-stratified
#   Metric 1 (same_cluster): null permuted within SPECIFICITY only.
#   Metric 2 (cross_donor):  null permuted within SPECIFICITY × DONOR.
#   Both: clone-weighted, same-clone neighbors excluded, candidates restricted
#   to the focal cell's cluster-specificity. kNN computed within each specificity.
#
#
library(Seurat); library(dplyr); library(tidyr); library(FNN); library(readr)

out_dir <- "exports/conga/conga_exports/tcrdist_phenotype_test"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj <- mdata_CMV
set.seed(1)
K        <- 100      # neighbors retained after same-clone + same-spec restriction
N_PERM   <- 10000
N_PCS    <- 30
OVERFETCH <- 1500    # raw neighbors fetched (within specificity) before exclusion

# ── 0. Per-cell table + cluster-level specificity (majority vote) ────────────
pcs <- Embeddings(obj, "pca")[, 1:N_PCS, drop = FALSE]
md  <- obj@meta.data
cells <- tibble(
  barcode  = colnames(obj),
  cluster  = as.character(md$TCR_cluster_et64),
  cloneId  = as.character(md$cloneId),
  donor    = as.character(md$donor_identity),
  spec_raw = as.character(md$specificity))
keep <- !is.na(cells$cluster) & !is.na(cells$cloneId) & !grepl("notCalled", cells$cloneId)
cells <- cells[keep, ]; pcs <- pcs[keep, , drop = FALSE]

# each cluster -> one specificity (majority of its cells); every cell inherits it
cluster_spec <- cells %>% dplyr::group_by(cluster) %>%
  dplyr::summarise(specificity = names(which.max(table(spec_raw))), .groups = "drop")
cells$spec <- cluster_spec$specificity[match(cells$cluster, cluster_spec$cluster)]
N <- nrow(cells)
message(sprintf("Cells: %d | clusters: %d | clones: %d | donors: %d | specificities: %s",
                N, dplyr::n_distinct(cells$cluster), dplyr::n_distinct(cells$cloneId),
                dplyr::n_distinct(cells$donor), paste(unique(cells$spec), collapse = "/")))

# ── 1. kNN within each specificity (so all candidates are same-spec) ──────────
knn_idx <- matrix(NA_integer_, nrow = N, ncol = OVERFETCH)  # global neighbor row indices
for (s in unique(cells$spec)) {
  ix  <- which(cells$spec == s)
  k_s <- min(OVERFETCH, length(ix) - 1)
  nn  <- get.knn(pcs[ix, , drop = FALSE], k = k_s)$nn.index   # indices into ix, distance-sorted
  knn_idx[ix, seq_len(k_s)] <- matrix(ix[nn], nrow = length(ix))
}

# ── 2. Per cell: K nearest surviving neighbors (drop same-clone), + x-donor flag ──
surv_idx  <- vector("list", N)
surv_xdon <- vector("list", N)
n_dropped <- 0L
for (r in seq_len(N)) {
  g <- knn_idx[r, ]; g <- g[!is.na(g)]
  g <- g[cells$cloneId[g] != cells$cloneId[r]]   # exclude same clone
  if (length(g) < K) n_dropped <- n_dropped + 1L
  g <- head(g, K)
  surv_idx[[r]]  <- g
  surv_xdon[[r]] <- cells$donor[g] != cells$donor[r]
}
message(sprintf("Cells with < K=%d same-spec non-same-clone neighbors: %d (%.1f%%)",
                K, n_dropped, 100 * n_dropped / N))

# ── 3. Per-cell statistic for a given cluster labeling ───────────────────────
compute_stats <- function(cluster_vec) {
  n_same <- numeric(N); n_xdon <- numeric(N)
  for (r in seq_len(N)) {
    g <- surv_idx[[r]]
    if (length(g) == 0) { n_same[r] <- NA; n_xdon[r] <- NA; next }
    sc <- cluster_vec[g] == cluster_vec[r]
    n_same[r] <- sum(sc)
    n_xdon[r] <- sum(sc & surv_xdon[[r]])
  }
  list(n_same = n_same, n_xdon = n_xdon)
}

# ── 4. Clone-weighted cluster aggregator ─────────────────────────────────────
agg_clone <- function(df, stat) {
  df %>% dplyr::group_by(cluster, cloneId) %>%
    dplyr::summarise(cm = mean(.data[[stat]], na.rm = TRUE), .groups = "drop") %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(score = mean(cm, na.rm = TRUE), .groups = "drop")
}

# ── 5. Permutation label generators ──────────────────────────────────────────
perm_spec <- function() {            # within specificity only
  out <- cells$cluster
  for (s in unique(cells$spec)) { ix <- which(cells$spec == s); out[ix] <- sample(cells$cluster[ix]) }
  out
}
perm_spec_donor <- function() {      # within specificity × donor
  out <- cells$cluster
  for (s in unique(cells$spec)) for (d in unique(cells$donor[cells$spec == s])) {
    ix <- which(cells$spec == s & cells$donor == d)
    if (length(ix) > 1) out[ix] <- sample(cells$cluster[ix])
  }
  out
}

# ── 6. Null runner ───────────────────────────────────────────────────────────
run_null <- function(stat, perm_fn) {
  obs_tbl  <- agg_clone(cells, stat) %>% dplyr::rename(obs_score = score)
  clusters <- obs_tbl$cluster
  null_mat <- matrix(NA_real_, nrow = N_PERM, ncol = length(clusters),
                     dimnames = list(NULL, clusters))
  for (p in seq_len(N_PERM)) {
    pl <- perm_fn(); st <- compute_stats(pl)
    tmp <- cells; tmp[[stat]] <- st[[stat]]; tmp$cluster <- pl
    a <- agg_clone(tmp, stat)
    null_mat[p, ] <- a$score[match(clusters, a$cluster)]
    if (p %% 1000 == 0) message(sprintf("    perm %d / %d", p, N_PERM))   # <- add this line
  }
  pvals <- sapply(seq_along(clusters), function(j)
    (1 + sum(null_mat[, j] >= obs_tbl$obs_score[j], na.rm = TRUE)) / (N_PERM + 1))
  obs_tbl %>% dplyr::mutate(
    null_mean = colMeans(null_mat, na.rm = TRUE),
    effect    = obs_score - null_mean,
    p_emp     = pvals,
    p_adj     = p.adjust(pvals, "BH"))           # BH within this metric
}

# ── 7. Observed stats + run both metrics with their designated nulls ─────────
obs <- compute_stats(cells$cluster); cells$n_same <- obs$n_same; cells$n_xdon <- obs$n_xdon
message("Running same_cluster (specificity-stratified null)...")
res_same <- run_null("n_same", perm_spec)
message("Running cross_donor (specificity + donor-stratified null)...")
res_xdon <- run_null("n_xdon", perm_spec_donor)

# ── 8. Cluster composition for the export ────────────────────────────────────
cluster_info <- cells %>% dplyr::group_by(cluster) %>%
  dplyr::summarise(specificity = dplyr::first(spec),
                   n_cells  = dplyr::n(),
                   n_clones = dplyr::n_distinct(cloneId),
                   n_donors = dplyr::n_distinct(donor),
                   donors   = paste(sort(unique(donor)), collapse = ","),
                   .groups = "drop")

final <- cluster_info %>%
  dplyr::left_join(res_same %>% dplyr::select(cluster,
        same_cluster_obs = obs_score, same_cluster_null_mean = null_mean,
        same_cluster_effect = effect, same_cluster_p_emp = p_emp, same_cluster_p_adj = p_adj),
      by = "cluster") %>%
  dplyr::left_join(res_xdon %>% dplyr::select(cluster,
        cross_donor_obs = obs_score, cross_donor_null_mean = null_mean,
        cross_donor_effect = effect, cross_donor_p_emp = p_emp, cross_donor_p_adj = p_adj),
      by = "cluster") %>%
  dplyr::arrange(dplyr::desc(same_cluster_effect))

# ── 9. Write CSV with an explanatory header block ────────────────────────────
csv_path <- file.path(out_dir, "tcrdist_phenotype_specstrat_results.csv")
header <- c(
"# TCR-cluster phenotypic coherence test (specificity-stratified). One row per cluster.",
"# Neighbor space: 30-PC integrated PCA. k = 100 nearest neighbors per cell, after excluding",
"#   same-cloneId neighbors and restricting candidates to the focal cell's cluster-specificity.",
"# Each cluster is assigned one specificity by majority vote of its cells; every cell in the",
"#   cluster is treated as that specificity, so cells are only ever compared within-specificity.",
"# Both metrics are clone-weighted: per-cell counts are averaged within each clone, then across",
"#   clones, so abundant clones cannot dominate a cluster's score.",
"#",
"# METRICS:",
"#   same_cluster = per cell, how many of its k neighbors are in the same TCR cluster.",
"#   cross_donor  = same count, but only neighbors from a DIFFERENT donor.",
"#",
"# COLUMN MEANINGS (prefix same_cluster_ or cross_donor_):",
"#   obs       = observed clone-weighted mean neighbor count for the cluster (the real data).",
"#   null_mean = mean of that score across 1000 permutations where cluster labels were shuffled",
"#               within the null's strata; i.e. the score expected if cluster membership carried",
"#               no phenotypic information.",
"#   effect    = obs - null_mean. How many more same-cluster (or cross-donor) neighbors each cell",
"#               has than chance. Larger = stronger phenotypic coherence.",
"#   p_emp     = empirical p-value: fraction of permutations with null score >= observed, with a",
"#               +1/+1 correction. The chance of coherence this strong if labels were random",
"#               within the strata.",
"#   p_adj     = p_emp after Benjamini-Hochberg FDR correction across clusters, separately per metric.",
"#",
"# NULL STRATA:",
"#   same_cluster p-values: cluster labels permuted within SPECIFICITY only (so significance means",
"#     'more coherent than a random grouping of same-specificity cells').",
"#   cross_donor  p-values: labels permuted within SPECIFICITY and DONOR jointly (so significance",
"#     means cross-donor coherence beyond same-specificity, same-donor-composition chance).",
"#   Single-donor clusters cannot have cross-donor neighbors; their cross_donor effect is ~0, p~1.")
con <- file(csv_path, "w"); writeLines(header, con); close(con)
write_csv(final, csv_path, append = TRUE, col_names = TRUE)

cat(sprintf("\nWrote %s\n", csv_path))
cat(sprintf("same_cluster hits (p_adj<0.05): %d | cross_donor hits: %d\n",
            sum(final$same_cluster_p_adj < 0.05, na.rm = TRUE),
            sum(final$cross_donor_p_adj  < 0.05, na.rm = TRUE)))
print(final %>% dplyr::filter(same_cluster_p_adj < 0.05) %>%
        dplyr::select(cluster, specificity, n_donors, same_cluster_effect,
                      same_cluster_p_adj, cross_donor_effect, cross_donor_p_adj), n = Inf)


#
#
# Cluster cards — kNN bias-test variant, TWO sets (all-cluster + cross-donor)
#   exports/TCR_clusters/cluster_cards/kNN_bias_test/all_cluster/
#   exports/TCR_clusters/cluster_cards/kNN_bias_test/cross_donor/
# Self-contained: all helpers namespaced kc_. Bottom-right panel = the chosen test.
#
#
library(Seurat); library(dplyr); library(tidyr); library(stringr)
library(ggplot2); library(cowplot); library(ggseqlogo); library(readr)

base_dir <- "exports/TCR_clusters/cluster_cards/kNN_bias_test"
res_csv  <- "exports/conga/conga_exports/tcrdist_phenotype_test/tcrdist_phenotype_specstrat_results.csv"

# ── Read results (skip the # header block) and build the two lookups ─────────
res <- read_csv(res_csv, comment = "#", show_col_types = FALSE)
lookup_all  <- res %>% dplyr::transmute(cluster = as.character(cluster),
                  effect = same_cluster_effect, p_adj = same_cluster_p_adj, n_donors)
lookup_xdon <- res %>% dplyr::transmute(cluster = as.character(cluster),
                  effect = cross_donor_effect,  p_adj = cross_donor_p_adj,  n_donors)

# ── kc_ helpers (palette, legend, summary, panels) ───────────────────────────
KC_QUAL_PALETTE <- c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e","#e6ab02",
  "#a6761d","#386cb0","#f0027f","#bf5b17","#4daf4a","#e41a1c","#984ea3","#ff7f00",
  "#a65628","#f781bf","#999999","#8dd3c7","#fb8072","#80b1d3")
kc_assign_palette <- function(groups) {
  n <- length(groups); cols <- KC_QUAL_PALETTE[seq_len(min(n, length(KC_QUAL_PALETTE)))]
  if (n > length(KC_QUAL_PALETTE)) cols <- colorRampPalette(KC_QUAL_PALETTE)(n)
  setNames(cols, groups)
}
KC_AA_CHEMISTRY <- data.frame(
  letter = c("G","S","T","Y","C","N","Q","K","R","H","D","E","P","A","W","F","L","I","M","V"),
  group  = c(rep("Polar",5),rep("Neutral",2),rep("Basic",3),rep("Acidic",2),rep("Hydrophobic",8)),
  col    = c(rep("#109648",5),rep("#5E239D",2),rep("#255C99",3),rep("#D62839",2),rep("#221E22",8)),
  stringsAsFactors = FALSE)
kc_export_aa_legend <- function(present_letters, out_file) {
  leg <- KC_AA_CHEMISTRY %>% dplyr::filter(letter %in% present_letters)
  if (nrow(leg) == 0) return(invisible(NULL))
  lv <- c("Polar","Neutral","Basic","Acidic","Hydrophobic")
  gs <- leg %>% dplyr::group_by(group) %>%
    dplyr::summarise(col = dplyr::first(col), aas = paste(base::sort(letter), collapse=", "),
                     .groups="drop") %>%
    dplyr::mutate(group = factor(group, levels = lv)) %>% dplyr::arrange(group) %>%
    dplyr::mutate(label = sprintf("%s (%s)", group, aas))
  cm <- setNames(gs$col, gs$label)
  d <- data.frame(x = seq_len(nrow(gs)), y = 1, label = factor(gs$label, levels = gs$label))
  lp <- ggplot(d, aes(x, y, fill = label)) + geom_tile() +
    scale_fill_manual(values = cm, name = "Amino acid chemistry") +
    guides(fill = guide_legend(override.aes = list(size = 4), ncol = 1)) +
    theme_void() + theme(legend.position = "right",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10), legend.key.size = unit(0.6, "cm"))
  ggsave(out_file, cowplot::ggdraw(cowplot::get_legend(lp)),
         width = 3.5, height = 3, dpi = 150, bg = "white")
}
kc_build_cluster_summary <- function(seurat_obj, meta_col, min_cells = 10) {
  md <- seurat_obj@meta.data %>% dplyr::filter(!is.na(.data[[meta_col]])) %>%
    dplyr::transmute(cluster = as.character(.data[[meta_col]]),
      clone_id = Global_Clone_ID, donor = paste0("D", str_extract(donor_identity, "\\d")),
      spec = specificity, cdr3 = t_cdr3s_aa, metacluster = metacluster)
  eligible <- md %>% dplyr::count(cluster, name = "total_cells") %>%
    dplyr::filter(total_cells > min_cells) %>% dplyr::arrange(dplyr::desc(total_cells))
  list(md = md, eligible = eligible)
}
kc_make_text_panel <- function(cluster_id, cmd) {
  lines <- c(sprintf("Cluster %s", cluster_id), sprintf("n = %d cells", nrow(cmd)), "")
  for (d in base::sort(base::unique(cmd$donor))) {
    dm <- cmd %>% dplyr::filter(donor == d)
    lines <- c(lines, sprintf("%s  (n = %d)", d, nrow(dm)))
    nm <- dm %>% dplyr::filter(!is.na(clone_id)) %>% dplyr::count(clone_id, spec, name = "n") %>%
      dplyr::mutate(k = as.integer(str_extract(clone_id, "\\d+$"))) %>% dplyr::arrange(spec, k)
    for (i in seq_len(nrow(nm)))
      lines <- c(lines, sprintf("  %s:  %d", str_remove(nm$clone_id[i], paste0("^", d, "_")), nm$n[i]))
    nmin <- dm %>% dplyr::filter(is.na(clone_id)) %>% nrow()
    if (nmin > 0) lines <- c(lines, sprintf("  minor:  %d", nmin)); lines <- c(lines, "")
  }
  ggplot() + annotate("text", x = 0, y = 1, label = paste(lines, collapse = "\n"),
      hjust = 0, vjust = 1, size = 3.2, family = "mono", lineheight = 1.4) +
    xlim(0, 1) + ylim(0, 1) + theme_void() + theme(plot.margin = margin(10, 0, 10, 10))
}
kc_make_umap_panel <- function(seurat_obj, cluster_id, meta_col, cmd, color_map) {
  ud <- as.data.frame(Embeddings(seurat_obj, "umap")); uc <- colnames(ud)
  mf <- seurat_obj@meta.data %>% dplyr::mutate(
    donor = paste0("D", str_extract(donor_identity, "\\d")),
    group_id = dplyr::case_when(!is.na(Global_Clone_ID) ~ as.character(Global_Clone_ID),
      !is.na(.data[[meta_col]]) & .data[[meta_col]] == as.character(cluster_id) ~ paste0(donor, "_minor-clones"),
      TRUE ~ NA_character_))
  inc <- !is.na(seurat_obj@meta.data[[meta_col]]) & seurat_obj@meta.data[[meta_col]] == as.character(cluster_id)
  ud$group_id <- NA_character_; ud$group_id[inc] <- mf$group_id[inc]
  ud$group_id <- factor(ud$group_id, levels = names(color_map))
  ggplot() +
    geom_point(data = ud[!inc, ], aes(.data[[uc[1]]], .data[[uc[2]]]), color = "lightgrey", size = 0.3) +
    geom_point(data = ud[inc, ], aes(.data[[uc[1]]], .data[[uc[2]]], color = group_id), size = 1.4) +
    scale_color_manual(values = color_map, name = "Group", drop = FALSE,
      guide = guide_legend(override.aes = list(size = 3), ncol = 1)) +
    theme_void() + theme(legend.text = element_text(size = 7),
      legend.title = element_text(size = 8, face = "bold"),
      legend.key.size = unit(0.4, "cm"), plot.margin = margin(5, 5, 5, 5))
}
kc_extract_cdr3 <- function(v, p) {
  s <- unlist(lapply(v, function(x) { if (is.na(x)) return(character(0))
    pr <- strsplit(x, ";")[[1]]; sub(paste0("^", p), "", pr[grepl(paste0("^", p), pr)]) }))
  s[nchar(s) > 0]
}
kc_filter_modal_length <- function(s) {
  if (length(s) == 0) return(character(0)); l <- nchar(s)
  m <- as.integer(names(base::sort(table(l), decreasing = TRUE))[1]); s[l == m]
}
kc_make_logo_panel <- function(v, weight_mode = "unique") {
  pull <- function(p) { s <- kc_extract_cdr3(v, p); if (weight_mode == "unique") s <- base::unique(s); kc_filter_modal_length(s) }
  lab <- if (weight_mode == "unique") "unique" else "cells"
  one <- function(s, t) { if (length(s) < 1)
    return(ggplot() + annotate("text", x=.5, y=.5, label=sprintf("%s: no\nsequences", t), size=3) + theme_void())
    ggseqlogo(s, method = if (length(s) == 1) "prob" else "bits") +
      ggtitle(sprintf("%s  (%d %s, len = %d)", t, length(s), lab, nchar(s[1]))) +
      theme(plot.title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 6), axis.title = element_text(size = 7), legend.position = "none") }
  cowplot::plot_grid(one(pull("TRA:"), "TRA"), one(pull("TRB:"), "TRB"), ncol = 1)
}
kc_make_knn_test_panel <- function(cluster_id, lookup, heading, test_type) {
  row <- lookup[lookup$cluster == as.character(cluster_id), ]
  if (nrow(row) == 0) { eff <- "Effect size: not tested"; pp <- "p (BH-adj): not tested" }
  else if (test_type == "cross_donor" && row$n_donors[1] == 1) { eff <- "n/a (single-donor cluster)"; pp <- "" }
  else { eff <- sprintf("Effect size: %.2f", row$effect[1])
    pp <- if (row$p_adj[1] < 1e-3) sprintf("p (BH-adj) = %.1e", row$p_adj[1]) else sprintf("p (BH-adj) = %.3f", row$p_adj[1]) }
  ggplot() +
    annotate("text", x = 0.5, y = 0.66, label = heading, fontface = "bold", size = 4.2) +
    annotate("text", x = 0.5, y = 0.46, label = eff, size = 4.0) +
    annotate("text", x = 0.5, y = 0.30, label = pp,  size = 4.0) +
    xlim(0, 1) + ylim(0, 1) + theme_void() + theme(plot.margin = margin(5, 5, 5, 5))
}
kc_save_card <- function(seurat_obj, cluster_id, meta_col, cmd, export_dir, lookup,
                         heading, test_type, weight_mode = "unique") {
  ng <- cmd %>% dplyr::filter(!is.na(clone_id)) %>% dplyr::pull(clone_id) %>% base::unique() %>% base::sort()
  mg <- cmd %>% dplyr::filter(is.na(clone_id)) %>% dplyr::pull(donor) %>% base::unique() %>% base::sort() %>% paste0("_minor-clones")
  cmap <- kc_assign_palette(c(ng, mg))
  top <- cowplot::plot_grid(kc_make_text_panel(cluster_id, cmd),
    kc_make_umap_panel(seurat_obj, cluster_id, meta_col, cmd, cmap), nrow = 1, rel_widths = c(0.3, 0.7))
  bot <- cowplot::plot_grid(kc_make_logo_panel(cmd$cdr3, weight_mode),
    kc_make_knn_test_panel(cluster_id, lookup, heading, test_type), nrow = 1, rel_widths = c(0.45, 0.55))
  ggsave(file.path(export_dir, sprintf("cluster_%s_card.png", cluster_id)),
    cowplot::plot_grid(top, bot, ncol = 1, rel_heights = c(0.55, 0.45)),
    width = 10, height = 10, dpi = 150, bg = "white")
}

# ── Driver: build both card sets ─────────────────────────────────────────────
res_sum  <- kc_build_cluster_summary(mdata_CMV, "TCR_cluster_et64", min_cells = 10)
md_all   <- res_sum$md; eligible <- res_sum$eligible
card_sets <- list(
  list(sub = "all_cluster", lookup = lookup_all,  heading = "All-cluster kNN bias test",  test_type = "same_cluster"),
  list(sub = "cross_donor", lookup = lookup_xdon, heading = "Cross-donor kNN bias test",  test_type = "cross_donor"))

for (cs in card_sets) {
  edir <- file.path(base_dir, cs$sub); dir.create(edir, recursive = TRUE, showWarnings = FALSE)
  present_aa <- character(0)
  message(sprintf("[%s] %d clusters > 10 cells", cs$sub, nrow(eligible)))
  for (i in seq_len(nrow(eligible))) {
    cl <- eligible$cluster[i]; cmd <- md_all %>% dplyr::filter(cluster == cl)
    tra <- kc_filter_modal_length(base::unique(kc_extract_cdr3(cmd$cdr3, "TRA:")))
    trb <- kc_filter_modal_length(base::unique(kc_extract_cdr3(cmd$cdr3, "TRB:")))
    if (length(tra)) present_aa <- base::union(present_aa, unlist(strsplit(tra, "")))
    if (length(trb)) present_aa <- base::union(present_aa, unlist(strsplit(trb, "")))
    kc_save_card(mdata_CMV, cl, "TCR_cluster_et64", cmd, edir, cs$lookup, cs$heading, cs$test_type)
  }
  kc_export_aa_legend(present_aa, file.path(edir, "aa_chemistry_legend.png"))
  message(sprintf("[%s] cards written to %s", cs$sub, edir))
}





#
#
# Constellation — Stage 0 (final): map each clustered cloneId to its file TCR
#   Node = cloneId. Cluster = TCR_cluster_et64 (object). Sequences = collaborator
#   file (alpha + beta), matched by beta CDR3, narrowed by alpha, ties by total_count.
#
#
library(dplyr); library(stringr); library(readr); library(tidyr)

conga_dir <- "exports/conga/mdata_CMV"
et64_file <- read_csv("cmv_d1-5_dataset_paired_clusters_et64_042826.csv",
                      show_col_types = FALSE)

count_cols <- c("IE1_Early","Bulk_Early","pp65_Early","Bulk_Late","IE1_Late","pp65_Late")
cl_file <- et64_file %>%
  dplyr::filter(!is.na(cluster), cluster != "") %>%
  dplyr::mutate(cluster = as.character(cluster),
                total_count = rowSums(across(all_of(count_cols)), na.rm = TRUE))

# ── 1. Clustered cloneIds with their object attributes + object CDR3 strings ──
node_base <- mdata_CMV@meta.data %>%
  dplyr::mutate(cloneId = as.character(cloneId)) %>%
  dplyr::filter(!is.na(TCR_cluster_et64), !grepl("notCalled", cloneId)) %>%
  dplyr::group_by(cloneId) %>%
  dplyr::summarise(
    cluster   = as.character(dplyr::first(TCR_cluster_et64)),
    n_cells   = dplyr::n(),
    donor     = names(which.max(table(donor_identity[!is.na(donor_identity)]))),
    spec      = names(which.max(table(specificity[!is.na(specificity)]))),
    is_global = any(!is.na(Global_Clone_ID)),
    global_id = { g <- Global_Clone_ID[!is.na(Global_Clone_ID)]
                  if (length(g)) names(which.max(table(g))) else NA_character_ },
    cdr3str   = dplyr::first(t_cdr3s_aa),
    .groups = "drop") %>%
  dplyr::mutate(obj_beta  = str_match(cdr3str, "TRB:([A-Z]+)")[,2],
                obj_alpha = str_match(cdr3str, "TRA:([A-Z]+)")[,2])

# ── 2. Match each cloneId to a file row (beta; narrow by alpha; tie by count) ─
match_file_row <- function(beta, alpha) {
  if (is.na(beta)) return(NA_integer_)
  hits <- which(cl_file$cdr3_b_aa == beta)
  if (length(hits) == 0) return(NA_integer_)
  if (!is.na(alpha)) {
    narrow <- hits[cl_file$cdr3_a_aa[hits] == alpha]
    if (length(narrow) > 0) hits <- narrow
  }
  hits[which.max(cl_file$total_count[hits])]
}
node_base$file_row <- mapply(match_file_row, node_base$obj_beta, node_base$obj_alpha)

# ── 3. Attach file sequences; report coverage + label agreement ──────────────
nodes <- node_base %>%
  dplyr::mutate(
    cdr3_a_aa = cl_file$cdr3_a_aa[file_row], v_a_gene = cl_file$v_a_gene[file_row],
    j_a_gene  = cl_file$j_a_gene[file_row],
    cdr3_b_aa = cl_file$cdr3_b_aa[file_row], v_b_gene = cl_file$v_b_gene[file_row],
    j_b_gene  = cl_file$j_b_gene[file_row],
    file_cluster = cl_file$cluster[file_row])

matched <- nodes %>% dplyr::filter(!is.na(file_row))
cat(sprintf("clustered cloneIds: %d | matched to file: %d (%.0f%%)\n",
            nrow(nodes), nrow(matched), 100*nrow(matched)/nrow(nodes)))
cat(sprintf("cluster-label agreement (object vs file) among matched: %d / %d (%.0f%%)\n",
            sum(matched$cluster == matched$file_cluster), nrow(matched),
            100*mean(matched$cluster == matched$file_cluster)))

# keep only matched nodes (need a sequence for tcrdist3); write tcrdist3 input
nodes <- matched
write_csv(nodes %>% dplyr::select(cloneId, cluster, cdr3_a_aa, cdr3_b_aa,
                                  v_a_gene, v_b_gene, j_a_gene, j_b_gene),
          file.path(conga_dir, "constellation_tcrdist3_input.csv"))
saveRDS(nodes, file.path(conga_dir, "constellation_nodes.rds"))
message("Wrote tcrdist3 input + node attributes")



#
#
# Constellation — Stage 1: build 600-node <=64 edge graph from the 84x84 tcrdist matrix
#   Expands unique-sequence distances back to cloneId nodes via the clone->seq map.
#   Validates edges fall within et64 clusters. Saves graph for the plot.
#
#
library(dplyr); library(readr); library(tidyr); library(igraph)

conga_dir      <- "exports/conga/mdata_CMV"
EDGE_THRESHOLD <- 64

D        <- as.matrix(read_csv(file.path(conga_dir, "constellation_tcrdist_matrix.csv"),
                               show_col_types = FALSE)); dimnames(D) <- NULL
row_ord  <- read_csv(file.path(conga_dir, "constellation_tcrdist_order.csv"), show_col_types = FALSE)
clone2seq<- read_csv(file.path(conga_dir, "constellation_clone_to_seq.csv"), show_col_types = FALSE)
nodes    <- readRDS(file.path(conga_dir, "constellation_nodes.rds"))

# ── 0. Reconcile id base: normalize every seq_id to 0-based ──────────────────
# tcrdist3 emitted 1-based clone_id; our clone2seq is 0-based. Shift to match.
row_ord <- row_ord %>% dplyr::mutate(seq_id0 = seq_id - min(seq_id))     # -> 0-based
clone2seq <- clone2seq %>% dplyr::mutate(seq_id0 = seq_id - min(seq_id)) # already 0-based; min=0
stopifnot(nrow(D) == nrow(row_ord))
stopifnot(setequal(row_ord$seq_id0, clone2seq$seq_id0))   # HARD CHECK: id sets must match
# matrix row k -> seq_id0 = row_ord$seq_id0[k]; build seq_id0 -> matrix-row lookup
seq_to_row <- setNames(row_ord$row_index + 1L, row_ord$seq_id0)  # +1: R 1-based rows
message(sprintf("Matrix: %d unique seqs | cloneIds mapped: %d | id sets match: TRUE",
                nrow(D), nrow(clone2seq)))

# ── 1. cloneId -> matrix row (many cloneIds may share a seq row) ─────────────
c2r <- clone2seq %>%
  dplyr::mutate(mrow = seq_to_row[as.character(seq_id0)]) %>%
  dplyr::select(cloneId, seq_id0, mrow)
node_tbl <- nodes %>% dplyr::inner_join(c2r, by = "cloneId")
message(sprintf("Nodes with sequence + attributes: %d", nrow(node_tbl)))

# ── 2. Intra-cluster edges at TCRdist <= threshold ───────────────────────────
# For each cluster, all cloneId pairs whose seq-rows are within threshold (incl. 0).
edges <- lapply(split(node_tbl, node_tbl$cluster), function(cl) {
  if (nrow(cl) < 2) return(NULL)
  cmb <- t(combn(nrow(cl), 2))
  d   <- D[cbind(cl$mrow[cmb[,1]], cl$mrow[cmb[,2]])]
  keep <- which(d <= EDGE_THRESHOLD)
  if (length(keep) == 0) return(NULL)
  tibble(cluster = cl$cluster[1],
         from = cl$cloneId[cmb[keep,1]], to = cl$cloneId[cmb[keep,2]], tcrdist = d[keep])
}) %>% bind_rows()

# ── 3. VALIDATION: cross-cluster leakage check + per-cluster connectivity ────
# (recompute ignoring cluster, to see if any <=64 edges would cross clusters)
all_pairs_leak <- {
  mr <- node_tbl$mrow; cl <- node_tbl$cluster
  cmb <- t(combn(nrow(node_tbl), 2))
  d <- D[cbind(mr[cmb[,1]], mr[cmb[,2]])]
  near <- which(d <= EDGE_THRESHOLD)
  sum(cl[cmb[near,1]] != cl[cmb[near,2]])
}
edge_summary <- node_tbl %>% dplyr::count(cluster, name = "n_clones") %>%
  dplyr::left_join(dplyr::count(edges, cluster, name = "n_edges"), by = "cluster") %>%
  dplyr::mutate(n_edges = tidyr::replace_na(n_edges, 0L)) %>%
  dplyr::arrange(dplyr::desc(n_clones))

cat(sprintf("\nIntra-cluster edges (<=%d): %d\n", EDGE_THRESHOLD, nrow(edges)))
cat(sprintf("Cross-cluster <=%d edges (should be ~0 if scale matches et64): %d\n",
            EDGE_THRESHOLD, all_pairs_leak))
print(edge_summary, n = Inf)

saveRDS(list(nodes = node_tbl, edges = edges, D = D, threshold = EDGE_THRESHOLD),
        file.path(conga_dir, "constellation_graph.rds"))




#
#
# Constellation plot — base version (specificity/donor/size/Global encodings)
#   Edges: within-cluster only. Convergent (identical-aa) clones = distance-0 clique;
#   backbone = MST on paired TCRdist. Layout: per-cluster force-directed, circles
#   packed largest-at-center. Two outputs: all clusters, and clusters with >=5 clones.
#
#
library(dplyr); library(igraph); library(ggraph); library(ggplot2)
library(packcircles); library(tidygraph); library(readr)

conga_dir <- "exports/conga/mdata_CMV"
out_dir   <- "exports/TCR_clusters/constellation"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g_data   <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
node_tbl <- g_data$nodes          # cloneId, cluster, n_cells, donor, spec, is_global, global_id, mrow
D        <- g_data$D

spec_colors  <- c("pp65" = "#1976D2", "IE1" = "#8B0000")
donor_shapes <- c("Donor 1" = 21, "Donor 2" = 22, "Donor 3" = 24, "Donor 4" = 25, "Donor 5" = 23)

# ── 1. Within-cluster edges: 0-distance cliques + MST backbone ───────────────
build_cluster_edges <- function(cl) {
  n <- nrow(cl); if (n < 2) return(NULL)
  sub <- D[cl$mrow, cl$mrow, drop = FALSE]
  # MST over ALL nodes (no zero-distance cliques). Shift weights by +1 so genuine
  # zero-distance pairs still register as edges, then subtract it back.
  gg    <- graph_from_adjacency_matrix(sub + 1, mode = "undirected", weighted = TRUE, diag = FALSE)
  mst_g <- igraph::mst(gg)
  if (igraph::ecount(mst_g) == 0) return(NULL)
  el <- igraph::as_edgelist(mst_g)
  tibble(from = cl$cloneId[as.integer(el[,1])],
         to   = cl$cloneId[as.integer(el[,2])],
         tcrdist = igraph::E(mst_g)$weight - 1)
}

# ── 2. Per-cluster force-directed layout, normalized into a unit disc ─────────
layout_one_cluster <- function(cl, edges) {
  n <- nrow(cl)
  if (n == 1) return(tibble(cloneId = cl$cloneId, x = 0, y = 0))
  gg <- graph_from_data_frame(
    if (is.null(edges) || nrow(edges) == 0)
      data.frame(from = cl$cloneId[1], to = cl$cloneId[1])[0,] else edges[,c("from","to")],
    vertices = cl$cloneId, directed = FALSE)
  set.seed(1)
  xy <- layout_with_fr(gg)                                   # force-directed
  xy <- scale(xy, center = TRUE, scale = FALSE)
  r  <- max(sqrt(rowSums(xy^2))); if (r == 0 || !is.finite(r)) r <- 1
  xy <- xy / r                                               # -> unit disc
  tibble(cloneId = igraph::V(gg)$name, x = xy[,1], y = xy[,2])
}

# ── 3. Assemble full graph: pack clusters, biggest at center ─────────────────
render_constellation <- function(keep_clusters, file_tag, title) {
  nt <- node_tbl %>% dplyr::filter(cluster %in% keep_clusters)
  # clones per cluster (drives packing radius) + cells per cluster, computed separately
  clones_per <- nt %>% dplyr::count(cluster, name = "n_clones")
  cells_per  <- nt %>% dplyr::group_by(cluster) %>%
    dplyr::summarise(cells = sum(n_cells), .groups = "drop")
  clu_sizes <- clones_per %>%
    dplyr::left_join(cells_per, by = "cluster") %>%
    dplyr::arrange(dplyr::desc(n_clones))
  # circle radius per cluster ~ sqrt(n_clones); pack with packcircles
  circ <- packcircles::circleProgressiveLayout(sqrt(clu_sizes$n_clones), sizetype = "radius")
  clu_sizes <- dplyr::bind_cols(clu_sizes, circ)             # x,y,radius per cluster
  # order so biggest is centered: progressive layout already spirals largest-first

  # per-cluster node coords, scaled to that cluster's packed circle
  all_nodes <- lapply(seq_len(nrow(clu_sizes)), function(i) {
    cid <- clu_sizes$cluster[i]
    cl  <- nt %>% dplyr::filter(cluster == cid)
    ed  <- build_cluster_edges(cl)
    lay <- layout_one_cluster(cl, ed)
    lay %>% dplyr::mutate(cluster = cid,
                          x = clu_sizes$x[i] + x * clu_sizes$radius[i] * 0.85,
                          y = clu_sizes$y[i] + y * clu_sizes$radius[i] * 0.85)
  }) %>% bind_rows() %>%
    dplyr::left_join(nt, by = c("cloneId","cluster"))

  # edges in packed coordinates
  all_edges <- lapply(clu_sizes$cluster, function(cid) {
    cl <- nt %>% dplyr::filter(cluster == cid); ed <- build_cluster_edges(cl)
    if (is.null(ed) || nrow(ed) == 0) return(NULL)
    ed %>% dplyr::mutate(cluster = cid)
  }) %>% bind_rows()
  pos <- all_nodes %>% dplyr::select(cloneId, x, y)
  all_edges <- all_edges %>%
    dplyr::left_join(pos, by = c("from"="cloneId")) %>% dplyr::rename(x = x, y = y) %>%
    dplyr::left_join(pos, by = c("to"="cloneId"), suffix = c("","end")) %>%
    dplyr::rename(xend = xend, yend = yend)

  # cluster labels at circle centers
  clab <- clu_sizes %>% dplyr::mutate(y = y + radius) %>% dplyr::select(cluster, x, y)

  p <- ggplot() +
    geom_segment(data = all_edges, aes(x = x, y = y, xend = xend, yend = yend),
                 color = "grey70", linewidth = 0.25, alpha = 0.7) +
    geom_point(data = all_nodes,
               aes(x, y, size = n_cells, shape = donor,
                   fill = ifelse(is_global, spec, NA), color = spec)) +
    scale_shape_manual(values = donor_shapes, name = "Donor") +
    scale_fill_manual(values = spec_colors, na.value = NA, name = "Specificity (Global)") +
    scale_color_manual(values = spec_colors, name = "Specificity") +
    scale_size_area(max_size = 9, name = "Cells in clone") +
    geom_text(data = clab, aes(x, y, label = cluster), size = 2.6, color = "grey30") +
    coord_equal() + theme_void() +
    ggtitle(title) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 4)),
           color = guide_legend(override.aes = list(shape = 21, size = 4)),
           shape = guide_legend(override.aes = list(size = 4)))

  ggsave(file.path(out_dir, sprintf("constellation_%s.png", file_tag)),
         p, width = 16, height = 14, dpi = 300, bg = "white")
  message(sprintf("[%s] %d clusters, %d nodes, %d edges -> %s",
                  file_tag, nrow(clu_sizes), nrow(all_nodes), nrow(all_edges),
                  file.path(out_dir, sprintf("constellation_%s.png", file_tag))))
}

# ── 4. Render both versions ──────────────────────────────────────────────────
all_clusters <- node_tbl %>% dplyr::distinct(cluster) %>% dplyr::pull(cluster)
big_clusters <- node_tbl %>% dplyr::count(cluster) %>% dplyr::filter(n >= 5) %>% dplyr::pull(cluster)

render_constellation(all_clusters, "all",        "et64 TCR constellations — all clusters")
render_constellation(big_clusters, "min5clones", "et64 TCR constellations — clusters with >=5 clones")



#
#
# Constellation — by sequence: log-sized nodes, two-tier edges, no cross-cluster spill
#   Node size: radius linear in log(cells), min->max across the whole plot.
#   Edges: MST per cluster guarantees every cluster is connected; edges that are ALSO
#   within EDGE_CUT TCRdist are drawn darker/thicker (a "these are close" tier).
#   Layout: log1p(TCRdist) via KK, repelled + hard-fit inside each packed cluster circle.
#   Output: exports/TCR_clusters/constellation/by_sequence_v2/
#
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr)

conga_dir <- "exports/conga/mdata_CMV"
out_dir   <- "exports/TCR_clusters/constellation/by_sequence_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g_data <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
node_tbl <- g_data$nodes; D <- g_data$D
donor_color_map <- c("D1"="#1b9e77","D2"="#d95f02","D3"="#7570b3","D4"="#e7298a","D5"="#66a61e")
spec_colors <- c("pp65"="#1976D2","IE1"="#8B0000")
node_tbl <- node_tbl %>% dplyr::mutate(donorS = paste0("D", str_extract(as.character(donor), "\\d")))

# ── Tunable settings ─────────────────────────────────────────────────────────
EDGE_CUT       <- 250                    # RAW TCRdist below which an edge is "close" (dark tier)
DIST_TRANSFORM <- function(d) log1p(d)     # near sequences closer; sqrt() = gentler
NODE_R_MIN     <- 0.35                    # smallest node radius (fewest cells)
NODE_R_MAX     <- 1.40                    # largest node radius (most cells)
CLUST_PAD      <- 3                       # spacing between cluster circles
GAP_MIN        <- 0.25                    # min gap between node edges within a cluster
EDGE_FAR_COL   <- "grey80"; EDGE_FAR_LW   <- 0.25   # MST edge, sequences far apart
EDGE_NEAR_COL  <- "grey35"; EDGE_NEAR_LW  <- 0.7    # MST edge that is ALSO within EDGE_CUT

# ── 1. Unique-sequence nodes + LOG-VALUE node radius ─────────────────────────
seq_nodes <- node_tbl %>%
  dplyr::group_by(cluster, mrow) %>%
  dplyr::summarise(total_cells=sum(n_cells), spec=names(which.max(table(spec))),
                   d1=sum(n_cells[donorS=="D1"]), d2=sum(n_cells[donorS=="D2"]),
                   d3=sum(n_cells[donorS=="D3"]), d4=sum(n_cells[donorS=="D4"]),
                   d5=sum(n_cells[donorS=="D5"]), .groups="drop") %>%
  dplyr::mutate(seq_node = paste0("s", dplyr::row_number()))

lc_min <- log(min(seq_nodes$total_cells)); lc_max <- log(max(seq_nodes$total_cells))
seq_nodes <- seq_nodes %>%
  dplyr::mutate(node_r = NODE_R_MIN + (NODE_R_MAX-NODE_R_MIN) *
                  (log(total_cells) - lc_min) / (lc_max - lc_min))
message(sprintf("Sequence nodes: %d | cells %d-%d -> radius %.2f-%.2f",
                nrow(seq_nodes), min(seq_nodes$total_cells), max(seq_nodes$total_cells),
                min(seq_nodes$node_r), max(seq_nodes$node_r)))

# ── 2. Local layout (log-distance KK) with guaranteed node separation ────────
layout_one <- function(sn) {
  n <- nrow(sn)
  if (n == 1) return(tibble(seq_node=sn$seq_node, x=0, y=0))
  sub   <- D[sn$mrow, sn$mrow, drop=FALSE]
  sub_t <- DIST_TRANSFORM(sub)
  if (n == 2) { xy <- matrix(c(0,0, sub_t[1,2],0), nrow=2, byrow=TRUE) }
  else {
    gg <- graph_from_adjacency_matrix(sub_t, mode="undirected", weighted=TRUE, diag=FALSE)
    set.seed(1); xy <- layout_with_kk(gg, weights = E(gg)$weight)
  }
  xy <- scale(xy, center=TRUE, scale=FALSE)
  tibble(seq_node=sn$seq_node, x=xy[,1], y=xy[,2])
}

layout_cluster_local <- function(sn) {
  lay <- layout_one(sn) %>% dplyr::left_join(sn, by="seq_node")
  if (nrow(lay) == 1)
    return(list(coords = lay %>% dplyr::transmute(seq_node, x=0, y=0, node_r), bound = lay$node_r[1]))
  sp <- max(sqrt(lay$x^2 + lay$y^2)); if (sp==0 || !is.finite(sp)) sp <- 1
  lay$x <- lay$x/sp; lay$y <- lay$y/sp
  set.seed(1)
  lay$x <- lay$x + rnorm(nrow(lay),0,0.01); lay$y <- lay$y + rnorm(nrow(lay),0,0.01)
  rep <- packcircles::circleRepelLayout(
    data.frame(x=lay$x, y=lay$y, radius=lay$node_r + GAP_MIN/2),
    xysizecols=1:3, sizetype="radius", maxiter=1000)$layout
  lay$x <- rep$x; lay$y <- rep$y
  cx <- mean(range(lay$x)); cy <- mean(range(lay$y)); lay$x <- lay$x - cx; lay$y <- lay$y - cy
  bound <- max(sqrt(lay$x^2 + lay$y^2) + lay$node_r)
  list(coords = lay %>% dplyr::select(seq_node, x, y, node_r), bound = bound)
}

# ── 3. Two-tier edges: MST (all nodes) + flag which MST edges are within EDGE_CUT ──
cluster_edges <- function(sn) {
  if (nrow(sn) < 2) return(NULL)
  sub <- D[sn$mrow, sn$mrow, drop=FALSE]
  gg  <- graph_from_adjacency_matrix(sub + 1, mode="undirected", weighted=TRUE, diag=FALSE)
  m   <- igraph::mst(gg)
  if (igraph::ecount(m) == 0) return(NULL)
  el <- igraph::as_edgelist(m)
  i1 <- as.integer(el[,1]); i2 <- as.integer(el[,2])
  raw_d <- sub[cbind(i1, i2)]                       # true TCRdist of each MST edge
  tibble(from = sn$seq_node[i1], to = sn$seq_node[i2],
         near = raw_d <= EDGE_CUT)                  # TRUE = close (dark tier)
}

# ── 4. Build: local layouts -> pack by true bounds -> place 1:1 ──────────────
build_geometry <- function(keep) {
  sn_all <- seq_nodes %>% dplyr::filter(cluster %in% keep)
  cl_ids <- sn_all %>% dplyr::count(cluster, name="n_seq") %>% dplyr::arrange(dplyr::desc(n_seq)) %>%
    dplyr::pull(cluster)
  locals <- lapply(cl_ids, function(cid) layout_cluster_local(sn_all %>% dplyr::filter(cluster==cid)))
  bounds <- sapply(locals, function(z) z$bound)
  circ <- packcircles::circleProgressiveLayout(bounds + CLUST_PAD, sizetype="radius")

  coords <- list(); edges_all <- list()
  for (i in seq_along(cl_ids)) {
    cid <- cl_ids[i]; lc <- locals[[i]]$coords
    lc <- lc %>% dplyr::mutate(x = circ$x[i] + x, y = circ$y[i] + y, cluster = cid)
    coords[[i]] <- lc
    ed <- cluster_edges(sn_all %>% dplyr::filter(cluster==cid))
    if (!is.null(ed)) { p <- lc %>% dplyr::select(seq_node,x,y)
      edges_all[[i]] <- ed %>% dplyr::left_join(p,by=c("from"="seq_node")) %>%
        dplyr::left_join(p,by=c("to"="seq_node"),suffix=c("",".e")) %>% dplyr::rename(xend=x.e,yend=y.e) }
  }
  nodes <- bind_rows(coords) %>% dplyr::left_join(sn_all %>% dplyr::select(-node_r), by=c("seq_node","cluster"))
  list(nodes = nodes, edges = bind_rows(edges_all))
}

# ── 5. Render (two edge tiers drawn separately) ──────────────────────────────
render <- function(keep, tag, title, mode) {
  geo <- build_geometry(keep); nodes <- geo$nodes; edges <- geo$edges
  e_far  <- edges %>% dplyr::filter(!near)
  e_near <- edges %>% dplyr::filter(near)

  p <- ggplot() +
    geom_segment(data=e_far,  aes(x,y,xend=xend,yend=yend), color=EDGE_FAR_COL,  linewidth=EDGE_FAR_LW) +
    geom_segment(data=e_near, aes(x,y,xend=xend,yend=yend), color=EDGE_NEAR_COL, linewidth=EDGE_NEAR_LW)
  if (mode=="donor") {
    p <- p +
      ggforce::geom_circle(data=nodes, aes(x0=x,y0=y,r=node_r), fill="white", color="grey30", linewidth=0.2) +
      scatterpie::geom_scatterpie(data=nodes, aes(x=x,y=y,r=node_r),
             cols=c("d1","d2","d3","d4","d5"), color=NA) +
      scale_fill_manual(values=setNames(unname(donor_color_map),c("d1","d2","d3","d4","d5")),
                        labels=c("D1","D2","D3","D4","D5"), name="Donor")
  } else {
    p <- p + ggforce::geom_circle(data=nodes, aes(x0=x,y0=y,r=node_r,fill=spec),
             color="grey30", linewidth=0.2) +
      scale_fill_manual(values=spec_colors, name="Specificity")
  }
  p <- p + coord_equal() + theme_void() + ggtitle(title) + theme(plot.title=element_text(size=14))
  ggsave(file.path(out_dir, sprintf("constellation_%s.png", tag)), p, width=16, height=14, dpi=300, bg="white")
  message(sprintf("[%s | %s] %d clusters, %d nodes, %d edges (%d near / %d far)",
                  tag, mode, dplyr::n_distinct(nodes$cluster), nrow(nodes), nrow(edges),
                  sum(edges$near), sum(!edges$near)))
}

all_cl <- seq_nodes %>% dplyr::distinct(cluster) %>% dplyr::pull(cluster)
big_cl <- node_tbl %>% dplyr::count(cluster) %>% dplyr::filter(n>=5) %>% dplyr::pull(cluster)
render(all_cl,"donor_all","et64 constellations (by sequence) — donor","donor")
render(all_cl,"spec_all","et64 constellations (by sequence) — specificity","spec")
render(big_cl,"donor_min5clones","et64 constellations (by sequence, >=5 clones) — donor","donor")
render(big_cl,"spec_min5clones","et64 constellations (by sequence, >=5 clones) — specificity","spec")



#
#
# Constellation — threshold-graph style (matches collaborator aesthetic)
#   Node = cloneId. Edges = ALL within-cluster pairs with TCRdist <= EDGE_CUT.
#   Larger nodes, tight packing. Two colorings: donor and specificity.
#   Output: exports/TCR_clusters/constellation/threshold_graph/
#
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles); library(stringr)

conga_dir <- "exports/conga/mdata_CMV"
out_dir   <- "exports/TCR_clusters/constellation/threshold_graph"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g_data   <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
node_tbl <- g_data$nodes; D <- g_data$D

EDGE_CUT <- 100     # within-cluster edge threshold on the alpha+beta scale — TUNE THIS
NODE_MAX <- 6       # max node size
donor_color_map <- c("D1"="#1b9e77","D2"="#d95f02","D3"="#7570b3",
                     "D4"="#e7298a","D5"="#66a61e","Other"="grey70")
spec_colors     <- c("pp65"="#1976D2","IE1"="#8B0000")
node_tbl <- node_tbl %>% dplyr::mutate(donorS = paste0("D", str_extract(as.character(donor), "\\d")))

# ── 1. Within-cluster threshold edges (ALL pairs <= EDGE_CUT) ─────────────────
cluster_edges <- function(cl) {
  if (nrow(cl) < 2) return(NULL)
  sub <- D[cl$mrow, cl$mrow, drop = FALSE]
  ij  <- which(upper.tri(sub) & sub <= EDGE_CUT, arr.ind = TRUE)
  if (nrow(ij) == 0) return(NULL)
  tibble(from = cl$cloneId[ij[,1]], to = cl$cloneId[ij[,2]])
}
layout_one <- function(cl, ed) {
  if (nrow(cl) == 1) return(tibble(cloneId = cl$cloneId, x=0, y=0))
  gg <- graph_from_data_frame(
    if (is.null(ed)) data.frame(from=cl$cloneId[1], to=cl$cloneId[1])[0,] else ed,
    vertices = cl$cloneId, directed = FALSE)
  set.seed(1); xy <- layout_with_fr(gg); xy <- scale(xy, center=TRUE, scale=FALSE)
  r <- max(sqrt(rowSums(xy^2))); if (r==0 || !is.finite(r)) r <- 1
  tibble(cloneId = igraph::V(gg)$name, x = xy[,1]/r, y = xy[,2]/r)
}

render <- function(keep, tag, title, mode) {
  nt <- node_tbl %>% dplyr::filter(cluster %in% keep)
  sizes <- nt %>% dplyr::count(cluster, name="n") %>% dplyr::arrange(dplyr::desc(n))
  circ  <- packcircles::circleProgressiveLayout(sqrt(sizes$n), sizetype="radius")
  sizes <- dplyr::bind_cols(sizes, circ)

  coords <- list(); edges_all <- list()
  for (i in seq_len(nrow(sizes))) {
    cid <- sizes$cluster[i]; cl <- nt %>% dplyr::filter(cluster==cid)
    ed  <- cluster_edges(cl); lay <- layout_one(cl, ed)
    lay <- lay %>% dplyr::mutate(x = sizes$x[i] + x*sizes$radius[i]*0.9,
                                 y = sizes$y[i] + y*sizes$radius[i]*0.9)
    coords[[i]] <- lay
    if (!is.null(ed)) {
      p <- lay %>% dplyr::select(cloneId, x, y)
      edges_all[[i]] <- ed %>%
        dplyr::left_join(p, by=c("from"="cloneId")) %>%
        dplyr::left_join(p, by=c("to"="cloneId"), suffix=c("",".e")) %>%
        dplyr::rename(xend=x.e, yend=y.e)
    }
  }
  nodes <- bind_rows(coords) %>% dplyr::left_join(nt, by="cloneId")
  edges <- bind_rows(edges_all)

  aes_fill <- if (mode=="donor") quote(donorS) else quote(spec)
  pal      <- if (mode=="donor") donor_color_map else spec_colors
  leg      <- if (mode=="donor") "Donor" else "Specificity"

  p <- ggplot() +
    geom_segment(data=edges, aes(x,y,xend=xend,yend=yend),
                 color="grey55", linewidth=0.4, alpha=0.8) +
    geom_point(data=nodes, aes(x, y, size=n_cells, fill=!!aes_fill),
               shape=21, color="grey25", stroke=0.3) +
    scale_fill_manual(values=pal, name=leg) +
    scale_size_area(max_size=NODE_MAX, name="Cells in clone") +
    coord_equal() + theme_void() + ggtitle(title) +
    guides(fill = guide_legend(override.aes = list(size=5)))

  ggsave(file.path(out_dir, sprintf("constellation_thr_%s.png", tag)),
         p, width=15, height=13, dpi=300, bg="white")
  message(sprintf("[%s] %d clusters, %d nodes, %d edges (cut=%d)",
                  tag, nrow(sizes), nrow(nodes), nrow(edges), EDGE_CUT))
}

all_cl <- node_tbl %>% dplyr::distinct(cluster) %>% dplyr::pull(cluster)
big_cl <- node_tbl %>% dplyr::count(cluster) %>% dplyr::filter(n>=5) %>% dplyr::pull(cluster)

render(all_cl, "donor_all",       "et64 TCR constellations — donor",       "donor")
render(all_cl, "spec_all",        "et64 TCR constellations — specificity", "spec")
render(big_cl, "donor_min5clones","et64 TCR constellations (>=5 clones) — donor",       "donor")
render(big_cl, "spec_min5clones", "et64 TCR constellations (>=5 clones) — specificity", "spec")



#
#
# Constellation — split by dominant specificity, donor overlay, with cluster labels
#   Node = unique sequence; size = log(cells); two-tier MST edges.
#   Clusters split into pp65-dominant vs IE1-dominant; cutoffs: all, and >=2 sequences.
#   Each cluster labeled "cluster X" floating above its topmost node.
#   Files: {cutoff}_{specificity}_{overlay}.png in .../constellation/by_specificity/
#
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr)

conga_dir <- "exports/conga/mdata_CMV"
out_dir   <- "exports/TCR_clusters/constellation/by_specificity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g_data <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
node_tbl <- g_data$nodes; D <- g_data$D
donor_color_map <- c("D1"="#1b9e77","D2"="#d95f02","D3"="#7570b3","D4"="#e7298a","D5"="#66a61e")
node_tbl <- node_tbl %>% dplyr::mutate(donorS = paste0("D", str_extract(as.character(donor), "\\d")))

# ── Tunable settings ─────────────────────────────────────────────────────────
EDGE_CUT       <- 250
DIST_TRANSFORM <- function(d) log1p(d)
NODE_R_MIN     <- 0.35
NODE_R_MAX     <- 1.40
CLUST_PAD      <- 3
GAP_MIN        <- 0.25
EDGE_FAR_COL   <- "grey80"; EDGE_FAR_LW   <- 0.25
EDGE_NEAR_COL  <- "grey35"; EDGE_NEAR_LW  <- 0.7
LABEL_OFFSET   <- 0.8      # vertical gap between a cluster's top node and its label
LABEL_SIZE     <- 3.5      # cluster label text size

# ── 1. Unique-sequence nodes + log-value radius + per-cluster dominant spec ──
seq_nodes <- node_tbl %>%
  dplyr::group_by(cluster, mrow) %>%
  dplyr::summarise(total_cells=sum(n_cells), spec=names(which.max(table(spec))),
                   d1=sum(n_cells[donorS=="D1"]), d2=sum(n_cells[donorS=="D2"]),
                   d3=sum(n_cells[donorS=="D3"]), d4=sum(n_cells[donorS=="D4"]),
                   d5=sum(n_cells[donorS=="D5"]), .groups="drop") %>%
  dplyr::mutate(seq_node = paste0("s", dplyr::row_number()))

# radius uses the GLOBAL cell range so sizes are comparable across all four plots
lc_min <- log(min(seq_nodes$total_cells)); lc_max <- log(max(seq_nodes$total_cells))
seq_nodes <- seq_nodes %>%
  dplyr::mutate(node_r = NODE_R_MIN + (NODE_R_MAX-NODE_R_MIN) *
                  (log(total_cells) - lc_min) / (lc_max - lc_min))

# dominant specificity per cluster (cell-weighted majority) + unique-seq count
cluster_meta <- seq_nodes %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n_seq = dplyr::n(),
                   dom_spec = { w <- tapply(total_cells, spec, sum)
                                names(w)[which.max(w)] }, .groups="drop")
message(sprintf("Clusters: %d | pp65-dom %d, IE1-dom %d | with >=2 seq: %d",
                nrow(cluster_meta), sum(cluster_meta$dom_spec=="pp65"),
                sum(cluster_meta$dom_spec=="IE1"), sum(cluster_meta$n_seq>=2)))

# ── 2. Layout + two-tier MST edges ───────────────────────────────────────────
layout_one <- function(sn) {
  n <- nrow(sn)
  if (n == 1) return(tibble(seq_node=sn$seq_node, x=0, y=0))
  sub   <- D[sn$mrow, sn$mrow, drop=FALSE]; sub_t <- DIST_TRANSFORM(sub)
  if (n == 2) { xy <- matrix(c(0,0, sub_t[1,2],0), nrow=2, byrow=TRUE) }
  else { gg <- graph_from_adjacency_matrix(sub_t, mode="undirected", weighted=TRUE, diag=FALSE)
    set.seed(1); xy <- layout_with_kk(gg, weights = E(gg)$weight) }
  xy <- scale(xy, center=TRUE, scale=FALSE)
  tibble(seq_node=sn$seq_node, x=xy[,1], y=xy[,2])
}
layout_cluster_local <- function(sn) {
  lay <- layout_one(sn) %>% dplyr::left_join(sn, by="seq_node")
  if (nrow(lay) == 1)
    return(list(coords = lay %>% dplyr::transmute(seq_node, x=0, y=0, node_r), bound = lay$node_r[1]))
  sp <- max(sqrt(lay$x^2 + lay$y^2)); if (sp==0 || !is.finite(sp)) sp <- 1
  lay$x <- lay$x/sp; lay$y <- lay$y/sp
  set.seed(1); lay$x <- lay$x + rnorm(nrow(lay),0,0.01); lay$y <- lay$y + rnorm(nrow(lay),0,0.01)
  rep <- packcircles::circleRepelLayout(
    data.frame(x=lay$x, y=lay$y, radius=lay$node_r + GAP_MIN/2),
    xysizecols=1:3, sizetype="radius", maxiter=1000)$layout
  lay$x <- rep$x; lay$y <- rep$y
  cx <- mean(range(lay$x)); cy <- mean(range(lay$y)); lay$x <- lay$x - cx; lay$y <- lay$y - cy
  bound <- max(sqrt(lay$x^2 + lay$y^2) + lay$node_r)
  list(coords = lay %>% dplyr::select(seq_node, x, y, node_r), bound = bound)
}
cluster_edges <- function(sn) {
  if (nrow(sn) < 2) return(NULL)
  sub <- D[sn$mrow, sn$mrow, drop=FALSE]
  gg  <- graph_from_adjacency_matrix(sub + 1, mode="undirected", weighted=TRUE, diag=FALSE)
  m   <- igraph::mst(gg); if (igraph::ecount(m) == 0) return(NULL)
  el <- igraph::as_edgelist(m); i1 <- as.integer(el[,1]); i2 <- as.integer(el[,2])
  tibble(from = sn$seq_node[i1], to = sn$seq_node[i2], near = sub[cbind(i1,i2)] <= EDGE_CUT)
}

# ── 3. Build geometry (nodes, edges, labels) for a set of clusters ───────────
build_geometry <- function(keep_clusters) {
  sn_all <- seq_nodes %>% dplyr::filter(cluster %in% keep_clusters)
  cl_ids <- sn_all %>% dplyr::count(cluster, name="n_seq") %>% dplyr::arrange(dplyr::desc(n_seq)) %>%
    dplyr::pull(cluster)
  locals <- lapply(cl_ids, function(cid) layout_cluster_local(sn_all %>% dplyr::filter(cluster==cid)))
  bounds <- sapply(locals, function(z) z$bound)
  circ <- packcircles::circleProgressiveLayout(bounds + CLUST_PAD, sizetype="radius")
  coords <- list(); edges_all <- list(); labels <- list()
  for (i in seq_along(cl_ids)) {
    cid <- cl_ids[i]; lc <- locals[[i]]$coords %>%
      dplyr::mutate(x = circ$x[i] + x, y = circ$y[i] + y, cluster = cid)
    coords[[i]] <- lc
    top_y <- max(lc$y + lc$node_r)                       # cluster's actual top edge
    labels[[i]] <- tibble(cluster = cid, x = circ$x[i], y = top_y + LABEL_OFFSET)
    ed <- cluster_edges(sn_all %>% dplyr::filter(cluster==cid))
    if (!is.null(ed)) { p <- lc %>% dplyr::select(seq_node,x,y)
      edges_all[[i]] <- ed %>% dplyr::left_join(p,by=c("from"="seq_node")) %>%
        dplyr::left_join(p,by=c("to"="seq_node"),suffix=c("",".e")) %>% dplyr::rename(xend=x.e,yend=y.e) }
  }
  nodes <- bind_rows(coords) %>% dplyr::left_join(sn_all %>% dplyr::select(-node_r), by=c("seq_node","cluster"))
  list(nodes = nodes, edges = bind_rows(edges_all), labels = bind_rows(labels))
}

# ── 4. Render one figure (donor overlay + cluster labels) ────────────────────
render <- function(keep_clusters, tag, title) {
  if (length(keep_clusters) == 0) { message(sprintf("[%s] no clusters, skipped", tag)); return(invisible()) }
  geo <- build_geometry(keep_clusters); nodes <- geo$nodes; edges <- geo$edges; labels <- geo$labels
  e_far  <- edges %>% dplyr::filter(!near); e_near <- edges %>% dplyr::filter(near)
  p <- ggplot() +
    geom_segment(data=e_far,  aes(x,y,xend=xend,yend=yend), color=EDGE_FAR_COL,  linewidth=EDGE_FAR_LW) +
    geom_segment(data=e_near, aes(x,y,xend=xend,yend=yend), color=EDGE_NEAR_COL, linewidth=EDGE_NEAR_LW) +
    ggforce::geom_circle(data=nodes, aes(x0=x,y0=y,r=node_r), fill="white", color="grey30", linewidth=0.2) +
    scatterpie::geom_scatterpie(data=nodes, aes(x=x,y=y,r=node_r),
           cols=c("d1","d2","d3","d4","d5"), color=NA) +
    geom_text(data=labels, aes(x=x, y=y, label=paste("cluster", cluster)),
              size=LABEL_SIZE, color="grey20") +
    scale_fill_manual(values=setNames(unname(donor_color_map),c("d1","d2","d3","d4","d5")),
                      labels=c("D1","D2","D3","D4","D5"), name="Donor") +
    coord_equal() + theme_void() + ggtitle(title) + theme(plot.title=element_text(size=14))
  ggsave(file.path(out_dir, sprintf("%s.png", tag)), p, width=15, height=13, dpi=300, bg="white")
  message(sprintf("[%s] %d clusters, %d nodes, %d edges (%d near/%d far)",
                  tag, dplyr::n_distinct(nodes$cluster), nrow(nodes), nrow(edges),
                  sum(edges$near), sum(!edges$near)))
}

# ── 5. Cluster sets: {all, >=2 seq} x {pp65, IE1} ────────────────────────────
pp65_all  <- cluster_meta %>% dplyr::filter(dom_spec=="pp65") %>% dplyr::pull(cluster)
ie1_all   <- cluster_meta %>% dplyr::filter(dom_spec=="IE1")  %>% dplyr::pull(cluster)
pp65_min2 <- cluster_meta %>% dplyr::filter(dom_spec=="pp65", n_seq>=2) %>% dplyr::pull(cluster)
ie1_min2  <- cluster_meta %>% dplyr::filter(dom_spec=="IE1",  n_seq>=2) %>% dplyr::pull(cluster)

render(pp65_all,  "all_pp65_donor",     "et64 constellations — pp65-dominant clusters (all) — donor")
render(ie1_all,   "all_IE1_donor",      "et64 constellations — IE1-dominant clusters (all) — donor")
render(pp65_min2, "min2seq_pp65_donor", "et64 constellations — pp65-dominant clusters (>=2 sequences) — donor")
render(ie1_min2,  "min2seq_IE1_donor",  "et64 constellations — IE1-dominant clusters (>=2 sequences) — donor")



#
#
# Constellation — continuous overlays: kNN bias effect + metacluster % (per node)
#   Reuses the by_specificity layout (log-sized nodes, two-tier MST edges, labels).
#   kNN: same_cluster_effect per cluster, blue-white-red diverging (white=0).
#   Metacluster: per-node % of cells in each phenotype, sequential white->color.
#   Files: {cutoff}_{specificity}_{overlay}.png in .../constellation/by_specificity/
#
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr); library(readr); library(tidyr)

conga_dir <- "exports/conga/mdata_CMV"
out_dir   <- "exports/TCR_clusters/constellation/by_specificity"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g_data <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
node_tbl <- g_data$nodes; D <- g_data$D
node_tbl <- node_tbl %>% dplyr::mutate(donorS = paste0("D", str_extract(as.character(donor), "\\d")))

# ── Tunable settings ─────────────────────────────────────────────────────────
EDGE_CUT       <- 250
DIST_TRANSFORM <- function(d) log1p(d)
NODE_R_MIN     <- 0.35
NODE_R_MAX     <- 1.40
CLUST_PAD      <- 3
GAP_MIN        <- 0.25
EDGE_FAR_COL   <- "grey80"; EDGE_FAR_LW   <- 0.25
EDGE_NEAR_COL  <- "grey35"; EDGE_NEAR_LW  <- 0.7
LABEL_OFFSET   <- 0.8
LABEL_SIZE     <- 3.5
NA_FILL        <- "grey85"   # nodes with no value (e.g. cluster not in kNN results)
META_LEVELS    <- c("Tcm","Tem","T-CMC1","NK-like","ISG+")

# ── 1. Per-cloneId metacluster cell counts (from the object) ─────────────────
# tally each cloneId's cells across metaclusters, keep clustered non-notCalled clones
meta_counts <- mdata_CMV@meta.data %>%
  dplyr::mutate(cloneId = as.character(cloneId)) %>%
  dplyr::filter(!is.na(TCR_cluster_et64), !grepl("notCalled", cloneId), !is.na(metacluster)) %>%
  dplyr::count(cloneId, metacluster, name="n") %>%
  tidyr::pivot_wider(names_from=metacluster, values_from=n, values_fill=0)
# ensure all 5 columns exist even if some absent
for (m in META_LEVELS) if (!m %in% names(meta_counts)) meta_counts[[m]] <- 0

# ── 2. Unique-sequence nodes: cells, donor split, dominant spec, metacluster % ──
# attach per-cloneId metacluster counts to node_tbl, then aggregate to sequence nodes
node_meta <- node_tbl %>% dplyr::left_join(meta_counts, by="cloneId")

seq_nodes <- node_meta %>%
  dplyr::group_by(cluster, mrow) %>%
  dplyr::summarise(
    total_cells=sum(n_cells), spec=names(which.max(table(spec))),
    d1=sum(n_cells[donorS=="D1"]), d2=sum(n_cells[donorS=="D2"]),
    d3=sum(n_cells[donorS=="D3"]), d4=sum(n_cells[donorS=="D4"]),
    d5=sum(n_cells[donorS=="D5"]),
    Tcm=sum(Tcm), Tem=sum(Tem), `T-CMC1`=sum(`T-CMC1`),
    `NK-like`=sum(`NK-like`), `ISG+`=sum(`ISG+`),
    .groups="drop") %>%
  dplyr::mutate(seq_node = paste0("s", dplyr::row_number()),
                meta_total = Tcm+Tem+`T-CMC1`+`NK-like`+`ISG+`)
# per-node metacluster PERCENT (of that node's cells with a metacluster label)
for (m in META_LEVELS) {
  pct_col <- paste0("pct_", m)
  seq_nodes[[pct_col]] <- ifelse(seq_nodes$meta_total>0, 100*seq_nodes[[m]]/seq_nodes$meta_total, NA_real_)
}

# log-value radius (global range) + dominant spec per cluster
lc_min <- log(min(seq_nodes$total_cells)); lc_max <- log(max(seq_nodes$total_cells))
seq_nodes <- seq_nodes %>%
  dplyr::mutate(node_r = NODE_R_MIN + (NODE_R_MAX-NODE_R_MIN)*(log(total_cells)-lc_min)/(lc_max-lc_min))
cluster_meta <- seq_nodes %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n_seq=dplyr::n(),
                   dom_spec={ w <- tapply(total_cells, spec, sum); names(w)[which.max(w)] }, .groups="drop")

# ── 3. kNN same_cluster_effect per cluster ───────────────────────────────────
knn <- read_csv(file.path(conga_dir, "..", "conga_exports", "tcrdist_phenotype_test",
                          "tcrdist_phenotype_specstrat_results.csv"), comment="#", show_col_types=FALSE) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  dplyr::select(cluster, knn_effect = same_cluster_effect)
seq_nodes <- seq_nodes %>% dplyr::mutate(cluster = as.character(cluster)) %>%
  dplyr::left_join(knn, by="cluster")
message(sprintf("Clusters: %d | pp65 %d, IE1 %d | with kNN score: %d",
                nrow(cluster_meta), sum(cluster_meta$dom_spec=="pp65"),
                sum(cluster_meta$dom_spec=="IE1"), sum(!is.na(unique(seq_nodes$knn_effect)))))

# ── 4. Layout + two-tier MST edges ───────────────────────────────────────────
layout_one <- function(sn) {
  n <- nrow(sn); if (n==1) return(tibble(seq_node=sn$seq_node, x=0, y=0))
  sub <- D[sn$mrow, sn$mrow, drop=FALSE]; sub_t <- DIST_TRANSFORM(sub)
  if (n==2) { xy <- matrix(c(0,0, sub_t[1,2],0), nrow=2, byrow=TRUE) }
  else { gg <- graph_from_adjacency_matrix(sub_t, mode="undirected", weighted=TRUE, diag=FALSE)
    set.seed(1); xy <- layout_with_kk(gg, weights=E(gg)$weight) }
  xy <- scale(xy, center=TRUE, scale=FALSE); tibble(seq_node=sn$seq_node, x=xy[,1], y=xy[,2])
}
layout_cluster_local <- function(sn) {
  lay <- layout_one(sn) %>% dplyr::left_join(sn, by="seq_node")
  if (nrow(lay)==1) return(list(coords=lay %>% dplyr::transmute(seq_node,x=0,y=0,node_r), bound=lay$node_r[1]))
  sp <- max(sqrt(lay$x^2+lay$y^2)); if (sp==0||!is.finite(sp)) sp <- 1
  lay$x <- lay$x/sp; lay$y <- lay$y/sp
  set.seed(1); lay$x <- lay$x+rnorm(nrow(lay),0,0.01); lay$y <- lay$y+rnorm(nrow(lay),0,0.01)
  rep <- packcircles::circleRepelLayout(data.frame(x=lay$x,y=lay$y,radius=lay$node_r+GAP_MIN/2),
    xysizecols=1:3, sizetype="radius", maxiter=1000)$layout
  lay$x <- rep$x; lay$y <- rep$y
  cx <- mean(range(lay$x)); cy <- mean(range(lay$y)); lay$x <- lay$x-cx; lay$y <- lay$y-cy
  list(coords=lay %>% dplyr::select(seq_node,x,y,node_r), bound=max(sqrt(lay$x^2+lay$y^2)+lay$node_r))
}
cluster_edges <- function(sn) {
  if (nrow(sn)<2) return(NULL)
  sub <- D[sn$mrow, sn$mrow, drop=FALSE]
  gg <- graph_from_adjacency_matrix(sub+1, mode="undirected", weighted=TRUE, diag=FALSE)
  m <- igraph::mst(gg); if (igraph::ecount(m)==0) return(NULL)
  el <- igraph::as_edgelist(m); i1 <- as.integer(el[,1]); i2 <- as.integer(el[,2])
  tibble(from=sn$seq_node[i1], to=sn$seq_node[i2], near=sub[cbind(i1,i2)] <= EDGE_CUT)
}

# ── 5. Build geometry (nodes carry all overlay columns) ──────────────────────
build_geometry <- function(keep_clusters) {
  sn_all <- seq_nodes %>% dplyr::filter(cluster %in% keep_clusters)
  cl_ids <- sn_all %>% dplyr::count(cluster, name="n_seq") %>% dplyr::arrange(dplyr::desc(n_seq)) %>% dplyr::pull(cluster)
  locals <- lapply(cl_ids, function(cid) layout_cluster_local(sn_all %>% dplyr::filter(cluster==cid)))
  bounds <- sapply(locals, function(z) z$bound)
  circ <- packcircles::circleProgressiveLayout(bounds + CLUST_PAD, sizetype="radius")
  coords <- list(); edges_all <- list(); labels <- list()
  for (i in seq_along(cl_ids)) {
    cid <- cl_ids[i]; lc <- locals[[i]]$coords %>% dplyr::mutate(x=circ$x[i]+x, y=circ$y[i]+y, cluster=cid)
    coords[[i]] <- lc
    labels[[i]] <- tibble(cluster=cid, x=circ$x[i], y=max(lc$y+lc$node_r)+LABEL_OFFSET)
    ed <- cluster_edges(sn_all %>% dplyr::filter(cluster==cid))
    if (!is.null(ed)) { p <- lc %>% dplyr::select(seq_node,x,y)
      edges_all[[i]] <- ed %>% dplyr::left_join(p,by=c("from"="seq_node")) %>%
        dplyr::left_join(p,by=c("to"="seq_node"),suffix=c("",".e")) %>% dplyr::rename(xend=x.e,yend=y.e) }
  }
  nodes <- bind_rows(coords) %>% dplyr::left_join(sn_all %>% dplyr::select(-node_r), by=c("seq_node","cluster"))
  list(nodes=nodes, edges=bind_rows(edges_all), labels=bind_rows(labels))
}

# ── 6. Generalized continuous-overlay renderer ───────────────────────────────
#   value_col: node column to color by. scale_type: "diverging" (blue-white-red, 0-centered)
#   or "sequential" (white -> high_col). 
render_continuous <- function(keep_clusters, tag, title, value_col, scale_type, high_col="#8B0000") {
  if (length(keep_clusters)==0) { message(sprintf("[%s] no clusters, skipped", tag)); return(invisible()) }
  geo <- build_geometry(keep_clusters); nodes <- geo$nodes; edges <- geo$edges; labels <- geo$labels
  nodes$val <- nodes[[value_col]]
  e_far <- edges %>% dplyr::filter(!near); e_near <- edges %>% dplyr::filter(near)
  p <- ggplot() +
    geom_segment(data=e_far,  aes(x,y,xend=xend,yend=yend), color=EDGE_FAR_COL,  linewidth=EDGE_FAR_LW) +
    geom_segment(data=e_near, aes(x,y,xend=xend,yend=yend), color=EDGE_NEAR_COL, linewidth=EDGE_NEAR_LW) +
    ggforce::geom_circle(data=nodes, aes(x0=x,y0=y,r=node_r,fill=val), color="grey30", linewidth=0.2)
  if (scale_type=="diverging") {
    lim <- max(abs(nodes$val), na.rm=TRUE)
    p <- p + scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
                                  limits=c(-lim, lim), na.value=NA_FILL, name=value_col)
  } else {
    p <- p + scale_fill_gradient(low="white", high=high_col, na.value=NA_FILL, name=value_col)
  }
  p <- p +
    geom_text(data=labels, aes(x=x, y=y, label=paste("cluster", cluster)), size=LABEL_SIZE, color="grey20") +
    coord_equal() + theme_void() + ggtitle(title) + theme(plot.title=element_text(size=14))
  ggsave(file.path(out_dir, sprintf("%s.png", tag)), p, width=15, height=13, dpi=300, bg="white")
  message(sprintf("[%s] %d clusters, %d nodes", tag, dplyr::n_distinct(nodes$cluster), nrow(nodes)))
}

# ── 7. Cluster sets and render all overlays ──────────────────────────────────
sets <- list(
  list(cl = cluster_meta %>% dplyr::filter(dom_spec=="pp65") %>% dplyr::pull(cluster),               cut="all",     spec="pp65"),
  list(cl = cluster_meta %>% dplyr::filter(dom_spec=="IE1")  %>% dplyr::pull(cluster),               cut="all",     spec="IE1"),
  list(cl = cluster_meta %>% dplyr::filter(dom_spec=="pp65", n_seq>=2) %>% dplyr::pull(cluster),     cut="min2seq", spec="pp65"),
  list(cl = cluster_meta %>% dplyr::filter(dom_spec=="IE1",  n_seq>=2) %>% dplyr::pull(cluster),     cut="min2seq", spec="IE1")
)
meta_palette <- c(Tcm="#E8A0A0", Tem="#7B8B3A", `T-CMC1`="#1B9E8A", `NK-like`="#888888", `ISG+`="#C050C0")

for (s in sets) {
  base <- sprintf("%s_%s", s$cut, s$spec)
  # kNN bias (diverging, 0-centered)
  render_continuous(s$cl, sprintf("%s_kNNbias", base),
                    sprintf("et64 — %s-dominant (%s) — kNN same-cluster effect", s$spec, s$cut),
                    "knn_effect", "diverging")
  # 5 metacluster % (sequential)
  for (m in META_LEVELS) {
    safe_m <- gsub("[^A-Za-z0-9]", "", m)   # file-safe (T-CMC1 -> TCMC1, ISG+ -> ISG)
    render_continuous(s$cl, sprintf("%s_meta%s", base, safe_m),
                      sprintf("et64 — %s-dominant (%s) — %% %s", s$spec, s$cut, m),
                      paste0("pct_", m), "sequential", high_col = meta_palette[[m]])
  }
}


#
#
# Audit of constellation clusters
#
#

library(dplyr); library(stringr); library(tidyr)

# Pull every cell in a plotted cluster, parse its alpha/beta CDR3 from t_cdr3s_aa
aud <- mdata_CMV@meta.data %>%
  dplyr::filter(!is.na(TCR_cluster_et64)) %>%
  dplyr::transmute(
    cluster   = as.character(TCR_cluster_et64),
    cloneId   = as.character(cloneId),
    on_clone  = !is.na(cloneId) & !grepl("notCalled", cloneId),   # real clone vs notCalled
    tra       = str_match(t_cdr3s_aa, "TRA:([^;]+)")[, 2],
    trb       = str_match(t_cdr3s_aa, "TRB:([^;]+)")[, 2])

# ---- Per-cluster diversity, split by clone standing ----
summ <- aud %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    n_cells        = dplyr::n(),
    n_clones       = dplyr::n_distinct(cloneId[on_clone]),
    # ALPHA
    a_all          = dplyr::n_distinct(tra[!is.na(tra)]),                 # any cell
    a_on_clone     = dplyr::n_distinct(tra[!is.na(tra) & on_clone]),      # on a real clone
    a_cells_w_alpha= sum(!is.na(tra)),
    # BETA
    b_all          = dplyr::n_distinct(trb[!is.na(trb)]),
    b_on_clone     = dplyr::n_distinct(trb[!is.na(trb) & on_clone]),
    # PAIRED combos (what the constellation node ~ is): distinct alpha+beta on real clones
    pair_on_clone  = dplyr::n_distinct(paste(tra, trb)[!is.na(tra) & !is.na(trb) & on_clone]),
    pair_all       = dplyr::n_distinct(paste(tra, trb)[!is.na(tra) & !is.na(trb)]),
    .groups = "drop")

# ---- How many of the "extra" sequences are singletons (1 cell)? ----
# alpha strings appearing on only 1 cell in the cluster
a_singletons <- aud %>%
  dplyr::filter(!is.na(tra)) %>%
  dplyr::count(cluster, tra, name = "cells") %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(a_singleton_seqs = sum(cells == 1), .groups = "drop")

out <- summ %>%
  dplyr::left_join(a_singletons, by = "cluster") %>%
  dplyr::mutate(a_singleton_seqs = tidyr::replace_na(a_singleton_seqs, 0L)) %>%
  dplyr::arrange(dplyr::desc(n_cells))

print(out, n = Inf)

# Quick headline: where does card-style alpha diversity exceed on-clone paired nodes?
cat("\n--- Clusters where raw alpha diversity > paired-on-clone nodes ---\n")
out %>%
  dplyr::mutate(inflation = a_all - pair_on_clone) %>%
  dplyr::filter(inflation > 0) %>%
  dplyr::arrange(dplyr::desc(inflation)) %>%
  dplyr::select(cluster, n_cells, n_clones,
                a_all, a_on_clone, a_singleton_seqs,
                b_all, b_on_clone, pair_on_clone, pair_all) %>%
  print(n = Inf)




#
#
# Constellation — continuous overlays, GLOBAL-CLONES-ONLY subset
#   Identical to the by_specificity continuous-overlay block, EXCEPT:
#     (1) node_tbl is downsampled to clones with a non-NA Global_Clone_ID,
#     (2) output goes to constellation/only_global_clones/.
#   Everything else (layout, packing, two-tier MST edges, kNN + metacluster
#   overlays, cutoffs, palettes) is unchanged.
#
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr); library(readr); library(tidyr)

conga_dir <- "exports/conga/mdata_CMV"
out_dir   <- "exports/TCR_clusters/constellation/only_global_clones"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g_data <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
node_tbl <- g_data$nodes; D <- g_data$D

# ── GLOBAL-ONLY DOWNSAMPLE ───────────────────────────────────────────────────
# Keep only nodes whose clone carries a Global_Clone_ID. D is unchanged (mrow
# still indexes the full matrix), so all distances/edges remain correct.
n_before <- nrow(node_tbl)
node_tbl <- node_tbl %>% dplyr::filter(is_global)
message(sprintf("Global-only downsample: %d -> %d nodes | clusters %d -> %d",
                n_before, nrow(node_tbl),
                dplyr::n_distinct(g_data$nodes$cluster), dplyr::n_distinct(node_tbl$cluster)))

node_tbl <- node_tbl %>% dplyr::mutate(donorS = paste0("D", str_extract(as.character(donor), "\\d")))

# ── Tunable settings (unchanged) ─────────────────────────────────────────────
EDGE_CUT       <- 250
DIST_TRANSFORM <- function(d) log1p(d)
NODE_R_MIN     <- 0.35
NODE_R_MAX     <- 1.40
CLUST_PAD      <- 3
GAP_MIN        <- 0.25
EDGE_FAR_COL   <- "grey80"; EDGE_FAR_LW   <- 0.25
EDGE_NEAR_COL  <- "grey35"; EDGE_NEAR_LW  <- 0.7
LABEL_OFFSET   <- 0.8
LABEL_SIZE     <- 3.5
NA_FILL        <- "grey85"
META_LEVELS    <- c("Tcm","Tem","T-CMC1","NK-like","ISG+")

# ── 1. Per-cloneId metacluster cell counts (unchanged) ───────────────────────
meta_counts <- mdata_CMV@meta.data %>%
  dplyr::mutate(cloneId = as.character(cloneId)) %>%
  dplyr::filter(!is.na(TCR_cluster_et64), !grepl("notCalled", cloneId), !is.na(metacluster)) %>%
  dplyr::count(cloneId, metacluster, name="n") %>%
  tidyr::pivot_wider(names_from=metacluster, values_from=n, values_fill=0)
for (m in META_LEVELS) if (!m %in% names(meta_counts)) meta_counts[[m]] <- 0

# ── 2. Unique-sequence nodes (unchanged, but built on the global-only node_tbl) ──
node_meta <- node_tbl %>% dplyr::left_join(meta_counts, by="cloneId")

seq_nodes <- node_meta %>%
  dplyr::group_by(cluster, mrow) %>%
  dplyr::summarise(
    total_cells=sum(n_cells), spec=names(which.max(table(spec))),
    d1=sum(n_cells[donorS=="D1"]), d2=sum(n_cells[donorS=="D2"]),
    d3=sum(n_cells[donorS=="D3"]), d4=sum(n_cells[donorS=="D4"]),
    d5=sum(n_cells[donorS=="D5"]),
    Tcm=sum(Tcm), Tem=sum(Tem), `T-CMC1`=sum(`T-CMC1`),
    `NK-like`=sum(`NK-like`), `ISG+`=sum(`ISG+`),
    .groups="drop") %>%
  dplyr::mutate(seq_node = paste0("s", dplyr::row_number()),
                meta_total = Tcm+Tem+`T-CMC1`+`NK-like`+`ISG+`)
for (m in META_LEVELS) {
  pct_col <- paste0("pct_", m)
  seq_nodes[[pct_col]] <- ifelse(seq_nodes$meta_total>0, 100*seq_nodes[[m]]/seq_nodes$meta_total, NA_real_)
}

# radius uses the GLOBAL-ONLY cell range (see note below re: comparability)
lc_min <- log(min(seq_nodes$total_cells)); lc_max <- log(max(seq_nodes$total_cells))
seq_nodes <- seq_nodes %>%
  dplyr::mutate(node_r = NODE_R_MIN + (NODE_R_MAX-NODE_R_MIN)*(log(total_cells)-lc_min)/(lc_max-lc_min))
cluster_meta <- seq_nodes %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n_seq=dplyr::n(),
                   dom_spec={ w <- tapply(total_cells, spec, sum); names(w)[which.max(w)] }, .groups="drop")

# ── 3. kNN same_cluster_effect per cluster (unchanged) ───────────────────────
knn <- read_csv(file.path(conga_dir, "..", "conga_exports", "tcrdist_phenotype_test",
                          "tcrdist_phenotype_specstrat_results.csv"), comment="#", show_col_types=FALSE) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  dplyr::select(cluster, knn_effect = same_cluster_effect)
seq_nodes <- seq_nodes %>% dplyr::mutate(cluster = as.character(cluster)) %>%
  dplyr::left_join(knn, by="cluster")
message(sprintf("Clusters: %d | pp65 %d, IE1 %d | with kNN score: %d",
                nrow(cluster_meta), sum(cluster_meta$dom_spec=="pp65"),
                sum(cluster_meta$dom_spec=="IE1"), sum(!is.na(unique(seq_nodes$knn_effect)))))

# ── 4. Layout + two-tier MST edges (unchanged) ───────────────────────────────
layout_one <- function(sn) {
  n <- nrow(sn); if (n==1) return(tibble(seq_node=sn$seq_node, x=0, y=0))
  sub <- D[sn$mrow, sn$mrow, drop=FALSE]; sub_t <- DIST_TRANSFORM(sub)
  if (n==2) { xy <- matrix(c(0,0, sub_t[1,2],0), nrow=2, byrow=TRUE) }
  else { gg <- graph_from_adjacency_matrix(sub_t, mode="undirected", weighted=TRUE, diag=FALSE)
    set.seed(1); xy <- layout_with_kk(gg, weights=E(gg)$weight) }
  xy <- scale(xy, center=TRUE, scale=FALSE); tibble(seq_node=sn$seq_node, x=xy[,1], y=xy[,2])
}
layout_cluster_local <- function(sn) {
  lay <- layout_one(sn) %>% dplyr::left_join(sn, by="seq_node")
  if (nrow(lay)==1) return(list(coords=lay %>% dplyr::transmute(seq_node,x=0,y=0,node_r), bound=lay$node_r[1]))
  sp <- max(sqrt(lay$x^2+lay$y^2)); if (sp==0||!is.finite(sp)) sp <- 1
  lay$x <- lay$x/sp; lay$y <- lay$y/sp
  set.seed(1); lay$x <- lay$x+rnorm(nrow(lay),0,0.01); lay$y <- lay$y+rnorm(nrow(lay),0,0.01)
  rep <- packcircles::circleRepelLayout(data.frame(x=lay$x,y=lay$y,radius=lay$node_r+GAP_MIN/2),
    xysizecols=1:3, sizetype="radius", maxiter=1000)$layout
  lay$x <- rep$x; lay$y <- rep$y
  cx <- mean(range(lay$x)); cy <- mean(range(lay$y)); lay$x <- lay$x-cx; lay$y <- lay$y-cy
  list(coords=lay %>% dplyr::select(seq_node,x,y,node_r), bound=max(sqrt(lay$x^2+lay$y^2)+lay$node_r))
}
cluster_edges <- function(sn) {
  if (nrow(sn)<2) return(NULL)
  sub <- D[sn$mrow, sn$mrow, drop=FALSE]
  gg <- graph_from_adjacency_matrix(sub+1, mode="undirected", weighted=TRUE, diag=FALSE)
  m <- igraph::mst(gg); if (igraph::ecount(m)==0) return(NULL)
  el <- igraph::as_edgelist(m); i1 <- as.integer(el[,1]); i2 <- as.integer(el[,2])
  tibble(from=sn$seq_node[i1], to=sn$seq_node[i2], near=sub[cbind(i1,i2)] <= EDGE_CUT)
}

# ── 5. Build geometry (unchanged) ────────────────────────────────────────────
build_geometry <- function(keep_clusters) {
  sn_all <- seq_nodes %>% dplyr::filter(cluster %in% keep_clusters)
  cl_ids <- sn_all %>% dplyr::count(cluster, name="n_seq") %>% dplyr::arrange(dplyr::desc(n_seq)) %>% dplyr::pull(cluster)
  locals <- lapply(cl_ids, function(cid) layout_cluster_local(sn_all %>% dplyr::filter(cluster==cid)))
  bounds <- sapply(locals, function(z) z$bound)
  circ <- packcircles::circleProgressiveLayout(bounds + CLUST_PAD, sizetype="radius")
  coords <- list(); edges_all <- list(); labels <- list()
  for (i in seq_along(cl_ids)) {
    cid <- cl_ids[i]; lc <- locals[[i]]$coords %>% dplyr::mutate(x=circ$x[i]+x, y=circ$y[i]+y, cluster=cid)
    coords[[i]] <- lc
    labels[[i]] <- tibble(cluster=cid, x=circ$x[i], y=max(lc$y+lc$node_r)+LABEL_OFFSET)
    ed <- cluster_edges(sn_all %>% dplyr::filter(cluster==cid))
    if (!is.null(ed)) { p <- lc %>% dplyr::select(seq_node,x,y)
      edges_all[[i]] <- ed %>% dplyr::left_join(p,by=c("from"="seq_node")) %>%
        dplyr::left_join(p,by=c("to"="seq_node"),suffix=c("",".e")) %>% dplyr::rename(xend=x.e,yend=y.e) }
  }
  nodes <- bind_rows(coords) %>% dplyr::left_join(sn_all %>% dplyr::select(-node_r), by=c("seq_node","cluster"))
  list(nodes=nodes, edges=bind_rows(edges_all), labels=bind_rows(labels))
}

# ── 6. Generalized continuous-overlay renderer (unchanged) ───────────────────
render_continuous <- function(keep_clusters, tag, title, value_col, scale_type, high_col="#8B0000") {
  if (length(keep_clusters)==0) { message(sprintf("[%s] no clusters, skipped", tag)); return(invisible()) }
  geo <- build_geometry(keep_clusters); nodes <- geo$nodes; edges <- geo$edges; labels <- geo$labels
  nodes$val <- nodes[[value_col]]
  e_far <- edges %>% dplyr::filter(!near); e_near <- edges %>% dplyr::filter(near)
  p <- ggplot() +
    geom_segment(data=e_far,  aes(x,y,xend=xend,yend=yend), color=EDGE_FAR_COL,  linewidth=EDGE_FAR_LW) +
    geom_segment(data=e_near, aes(x,y,xend=xend,yend=yend), color=EDGE_NEAR_COL, linewidth=EDGE_NEAR_LW) +
    ggforce::geom_circle(data=nodes, aes(x0=x,y0=y,r=node_r,fill=val), color="grey30", linewidth=0.2)
  if (scale_type=="diverging") {
    lim <- max(abs(nodes$val), na.rm=TRUE)
    p <- p + scale_fill_gradient2(low="#2166AC", mid="white", high="#B2182B", midpoint=0,
                                  limits=c(-lim, lim), na.value=NA_FILL, name=value_col)
  } else {
    p <- p + scale_fill_gradient(low="white", high=high_col, na.value=NA_FILL, name=value_col)
  }
  p <- p +
    geom_text(data=labels, aes(x=x, y=y, label=paste("cluster", cluster)), size=LABEL_SIZE, color="grey20") +
    coord_equal() + theme_void() + ggtitle(title) + theme(plot.title=element_text(size=14))
  ggsave(file.path(out_dir, sprintf("%s.png", tag)), p, width=15, height=13, dpi=300, bg="white")
  message(sprintf("[%s] %d clusters, %d nodes", tag, dplyr::n_distinct(nodes$cluster), nrow(nodes)))
}

# ── 7. Cluster sets and render all overlays (unchanged, incl. donor overlay) ──
sets <- list(
  list(cl = cluster_meta %>% dplyr::filter(dom_spec=="pp65") %>% dplyr::pull(cluster),           cut="all",     spec="pp65"),
  list(cl = cluster_meta %>% dplyr::filter(dom_spec=="IE1")  %>% dplyr::pull(cluster),           cut="all",     spec="IE1"),
  list(cl = cluster_meta %>% dplyr::filter(dom_spec=="pp65", n_seq>=2) %>% dplyr::pull(cluster), cut="min2seq", spec="pp65"),
  list(cl = cluster_meta %>% dplyr::filter(dom_spec=="IE1",  n_seq>=2) %>% dplyr::pull(cluster), cut="min2seq", spec="IE1")
)
meta_palette <- c(Tcm="#E8A0A0", Tem="#7B8B3A", `T-CMC1`="#1B9E8A", `NK-like`="#888888", `ISG+`="#C050C0")

for (s in sets) {
  base <- sprintf("%s_%s", s$cut, s$spec)
  render_continuous(s$cl, sprintf("%s_kNNbias", base),
                    sprintf("et64 (Global clones) — %s-dominant (%s) — kNN same-cluster effect", s$spec, s$cut),
                    "knn_effect", "diverging")
  for (m in META_LEVELS) {
    safe_m <- gsub("[^A-Za-z0-9]", "", m)
    render_continuous(s$cl, sprintf("%s_meta%s", base, safe_m),
                      sprintf("et64 (Global clones) — %s-dominant (%s) — %% %s", s$spec, s$cut, m),
                      paste0("pct_", m), "sequential", high_col = meta_palette[[m]])
  }
}



#
#
# Constellation — donor-pie overlay, GLOBAL-CLONES-ONLY subset
#   Identical to the by_specificity donor block, EXCEPT:
#     (1) node_tbl downsampled to clones with a non-NA Global_Clone_ID,
#     (2) output goes to constellation/only_global_clones/.
#   Node = unique sequence; size = log(cells); two-tier MST edges; donor pie wedges.
#
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr)

conga_dir <- "exports/conga/mdata_CMV"
out_dir   <- "exports/TCR_clusters/constellation/only_global_clones"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

g_data <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
node_tbl <- g_data$nodes; D <- g_data$D
donor_color_map <- c("D1"="#1b9e77","D2"="#d95f02","D3"="#7570b3","D4"="#e7298a","D5"="#66a61e")

# ── GLOBAL-ONLY DOWNSAMPLE ───────────────────────────────────────────────────
n_before <- nrow(node_tbl)
node_tbl <- node_tbl %>% dplyr::filter(is_global)
message(sprintf("Global-only downsample: %d -> %d nodes | clusters %d -> %d",
                n_before, nrow(node_tbl),
                dplyr::n_distinct(g_data$nodes$cluster), dplyr::n_distinct(node_tbl$cluster)))

node_tbl <- node_tbl %>% dplyr::mutate(donorS = paste0("D", str_extract(as.character(donor), "\\d")))

# ── Tunable settings (unchanged) ─────────────────────────────────────────────
EDGE_CUT       <- 250
DIST_TRANSFORM <- function(d) log1p(d)
NODE_R_MIN     <- 0.35
NODE_R_MAX     <- 1.40
CLUST_PAD      <- 3
GAP_MIN        <- 0.25
EDGE_FAR_COL   <- "grey80"; EDGE_FAR_LW   <- 0.25
EDGE_NEAR_COL  <- "grey35"; EDGE_NEAR_LW  <- 0.7
LABEL_OFFSET   <- 0.8
LABEL_SIZE     <- 3.5

# ── 1. Unique-sequence nodes + log-value radius + per-cluster dominant spec ──
seq_nodes <- node_tbl %>%
  dplyr::group_by(cluster, mrow) %>%
  dplyr::summarise(total_cells=sum(n_cells), spec=names(which.max(table(spec))),
                   d1=sum(n_cells[donorS=="D1"]), d2=sum(n_cells[donorS=="D2"]),
                   d3=sum(n_cells[donorS=="D3"]), d4=sum(n_cells[donorS=="D4"]),
                   d5=sum(n_cells[donorS=="D5"]), .groups="drop") %>%
  dplyr::mutate(seq_node = paste0("s", dplyr::row_number()))

# radius uses the GLOBAL-ONLY cell range (see note below re: comparability)
lc_min <- log(min(seq_nodes$total_cells)); lc_max <- log(max(seq_nodes$total_cells))
seq_nodes <- seq_nodes %>%
  dplyr::mutate(node_r = NODE_R_MIN + (NODE_R_MAX-NODE_R_MIN) *
                  (log(total_cells) - lc_min) / (lc_max - lc_min))

cluster_meta <- seq_nodes %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(n_seq = dplyr::n(),
                   dom_spec = { w <- tapply(total_cells, spec, sum)
                                names(w)[which.max(w)] }, .groups="drop")
message(sprintf("Clusters: %d | pp65-dom %d, IE1-dom %d | with >=2 seq: %d",
                nrow(cluster_meta), sum(cluster_meta$dom_spec=="pp65"),
                sum(cluster_meta$dom_spec=="IE1"), sum(cluster_meta$n_seq>=2)))

# ── 2. Layout + two-tier MST edges (unchanged) ───────────────────────────────
layout_one <- function(sn) {
  n <- nrow(sn)
  if (n == 1) return(tibble(seq_node=sn$seq_node, x=0, y=0))
  sub   <- D[sn$mrow, sn$mrow, drop=FALSE]; sub_t <- DIST_TRANSFORM(sub)
  if (n == 2) { xy <- matrix(c(0,0, sub_t[1,2],0), nrow=2, byrow=TRUE) }
  else { gg <- graph_from_adjacency_matrix(sub_t, mode="undirected", weighted=TRUE, diag=FALSE)
    set.seed(1); xy <- layout_with_kk(gg, weights = E(gg)$weight) }
  xy <- scale(xy, center=TRUE, scale=FALSE)
  tibble(seq_node=sn$seq_node, x=xy[,1], y=xy[,2])
}
layout_cluster_local <- function(sn) {
  lay <- layout_one(sn) %>% dplyr::left_join(sn, by="seq_node")
  if (nrow(lay) == 1)
    return(list(coords = lay %>% dplyr::transmute(seq_node, x=0, y=0, node_r), bound = lay$node_r[1]))
  sp <- max(sqrt(lay$x^2 + lay$y^2)); if (sp==0 || !is.finite(sp)) sp <- 1
  lay$x <- lay$x/sp; lay$y <- lay$y/sp
  set.seed(1); lay$x <- lay$x + rnorm(nrow(lay),0,0.01); lay$y <- lay$y + rnorm(nrow(lay),0,0.01)
  rep <- packcircles::circleRepelLayout(
    data.frame(x=lay$x, y=lay$y, radius=lay$node_r + GAP_MIN/2),
    xysizecols=1:3, sizetype="radius", maxiter=1000)$layout
  lay$x <- rep$x; lay$y <- rep$y
  cx <- mean(range(lay$x)); cy <- mean(range(lay$y)); lay$x <- lay$x - cx; lay$y <- lay$y - cy
  bound <- max(sqrt(lay$x^2 + lay$y^2) + lay$node_r)
  list(coords = lay %>% dplyr::select(seq_node, x, y, node_r), bound = bound)
}
cluster_edges <- function(sn) {
  if (nrow(sn) < 2) return(NULL)
  sub <- D[sn$mrow, sn$mrow, drop=FALSE]
  gg  <- graph_from_adjacency_matrix(sub + 1, mode="undirected", weighted=TRUE, diag=FALSE)
  m   <- igraph::mst(gg); if (igraph::ecount(m) == 0) return(NULL)
  el <- igraph::as_edgelist(m); i1 <- as.integer(el[,1]); i2 <- as.integer(el[,2])
  tibble(from = sn$seq_node[i1], to = sn$seq_node[i2], near = sub[cbind(i1,i2)] <= EDGE_CUT)
}

# ── 3. Build geometry (unchanged) ────────────────────────────────────────────
build_geometry <- function(keep_clusters) {
  sn_all <- seq_nodes %>% dplyr::filter(cluster %in% keep_clusters)
  cl_ids <- sn_all %>% dplyr::count(cluster, name="n_seq") %>% dplyr::arrange(dplyr::desc(n_seq)) %>%
    dplyr::pull(cluster)
  locals <- lapply(cl_ids, function(cid) layout_cluster_local(sn_all %>% dplyr::filter(cluster==cid)))
  bounds <- sapply(locals, function(z) z$bound)
  circ <- packcircles::circleProgressiveLayout(bounds + CLUST_PAD, sizetype="radius")
  coords <- list(); edges_all <- list(); labels <- list()
  for (i in seq_along(cl_ids)) {
    cid <- cl_ids[i]; lc <- locals[[i]]$coords %>%
      dplyr::mutate(x = circ$x[i] + x, y = circ$y[i] + y, cluster = cid)
    coords[[i]] <- lc
    top_y <- max(lc$y + lc$node_r)
    labels[[i]] <- tibble(cluster = cid, x = circ$x[i], y = top_y + LABEL_OFFSET)
    ed <- cluster_edges(sn_all %>% dplyr::filter(cluster==cid))
    if (!is.null(ed)) { p <- lc %>% dplyr::select(seq_node,x,y)
      edges_all[[i]] <- ed %>% dplyr::left_join(p,by=c("from"="seq_node")) %>%
        dplyr::left_join(p,by=c("to"="seq_node"),suffix=c("",".e")) %>% dplyr::rename(xend=x.e,yend=y.e) }
  }
  nodes <- bind_rows(coords) %>% dplyr::left_join(sn_all %>% dplyr::select(-node_r), by=c("seq_node","cluster"))
  list(nodes = nodes, edges = bind_rows(edges_all), labels = bind_rows(labels))
}

# ── 4. Render one figure (donor overlay + cluster labels) (unchanged) ────────
render <- function(keep_clusters, tag, title) {
  if (length(keep_clusters) == 0) { message(sprintf("[%s] no clusters, skipped", tag)); return(invisible()) }
  geo <- build_geometry(keep_clusters); nodes <- geo$nodes; edges <- geo$edges; labels <- geo$labels
  e_far  <- edges %>% dplyr::filter(!near); e_near <- edges %>% dplyr::filter(near)
  p <- ggplot() +
    geom_segment(data=e_far,  aes(x,y,xend=xend,yend=yend), color=EDGE_FAR_COL,  linewidth=EDGE_FAR_LW) +
    geom_segment(data=e_near, aes(x,y,xend=xend,yend=yend), color=EDGE_NEAR_COL, linewidth=EDGE_NEAR_LW) +
    ggforce::geom_circle(data=nodes, aes(x0=x,y0=y,r=node_r), fill="white", color="grey30", linewidth=0.2) +
    scatterpie::geom_scatterpie(data=nodes, aes(x=x,y=y,r=node_r),
           cols=c("d1","d2","d3","d4","d5"), color=NA) +
    geom_text(data=labels, aes(x=x, y=y, label=paste("cluster", cluster)),
              size=LABEL_SIZE, color="grey20") +
    scale_fill_manual(values=setNames(unname(donor_color_map),c("d1","d2","d3","d4","d5")),
                      labels=c("D1","D2","D3","D4","D5"), name="Donor") +
    coord_equal() + theme_void() + ggtitle(title) + theme(plot.title=element_text(size=14))
  ggsave(file.path(out_dir, sprintf("%s.png", tag)), p, width=15, height=13, dpi=300, bg="white")
  message(sprintf("[%s] %d clusters, %d nodes, %d edges (%d near/%d far)",
                  tag, dplyr::n_distinct(nodes$cluster), nrow(nodes), nrow(edges),
                  sum(edges$near), sum(!edges$near)))
}

# ── 5. Cluster sets: {all, >=2 seq} x {pp65, IE1} (unchanged) ────────────────
pp65_all  <- cluster_meta %>% dplyr::filter(dom_spec=="pp65") %>% dplyr::pull(cluster)
ie1_all   <- cluster_meta %>% dplyr::filter(dom_spec=="IE1")  %>% dplyr::pull(cluster)
pp65_min2 <- cluster_meta %>% dplyr::filter(dom_spec=="pp65", n_seq>=2) %>% dplyr::pull(cluster)
ie1_min2  <- cluster_meta %>% dplyr::filter(dom_spec=="IE1",  n_seq>=2) %>% dplyr::pull(cluster)

render(pp65_all,  "all_pp65_donor",     "et64 (Global clones) — pp65-dominant clusters (all) — donor")
render(ie1_all,   "all_IE1_donor",      "et64 (Global clones) — IE1-dominant clusters (all) — donor")
render(pp65_min2, "min2seq_pp65_donor", "et64 (Global clones) — pp65-dominant clusters (>=2 sequences) — donor")
render(ie1_min2,  "min2seq_IE1_donor",  "et64 (Global clones) — IE1-dominant clusters (>=2 sequences) — donor")



#
#
# Correcting mislabeled CITEseq markers
#
#

library(Seurat)
library(Matrix)

pairs <- list(
  CCR7    = c("anti-human CCR7 (custom)",  "CCR7.1"),
  CXCR3   = c("anti-human CD183 (CXCR3)",  "CXCR3.1"),
  CD226   = c("anti-human CD226 (DNAM-1)", "CD226.1"),
  `PD-1`  = c("anti-human CD279 (PD-1)",   "PD-1"),
  `Tim-3` = c("anti-human CD366 (Tim-3)",  "Tim-3"),
  CD45RA  = c("anti-human CD45RA",         "CD45RA"),
  CD57    = c("anti-human CD57",           "CD57"),
  CD95    = c("anti-human CD95 (Fas)",     "CD95"),
  KLRG1   = c("anti-human KLRG1 (custom)", "KLRG1.1"),
  TIGIT   = c("anti-human TIGIT (VSTM3)",  "TIGIT.1")
)

# ---- mdata_CMV: collapse variant pairs into one canonical row ----
cts <- GetAssayData(mdata_CMV, assay = "FB", layer = "counts")
stopifnot(all(unlist(pairs, use.names = FALSE) %in% rownames(cts)))

canon     <- names(pairs)
collapsed <- matrix(0, nrow = length(canon), ncol = ncol(cts),
                    dimnames = list(canon, colnames(cts)))
for (m in canon) collapsed[m, ] <- cts[pairs[[m]][1], ] + cts[pairs[[m]][2], ]
collapsed <- as(collapsed, "CsparseMatrix")

mdata_CMV[["FB"]] <- CreateAssayObject(counts = collapsed)
DefaultAssay(mdata_CMV) <- "FB"
mdata_CMV <- NormalizeData(mdata_CMV, assay = "FB",
                           normalization.method = "CLR", margin = 2)
mdata_CMV <- ScaleData(mdata_CMV, assay = "FB", verbose = FALSE)

# ---- D0_data_CMV: rename only ----
rename_map <- c(
  "anti-human CD183 (CXCR3)"  = "CXCR3",
  "anti-human CD226 (DNAM-1)" = "CD226",
  "anti-human CD279 (PD-1)"   = "PD-1",
  "anti-human CD366 (Tim-3)"  = "Tim-3",
  "anti-human CD45RA"         = "CD45RA",
  "anti-human CD57"           = "CD57",
  "anti-human CD95 (Fas)"     = "CD95",
  "anti-human TIGIT (VSTM3)"  = "TIGIT"
)

cts_d0 <- GetAssayData(D0_data_CMV, assay = "FB", layer = "counts")
rn <- rownames(cts_d0); hit <- rn %in% names(rename_map)
rn[hit] <- rename_map[rn[hit]]; rownames(cts_d0) <- rn

D0_data_CMV[["FB"]] <- CreateAssayObject(counts = cts_d0)
DefaultAssay(D0_data_CMV) <- "FB"
D0_data_CMV <- NormalizeData(D0_data_CMV, assay = "FB",
                             normalization.method = "CLR", margin = 2)
D0_data_CMV <- ScaleData(D0_data_CMV, assay = "FB", verbose = FALSE)

# ---- save ----
dir.create("Donor_0", showWarnings = FALSE, recursive = TRUE)
saveRDS(mdata_CMV,   file = "final_mdata_CMV.rds")
saveRDS(D0_data_CMV, file = "Donor_0/D0_data_CMV.rds")


library(ggplot2)

# Assumption: D0_data_CMV is donor 0; mdata_CMV is the donors 1–5 merge.
# Adjust out_dir/prefix below if that mapping is wrong.
export_fb_umaps <- function(obj, out_dir, prefix, reduction = "umap") {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  DefaultAssay(obj) <- "FB"                      # ensures FB is plotted, not same-named RNA genes
  feats <- rownames(obj[["FB"]])
  for (f in feats) {
    fname <- gsub("[^A-Za-z0-9]+", "_", f)        # filename-safe (feature name itself unchanged)
    ggsave(
      filename = file.path(out_dir, sprintf("%s_UMAP_%s.png", prefix, fname)),
      plot = FeaturePlot(obj, features = f, reduction = reduction,
                         cols = c("lightgrey", "darkred")),
      width = 6, height = 5, dpi = 300
    )
  }
}

export_fb_umaps(D0_data_CMV, "corrected_FB_exports/D0",   "D0_CMV")
export_fb_umaps(mdata_CMV,   "corrected_FB_exports/D1-5", "D1-5_CMV")


library(ggplot2)

# ---- D0: all FB markers x 8 seurat_clusters ----
d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

p_d0 <- DotPlot(D0_data_CMV, features = d0_markers,
                group.by = "seurat_clusters", assay = "FB")
plot_data <- p_d0$data

# clusters numeric-sorted on x (robust to 0–7 vs 1–8); markers on y
clus_order <- as.character(sort(as.numeric(as.character(unique(plot_data$id)))))
plot_data$id <- factor(plot_data$id, levels = clus_order)
plot_data$features.plot <- factor(plot_data$features.plot, levels = rev(d0_markers))

p <- ggplot(plot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title  = element_text(size = 12)) +
  labs(x = "Seurat cluster", y = "CITEseq marker")

dir.create("corrected_FB_exports/D0", showWarnings = FALSE, recursive = TRUE)
ggsave("corrected_FB_exports/D0/D0_CMV_FB_bubbleplot_by_seurat_cluster_short.png",
       p, width = 7, height = 4, dpi = 300)

# 
# ---- D0: FB markers x specificity (pp65 vs IE1) ----

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

spec_levels <- c("pp65", "IE1")

# keep only cells with a defined specificity
keep_cells <- colnames(D0_data_CMV)[
  !is.na(D0_data_CMV$specificity) & D0_data_CMV$specificity %in% spec_levels
]
d0_spec <- subset(D0_data_CMV, cells = keep_cells)
d0_spec$specificity <- factor(as.character(d0_spec$specificity), levels = spec_levels)

message("Cells per specificity:")
print(table(d0_spec$specificity))

# scale = FALSE -> avg.exp.scaled is log1p(avg.exp), not a 2-group z-score
p_spec <- DotPlot(d0_spec, features = d0_markers,
                  group.by = "specificity", assay = "FB", scale = FALSE)
plot_data <- p_spec$data

# markers along x in the order given; specificity as 2 rows, pp65 on top
plot_data$features.plot <- factor(plot_data$features.plot, levels = d0_markers)
plot_data$id <- factor(as.character(plot_data$id), levels = rev(spec_levels))

# y-axis label colors, bottom -> top to match factor level order
spec_cols <- c(pp65 = "#1976D2", IE1 = "#8B0000")
y_lab_cols <- unname(spec_cols[levels(plot_data$id)])

p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 8), name = "% cells\nexpressing") +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red"),
                        name = "Mean CLR\nexpression\n(log1p)") +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y  = element_text(size = 11, face = "bold", color = y_lab_cols),
        axis.title   = element_text(size = 12),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()) +
  labs(x = "CITEseq marker", y = "Specificity",
       title = "D0 CMV-specific CD8 T cells: CITEseq markers by specificity")

dir.create("corrected_FB_exports/D0", showWarnings = FALSE, recursive = TRUE)
ggsave("corrected_FB_exports/D0/D0_CMV_FB_bubbleplot_by_specificity.png",
       p, width = 7.5, height = 3.2, dpi = 300)


# ---- D0: FB markers x major clones (IE1 1-7, pp65 1-20) ----

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

# numeric ordering, not lexicographic (clone_2 before clone_10)
ie1_clones  <- paste0("D0_IE1_clone_",  1:7)
pp65_clones <- paste0("D0_pp65_clone_", 1:20)
clone_levels_all <- c(ie1_clones, pp65_clones)

obs_clones <- unique(as.character(D0_data_CMV$Global_Clone_ID))
missing    <- setdiff(clone_levels_all, obs_clones)
if (length(missing)) message("Not found in object: ", paste(missing, collapse = ", "))
clone_levels <- clone_levels_all[clone_levels_all %in% obs_clones]

keep_cells <- colnames(D0_data_CMV)[
  !is.na(D0_data_CMV$Global_Clone_ID) &
    as.character(D0_data_CMV$Global_Clone_ID) %in% clone_levels
]
d0_clones <- subset(D0_data_CMV, cells = keep_cells)
d0_clones$Global_Clone_ID <- factor(as.character(d0_clones$Global_Clone_ID),
                                    levels = clone_levels)

clone_n <- table(d0_clones$Global_Clone_ID)
message("Cells per clone:"); print(clone_n)
message("Total cells retained: ", ncol(d0_clones))

p_clone <- DotPlot(d0_clones, features = d0_markers,
                   group.by = "Global_Clone_ID", assay = "FB", scale = TRUE)
plot_data <- p_clone$data

# markers along x; clones as rows, IE1 block on top, clone_1 at top of each block
plot_data$features.plot <- factor(plot_data$features.plot, levels = d0_markers)
plot_data$id <- factor(as.character(plot_data$id), levels = rev(clone_levels))

# y labels: color by specificity, append n
y_levels  <- levels(plot_data$id)
y_lab_col <- ifelse(grepl("_pp65_", y_levels), "#1976D2", "#8B0000")
y_lab_txt <- setNames(
  paste0(sub("^D0_", "", y_levels), "  (n=", as.integer(clone_n[y_levels]), ")"),
  y_levels
)

# divider between the pp65 block and the IE1 block
split_at <- sum(grepl("_pp65_", y_levels)) + 0.5

p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_hline(yintercept = split_at, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6), name = "% cells\nexpressing") +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red"),
                        name = "Scaled mean\nexpression") +
  scale_y_discrete(labels = y_lab_txt) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y  = element_text(size = 9, face = "bold", color = y_lab_col),
        axis.title   = element_text(size = 12),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()) +
  labs(x = "CITEseq marker", y = "Clone (Global_Clone_ID)",
       title = "D0 CMV-specific CD8 T cells: CITEseq markers by major clone")

dir.create("corrected_FB_exports/D0", showWarnings = FALSE, recursive = TRUE)
ggsave("corrected_FB_exports/D0/D0_CMV_FB_bubbleplot_by_clone.png",
       p, width = 8, height = 8.5, dpi = 300)



# ---- mdata_CMV: all FB markers x 5 metaclusters (flipped style) ----
mdata_markers <- intersect(
  c("CD45RA", "CCR7", "CD95", "CD57", "KLRG1", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(mdata_CMV[["FB"]])
)

p_md <- DotPlot(mdata_CMV, features = mdata_markers,
                group.by = "metacluster", assay = "FB")
plot_data <- p_md$data

# same metacluster order you used in the CMV bubble plots
metacluster_order <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")
plot_data$id <- factor(plot_data$id, levels = rev(metacluster_order))
plot_data$features.plot <- factor(plot_data$features.plot, levels = mdata_markers)

p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        axis.title  = element_blank())

dir.create("corrected_FB_exports/D1-5", showWarnings = FALSE, recursive = TRUE)
ggsave("corrected_FB_exports/D1-5/D1-5_CMV_FB_bubbleplot_by_metacluster.png",
       p, width = 7, height = 2.5, dpi = 300)




#
#
# Global-clone-group phenotypic coherence test (Dan's 6 curated pp65 groups)
#   Mirrors the et64 specificity-stratified test, but the "cluster" is one of six
#   hand-curated groups of 2-3 pp65 Global clones. For each cell in a group, we ask
#   how many of its k neighbors belong to the OTHER clone(s) in its group.
#   - candidates restricted to pp65 cells (all groups are pp65)
#   - same-cloneId neighbors excluded
#   - clone-weighted; null shuffles group labels across ALL pp65 cells (matches et64)
#   - each group tested self-contained (a group's target = its other member clones)
#
#
library(Seurat); library(dplyr); library(tidyr); library(FNN); library(readr)

out_dir <- "exports/Dan's_clone_groups_kNN_bias_analysis"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
obj <- mdata_CMV
set.seed(1)
K         <- 100
N_PERM    <- 10000
N_PCS     <- 30
OVERFETCH <- 1500

# ── 0. The six curated groups (edit here if labels change) ────────────────────
clone_groups <- list(
  F1 = c("D2_pp65-clone-5", "D3_pp65-clone-5", "D4_pp65-clone-2"),
  F2 = c("D1_pp65-clone-6", "D3_pp65-clone-1", "D3_pp65-clone-2"),
  F3 = c("D3_pp65-clone-3", "D4_pp65-clone-8"),
  F4 = c("D3_pp65-clone-8", "D4_pp65-clone-5"),
  F5 = c("D3_pp65-clone-6", "D3_pp65-clone-7"),
  F6 = c("D1_pp65-clone-1", "D2_pp65-clone-2")
)

# ── 1. Per-cell table over ALL pp65 cells (the candidate universe) ────────────
# All group clones are pp65, so the neighbor pool is pp65 cells (mirrors the
# et64 test, where candidates were restricted to the cluster's specificity).
pcs <- Embeddings(obj, "pca")[, 1:N_PCS, drop = FALSE]
md  <- obj@meta.data
cells <- tibble(
  barcode   = colnames(obj),
  cloneId   = as.character(md$cloneId),
  donor     = as.character(md$donor_identity),
  spec_raw  = as.character(md$specificity),
  global_id = as.character(md$Global_Clone_ID))
keep <- cells$spec_raw == "pp65" & !is.na(cells$spec_raw) &
        !is.na(cells$cloneId) & !grepl("notCalled", cells$cloneId)
cells <- cells[keep, ]; pcs <- pcs[keep, , drop = FALSE]
N <- nrow(cells)
message(sprintf("pp65 candidate cells: %d | clones: %d | donors: %d",
                N, dplyr::n_distinct(cells$cloneId), dplyr::n_distinct(cells$donor)))

# map each cell to which group(s) it belongs to via its Global_Clone_ID
clone_to_group <- stack(clone_groups) %>%           # values, ind(=group)
  dplyr::transmute(global_id = as.character(values), group = as.character(ind))
cells <- cells %>% dplyr::left_join(clone_to_group, by = "global_id")
grp_counts <- cells %>% dplyr::filter(!is.na(group)) %>%
  dplyr::count(group, global_id, name = "n_cells")
message("Cells per group clone:")
print(grp_counts, n = Inf)

# ── 2. kNN among ALL pp65 cells (single specificity, so one pool) ─────────────
k_s <- min(OVERFETCH, N - 1)
nn  <- get.knn(pcs, k = k_s)$nn.index          # row indices into `cells`, distance-sorted

# ── 3. Per cell: K nearest surviving neighbors (drop same-clone) ─────────────
surv_idx  <- vector("list", N)
n_dropped <- 0L
for (r in seq_len(N)) {
  g <- nn[r, ]
  g <- g[cells$cloneId[g] != cells$cloneId[r]]   # exclude same clone
  if (length(g) < K) n_dropped <- n_dropped + 1L
  surv_idx[[r]] <- head(g, K)
}
message(sprintf("Cells with < K=%d same-spec non-same-clone neighbors: %d (%.1f%%)",
                K, n_dropped, 100 * n_dropped / N))

# ── 4. Per-group self-contained statistic ────────────────────────────────────
# For a given group G and a group-label vector `lab` (each cell's group or NA):
#   focal cells = cells whose lab == G; for each, count k-neighbors also labeled G.
# Because same-clone neighbors are already excluded, a cell in a 1-clone-represented
# group scores 0 by construction (no other clone present) — reported, not hidden.
count_same_group <- function(lab, G) {
  focal <- which(lab == G & !is.na(lab))
  sc <- rep(NA_real_, N)
  for (r in focal) {
    g <- surv_idx[[r]]
    sc[r] <- if (length(g) == 0) NA_real_ else sum(lab[g] == G & !is.na(lab[g]))
  }
  sc
}

# clone-weighted group score: mean within clone, then mean across the group's clones
agg_group_score <- function(sc, lab, G) {
  idx <- which(lab == G & !is.na(lab))
  if (length(idx) == 0) return(NA_real_)
  tibble(cloneId = cells$cloneId[idx], v = sc[idx]) %>%
    dplyr::group_by(cloneId) %>% dplyr::summarise(cm = mean(v, na.rm = TRUE), .groups = "drop") %>%
    dplyr::summarise(s = mean(cm, na.rm = TRUE)) %>% dplyr::pull(s)
}

# ── 5. Observed group scores ─────────────────────────────────────────────────
obs_lab <- cells$group
groups  <- names(clone_groups)
obs_score <- sapply(groups, function(G) agg_group_score(count_same_group(obs_lab, G), obs_lab, G))

# ── 6. Null: shuffle group labels across ALL pp65 cells (matches et64) ────────
# Preserve the number of cells carrying each group label; assign them to random
# pp65 cells. This asks: are the curated groups more coherent than an equally-sized
# random grouping of pp65 cells?
lab_nonNA <- which(!is.na(obs_lab))
n_labeled <- length(lab_nonNA)
lab_values <- obs_lab[lab_nonNA]     # the multiset of group labels to redistribute

null_mat <- matrix(NA_real_, nrow = N_PERM, ncol = length(groups),
                   dimnames = list(NULL, groups))
for (p in seq_len(N_PERM)) {
  pl <- rep(NA_character_, N)
  pl[sample.int(N, n_labeled)] <- sample(lab_values)   # scatter labels over all pp65 cells
  for (j in seq_along(groups)) {
    G <- groups[j]
    null_mat[p, j] <- agg_group_score(count_same_group(pl, G), pl, G)
  }
  if (p %% 1000 == 0) message(sprintf("    perm %d / %d", p, N_PERM))
}

# ── 7. Assemble results ──────────────────────────────────────────────────────
pvals <- sapply(seq_along(groups), function(j)
  (1 + sum(null_mat[, j] >= obs_score[j], na.rm = TRUE)) / (N_PERM + 1))

group_info <- cells %>% dplyr::filter(!is.na(group)) %>%
  dplyr::group_by(group) %>%
  dplyr::summarise(clones   = paste(sort(unique(global_id)), collapse = ","),
                   n_clones = dplyr::n_distinct(global_id),
                   n_cells  = dplyr::n(),
                   n_donors = dplyr::n_distinct(donor),
                   donors   = paste(sort(unique(donor)), collapse = ","),
                   .groups = "drop")

final <- tibble(group = groups,
                same_group_obs       = obs_score,
                same_group_null_mean = colMeans(null_mat, na.rm = TRUE),
                same_group_effect    = obs_score - colMeans(null_mat, na.rm = TRUE),
                same_group_p_emp     = pvals,
                same_group_p_adj     = p.adjust(pvals, "BH")) %>%
  dplyr::left_join(group_info, by = "group") %>%
  dplyr::select(group, clones, n_clones, n_cells, n_donors, donors,
                same_group_obs, same_group_null_mean, same_group_effect,
                same_group_p_emp, same_group_p_adj) %>%
  dplyr::arrange(dplyr::desc(same_group_effect))

# ── 8. Write CSV with explanatory header ─────────────────────────────────────
csv_path <- file.path(out_dir, "clone_group_kNN_bias_results.csv")
header <- c(
"# Global-clone-group phenotypic coherence test. One row per curated group (F1-F6).",
"# Each group is 2-3 pp65-specific Global clones (Dan's TCR families). Neighbor space:",
"#   30-PC integrated PCA. k = 100 nearest neighbors per cell, after excluding same-cloneId",
"#   neighbors and restricting the candidate pool to pp65 cells (all group clones are pp65).",
"# Each group is tested self-contained: for a cell in group G, the statistic counts how many",
"#   of its k neighbors belong to G (i.e. to the OTHER member clone(s) of G, since same-clone",
"#   neighbors are excluded). Scores are clone-weighted: averaged within each clone, then",
"#   across the group's clones, so an abundant clone cannot dominate the group score.",
"#",
"# COLUMNS (prefix same_group_):",
"#   obs       = observed clone-weighted mean same-group neighbor count for the group.",
"#   null_mean = mean of that score across 10000 permutations where group labels were randomly",
"#               scattered across ALL pp65 cells (preserving how many cells carry each label);",
"#               the score expected if these clones grouped together no more than chance.",
"#   effect    = obs - null_mean. Extra same-group neighbors per cell beyond chance. Note the",
"#               target is tiny (1-2 other clones), so effects are small in absolute count even",
"#               when highly significant; read significance from p_adj and the obs/null ratio.",
"#   p_emp     = empirical p-value: fraction of permutations with null >= observed (+1/+1).",
"#   p_adj     = p_emp after Benjamini-Hochberg FDR correction across the six groups.",
"#",
"# clones/n_clones/n_cells/n_donors/donors describe the group's composition in the object.")
con <- file(csv_path, "w"); writeLines(header, con); close(con)
write_csv(final, csv_path, append = TRUE, col_names = TRUE)

cat(sprintf("\nWrote %s\n", csv_path))
cat(sprintf("Significant groups (p_adj < 0.05): %d / %d\n",
            sum(final$same_group_p_adj < 0.05, na.rm = TRUE), nrow(final)))
print(final %>% dplyr::select(group, clones, n_cells, n_donors,
        same_group_obs, same_group_null_mean, same_group_effect, same_group_p_adj), n = Inf)




#
#
# Per-cluster TCR summary tables (one PNG per TCRdist cluster)
#   Title row = cluster name. Then, in this order:
#     - one row per Global clone in the cluster (grouped by Global_Clone_ID)
#     - one row per unique (TRB, TRA-set, donor) among MINOR (non-Global) cells
#   Columns: Donor | Clone | Number of Cells | TRB cdr3 | TRA cdr3(s)
#   Sequences parsed from the object's t_cdr3s_aa (per-cell). notCalled excluded.
#   Output: exports/TCR_clusters/cluster_tables/
#
#
library(dplyr); library(stringr); library(tidyr); library(gt); library(purrr)

out_dir <- "exports/TCR_clusters/cluster_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 0. Helpers ───────────────────────────────────────────────────────────────
# extract ALL TRA / TRB CDR3s from a t_cdr3s_aa string (handles dual-alpha),
# returned as a single "; "-joined, de-duplicated, sorted string ("" if none).
extract_chain <- function(s, prefix) {
  if (is.na(s)) return(NA_character_)
  parts <- str_split(s, ";")[[1]]
  hit   <- parts[str_detect(parts, paste0("^", prefix, ":"))]
  if (length(hit) == 0) return(NA_character_)
  seqs  <- str_replace(hit, paste0("^", prefix, ":"), "")
  seqs  <- seqs[nchar(seqs) > 0]
  if (length(seqs) == 0) return(NA_character_)
  paste(sort(unique(seqs)), collapse = "; ")
}

# format a Global_Clone_ID like "D1_pp65-clone-3" -> "pp65 - clone 3"
fmt_clone <- function(gid) {
  spec <- str_match(gid, "_(pp65|IE1)-")[, 2]
  num  <- str_match(gid, "clone-(\\d+)")[, 2]
  ifelse(is.na(spec) | is.na(num), gid, sprintf("%s - clone %s", spec, num))
}

# short donor label "Donor 2" / "D2" -> "D2"
short_donor <- function(d) paste0("D", str_extract(as.character(d), "\\d"))

# ── 1. Per-cell working table (clustered, non-notCalled) ─────────────────────
md <- mdata_CMV@meta.data
cell_tbl <- tibble(
  cluster   = as.character(md$TCR_cluster_et64),
  cloneId   = as.character(md$cloneId),
  donor     = short_donor(md$donor_identity),
  global_id = as.character(md$Global_Clone_ID),
  cdr3str   = as.character(md$t_cdr3s_aa)) %>%
  dplyr::filter(!is.na(cluster), !grepl("notCalled", cloneId)) %>%
  dplyr::mutate(trb = vapply(cdr3str, extract_chain, character(1), prefix = "TRB"),
                tra = vapply(cdr3str, extract_chain, character(1), prefix = "TRA"),
                is_global = !is.na(global_id))

# ── 2. GLOBAL-clone rows: one per Global_Clone_ID ────────────────────────────
# dominant TRB / TRA-set across the clone's cells; total cells; its donor.
mode_or_join <- function(x) {                    # most common non-NA value
  x <- x[!is.na(x)]; if (!length(x)) return(NA_character_)
  names(which.max(table(x)))
}
global_rows <- cell_tbl %>%
  dplyr::filter(is_global) %>%
  dplyr::group_by(cluster, global_id) %>%
  dplyr::summarise(
    donor    = mode_or_join(donor),
    clone    = fmt_clone(dplyr::first(global_id)),
    n_cells  = dplyr::n(),
    trb      = mode_or_join(trb),
    tra      = mode_or_join(tra),
    .groups  = "drop") %>%
  dplyr::mutate(kind = "global")

# ── 3. MINOR rows: one per unique (TRB, TRA, donor) among non-Global cells ───
minor_rows <- cell_tbl %>%
  dplyr::filter(!is_global) %>%
  dplyr::group_by(cluster, donor, trb, tra) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::mutate(clone = "minor", global_id = NA_character_, kind = "minor")

# ── 4. Combine; order Global first (by clone number) then minor (by cells) ────
tbl_all <- dplyr::bind_rows(global_rows, minor_rows) %>%
  dplyr::mutate(
    trb = tidyr::replace_na(trb, ""),
    tra = tidyr::replace_na(tra, ""),
    clone_num = suppressWarnings(as.integer(str_match(global_id, "clone-(\\d+)")[, 2]))) %>%
  dplyr::arrange(cluster,
                 dplyr::desc(kind == "global"),      # global block first
                 clone_num,                          # global sorted by clone number
                 dplyr::desc(n_cells)) %>%           # minor sorted by cell count
  dplyr::select(cluster, kind, Donor = donor, Clone = clone,
                `Number of Cells` = n_cells, `TRB cdr3` = trb, `TRA cdr3(s)` = tra)

# ── 5. Render one PNG table per cluster (browserless: gridExtra/gtable) ───────
library(gridExtra); library(grid); library(gtable)

clusters <- tbl_all %>% dplyr::distinct(cluster) %>% dplyr::pull(cluster)
clusters <- clusters[order(suppressWarnings(as.integer(clusters)), clusters)]

render_cluster_table <- function(cl) {
  d <- tbl_all %>% dplyr::filter(cluster == cl)
  n_global <- sum(d$kind == "global"); n_minor <- sum(d$kind == "minor")
  is_glob  <- d$kind == "global"
  d_show   <- d %>% dplyr::select(Donor, Clone, `Number of Cells`,
                                  `TRB cdr3`, `TRA cdr3(s)`)

  # base theme: left-aligned, small; CDR3 cols will be overridden to monospace
  th <- ttheme_minimal(
    base_size = 9,
    core = list(
      fg_params = list(hjust = 0, x = 0.02, fontsize = 9),
      bg_params = list(fill = ifelse(is_glob, "#EAF2FB", "white"),
                       col = "grey90", lwd = 0.5)),
    colhead = list(
      fg_params = list(fontface = "bold", hjust = 0, x = 0.02, fontsize = 9),
      bg_params = list(fill = "grey95", col = "grey80", lwd = 0.5)))

  tg <- tableGrob(d_show, rows = NULL, theme = th)

  # monospace font for the two CDR3 columns (cols 4 and 5), body cells only
  mono_cols <- which(colnames(d_show) %in% c("TRB cdr3", "TRA cdr3(s)"))
  for (cc in mono_cols) {
    cells <- which(grepl(sprintf("core-fg", cc), tg$layout$name) &
                   tg$layout$l == cc)
    for (k in cells) tg$grobs[[k]]$gp <- gpar(fontfamily = "mono", fontsize = 8.5)
  }
  # bold the Global clone rows' first two columns for emphasis
  for (rr in which(is_glob)) {
    for (cc in 1:2) {
      k <- which(tg$layout$name == "core-fg" & tg$layout$t == rr + 1 & tg$layout$l == cc)
      if (length(k)) tg$grobs[[k]]$gp <- gpar(fontface = "bold", fontsize = 9)
    }
  }

  # title + footnote as stacked grobs
  title <- textGrob(sprintf("Cluster %s", cl),
                    gp = gpar(fontface = "bold", fontsize = 14), hjust = 0, x = 0.02)
  foot  <- textGrob(sprintf(
      "%d Global clone%s, %d unique minor TRB+TRA combination%s.  Sequences from t_cdr3s_aa; blank TRA = no alpha recovered.",
      n_global, ifelse(n_global == 1, "", "s"),
      n_minor,  ifelse(n_minor == 1, "", "s")),
      gp = gpar(fontsize = 7.5, col = "grey40", fontface = "italic"), hjust = 0, x = 0.02)

  body <- gtable_add_rows(tg, unit(1.2, "lines"), pos = 0)             # title space
  body <- gtable_add_grob(body, title, t = 1, l = 1, r = ncol(tg))
  body <- gtable_add_rows(body, unit(1.4, "lines"), pos = -1)          # footnote space
  body <- gtable_add_grob(body, foot, t = nrow(body), l = 1, r = ncol(tg))

  # size the canvas to the content
  w_in <- 11
  h_in <- 0.9 + 0.24 * (n_global + n_minor + 1)   # header+title+footnote overhead + rows
  fn <- file.path(out_dir, sprintf("cluster_%s_table.png", cl))
  ggplot2::ggsave(fn, body, width = w_in, height = h_in, dpi = 200,
                  limitsize = FALSE, bg = "white")
  message(sprintf("cluster %s: %d global + %d minor rows -> %s", cl, n_global, n_minor, fn))
}

invisible(lapply(clusters, render_cluster_table))
message(sprintf("\nDone: %d cluster tables in %s", length(clusters), out_dir))


#
#
# Per-cluster TCR summary tables — FILTERED version
#   Changes vs. first version:
#     - Donor shown as "Donor 3" (full word)
#     - Number of Cells centered
#     - Global clones show ALL unique TRB / TRA cdr3s across the clone (not dominant),
#       so a clone with >1 TRB or >2 TRAs is made visible
#     - Minor rows kept only if >= 5 cells
#     - Filename prefixed with dominant specificity (pp65_ / IE1_)
#   Output: exports/TCR_clusters/cluster_tables/cluster_tables_filtered/
#
#
library(dplyr); library(stringr); library(tidyr); library(purrr)
library(gridExtra); library(grid); library(gtable)

out_dir <- "exports/TCR_clusters/cluster_tables/cluster_tables_filtered"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
MIN_MINOR_CELLS <- 5

# ── 0. Helpers ───────────────────────────────────────────────────────────────
# ALL distinct non-NA CDR3s of one chain across a set of t_cdr3s_aa strings,
# "; "-joined (handles dual-alpha; "" if none present).
extract_all_chain <- function(strings, prefix) {
  seqs <- unlist(lapply(strings, function(s) {
    if (is.na(s)) return(character(0))
    parts <- str_split(s, ";")[[1]]
    hit   <- parts[str_detect(parts, paste0("^", prefix, ":"))]
    str_replace(hit, paste0("^", prefix, ":"), "")
  }))
  seqs <- seqs[nchar(seqs) > 0]
  if (length(seqs) == 0) return("")
  paste(sort(unique(seqs)), collapse = "; ")
}
# per-cell single-chain extractor (for minor-row keying)
extract_cell_chain <- function(s, prefix) {
  if (is.na(s)) return(NA_character_)
  parts <- str_split(s, ";")[[1]]
  hit   <- parts[str_detect(parts, paste0("^", prefix, ":"))]
  seqs  <- str_replace(hit, paste0("^", prefix, ":"), "")
  seqs  <- seqs[nchar(seqs) > 0]
  if (length(seqs) == 0) return(NA_character_)
  paste(sort(unique(seqs)), collapse = "; ")   # keep dual-alpha together per cell
}
fmt_clone <- function(gid) {
  spec <- str_match(gid, "_(pp65|IE1)-")[, 2]
  num  <- str_match(gid, "clone-(\\d+)")[, 2]
  ifelse(is.na(spec) | is.na(num), gid, sprintf("%s - clone %s", spec, num))
}
full_donor <- function(d) paste("Donor", str_extract(as.character(d), "\\d"))  # "Donor 3"

# ── 1. Per-cell working table (clustered, non-notCalled) ─────────────────────
md <- mdata_CMV@meta.data
cell_tbl <- tibble(
  cluster   = as.character(md$TCR_cluster_et64),
  cloneId   = as.character(md$cloneId),
  donor     = full_donor(md$donor_identity),
  spec      = as.character(md$specificity),
  global_id = as.character(md$Global_Clone_ID),
  cdr3str   = as.character(md$t_cdr3s_aa)) %>%
  dplyr::filter(!is.na(cluster), !grepl("notCalled", cloneId)) %>%
  dplyr::mutate(is_global = !is.na(global_id))

# ── 2. GLOBAL-clone rows: ALL unique TRB / TRA across the clone ──────────────
global_rows <- cell_tbl %>%
  dplyr::filter(is_global) %>%
  dplyr::group_by(cluster, global_id) %>%
  dplyr::summarise(
    donor   = full_donor(str_extract(dplyr::first(donor), "\\d")),  # clone's donor
    clone   = fmt_clone(dplyr::first(global_id)),
    n_cells = dplyr::n(),
    trb     = extract_all_chain(cdr3str, "TRB"),
    tra     = extract_all_chain(cdr3str, "TRA"),
    .groups = "drop") %>%
  dplyr::mutate(kind = "global")

# ── 3. MINOR rows: unique (TRB, TRA, donor); keep only >= MIN_MINOR_CELLS ────
minor_rows <- cell_tbl %>%
  dplyr::filter(!is_global) %>%
  dplyr::mutate(trb = vapply(cdr3str, extract_cell_chain, character(1), prefix = "TRB"),
                tra = vapply(cdr3str, extract_cell_chain, character(1), prefix = "TRA")) %>%
  dplyr::group_by(cluster, donor, trb, tra) %>%
  dplyr::summarise(n_cells = dplyr::n(), .groups = "drop") %>%
  dplyr::filter(n_cells >= MIN_MINOR_CELLS) %>%
  dplyr::mutate(clone = "minor", global_id = NA_character_, kind = "minor")

# ── 4. Combine; Global first (by clone number), then minor (by cells) ────────
tbl_all <- dplyr::bind_rows(global_rows, minor_rows) %>%
  dplyr::mutate(
    trb = tidyr::replace_na(trb, ""),
    tra = tidyr::replace_na(tra, ""),
    clone_num = suppressWarnings(as.integer(str_match(global_id, "clone-(\\d+)")[, 2]))) %>%
  dplyr::arrange(cluster, dplyr::desc(kind == "global"), clone_num, dplyr::desc(n_cells)) %>%
  dplyr::select(cluster, kind, Donor = donor, Clone = clone,
                `Number of Cells` = n_cells, `TRB cdr3` = trb, `TRA cdr3(s)` = tra)

# ── 4b. Dominant specificity per cluster (cell-weighted) for filename prefix ──
cluster_spec <- cell_tbl %>%
  dplyr::filter(!is.na(spec)) %>%
  dplyr::count(cluster, spec, name = "n") %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(dom_spec = spec[which.max(n)], .groups = "drop")

# ── 5. Render one PNG per cluster (browserless) ──────────────────────────────
clusters <- tbl_all %>% dplyr::distinct(cluster) %>% dplyr::pull(cluster)
clusters <- clusters[order(suppressWarnings(as.integer(clusters)), clusters)]
skipped <- character(0)

render_cluster_table <- function(cl) {
  d <- tbl_all %>% dplyr::filter(cluster == cl)
  if (nrow(d) == 0) { skipped <<- c(skipped, cl); return(invisible()) }  # nothing survives filter
  n_global <- sum(d$kind == "global"); n_minor <- sum(d$kind == "minor")
  is_glob  <- d$kind == "global"
  d_show   <- d %>% dplyr::select(Donor, Clone, `Number of Cells`, `TRB cdr3`, `TRA cdr3(s)`)

  # per-column horizontal justification: center "Number of Cells" (col 3), rest left
  ncol_show <- ncol(d_show)
  hjust_vec <- c(0, 0, 0.5, 0, 0); x_vec <- c(0.02, 0.02, 0.5, 0.02, 0.02)

  th <- ttheme_minimal(
    base_size = 9,
    core = list(
      fg_params = list(hjust = rep(hjust_vec, each = nrow(d_show)),
                       x     = rep(x_vec,     each = nrow(d_show)), fontsize = 9),
      bg_params = list(fill = ifelse(is_glob, "#EAF2FB", "white"),
                       col = "grey90", lwd = 0.5)),
    colhead = list(
      fg_params = list(fontface = "bold",
                       hjust = rep(hjust_vec, each = 1), x = rep(x_vec, each = 1),
                       fontsize = 9),
      bg_params = list(fill = "grey95", col = "grey80", lwd = 0.5)))

  tg <- tableGrob(d_show, rows = NULL, theme = th)

  # monospace CDR3 columns (4,5), body only
  for (cc in c(4, 5)) {
    k <- which(tg$layout$name == "core-fg" & tg$layout$l == cc)
    for (kk in k) tg$grobs[[kk]]$gp <- gpar(fontfamily = "mono", fontsize = 8.5)
  }
  # bold Global rows' first two columns
  for (rr in which(is_glob)) for (cc in 1:2) {
    k <- which(tg$layout$name == "core-fg" & tg$layout$t == rr + 1 & tg$layout$l == cc)
    if (length(k)) tg$grobs[[k]]$gp <- gpar(fontface = "bold", fontsize = 9)
  }

  title <- textGrob(sprintf("Cluster %s", cl),
                    gp = gpar(fontface = "bold", fontsize = 14), hjust = 0, x = 0.02)
  foot  <- textGrob(sprintf(
      "%d Global clone%s, %d minor TRB+TRA combination%s (>= %d cells).  Sequences from t_cdr3s_aa; blank TRA = no alpha recovered.",
      n_global, ifelse(n_global == 1, "", "s"),
      n_minor,  ifelse(n_minor == 1, "", "s"), MIN_MINOR_CELLS),
      gp = gpar(fontsize = 7.5, col = "grey40", fontface = "italic"), hjust = 0, x = 0.02)

  body <- gtable_add_rows(tg, unit(1.2, "lines"), pos = 0)
  body <- gtable_add_grob(body, title, t = 1, l = 1, r = ncol(tg))
  body <- gtable_add_rows(body, unit(1.4, "lines"), pos = -1)
  body <- gtable_add_grob(body, foot, t = nrow(body), l = 1, r = ncol(tg))

  dom <- cluster_spec$dom_spec[match(cl, cluster_spec$cluster)]
  dom <- ifelse(is.na(dom), "NA", dom)
  w_in <- 11; h_in <- 0.9 + 0.24 * (n_global + n_minor + 1)
  fn <- file.path(out_dir, sprintf("%s_cluster_%s_table.png", dom, cl))
  ggplot2::ggsave(fn, body, width = w_in, height = h_in, dpi = 200, limitsize = FALSE, bg = "white")
  message(sprintf("[%s] cluster %s: %d global + %d minor -> %s", dom, cl, n_global, n_minor, basename(fn)))
}

invisible(lapply(clusters, render_cluster_table))
message(sprintf("\nDone: %d tables in %s", length(clusters) - length(skipped), out_dir))
if (length(skipped))
  message(sprintf("Skipped (no Global clones and no minor rows >= %d cells): %s",
                  MIN_MINOR_CELLS, paste(sort(skipped), collapse = ", ")))



#
#
# Cluster cards v2 — filtered table (top) + contour UMAP (left) + metacluster bars (right)
#   Table logic MATCHES the standalone filtered tables:
#     - Global rows: 1 TRB (most cells) + top-2 TRA
#     - Minor rows: cloneId-collapsed, shared-TRB grouping w/ fragment absorption,
#       merged across donors; compact donor label "Donor 2, 3, 4"
#     - Row filter: keep rows (major AND minor) with >= MIN_CELLS cells
#   Cards only for clusters with >= 2 unique TCR sequences (constellation min2seq set).
#   UMAP: per-clone KDE contours (bins=10) over grey full-object backdrop, square panel.
#   Bars: per-clone metacluster proportion (0-100%). Labels always include donor.
#   Metacluster fill uses the project-standard palette (matches the metacluster UMAP
#   and the corrected_colors stacked bars), not a card-local one.
#   Output: exports/TCR_clusters/cluster_cards_v2/
#
#
library(Seurat); library(dplyr); library(stringr); library(tidyr); library(purrr)
library(ggplot2); library(gridExtra); library(grid); library(gtable); library(cowplot)

out_dir <- "exports/TCR_clusters/cluster_cards_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
MIN_CELLS   <- 5
META_LEVELS <- c("Tcm","Tem","T-CMC1","NK-like","ISG+")

# Canonical metacluster palette — identical to scales::hue_pal()(5) applied to
# META_LEVELS, i.e. the colors DimPlot(group.by = "metacluster") assigns.
meta_palette <- c(
  "Tcm"     = "#F8766D",
  "Tem"     = "#A3A500",
  "T-CMC1"  = "#00BF7D",
  "NK-like" = "#00B0F6",
  "ISG+"    = "#E76BF3"
)

conga_dir   <- "exports/conga/mdata_CMV"

# ── Helpers ──────────────────────────────────────────────────────────────────
pull_chain_reads <- function(strings, prefix) {
  unlist(lapply(strings, function(s) {
    if (is.na(s)) return(character(0))
    parts <- str_split(s, ";")[[1]]; hit <- parts[str_detect(parts, paste0("^", prefix, ":"))]
    seqs  <- str_replace(hit, paste0("^", prefix, ":"), ""); seqs[nchar(seqs) > 0] }))
}
top_trb <- function(strings) { r <- pull_chain_reads(strings, "TRB"); if (!length(r)) return("")
  names(which.max(table(r))) }
top_tra <- function(strings, n = 2) { r <- pull_chain_reads(strings, "TRA"); if (!length(r)) return("")
  tb <- sort(table(r), decreasing = TRUE); paste(names(tb)[seq_len(min(n, length(tb)))], collapse = "; ") }
fmt_clone  <- function(gid) { spec <- str_match(gid, "_(pp65|IE1)-")[,2]; num <- str_match(gid, "clone-(\\d+)")[,2]
  ifelse(is.na(spec)|is.na(num), gid, sprintf("%s - clone %s", spec, num)) }
full_donor <- function(d) paste("Donor", str_extract(as.character(d), "\\d"))

clone_labels <- function(major_ids) {
  info <- tibble(global_id = major_ids,
                 base = fmt_clone(major_ids),
                 num  = suppressWarnings(as.integer(str_match(major_ids, "clone-(\\d+)")[,2])),
                 don  = str_match(major_ids, "^D(\\d)")[,2]) %>%
    dplyr::mutate(label = sprintf("%s (D%s)", base, don)) %>%
    dplyr::arrange(num, don)
  setNames(info$label, info$global_id)
}

# ── 0. Per-cell table + qualifying clusters ──────────────────────────────────
md <- mdata_CMV@meta.data
cell_tbl <- tibble(
  cluster=as.character(md$TCR_cluster_et64), cloneId=as.character(md$cloneId),
  donor=full_donor(md$donor_identity), spec=as.character(md$specificity),
  metacluster=as.character(md$metacluster),
  global_id=as.character(md$Global_Clone_ID), cdr3str=as.character(md$t_cdr3s_aa)) %>%
  dplyr::filter(!is.na(cluster), !grepl("notCalled", cloneId)) %>%
  dplyr::mutate(is_global = !is.na(global_id))

g_data <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
cluster_nseq <- g_data$nodes %>% dplyr::distinct(cluster, mrow) %>% dplyr::count(cluster, name="n_seq")
keep_clusters <- cluster_nseq %>% dplyr::filter(n_seq >= 2) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>% dplyr::pull(cluster)
message(sprintf("Clusters with >=2 unique sequences: %d", length(keep_clusters)))

cluster_spec <- cell_tbl %>% dplyr::filter(!is.na(spec)) %>%
  dplyr::count(cluster, spec, name="n") %>% dplyr::group_by(cluster) %>%
  dplyr::summarise(dom_spec = spec[which.max(n)], .groups="drop")

# ── 1. Table df for ONE cluster (matches standalone filtered logic) ──────────
build_cluster_table_df <- function(cl) {
  cc <- cell_tbl %>% dplyr::filter(cluster == cl)

  global_rows <- cc %>% dplyr::filter(is_global) %>%
    dplyr::group_by(global_id) %>%
    dplyr::summarise(donor=full_donor(str_extract(dplyr::first(donor),"\\d")),
                     clone=fmt_clone(dplyr::first(global_id)), n_cells=dplyr::n(),
                     trb=top_trb(cdr3str), tra=top_tra(cdr3str, 2), .groups="drop") %>%
    dplyr::mutate(kind="global",
                  clone_num=suppressWarnings(as.integer(str_match(global_id,"clone-(\\d+)")[,2])))

  mc <- cc %>% dplyr::filter(!is_global) %>%
    dplyr::group_by(cloneId) %>%
    dplyr::summarise(donor=full_donor(str_extract(dplyr::first(donor),"\\d")),
                     n_cells=dplyr::n(), trb=top_trb(cdr3str), tra=top_tra(cdr3str, 1),
                     .groups="drop") %>%
    dplyr::mutate(trb = dplyr::na_if(trb, ""), tra = dplyr::na_if(tra, ""))

  minor_rows <- tibble()
  if (nrow(mc) > 0) {
    has_trb <- mc %>% dplyr::filter(!is.na(trb)); no_trb <- mc %>% dplyr::filter(is.na(trb))
    trb_rows <- list()
    for (b in unique(has_trb$trb)) {
      grp <- has_trb %>% dplyr::filter(trb == b)
      with_a <- grp %>% dplyr::filter(!is.na(tra)); only_b <- grp %>% dplyr::filter(is.na(tra))
      if (nrow(with_a) == 0) {
        trb_rows[[length(trb_rows)+1]] <- only_b %>% dplyr::mutate(row_trb=b, row_tra=NA_character_)
      } else {
        tra_sizes <- with_a %>% dplyr::group_by(tra) %>%
          dplyr::summarise(cells=sum(n_cells), .groups="drop") %>% dplyr::arrange(dplyr::desc(cells))
        dom_a <- tra_sizes$tra[1]
        trb_rows[[length(trb_rows)+1]] <- dplyr::bind_rows(
          with_a %>% dplyr::mutate(row_trb=b, row_tra=tra),
          only_b %>% dplyr::mutate(row_trb=b, row_tra=dom_a))
      }
    }
    trb_assigned <- if (length(trb_rows) == 0)
      tibble(cloneId=character(), donor=character(), n_cells=integer(),
             trb=character(), tra=character(), row_trb=character(), row_tra=character())
      else dplyr::bind_rows(trb_rows)
    existing_tra <- unique(stats::na.omit(trb_assigned$row_tra))
    ta_match <- no_trb %>% dplyr::filter(!is.na(tra) & tra %in% existing_tra) %>%
      dplyr::mutate(row_trb = vapply(tra, function(a) {
          cand <- trb_assigned %>% dplyr::filter(row_tra == a) %>%
            dplyr::group_by(row_trb) %>% dplyr::summarise(c=sum(n_cells), .groups="drop") %>%
            dplyr::arrange(dplyr::desc(c))
          if (nrow(cand) == 0) NA_character_ else as.character(cand$row_trb[1]) }, character(1)),
        row_tra = tra)
    ta_own <- no_trb %>% dplyr::filter(is.na(tra) | !(tra %in% existing_tra)) %>%
      dplyr::mutate(row_trb=NA_character_, row_tra=tra)
    assigned <- dplyr::bind_rows(trb_assigned, ta_match, ta_own)
    minor_rows <- assigned %>%
      dplyr::group_by(row_trb, row_tra) %>%
      dplyr::summarise(n_cells=sum(n_cells), n_clones=dplyr::n_distinct(cloneId),
                       donor=paste("Donor", paste(sort(as.integer(str_extract(unique(donor), "\\d"))),
                                                   collapse=", ")), .groups="drop") %>%
      dplyr::transmute(donor,
        clone = dplyr::if_else(n_clones==1L, "minor clone", sprintf("%d minor clones", n_clones)),
        n_cells, trb=tidyr::replace_na(row_trb,""), tra=tidyr::replace_na(row_tra,""),
        kind="minor", clone_num=NA_integer_)
  }
  if (nrow(minor_rows) == 0)
    minor_rows <- tibble(donor=character(), clone=character(), n_cells=integer(),
                         trb=character(), tra=character(), kind=character(), clone_num=integer())
  if (nrow(global_rows) == 0)
    global_rows <- tibble(global_id=character(), donor=character(), clone=character(),
                          n_cells=integer(), trb=character(), tra=character(),
                          kind=character(), clone_num=integer())

  dplyr::bind_rows(global_rows, minor_rows) %>%
    dplyr::mutate(trb=tidyr::replace_na(trb,""), tra=tidyr::replace_na(tra,"")) %>%
    dplyr::filter(n_cells >= MIN_CELLS) %>%                       # ROW filter (major + minor)
    dplyr::arrange(dplyr::desc(kind=="global"), clone_num, dplyr::desc(n_cells)) %>%
    dplyr::select(kind, Donor=donor, Clone=clone, `Number of Cells`=n_cells,
                  `TRB cdr3`=trb, `TRA cdr3(s)`=tra)
}

# ── 2. Styled tableGrob ──────────────────────────────────────────────────────
make_table_grob <- function(d, cl) {
  is_glob <- d$kind == "global"
  d_show  <- d %>% dplyr::select(Donor, Clone, `Number of Cells`, `TRB cdr3`, `TRA cdr3(s)`)
  hjust_vec <- c(0,0,0.5,0,0); x_vec <- c(0.02,0.02,0.5,0.02,0.02)
  th <- ttheme_minimal(base_size=9,
    core=list(fg_params=list(hjust=rep(hjust_vec,each=nrow(d_show)),
                             x=rep(x_vec,each=nrow(d_show)), fontsize=9),
              bg_params=list(fill=ifelse(is_glob,"#EAF2FB","white"), col="grey90", lwd=0.5)),
    colhead=list(fg_params=list(fontface="bold", hjust=hjust_vec, x=x_vec, fontsize=9),
                 bg_params=list(fill="grey95", col="grey80", lwd=0.5)))
  tg <- tableGrob(d_show, rows=NULL, theme=th)
  for (ccx in c(4,5)) { k <- which(tg$layout$name=="core-fg" & tg$layout$l==ccx)
    for (kk in k) tg$grobs[[kk]]$gp <- gpar(fontfamily="mono", fontsize=8.5) }
  for (rr in which(is_glob)) for (ccx in 1:2) {
    k <- which(tg$layout$name=="core-fg" & tg$layout$t==rr+1 & tg$layout$l==ccx)
    if (length(k)) tg$grobs[[k]]$gp <- gpar(fontface="bold", fontsize=9) }
  title <- textGrob(sprintf("Cluster %s", cl), gp=gpar(fontface="bold", fontsize=14), hjust=0, x=0.02)
  body <- gtable_add_rows(tg, unit(1.4,"lines"), pos=0)
  body <- gtable_add_grob(body, title, t=1, l=1, r=ncol(tg))
  list(grob = body, n_rows = nrow(d_show))
}

# ── 3. Contour UMAP ──────────────────────────────────────────────────────────
umap_all <- as.data.frame(Embeddings(mdata_CMV, "umap"))
colnames(umap_all)[1:2] <- c("UMAP_1","UMAP_2")
umap_all$cloneId   <- as.character(md$cloneId)
umap_all$global_id <- as.character(md$Global_Clone_ID)

make_umap <- function(cl, major_ids) {
  lab <- clone_labels(major_ids)
  clone_cells <- umap_all %>% dplyr::filter(global_id %in% major_ids) %>%
    dplyr::mutate(clone = factor(lab[global_id], levels = unname(lab)))
  ggplot() +
    geom_point(data = umap_all, aes(UMAP_1, UMAP_2), color="grey70", size=0.35, stroke=0) +
    geom_density_2d(data = clone_cells, aes(UMAP_1, UMAP_2, color=clone),
                    linewidth=0.4, bins=10) +
    scale_color_brewer(palette="Set1", name="Major clone", drop=FALSE) +
    theme_classic(base_size=11) +
    theme(legend.position="right", legend.key.size=unit(0.8,"lines"), aspect.ratio = 1) +
    labs(title=NULL, x="UMAP 1", y="UMAP 2") +
    guides(color=guide_legend(override.aes=list(linewidth=1.2)))
}

# ── 4. Metacluster proportion bars ───────────────────────────────────────────
make_bars <- function(cl, major_ids) {
  lab <- clone_labels(major_ids)
  bd <- cell_tbl %>% dplyr::filter(cluster==cl, global_id %in% major_ids, !is.na(metacluster)) %>%
    dplyr::mutate(clone = factor(lab[global_id], levels = unname(lab)),
                  metacluster = factor(metacluster, levels=META_LEVELS)) %>%
    dplyr::count(clone, metacluster, name="n") %>%
    dplyr::group_by(clone) %>% dplyr::mutate(prop = 100*n/sum(n)) %>% dplyr::ungroup()
  ggplot(bd, aes(clone, prop, fill=metacluster)) +
    geom_col(width=0.75, color="white", linewidth=0.2) +
    scale_fill_manual(values=meta_palette, name="Metacluster", drop=FALSE) +
    scale_y_continuous(expand=expansion(mult=c(0,0.02))) +
    labs(title=NULL, x=NULL, y="% of cells") +
    theme_classic(base_size=11) +
    theme(axis.text.x = element_text(angle=30, hjust=1), legend.key.size=unit(0.8,"lines"))
}

# ── 5. Compose one card ──────────────────────────────────────────────────────
render_card <- function(cl) {
  d <- build_cluster_table_df(cl)
  major_ids <- cell_tbl %>% dplyr::filter(cluster==cl, is_global) %>%
    dplyr::count(global_id, name="n") %>% dplyr::filter(n >= MIN_CELLS) %>% dplyr::pull(global_id)
  if (length(major_ids) == 0) { message(sprintf("cluster %s: no major clones >= %d cells, skipped", cl, MIN_CELLS)); return(invisible()) }
  if (nrow(d) == 0)            { message(sprintf("cluster %s: empty table, skipped", cl)); return(invisible()) }

  tg <- make_table_grob(d, cl)
  tbl_h_in <- grid::convertHeight(sum(tg$grob$heights), "in", valueOnly = TRUE) + 0.15
  umap_panel <- make_umap(cl, major_ids); bar_panel <- make_bars(cl, major_ids)
  bottom_h_in <- 6.0
  bottom <- cowplot::plot_grid(umap_panel, bar_panel, ncol = 2, rel_widths = c(0.65, 0.35))
  card <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_grob(tg$grob),
                             bottom, ncol = 1, rel_heights = c(tbl_h_in, bottom_h_in))
  dom <- cluster_spec$dom_spec[match(cl, cluster_spec$cluster)]; dom <- ifelse(is.na(dom),"NA",dom)
  total_h <- tbl_h_in + bottom_h_in + 0.2
  fn <- file.path(out_dir, sprintf("%s_cluster_%s_card.png", dom, cl))
  ggsave(fn, card, width = 13, height = total_h, dpi = 200, limitsize = FALSE, bg="white")
  message(sprintf("[%s] cluster %s: %d rows, %d major clones -> %s",
                  dom, cl, tg$n_rows, length(major_ids), basename(fn)))
}

# ── 6. Render all qualifying clusters ────────────────────────────────────────
kc <- keep_clusters[order(suppressWarnings(as.integer(keep_clusters)), keep_clusters)]
invisible(lapply(kc, render_card))
message(sprintf("\nDone: cards in %s", out_dir))



#
#
# Per-cluster TCR summary tables — cloneId-collapsed minor rows
#   Global clones: one row each; TRB = most-cells beta; TRA = top-2 most-cells alphas.
#   Minor cells: collapse to cloneId first (each cloneId -> its most-cells TRB + top-2 TRA),
#     then group cloneIds sharing the same (TRB, TRA-set) representative into one row,
#     MERGING ACROSS DONORS. Label "minor clone" (1 cloneId) or "N minor clones" (merged);
#     Donor column lists all contributing donors. notCalled excluded. All clusters, no >=5 filter.
#   Sequences from t_cdr3s_aa. Output: exports/TCR_clusters/cluster_tables/
#
#
library(dplyr); library(stringr); library(tidyr); library(purrr)
library(gridExtra); library(grid); library(gtable)

out_dir <- "exports/TCR_clusters/cluster_tables"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helpers ──────────────────────────────────────────────────────────────────
# ALL per-cell TRB/TRA occurrences across a set of cells' t_cdr3s_aa, as a long
# vector (one entry per read; dual-alpha cells contribute 2 TRA entries).
pull_chain_reads <- function(strings, prefix) {
  unlist(lapply(strings, function(s) {
    if (is.na(s)) return(character(0))
    parts <- str_split(s, ";")[[1]]; hit <- parts[str_detect(parts, paste0("^", prefix, ":"))]
    seqs  <- str_replace(hit, paste0("^", prefix, ":"), ""); seqs[nchar(seqs) > 0] }))
}
# most-frequent TRB (by cell) across a set of cells
top_trb <- function(strings) {
  r <- pull_chain_reads(strings, "TRB"); if (!length(r)) return("")
  names(which.max(table(r)))
}
# top-N most-frequent TRAs (by read count), "; "-joined
top_tra <- function(strings, n = 2) {
  r <- pull_chain_reads(strings, "TRA"); if (!length(r)) return("")
  tb <- sort(table(r), decreasing = TRUE)
  paste(names(tb)[seq_len(min(n, length(tb)))], collapse = "; ")
}
fmt_clone  <- function(gid) { spec <- str_match(gid, "_(pp65|IE1)-")[,2]; num <- str_match(gid, "clone-(\\d+)")[,2]
  ifelse(is.na(spec)|is.na(num), gid, sprintf("%s - clone %s", spec, num)) }
full_donor <- function(d) paste("Donor", str_extract(as.character(d), "\\d"))

# ── 1. Per-cell working table ────────────────────────────────────────────────
md <- mdata_CMV@meta.data
cell_tbl <- tibble(
  cluster=as.character(md$TCR_cluster_et64), cloneId=as.character(md$cloneId),
  donor=full_donor(md$donor_identity), spec=as.character(md$specificity),
  global_id=as.character(md$Global_Clone_ID), cdr3str=as.character(md$t_cdr3s_aa)) %>%
  dplyr::filter(!is.na(cluster), !grepl("notCalled", cloneId)) %>%
  dplyr::mutate(is_global = !is.na(global_id))

# dominant specificity per cluster (filename prefix)
cluster_spec <- cell_tbl %>% dplyr::filter(!is.na(spec)) %>%
  dplyr::count(cluster, spec, name="n") %>% dplyr::group_by(cluster) %>%
  dplyr::summarise(dom_spec = spec[which.max(n)], .groups="drop")

# ── 2. Build table df for ONE cluster ────────────────────────────────────────
#   Minor rows built by shared-TRB grouping with fragment absorption:
#     - each distinct dominant TRA within a TRB group -> its own (TRB,TRA) row
#     - TRB-only cloneIds -> dominant TRA row of their TRB group (or a blank-TRA
#       row if the group has no TRA at all)
#     - TRA-only cloneIds -> absorbed into a matching-TRA row if one exists,
#       else grouped into their own TRA-only (blank-TRB) row
#   Each cloneId is assigned to exactly one row (its dominant TRA) so cells aren't
#   double counted; rows merge cloneIds across donors.
build_cluster_table_df <- function(cl) {
  cc <- cell_tbl %>% dplyr::filter(cluster == cl)

  # -- GLOBAL rows: one per Global clone; 1 TRB (most cells) + top-2 TRA --
  global_rows <- cc %>% dplyr::filter(is_global) %>%
    dplyr::group_by(global_id) %>%
    dplyr::summarise(donor=full_donor(str_extract(dplyr::first(donor),"\\d")),
                     clone=fmt_clone(dplyr::first(global_id)), n_cells=dplyr::n(),
                     trb=top_trb(cdr3str), tra=top_tra(cdr3str, 2),
                     .groups="drop") %>%
    dplyr::mutate(kind="global",
                  clone_num=suppressWarnings(as.integer(str_match(global_id,"clone-(\\d+)")[,2])))

  # -- MINOR: per-cloneId representative (dominant TRB, dominant TRA, cells, donor) --
  mc <- cc %>% dplyr::filter(!is_global) %>%
    dplyr::group_by(cloneId) %>%
    dplyr::summarise(donor   = full_donor(str_extract(dplyr::first(donor),"\\d")),
                     n_cells = dplyr::n(),
                     trb     = top_trb(cdr3str),      # "" if none
                     tra     = top_tra(cdr3str, 1),   # dominant single TRA for keying ("" if none)
                     .groups = "drop") %>%
    dplyr::mutate(trb = dplyr::na_if(trb, ""), tra = dplyr::na_if(tra, ""))

  minor_rows <- tibble()
  if (nrow(mc) > 0) {
    has_trb <- mc %>% dplyr::filter(!is.na(trb))
    no_trb  <- mc %>% dplyr::filter(is.na(trb))    # TRA-only (or nothing) cloneIds

    # ---- rows from TRB groups ----
    trb_rows <- list()
    for (b in unique(has_trb$trb)) {
      grp    <- has_trb %>% dplyr::filter(trb == b)
      with_a <- grp %>% dplyr::filter(!is.na(tra))
      only_b <- grp %>% dplyr::filter(is.na(tra))
      if (nrow(with_a) == 0) {
        # group has TRB only, no alpha anywhere -> single blank-TRA row
        trb_rows[[length(trb_rows)+1]] <- only_b %>%
          dplyr::mutate(row_trb = b, row_tra = NA_character_)
      } else {
        # one row per distinct dominant TRA; TRB-only cloneIds go to dominant TRA
        tra_sizes <- with_a %>% dplyr::group_by(tra) %>%
          dplyr::summarise(cells=sum(n_cells), .groups="drop") %>%
          dplyr::arrange(dplyr::desc(cells))
        dom_a   <- tra_sizes$tra[1]
        with_a2 <- with_a %>% dplyr::mutate(row_trb = b, row_tra = tra)
        only_b2 <- only_b %>% dplyr::mutate(row_trb = b, row_tra = dom_a)   # absorb into dominant
        trb_rows[[length(trb_rows)+1]] <- dplyr::bind_rows(with_a2, only_b2)
      }
    }
    # guard: cluster whose minor cells are entirely TRA-only -> empty trb_rows
    trb_assigned <- if (length(trb_rows) == 0)
      tibble(cloneId=character(), donor=character(), n_cells=integer(),
             trb=character(), tra=character(), row_trb=character(), row_tra=character())
      else dplyr::bind_rows(trb_rows)

    # ---- TRA-only cloneIds: absorb into a matching-TRA row if one exists ----
    existing_tra <- unique(stats::na.omit(trb_assigned$row_tra))
    ta_match <- no_trb %>% dplyr::filter(!is.na(tra) & tra %in% existing_tra) %>%
      dplyr::mutate(row_trb = vapply(tra, function(a) {
          cand <- trb_assigned %>% dplyr::filter(row_tra == a) %>%
            dplyr::group_by(row_trb) %>% dplyr::summarise(c=sum(n_cells), .groups="drop") %>%
            dplyr::arrange(dplyr::desc(c))
          if (nrow(cand) == 0) NA_character_ else as.character(cand$row_trb[1])
        }, character(1)),
        row_tra = tra)
    ta_own <- no_trb %>% dplyr::filter(is.na(tra) | !(tra %in% existing_tra)) %>%
      dplyr::mutate(row_trb = NA_character_, row_tra = tra)   # own TRA-only rows (blank TRB)

    assigned <- dplyr::bind_rows(trb_assigned, ta_match, ta_own)

    # ---- aggregate rows: merge cloneIds across donors ----
    minor_rows <- assigned %>%
      dplyr::group_by(row_trb, row_tra) %>%
      dplyr::summarise(n_cells  = sum(n_cells),
                       n_clones = dplyr::n_distinct(cloneId),
                       donor    = paste(sort(unique(donor)), collapse = ", "),
                       .groups="drop") %>%
      dplyr::transmute(
        donor,
        clone     = dplyr::if_else(n_clones == 1L, "minor clone",
                                   sprintf("%d minor clones", n_clones)),
        n_cells,
        trb       = tidyr::replace_na(row_trb, ""),
        tra       = tidyr::replace_na(row_tra, ""),
        kind      = "minor",
        clone_num = NA_integer_)
  }

  # empty-safe typing so bind_rows never sees a logical/list column
  if (nrow(minor_rows) == 0)
    minor_rows <- tibble(donor=character(), clone=character(), n_cells=integer(),
                         trb=character(), tra=character(), kind=character(), clone_num=integer())
  if (nrow(global_rows) == 0)
    global_rows <- tibble(global_id=character(), donor=character(), clone=character(),
                          n_cells=integer(), trb=character(), tra=character(),
                          kind=character(), clone_num=integer())

  dplyr::bind_rows(global_rows, minor_rows) %>%
    dplyr::mutate(trb=tidyr::replace_na(trb,""), tra=tidyr::replace_na(tra,"")) %>%
    dplyr::arrange(dplyr::desc(kind=="global"), clone_num, dplyr::desc(n_cells)) %>%
    dplyr::select(kind, Donor=donor, Clone=clone, `Number of Cells`=n_cells,
                  `TRB cdr3`=trb, `TRA cdr3(s)`=tra)
}

# ── 3. Styled tableGrob ──────────────────────────────────────────────────────
make_table_grob <- function(d, cl) {
  is_glob <- d$kind == "global"
  d_show  <- d %>% dplyr::select(Donor, Clone, `Number of Cells`, `TRB cdr3`, `TRA cdr3(s)`)
  hjust_vec <- c(0,0,0.5,0,0); x_vec <- c(0.02,0.02,0.5,0.02,0.02)
  th <- ttheme_minimal(base_size=9,
    core=list(fg_params=list(hjust=rep(hjust_vec,each=nrow(d_show)),
                             x=rep(x_vec,each=nrow(d_show)), fontsize=9),
              bg_params=list(fill=ifelse(is_glob,"#EAF2FB","white"), col="grey90", lwd=0.5)),
    colhead=list(fg_params=list(fontface="bold", hjust=hjust_vec, x=x_vec, fontsize=9),
                 bg_params=list(fill="grey95", col="grey80", lwd=0.5)))
  tg <- tableGrob(d_show, rows=NULL, theme=th)
  for (ccx in c(4,5)) { k <- which(tg$layout$name=="core-fg" & tg$layout$l==ccx)
    for (kk in k) tg$grobs[[kk]]$gp <- gpar(fontfamily="mono", fontsize=8.5) }
  for (rr in which(is_glob)) for (ccx in 1:2) {
    k <- which(tg$layout$name=="core-fg" & tg$layout$t==rr+1 & tg$layout$l==ccx)
    if (length(k)) tg$grobs[[k]]$gp <- gpar(fontface="bold", fontsize=9) }
  title <- textGrob(sprintf("Cluster %s", cl), gp=gpar(fontface="bold", fontsize=14), hjust=0, x=0.02)
  body <- gtable_add_rows(tg, unit(1.4,"lines"), pos=0)
  body <- gtable_add_grob(body, title, t=1, l=1, r=ncol(tg))
  list(grob = body, n_rows = nrow(d_show))
}

# ── 4. Render one PNG per cluster (height measured, no clipping) ─────────────
render_table <- function(cl) {
  d <- build_cluster_table_df(cl)
  if (nrow(d) == 0) { message(sprintf("cluster %s: empty, skipped", cl)); return(invisible()) }
  tg <- make_table_grob(d, cl)
  tbl_h_in <- grid::convertHeight(sum(tg$grob$heights), "in", valueOnly = TRUE) + 0.3
  dom <- cluster_spec$dom_spec[match(cl, cluster_spec$cluster)]; dom <- ifelse(is.na(dom),"NA",dom)
  fn <- file.path(out_dir, sprintf("%s_cluster_%s_table.png", dom, cl))
  ggplot2::ggsave(fn, tg$grob, width = 13, height = tbl_h_in, dpi = 200, limitsize = FALSE, bg = "white")
  message(sprintf("[%s] cluster %s: %d rows (%d global) -> %s",
                  dom, cl, tg$n_rows, sum(d$kind=="global"), basename(fn)))
}

# ── 5. All clusters ──────────────────────────────────────────────────────────
clusters <- cell_tbl %>% dplyr::distinct(cluster) %>% dplyr::pull(cluster)
clusters <- clusters[order(suppressWarnings(as.integer(clusters)), clusters)]
invisible(lapply(clusters, render_table))
message(sprintf("\nDone: %d tables in %s", length(clusters), out_dir))




#
#
# Per-cluster TCR summary tables — FILTERED
#   Row filter: keep rows (major AND minor) with >= MIN_CELLS cells.
#   Cluster filter: render only if, after row filtering, the cluster has
#     >= 1 major (Global) row AND >= 2 rows total.
#   Minor rows: cloneId-collapsed, shared-TRB grouping w/ fragment absorption,
#     merged across donors (as in the current logic). Global: 1 TRB + top-2 TRA.
#   Output: exports/TCR_clusters/cluster_tables/
#
#
library(dplyr); library(stringr); library(tidyr); library(purrr)
library(gridExtra); library(grid); library(gtable)

out_dir <- "exports/TCR_clusters/cluster_tables/cluster_tables_filtered"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
MIN_CELLS <- 5

# ── Helpers ──────────────────────────────────────────────────────────────────
pull_chain_reads <- function(strings, prefix) {
  unlist(lapply(strings, function(s) {
    if (is.na(s)) return(character(0))
    parts <- str_split(s, ";")[[1]]; hit <- parts[str_detect(parts, paste0("^", prefix, ":"))]
    seqs  <- str_replace(hit, paste0("^", prefix, ":"), ""); seqs[nchar(seqs) > 0] }))
}
top_trb <- function(strings) { r <- pull_chain_reads(strings, "TRB"); if (!length(r)) return("")
  names(which.max(table(r))) }
top_tra <- function(strings, n = 2) { r <- pull_chain_reads(strings, "TRA"); if (!length(r)) return("")
  tb <- sort(table(r), decreasing = TRUE); paste(names(tb)[seq_len(min(n, length(tb)))], collapse = "; ") }
fmt_clone  <- function(gid) { spec <- str_match(gid, "_(pp65|IE1)-")[,2]; num <- str_match(gid, "clone-(\\d+)")[,2]
  ifelse(is.na(spec)|is.na(num), gid, sprintf("%s - clone %s", spec, num)) }
full_donor <- function(d) paste("Donor", str_extract(as.character(d), "\\d"))

# ── 1. Per-cell working table ────────────────────────────────────────────────
md <- mdata_CMV@meta.data
cell_tbl <- tibble(
  cluster=as.character(md$TCR_cluster_et64), cloneId=as.character(md$cloneId),
  donor=full_donor(md$donor_identity), spec=as.character(md$specificity),
  global_id=as.character(md$Global_Clone_ID), cdr3str=as.character(md$t_cdr3s_aa)) %>%
  dplyr::filter(!is.na(cluster), !grepl("notCalled", cloneId)) %>%
  dplyr::mutate(is_global = !is.na(global_id))

cluster_spec <- cell_tbl %>% dplyr::filter(!is.na(spec)) %>%
  dplyr::count(cluster, spec, name="n") %>% dplyr::group_by(cluster) %>%
  dplyr::summarise(dom_spec = spec[which.max(n)], .groups="drop")

# ── 2. Build table df for ONE cluster (cloneId-collapsed minor rows) ─────────
build_cluster_table_df <- function(cl) {
  cc <- cell_tbl %>% dplyr::filter(cluster == cl)

  global_rows <- cc %>% dplyr::filter(is_global) %>%
    dplyr::group_by(global_id) %>%
    dplyr::summarise(donor=full_donor(str_extract(dplyr::first(donor),"\\d")),
                     clone=fmt_clone(dplyr::first(global_id)), n_cells=dplyr::n(),
                     trb=top_trb(cdr3str), tra=top_tra(cdr3str, 2),
                     .groups="drop") %>%
    dplyr::mutate(kind="global",
                  clone_num=suppressWarnings(as.integer(str_match(global_id,"clone-(\\d+)")[,2])))

  mc <- cc %>% dplyr::filter(!is_global) %>%
    dplyr::group_by(cloneId) %>%
    dplyr::summarise(donor   = full_donor(str_extract(dplyr::first(donor),"\\d")),
                     n_cells = dplyr::n(),
                     trb     = top_trb(cdr3str),
                     tra     = top_tra(cdr3str, 1),
                     .groups = "drop") %>%
    dplyr::mutate(trb = dplyr::na_if(trb, ""), tra = dplyr::na_if(tra, ""))

  minor_rows <- tibble()
  if (nrow(mc) > 0) {
    has_trb <- mc %>% dplyr::filter(!is.na(trb))
    no_trb  <- mc %>% dplyr::filter(is.na(trb))
    trb_rows <- list()
    for (b in unique(has_trb$trb)) {
      grp    <- has_trb %>% dplyr::filter(trb == b)
      with_a <- grp %>% dplyr::filter(!is.na(tra))
      only_b <- grp %>% dplyr::filter(is.na(tra))
      if (nrow(with_a) == 0) {
        trb_rows[[length(trb_rows)+1]] <- only_b %>% dplyr::mutate(row_trb=b, row_tra=NA_character_)
      } else {
        tra_sizes <- with_a %>% dplyr::group_by(tra) %>%
          dplyr::summarise(cells=sum(n_cells), .groups="drop") %>% dplyr::arrange(dplyr::desc(cells))
        dom_a   <- tra_sizes$tra[1]
        with_a2 <- with_a %>% dplyr::mutate(row_trb=b, row_tra=tra)
        only_b2 <- only_b %>% dplyr::mutate(row_trb=b, row_tra=dom_a)
        trb_rows[[length(trb_rows)+1]] <- dplyr::bind_rows(with_a2, only_b2)
      }
    }
    trb_assigned <- if (length(trb_rows) == 0)
      tibble(cloneId=character(), donor=character(), n_cells=integer(),
             trb=character(), tra=character(), row_trb=character(), row_tra=character())
      else dplyr::bind_rows(trb_rows)

    existing_tra <- unique(stats::na.omit(trb_assigned$row_tra))
    ta_match <- no_trb %>% dplyr::filter(!is.na(tra) & tra %in% existing_tra) %>%
      dplyr::mutate(row_trb = vapply(tra, function(a) {
          cand <- trb_assigned %>% dplyr::filter(row_tra == a) %>%
            dplyr::group_by(row_trb) %>% dplyr::summarise(c=sum(n_cells), .groups="drop") %>%
            dplyr::arrange(dplyr::desc(c))
          if (nrow(cand) == 0) NA_character_ else as.character(cand$row_trb[1])
        }, character(1)),
        row_tra = tra)
    ta_own <- no_trb %>% dplyr::filter(is.na(tra) | !(tra %in% existing_tra)) %>%
      dplyr::mutate(row_trb = NA_character_, row_tra = tra)

    assigned <- dplyr::bind_rows(trb_assigned, ta_match, ta_own)
    minor_rows <- assigned %>%
      dplyr::group_by(row_trb, row_tra) %>%
      dplyr::summarise(n_cells=sum(n_cells), n_clones=dplyr::n_distinct(cloneId),
                       donor=paste("Donor", paste(sort(as.integer(str_extract(unique(donor), "\\d"))),
                                                  collapse=", ")), .groups="drop") %>%
      dplyr::transmute(donor,
        clone = dplyr::if_else(n_clones==1L, "minor clone", sprintf("%d minor clones", n_clones)),
        n_cells, trb=tidyr::replace_na(row_trb,""), tra=tidyr::replace_na(row_tra,""),
        kind="minor", clone_num=NA_integer_)
  }
  if (nrow(minor_rows) == 0)
    minor_rows <- tibble(donor=character(), clone=character(), n_cells=integer(),
                         trb=character(), tra=character(), kind=character(), clone_num=integer())
  if (nrow(global_rows) == 0)
    global_rows <- tibble(global_id=character(), donor=character(), clone=character(),
                          n_cells=integer(), trb=character(), tra=character(),
                          kind=character(), clone_num=integer())

  dplyr::bind_rows(global_rows, minor_rows) %>%
    dplyr::mutate(trb=tidyr::replace_na(trb,""), tra=tidyr::replace_na(tra,"")) %>%
    dplyr::arrange(dplyr::desc(kind=="global"), clone_num, dplyr::desc(n_cells)) %>%
    dplyr::select(kind, Donor=donor, Clone=clone, `Number of Cells`=n_cells,
                  `TRB cdr3`=trb, `TRA cdr3(s)`=tra)
}

# ── 3. Styled tableGrob ──────────────────────────────────────────────────────
make_table_grob <- function(d, cl) {
  is_glob <- d$kind == "global"
  d_show  <- d %>% dplyr::select(Donor, Clone, `Number of Cells`, `TRB cdr3`, `TRA cdr3(s)`)
  hjust_vec <- c(0,0,0.5,0,0); x_vec <- c(0.02,0.02,0.5,0.02,0.02)
  th <- ttheme_minimal(base_size=9,
    core=list(fg_params=list(hjust=rep(hjust_vec,each=nrow(d_show)),
                             x=rep(x_vec,each=nrow(d_show)), fontsize=9),
              bg_params=list(fill=ifelse(is_glob,"#EAF2FB","white"), col="grey90", lwd=0.5)),
    colhead=list(fg_params=list(fontface="bold", hjust=hjust_vec, x=x_vec, fontsize=9),
                 bg_params=list(fill="grey95", col="grey80", lwd=0.5)))
  tg <- tableGrob(d_show, rows=NULL, theme=th)
  for (ccx in c(4,5)) { k <- which(tg$layout$name=="core-fg" & tg$layout$l==ccx)
    for (kk in k) tg$grobs[[kk]]$gp <- gpar(fontfamily="mono", fontsize=8.5) }
  for (rr in which(is_glob)) for (ccx in 1:2) {
    k <- which(tg$layout$name=="core-fg" & tg$layout$t==rr+1 & tg$layout$l==ccx)
    if (length(k)) tg$grobs[[k]]$gp <- gpar(fontface="bold", fontsize=9) }
  title <- textGrob(sprintf("Cluster %s", cl), gp=gpar(fontface="bold", fontsize=14), hjust=0, x=0.02)
  body <- gtable_add_rows(tg, unit(1.4,"lines"), pos=0)
  body <- gtable_add_grob(body, title, t=1, l=1, r=ncol(tg))
  list(grob = body, n_rows = nrow(d_show))
}

# ── 4. Render one PNG per cluster (with row + cluster filters) ───────────────
render_table <- function(cl) {
  d <- build_cluster_table_df(cl)
  d <- d %>% dplyr::filter(`Number of Cells` >= MIN_CELLS)     # ROW filter (major + minor)
  n_major <- sum(d$kind == "global")
  if (n_major < 1 || nrow(d) < 2) {                            # CLUSTER filter
    message(sprintf("cluster %s: %d rows (%d major) after >=%d filter — skipped",
                    cl, nrow(d), n_major, MIN_CELLS)); return(invisible()) }
  tg <- make_table_grob(d, cl)
  tbl_h_in <- grid::convertHeight(sum(tg$grob$heights), "in", valueOnly = TRUE) + 0.3
  dom <- cluster_spec$dom_spec[match(cl, cluster_spec$cluster)]; dom <- ifelse(is.na(dom),"NA",dom)
  fn <- file.path(out_dir, sprintf("%s_cluster_%s_table.png", dom, cl))
  ggplot2::ggsave(fn, tg$grob, width = 13, height = tbl_h_in, dpi = 200, limitsize = FALSE, bg = "white")
  message(sprintf("[%s] cluster %s: %d rows (%d major) -> %s", dom, cl, nrow(d), n_major, basename(fn)))
}

# ── 5. All clusters ──────────────────────────────────────────────────────────
clusters <- cell_tbl %>% dplyr::distinct(cluster) %>% dplyr::pull(cluster)
clusters <- clusters[order(suppressWarnings(as.integer(clusters)), clusters)]
invisible(lapply(clusters, render_table))
message(sprintf("\nDone: tables in %s", out_dir))



#
#
# Cluster cards v2 — filtered table (top) + contour UMAP (left) + metacluster bars (right)
#   Table logic MATCHES the standalone filtered tables:
#     - Global rows: 1 TRB (most cells) + top-2 TRA
#     - Minor rows: cloneId-collapsed, shared-TRB grouping w/ fragment absorption,
#       merged across donors; compact donor label "Donor 2, 3, 4"
#     - Row filter: keep rows (major AND minor) with >= MIN_CELLS cells
#   Cards only for clusters with >= 2 unique TCR sequences (constellation min2seq set).
#   UMAP: per-clone KDE contours (bins=10) over grey full-object backdrop, square panel.
#   Bars: per-clone metacluster proportion (0-100%). Labels always include donor.
#   Output: exports/TCR_clusters/cluster_cards_v2/
#
#
library(Seurat); library(dplyr); library(stringr); library(tidyr); library(purrr)
library(ggplot2); library(gridExtra); library(grid); library(gtable); library(cowplot)

out_dir <- "exports/TCR_clusters/cluster_cards_v2"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
MIN_CELLS   <- 5
META_LEVELS <- c("Tcm","Tem","T-CMC1","NK-like","ISG+")
meta_palette <- c(Tcm="#E8A0A0", Tem="#7B8B3A", `T-CMC1`="#1B9E8A", `NK-like`="#888888", `ISG+`="#C050C0")
conga_dir   <- "exports/conga/mdata_CMV"

# ── Helpers ──────────────────────────────────────────────────────────────────
pull_chain_reads <- function(strings, prefix) {
  unlist(lapply(strings, function(s) {
    if (is.na(s)) return(character(0))
    parts <- str_split(s, ";")[[1]]; hit <- parts[str_detect(parts, paste0("^", prefix, ":"))]
    seqs  <- str_replace(hit, paste0("^", prefix, ":"), ""); seqs[nchar(seqs) > 0] }))
}
top_trb <- function(strings) { r <- pull_chain_reads(strings, "TRB"); if (!length(r)) return("")
  names(which.max(table(r))) }
top_tra <- function(strings, n = 2) { r <- pull_chain_reads(strings, "TRA"); if (!length(r)) return("")
  tb <- sort(table(r), decreasing = TRUE); paste(names(tb)[seq_len(min(n, length(tb)))], collapse = "; ") }
fmt_clone  <- function(gid) { spec <- str_match(gid, "_(pp65|IE1)-")[,2]; num <- str_match(gid, "clone-(\\d+)")[,2]
  ifelse(is.na(spec)|is.na(num), gid, sprintf("%s - clone %s", spec, num)) }
full_donor <- function(d) paste("Donor", str_extract(as.character(d), "\\d"))

clone_labels <- function(major_ids) {
  info <- tibble(global_id = major_ids,
                 base = fmt_clone(major_ids),
                 num  = suppressWarnings(as.integer(str_match(major_ids, "clone-(\\d+)")[,2])),
                 don  = str_match(major_ids, "^D(\\d)")[,2]) %>%
    dplyr::mutate(label = sprintf("%s (D%s)", base, don)) %>%
    dplyr::arrange(num, don)
  setNames(info$label, info$global_id)
}

# ── 0. Per-cell table + qualifying clusters ──────────────────────────────────
md <- mdata_CMV@meta.data
cell_tbl <- tibble(
  cluster=as.character(md$TCR_cluster_et64), cloneId=as.character(md$cloneId),
  donor=full_donor(md$donor_identity), spec=as.character(md$specificity),
  metacluster=as.character(md$metacluster),
  global_id=as.character(md$Global_Clone_ID), cdr3str=as.character(md$t_cdr3s_aa)) %>%
  dplyr::filter(!is.na(cluster), !grepl("notCalled", cloneId)) %>%
  dplyr::mutate(is_global = !is.na(global_id))

g_data <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
cluster_nseq <- g_data$nodes %>% dplyr::distinct(cluster, mrow) %>% dplyr::count(cluster, name="n_seq")
keep_clusters <- cluster_nseq %>% dplyr::filter(n_seq >= 2) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>% dplyr::pull(cluster)
message(sprintf("Clusters with >=2 unique sequences: %d", length(keep_clusters)))

cluster_spec <- cell_tbl %>% dplyr::filter(!is.na(spec)) %>%
  dplyr::count(cluster, spec, name="n") %>% dplyr::group_by(cluster) %>%
  dplyr::summarise(dom_spec = spec[which.max(n)], .groups="drop")

# ── 1. Table df for ONE cluster (matches standalone filtered logic) ──────────
build_cluster_table_df <- function(cl) {
  cc <- cell_tbl %>% dplyr::filter(cluster == cl)

  global_rows <- cc %>% dplyr::filter(is_global) %>%
    dplyr::group_by(global_id) %>%
    dplyr::summarise(donor=full_donor(str_extract(dplyr::first(donor),"\\d")),
                     clone=fmt_clone(dplyr::first(global_id)), n_cells=dplyr::n(),
                     trb=top_trb(cdr3str), tra=top_tra(cdr3str, 2), .groups="drop") %>%
    dplyr::mutate(kind="global",
                  clone_num=suppressWarnings(as.integer(str_match(global_id,"clone-(\\d+)")[,2])))

  mc <- cc %>% dplyr::filter(!is_global) %>%
    dplyr::group_by(cloneId) %>%
    dplyr::summarise(donor=full_donor(str_extract(dplyr::first(donor),"\\d")),
                     n_cells=dplyr::n(), trb=top_trb(cdr3str), tra=top_tra(cdr3str, 1),
                     .groups="drop") %>%
    dplyr::mutate(trb = dplyr::na_if(trb, ""), tra = dplyr::na_if(tra, ""))

  minor_rows <- tibble()
  if (nrow(mc) > 0) {
    has_trb <- mc %>% dplyr::filter(!is.na(trb)); no_trb <- mc %>% dplyr::filter(is.na(trb))
    trb_rows <- list()
    for (b in unique(has_trb$trb)) {
      grp <- has_trb %>% dplyr::filter(trb == b)
      with_a <- grp %>% dplyr::filter(!is.na(tra)); only_b <- grp %>% dplyr::filter(is.na(tra))
      if (nrow(with_a) == 0) {
        trb_rows[[length(trb_rows)+1]] <- only_b %>% dplyr::mutate(row_trb=b, row_tra=NA_character_)
      } else {
        tra_sizes <- with_a %>% dplyr::group_by(tra) %>%
          dplyr::summarise(cells=sum(n_cells), .groups="drop") %>% dplyr::arrange(dplyr::desc(cells))
        dom_a <- tra_sizes$tra[1]
        trb_rows[[length(trb_rows)+1]] <- dplyr::bind_rows(
          with_a %>% dplyr::mutate(row_trb=b, row_tra=tra),
          only_b %>% dplyr::mutate(row_trb=b, row_tra=dom_a))
      }
    }
    trb_assigned <- if (length(trb_rows) == 0)
      tibble(cloneId=character(), donor=character(), n_cells=integer(),
             trb=character(), tra=character(), row_trb=character(), row_tra=character())
      else dplyr::bind_rows(trb_rows)
    existing_tra <- unique(stats::na.omit(trb_assigned$row_tra))
    ta_match <- no_trb %>% dplyr::filter(!is.na(tra) & tra %in% existing_tra) %>%
      dplyr::mutate(row_trb = vapply(tra, function(a) {
          cand <- trb_assigned %>% dplyr::filter(row_tra == a) %>%
            dplyr::group_by(row_trb) %>% dplyr::summarise(c=sum(n_cells), .groups="drop") %>%
            dplyr::arrange(dplyr::desc(c))
          if (nrow(cand) == 0) NA_character_ else as.character(cand$row_trb[1]) }, character(1)),
        row_tra = tra)
    ta_own <- no_trb %>% dplyr::filter(is.na(tra) | !(tra %in% existing_tra)) %>%
      dplyr::mutate(row_trb=NA_character_, row_tra=tra)
    assigned <- dplyr::bind_rows(trb_assigned, ta_match, ta_own)
    minor_rows <- assigned %>%
      dplyr::group_by(row_trb, row_tra) %>%
      dplyr::summarise(n_cells=sum(n_cells), n_clones=dplyr::n_distinct(cloneId),
                       donor=paste("Donor", paste(sort(as.integer(str_extract(unique(donor), "\\d"))),
                                                   collapse=", ")), .groups="drop") %>%
      dplyr::transmute(donor,
        clone = dplyr::if_else(n_clones==1L, "minor clone", sprintf("%d minor clones", n_clones)),
        n_cells, trb=tidyr::replace_na(row_trb,""), tra=tidyr::replace_na(row_tra,""),
        kind="minor", clone_num=NA_integer_)
  }
  if (nrow(minor_rows) == 0)
    minor_rows <- tibble(donor=character(), clone=character(), n_cells=integer(),
                         trb=character(), tra=character(), kind=character(), clone_num=integer())
  if (nrow(global_rows) == 0)
    global_rows <- tibble(global_id=character(), donor=character(), clone=character(),
                          n_cells=integer(), trb=character(), tra=character(),
                          kind=character(), clone_num=integer())

  dplyr::bind_rows(global_rows, minor_rows) %>%
    dplyr::mutate(trb=tidyr::replace_na(trb,""), tra=tidyr::replace_na(tra,"")) %>%
    dplyr::filter(n_cells >= MIN_CELLS) %>%                       # ROW filter (major + minor)
    dplyr::arrange(dplyr::desc(kind=="global"), clone_num, dplyr::desc(n_cells)) %>%
    dplyr::select(kind, Donor=donor, Clone=clone, `Number of Cells`=n_cells,
                  `TRB cdr3`=trb, `TRA cdr3(s)`=tra)
}

# ── 2. Styled tableGrob ──────────────────────────────────────────────────────
make_table_grob <- function(d, cl) {
  is_glob <- d$kind == "global"
  d_show  <- d %>% dplyr::select(Donor, Clone, `Number of Cells`, `TRB cdr3`, `TRA cdr3(s)`)
  hjust_vec <- c(0,0,0.5,0,0); x_vec <- c(0.02,0.02,0.5,0.02,0.02)
  th <- ttheme_minimal(base_size=9,
    core=list(fg_params=list(hjust=rep(hjust_vec,each=nrow(d_show)),
                             x=rep(x_vec,each=nrow(d_show)), fontsize=9),
              bg_params=list(fill=ifelse(is_glob,"#EAF2FB","white"), col="grey90", lwd=0.5)),
    colhead=list(fg_params=list(fontface="bold", hjust=hjust_vec, x=x_vec, fontsize=9),
                 bg_params=list(fill="grey95", col="grey80", lwd=0.5)))
  tg <- tableGrob(d_show, rows=NULL, theme=th)
  for (ccx in c(4,5)) { k <- which(tg$layout$name=="core-fg" & tg$layout$l==ccx)
    for (kk in k) tg$grobs[[kk]]$gp <- gpar(fontfamily="mono", fontsize=8.5) }
  for (rr in which(is_glob)) for (ccx in 1:2) {
    k <- which(tg$layout$name=="core-fg" & tg$layout$t==rr+1 & tg$layout$l==ccx)
    if (length(k)) tg$grobs[[k]]$gp <- gpar(fontface="bold", fontsize=9) }
  title <- textGrob(sprintf("Cluster %s", cl), gp=gpar(fontface="bold", fontsize=14), hjust=0, x=0.02)
  body <- gtable_add_rows(tg, unit(1.4,"lines"), pos=0)
  body <- gtable_add_grob(body, title, t=1, l=1, r=ncol(tg))
  list(grob = body, n_rows = nrow(d_show))
}

# ── 3. Contour UMAP ──────────────────────────────────────────────────────────
umap_all <- as.data.frame(Embeddings(mdata_CMV, "umap"))
colnames(umap_all)[1:2] <- c("UMAP_1","UMAP_2")
umap_all$cloneId   <- as.character(md$cloneId)
umap_all$global_id <- as.character(md$Global_Clone_ID)

make_umap <- function(cl, major_ids) {
  lab <- clone_labels(major_ids)
  clone_cells <- umap_all %>% dplyr::filter(global_id %in% major_ids) %>%
    dplyr::mutate(clone = factor(lab[global_id], levels = unname(lab)))
  ggplot() +
    geom_point(data = umap_all, aes(UMAP_1, UMAP_2), color="grey70", size=0.35, stroke=0) +
    geom_density_2d(data = clone_cells, aes(UMAP_1, UMAP_2, color=clone),
                    linewidth=0.4, bins=10) +
    scale_color_brewer(palette="Set1", name="Major clone", drop=FALSE) +
    theme_classic(base_size=11) +
    theme(legend.position="right", legend.key.size=unit(0.8,"lines"), aspect.ratio = 1) +
    labs(title=NULL, x="UMAP 1", y="UMAP 2") +
    guides(color=guide_legend(override.aes=list(linewidth=1.2)))
}

# ── 4. Metacluster proportion bars ───────────────────────────────────────────
make_bars <- function(cl, major_ids) {
  lab <- clone_labels(major_ids)
  bd <- cell_tbl %>% dplyr::filter(cluster==cl, global_id %in% major_ids, !is.na(metacluster)) %>%
    dplyr::mutate(clone = factor(lab[global_id], levels = unname(lab)),
                  metacluster = factor(metacluster, levels=META_LEVELS)) %>%
    dplyr::count(clone, metacluster, name="n") %>%
    dplyr::group_by(clone) %>% dplyr::mutate(prop = 100*n/sum(n)) %>% dplyr::ungroup()
  ggplot(bd, aes(clone, prop, fill=metacluster)) +
    geom_col(width=0.75, color="white", linewidth=0.2) +
    scale_fill_manual(values=meta_palette, name="Metacluster", drop=FALSE) +
    scale_y_continuous(expand=expansion(mult=c(0,0.02))) +
    labs(title=NULL, x=NULL, y="% of cells") +
    theme_classic(base_size=11) +
    theme(axis.text.x = element_text(angle=30, hjust=1), legend.key.size=unit(0.8,"lines"))
}

# ── 5. Compose one card ──────────────────────────────────────────────────────
render_card <- function(cl) {
  d <- build_cluster_table_df(cl)
  major_ids <- cell_tbl %>% dplyr::filter(cluster==cl, is_global) %>%
    dplyr::count(global_id, name="n") %>% dplyr::filter(n >= MIN_CELLS) %>% dplyr::pull(global_id)
  if (length(major_ids) == 0) { message(sprintf("cluster %s: no major clones >= %d cells, skipped", cl, MIN_CELLS)); return(invisible()) }
  if (nrow(d) == 0)            { message(sprintf("cluster %s: empty table, skipped", cl)); return(invisible()) }

  tg <- make_table_grob(d, cl)
  tbl_h_in <- grid::convertHeight(sum(tg$grob$heights), "in", valueOnly = TRUE) + 0.15
  umap_panel <- make_umap(cl, major_ids); bar_panel <- make_bars(cl, major_ids)
  bottom_h_in <- 6.0
  bottom <- cowplot::plot_grid(umap_panel, bar_panel, ncol = 2, rel_widths = c(0.65, 0.35))
  card <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_grob(tg$grob),
                             bottom, ncol = 1, rel_heights = c(tbl_h_in, bottom_h_in))
  dom <- cluster_spec$dom_spec[match(cl, cluster_spec$cluster)]; dom <- ifelse(is.na(dom),"NA",dom)
  total_h <- tbl_h_in + bottom_h_in + 0.2
  fn <- file.path(out_dir, sprintf("%s_cluster_%s_card.png", dom, cl))
  ggsave(fn, card, width = 13, height = total_h, dpi = 200, limitsize = FALSE, bg="white")
  message(sprintf("[%s] cluster %s: %d rows, %d major clones -> %s",
                  dom, cl, tg$n_rows, length(major_ids), basename(fn)))
}

# ── 6. Render all qualifying clusters ────────────────────────────────────────
kc <- keep_clusters[order(suppressWarnings(as.integer(keep_clusters)), keep_clusters)]
invisible(lapply(kc, render_card))
message(sprintf("\nDone: cards in %s", out_dir))






#
# Constellation v2 — Stage 1 (REVISED): object nodes matching cluster_cards_v2,
#   now carrying per-donor (d1..d5) and per-metacluster cell counts per node.
#   Node = table row (Global clone or shared-TRB minor group), >=5 cells.
#   Output sequences unchanged -> existing 32x32 tcrdist matrix stays valid.
#   Writes: constellation_v2_nodes.rds + constellation_v2_tcrdist3_input.csv
#
library(dplyr); library(stringr); library(tidyr); library(readr)

conga_dir <- "exports/conga/mdata_CMV"
v2_dir    <- file.path(conga_dir, "v2"); dir.create(v2_dir, recursive=TRUE, showWarnings=FALSE)
MIN_CELLS   <- 5
META_LEVELS <- c("Tcm","Tem","T-CMC1","NK-like","ISG+")

pull_reads <- function(strings, prefix) unlist(lapply(strings, function(s){
  if (is.na(s)) return(character(0)); parts <- str_split(s,";")[[1]]
  hit <- parts[str_detect(parts, paste0("^",prefix,":"))]
  seqs <- str_replace(hit, paste0("^",prefix,":"),""); seqs[nchar(seqs)>0] }))
top1 <- function(strings, prefix){ r <- pull_reads(strings, prefix); if(!length(r)) return(NA_character_); names(which.max(table(r))) }
donorS_of <- function(d) paste0("D", str_extract(as.character(d), "\\d"))

md <- mdata_CMV@meta.data
cell_tbl <- tibble(
  cluster=as.character(md$TCR_cluster_et64), cloneId=as.character(md$cloneId),
  donorS=donorS_of(md$donor_identity), spec=as.character(md$specificity),
  metacluster=as.character(md$metacluster),
  global_id=as.character(md$Global_Clone_ID), cdr3str=as.character(md$t_cdr3s_aa)) %>%
  dplyr::filter(!is.na(cluster), !grepl("notCalled", cloneId)) %>%
  dplyr::mutate(is_global = !is.na(global_id))

# candidate clusters: >=2 unique file sequences (same candidate set as cards)
g_old <- readRDS(file.path(conga_dir,"constellation_graph.rds"))
cand <- g_old$nodes %>% dplyr::distinct(cluster, mrow) %>% dplyr::count(cluster, name="n") %>%
  dplyr::filter(n>=2) %>% dplyr::mutate(cluster=as.character(cluster)) %>% dplyr::pull(cluster)

# per-set donor + metacluster count columns from a set of cells (tibble w/ donorS, metacluster)
donor_meta_counts <- function(sub){
  d <- c(D1=0L,D2=0L,D3=0L,D4=0L,D5=0L)
  td <- table(sub$donorS); d[names(td)] <- as.integer(td)
  m <- setNames(integer(length(META_LEVELS)), META_LEVELS)
  tm <- table(factor(sub$metacluster, levels=META_LEVELS)); m[names(tm)] <- as.integer(tm)
  c(d, m)
}

build_nodes_one <- function(cl){
  cc <- cell_tbl %>% dplyr::filter(cluster==cl)

  # ---- GLOBAL nodes: one per clone ----
  gl_ids <- cc %>% dplyr::filter(is_global) %>%
    dplyr::count(global_id, name="n_cells") %>% dplyr::filter(n_cells>=MIN_CELLS) %>%
    dplyr::pull(global_id)
  gl <- lapply(gl_ids, function(gid){
    sub <- cc %>% dplyr::filter(global_id==gid)
    dm <- donor_meta_counts(sub)
    tibble(cluster=cl, node_kind="global", node_id=gid,
           spec=names(which.max(table(sub$spec))), n_cells=nrow(sub),
           trb=top1(sub$cdr3str,"TRB"), tra=top1(sub$cdr3str,"TRA"),
           D1=dm["D1"],D2=dm["D2"],D3=dm["D3"],D4=dm["D4"],D5=dm["D5"],
           Tcm=dm["Tcm"],Tem=dm["Tem"],`T-CMC1`=dm["T-CMC1"],`NK-like`=dm["NK-like"],`ISG+`=dm["ISG+"]) })
  gl <- if(length(gl)) bind_rows(gl) else tibble()

  # ---- MINOR nodes: per-cloneId dominant TRB/TRA, then shared-(TRB,TRA) grouping ----
  mc <- cc %>% dplyr::filter(!is_global) %>% dplyr::group_by(cloneId) %>%
    dplyr::summarise(trb=top1(cdr3str,"TRB"), tra=top1(cdr3str,"TRA"), .groups="drop") %>%
    dplyr::mutate(gkey=paste(tidyr::replace_na(trb,"NA"), tidyr::replace_na(tra,"NA")))
  mn <- tibble()
  if (nrow(mc) > 0){
    # map each minor cell to its cloneId's gkey, so we can pull member cells per group
    minor_cells <- cc %>% dplyr::filter(!is_global) %>%
      dplyr::left_join(mc %>% dplyr::select(cloneId, gkey, gtrb=trb, gtra=tra), by="cloneId")
    grp_ids <- minor_cells %>% dplyr::count(gkey, name="n_cells") %>%
      dplyr::filter(n_cells>=MIN_CELLS) %>% dplyr::pull(gkey)
    mn <- lapply(seq_along(grp_ids), function(i){
      gk <- grp_ids[i]; sub <- minor_cells %>% dplyr::filter(gkey==gk)
      dm <- donor_meta_counts(sub)
      tibble(cluster=cl, node_kind="minor", node_id=paste0("minor_",i),
             spec=names(which.max(table(sub$spec))), n_cells=nrow(sub),
             trb=dplyr::first(sub$gtrb), tra=dplyr::first(sub$gtra),
             D1=dm["D1"],D2=dm["D2"],D3=dm["D3"],D4=dm["D4"],D5=dm["D5"],
             Tcm=dm["Tcm"],Tem=dm["Tem"],`T-CMC1`=dm["T-CMC1"],`NK-like`=dm["NK-like"],`ISG+`=dm["ISG+"]) }) %>%
      bind_rows()
  }
  dplyr::bind_rows(gl, mn)
}

nodes <- bind_rows(lapply(cand, build_nodes_one))

# card gate: keep clusters with >=1 major clone >=5 cells
keep_cl <- nodes %>% dplyr::group_by(cluster) %>%
  dplyr::summarise(has_major=any(node_kind=="global"), .groups="drop") %>%
  dplyr::filter(has_major) %>% dplyr::pull(cluster)
nodes <- nodes %>% dplyr::filter(cluster %in% keep_cl) %>%
  dplyr::mutate(paired = !is.na(trb) & !is.na(tra),
                total_cells = n_cells,
                uid = paste(cluster, node_id, sep="::"))

message(sprintf("v2 nodes: %d across %d clusters | paired: %d | unpaired: %d",
                nrow(nodes), dplyr::n_distinct(nodes$cluster), sum(nodes$paired), sum(!nodes$paired)))

saveRDS(nodes, file.path(v2_dir,"constellation_v2_nodes.rds"))
paired <- nodes %>% dplyr::filter(paired) %>% dplyr::distinct(trb, tra)
write_csv(paired %>% dplyr::transmute(cdr3_b_aa=trb, cdr3_a_aa=tra),
          file.path(v2_dir,"constellation_v2_tcrdist3_input.csv"))
message(sprintf("Unique paired sequences for tcrdist3: %d", nrow(paired)))



#
# Constellation v2 — Stage 3: map nodes -> matrix rows, save graph object.
#   Paired nodes matched to their tcrdist row by (trb,tra); unpaired flagged (mrow=NA).
#   Output list(nodes, D) mirrors constellation_graph.rds for the plotting code.
#
library(dplyr); library(readr)

conga_dir <- "exports/conga/mdata_CMV"; v2_dir <- file.path(conga_dir,"v2")

nodes <- readRDS(file.path(v2_dir,"constellation_v2_nodes.rds"))
D     <- as.matrix(read_csv(file.path(v2_dir,"constellation_v2_tcrdist_matrix.csv"),
                            col_names=FALSE, show_col_types=FALSE)); dimnames(D) <- NULL
ord   <- read_csv(file.path(v2_dir,"constellation_v2_tcrdist_order.csv"), show_col_types=FALSE)

# map paired nodes to matrix row (1-based) via (trb,tra); unpaired -> NA
seq_key  <- paste(ord$cdr3_b_aa, ord$cdr3_a_aa)
node_key <- paste(nodes$trb, nodes$tra)
nodes$mrow <- ifelse(nodes$paired, match(node_key, seq_key), NA_integer_)
stopifnot(all(!is.na(nodes$mrow[nodes$paired])))          # every paired node must map
stopifnot(nrow(D) == length(seq_key))

saveRDS(list(nodes=nodes, D=D), file.path(v2_dir,"constellation_v2_graph.rds"))
message(sprintf("v2 graph saved: %d nodes (%d mapped to %dx%d matrix, %d unpaired/disconnected)",
                nrow(nodes), sum(nodes$paired), nrow(D), ncol(D), sum(!nodes$paired)))


                
#
# constellation_v2 — rebuild graph with NODE = unique (TRB,TRA) combo (merges clones)
#   Matches the annotated plots' node logic: per-cloneId dominant TRB/TRA, merge on
#   exact (TRB,TRA) within a cluster, sum donor + metacluster counts, keep >=5 cells.
#   Writes constellation_v2_graph_merged.rds (same shape as the plotting code expects).
#
library(dplyr); library(stringr); library(tidyr); library(readr)

conga_dir <- "exports/conga/mdata_CMV"; v2_dir <- file.path(conga_dir,"v2")
MIN_CELLS <- 5; META_LEVELS <- c("Tcm","Tem","T-CMC1","NK-like","ISG+")

g <- readRDS(file.path(v2_dir,"constellation_v2_graph.rds")); D <- g$D
ord <- read_csv(file.path(v2_dir,"constellation_v2_tcrdist_order.csv"), show_col_types=FALSE)
seq_key_matrix <- paste(ord$cdr3_b_aa, ord$cdr3_a_aa)

pull_reads <- function(strings, prefix) unlist(lapply(strings, function(s){
  if (is.na(s)) return(character(0)); parts<-str_split(s,";")[[1]]
  hit<-parts[str_detect(parts,paste0("^",prefix,":"))]; seqs<-str_replace(hit,paste0("^",prefix,":"),"")
  seqs[nchar(seqs)>0] }))
top1 <- function(strings, prefix){ r<-pull_reads(strings,prefix); if(!length(r)) return(NA_character_); names(which.max(table(r))) }

md <- mdata_CMV@meta.data %>% dplyr::mutate(cloneId=as.character(cloneId),
        cluster=as.character(TCR_cluster_et64), global_id=as.character(Global_Clone_ID),
        donorS=paste0("D", str_extract(as.character(donor_identity),"\\d")),
        spec=as.character(specificity), metacluster=as.character(metacluster))

# which clusters were in the v2 set (same candidates as before)
v2_clusters <- unique(as.character(g$nodes$cluster))

build_merged_one <- function(cl){
  cc <- md %>% dplyr::filter(cluster==cl, !grepl("notCalled", cloneId))
  if (nrow(cc)==0) return(NULL)
  per_clone <- cc %>% dplyr::group_by(cloneId) %>%
    dplyr::summarise(trb=top1(t_cdr3s_aa,"TRB"), tra=top1(t_cdr3s_aa,"TRA"),
                     global_id=ifelse(any(!is.na(global_id)), global_id[!is.na(global_id)][1], NA_character_),
                     .groups="drop") %>%
    dplyr::mutate(nkey=paste(tidyr::replace_na(trb,"NA"), tidyr::replace_na(tra,"NA"), sep="||"))
  cell_nkey <- cc %>% dplyr::left_join(per_clone %>% dplyr::select(cloneId,nkey,trb,tra), by="cloneId")
  cell_nkey %>% dplyr::group_by(nkey) %>%
    dplyr::summarise(cluster=cl, trb=dplyr::first(trb), tra=dplyr::first(tra),
      total_cells=dplyr::n(), spec=names(which.max(table(spec[!is.na(spec)]))),
      n_major=dplyr::n_distinct(global_id[!is.na(global_id)]),
      n_minor=dplyr::n_distinct(cloneId[is.na(global_id)]),
      D1=sum(donorS=="D1"),D2=sum(donorS=="D2"),D3=sum(donorS=="D3"),D4=sum(donorS=="D4"),D5=sum(donorS=="D5"),
      Tcm=sum(metacluster=="Tcm",na.rm=TRUE), Tem=sum(metacluster=="Tem",na.rm=TRUE),
      `T-CMC1`=sum(metacluster=="T-CMC1",na.rm=TRUE), `NK-like`=sum(metacluster=="NK-like",na.rm=TRUE),
      `ISG+`=sum(metacluster=="ISG+",na.rm=TRUE), .groups="drop")
}

merged <- bind_rows(lapply(v2_clusters, build_merged_one)) %>%
  dplyr::filter(total_cells >= MIN_CELLS) %>%
  dplyr::mutate(paired = !is.na(trb) & !is.na(tra),
                total_cells = total_cells,
                uid = paste(cluster, nkey, sep="::"),
                mrow = ifelse(paired, match(paste(trb,tra), seq_key_matrix), NA_integer_))

message(sprintf("merged nodes: %d across %d clusters | paired %d, unpaired %d",
                nrow(merged), dplyr::n_distinct(merged$cluster), sum(merged$paired), sum(!merged$paired)))

saveRDS(list(nodes=merged, D=D), file.path(v2_dir,"constellation_v2_graph_merged.rds"))

#
# constellation_v2 — overview plots on MERGED (unique-sequence) nodes
#   Same overlays as before (donor, kNN, 5 metacluster %), now with merged nodes.
#   Unpaired nodes shown disconnected. Output: constellation_v2/
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr); library(readr); library(tidyr)

conga_dir <- "exports/conga/mdata_CMV"; v2_dir <- file.path(conga_dir,"v2")
out_dir   <- "exports/TCR_clusters/constellation_v2"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

g <- readRDS(file.path(v2_dir,"constellation_v2_graph_merged.rds"))
node_tbl <- g$nodes; D <- g$D
donor_color_map <- c("D1"="#1b9e77","D2"="#d95f02","D3"="#7570b3","D4"="#e7298a","D5"="#66a61e")

EDGE_CUT<-250; DIST_TRANSFORM<-function(d) log1p(d)
NODE_R_MIN<-0.35; NODE_R_MAX<-1.40; CLUST_PAD<-3; GAP_MIN<-0.25
EDGE_FAR_COL<-"grey80"; EDGE_FAR_LW<-0.25; EDGE_NEAR_COL<-"grey35"; EDGE_NEAR_LW<-0.7
LABEL_OFFSET<-0.8; LABEL_SIZE<-3.5; NA_FILL<-"grey85"
META_LEVELS<-c("Tcm","Tem","T-CMC1","NK-like","ISG+")
meta_palette<-c(Tcm="#E8A0A0",Tem="#7B8B3A",`T-CMC1`="#1B9E8A",`NK-like`="#888888",`ISG+`="#C050C0")

lc_min<-log(min(node_tbl$total_cells)); lc_max<-log(max(node_tbl$total_cells))
seq_nodes <- node_tbl %>% dplyr::mutate(
  node_r=NODE_R_MIN+(NODE_R_MAX-NODE_R_MIN)*(log(total_cells)-lc_min)/(lc_max-lc_min),
  d1=D1,d2=D2,d3=D3,d4=D4,d5=D5,
  meta_total=Tcm+Tem+`T-CMC1`+`NK-like`+`ISG+`, seq_node=uid)
for (m in META_LEVELS){ pc<-paste0("pct_",m)
  seq_nodes[[pc]]<-ifelse(seq_nodes$meta_total>0, 100*seq_nodes[[m]]/seq_nodes$meta_total, NA_real_) }

cluster_meta <- seq_nodes %>% dplyr::group_by(cluster) %>%
  dplyr::summarise(n_seq=dplyr::n(),
                   dom_spec={w<-tapply(total_cells,spec,sum); names(w)[which.max(w)]}, .groups="drop")
knn <- read_csv(file.path(conga_dir,"..","conga_exports","tcrdist_phenotype_test",
                          "tcrdist_phenotype_specstrat_results.csv"), comment="#", show_col_types=FALSE) %>%
  dplyr::mutate(cluster=as.character(cluster)) %>% dplyr::select(cluster, knn_effect=same_cluster_effect)
seq_nodes <- seq_nodes %>% dplyr::mutate(cluster=as.character(cluster)) %>% dplyr::left_join(knn, by="cluster")

# layout: paired KK; unpaired ringed; repel
layout_one <- function(sn){
  paired<-sn %>% dplyr::filter(paired); unp<-sn %>% dplyr::filter(!paired)
  coords<-tibble(seq_node=character(),x=numeric(),y=numeric())
  if (nrow(paired)>=2){ sub<-D[paired$mrow,paired$mrow,drop=FALSE]; sub_t<-DIST_TRANSFORM(sub)
    gg<-graph_from_adjacency_matrix(sub_t,mode="undirected",weighted=TRUE,diag=FALSE)
    set.seed(1); xy<-layout_with_kk(gg,weights=E(gg)$weight); xy<-scale(xy,TRUE,FALSE)
    coords<-tibble(seq_node=paired$seq_node,x=xy[,1],y=xy[,2])
  } else if (nrow(paired)==1){ coords<-tibble(seq_node=paired$seq_node,x=0,y=0) }
  if (nrow(unp)>0){ R0<-if(nrow(coords)) max(sqrt(coords$x^2+coords$y^2))+1 else 0
    a<-seq(0,2*pi,length.out=nrow(unp)+1)[-1]
    coords<-bind_rows(coords, tibble(seq_node=unp$seq_node,x=R0*cos(a),y=R0*sin(a))) }
  coords
}
layout_cluster_local <- function(sn){
  lay<-layout_one(sn) %>% dplyr::left_join(sn,by="seq_node")
  if (nrow(lay)==1) return(list(coords=lay %>% dplyr::transmute(seq_node,x=0,y=0,node_r),bound=lay$node_r[1]))
  sp<-max(sqrt(lay$x^2+lay$y^2)); if(sp==0||!is.finite(sp)) sp<-1
  lay$x<-lay$x/sp; lay$y<-lay$y/sp
  set.seed(1); lay$x<-lay$x+rnorm(nrow(lay),0,0.01); lay$y<-lay$y+rnorm(nrow(lay),0,0.01)
  rp<-packcircles::circleRepelLayout(data.frame(x=lay$x,y=lay$y,radius=lay$node_r+GAP_MIN/2),
    xysizecols=1:3,sizetype="radius",maxiter=1000)$layout
  lay$x<-rp$x; lay$y<-rp$y
  cx<-mean(range(lay$x)); cy<-mean(range(lay$y)); lay$x<-lay$x-cx; lay$y<-lay$y-cy
  list(coords=lay %>% dplyr::select(seq_node,x,y,node_r), bound=max(sqrt(lay$x^2+lay$y^2)+lay$node_r))
}
cluster_edges <- function(sn){
  paired<-sn %>% dplyr::filter(paired); if(nrow(paired)<2) return(NULL)
  sub<-D[paired$mrow,paired$mrow,drop=FALSE]
  gg<-graph_from_adjacency_matrix(sub+1,mode="undirected",weighted=TRUE,diag=FALSE)
  m<-igraph::mst(gg); if(igraph::ecount(m)==0) return(NULL)
  el<-igraph::as_edgelist(m); i1<-as.integer(el[,1]); i2<-as.integer(el[,2])
  tibble(from=paired$seq_node[i1],to=paired$seq_node[i2],near=sub[cbind(i1,i2)]<=EDGE_CUT)
}
build_geometry <- function(keep){
  sn_all<-seq_nodes %>% dplyr::filter(cluster %in% keep)
  cl_ids<-sn_all %>% dplyr::count(cluster,name="n") %>% dplyr::arrange(dplyr::desc(n)) %>% dplyr::pull(cluster)
  locals<-lapply(cl_ids, function(cid) layout_cluster_local(sn_all %>% dplyr::filter(cluster==cid)))
  bounds<-sapply(locals,function(z)z$bound)
  circ<-packcircles::circleProgressiveLayout(bounds+CLUST_PAD,sizetype="radius")
  coords<-list(); edges_all<-list(); labels<-list()
  for (i in seq_along(cl_ids)){ cid<-cl_ids[i]; lc<-locals[[i]]$coords %>% dplyr::mutate(x=circ$x[i]+x,y=circ$y[i]+y,cluster=cid)
    coords[[i]]<-lc; labels[[i]]<-tibble(cluster=cid,x=circ$x[i],y=max(lc$y+lc$node_r)+LABEL_OFFSET)
    ed<-cluster_edges(sn_all %>% dplyr::filter(cluster==cid))
    if(!is.null(ed)){ p<-lc %>% dplyr::select(seq_node,x,y)
      edges_all[[i]]<-ed %>% dplyr::left_join(p,by=c("from"="seq_node")) %>%
        dplyr::left_join(p,by=c("to"="seq_node"),suffix=c("",".e")) %>% dplyr::rename(xend=x.e,yend=y.e) } }
  nodes<-bind_rows(coords) %>% dplyr::left_join(sn_all %>% dplyr::select(-node_r),by=c("seq_node","cluster"))
  list(nodes=nodes,edges=bind_rows(edges_all),labels=bind_rows(labels))
}
render_donor <- function(keep,tag,title){ if(!length(keep)){message(sprintf("[%s] none",tag));return(invisible())}
  geo<-build_geometry(keep); nodes<-geo$nodes; edges<-geo$edges; labels<-geo$labels
  ef<-edges %>% dplyr::filter(!near); en<-edges %>% dplyr::filter(near)
  p<-ggplot()+
    geom_segment(data=ef,aes(x,y,xend=xend,yend=yend),color=EDGE_FAR_COL,linewidth=EDGE_FAR_LW)+
    geom_segment(data=en,aes(x,y,xend=xend,yend=yend),color=EDGE_NEAR_COL,linewidth=EDGE_NEAR_LW)+
    ggforce::geom_circle(data=nodes,aes(x0=x,y0=y,r=node_r),fill="white",color="grey30",linewidth=0.2)+
    scatterpie::geom_scatterpie(data=nodes,aes(x=x,y=y,r=node_r),cols=c("d1","d2","d3","d4","d5"),color=NA)+
    geom_text(data=labels,aes(x,y,label=paste("cluster",cluster)),size=LABEL_SIZE,color="grey20")+
    scale_fill_manual(values=setNames(unname(donor_color_map),c("d1","d2","d3","d4","d5")),
                      labels=c("D1","D2","D3","D4","D5"),name="Donor")+
    coord_equal()+theme_void()+ggtitle(title)+theme(plot.title=element_text(size=14))
  ggsave(file.path(out_dir,sprintf("%s.png",tag)),p,width=15,height=13,dpi=300,bg="white")
  message(sprintf("[%s] %d clusters, %d nodes",tag,dplyr::n_distinct(nodes$cluster),nrow(nodes))) }
render_continuous <- function(keep,tag,title,value_col,scale_type,high_col="#8B0000"){
  if(!length(keep)){message(sprintf("[%s] none",tag));return(invisible())}
  geo<-build_geometry(keep); nodes<-geo$nodes; edges<-geo$edges; labels<-geo$labels
  nodes$val<-nodes[[value_col]]; ef<-edges %>% dplyr::filter(!near); en<-edges %>% dplyr::filter(near)
  p<-ggplot()+
    geom_segment(data=ef,aes(x,y,xend=xend,yend=yend),color=EDGE_FAR_COL,linewidth=EDGE_FAR_LW)+
    geom_segment(data=en,aes(x,y,xend=xend,yend=yend),color=EDGE_NEAR_COL,linewidth=EDGE_NEAR_LW)+
    ggforce::geom_circle(data=nodes,aes(x0=x,y0=y,r=node_r,fill=val),color="grey30",linewidth=0.2)
  if(scale_type=="diverging"){lim<-max(abs(nodes$val),na.rm=TRUE)
    p<-p+scale_fill_gradient2(low="#2166AC",mid="white",high="#B2182B",midpoint=0,limits=c(-lim,lim),na.value=NA_FILL,name=value_col)
  } else p<-p+scale_fill_gradient(low="white",high=high_col,na.value=NA_FILL,name=value_col)
  p<-p+geom_text(data=labels,aes(x,y,label=paste("cluster",cluster)),size=LABEL_SIZE,color="grey20")+
    coord_equal()+theme_void()+ggtitle(title)+theme(plot.title=element_text(size=14))
  ggsave(file.path(out_dir,sprintf("%s.png",tag)),p,width=15,height=13,dpi=300,bg="white")
  message(sprintf("[%s] %d clusters, %d nodes",tag,dplyr::n_distinct(nodes$cluster),nrow(nodes))) }

sets <- list(
  list(cl=cluster_meta %>% dplyr::filter(dom_spec=="pp65") %>% dplyr::pull(cluster),          cut="all",spec="pp65"),
  list(cl=cluster_meta %>% dplyr::filter(dom_spec=="IE1")  %>% dplyr::pull(cluster),          cut="all",spec="IE1"),
  list(cl=cluster_meta %>% dplyr::filter(dom_spec=="pp65",n_seq>=2) %>% dplyr::pull(cluster), cut="min2seq",spec="pp65"),
  list(cl=cluster_meta %>% dplyr::filter(dom_spec=="IE1", n_seq>=2) %>% dplyr::pull(cluster), cut="min2seq",spec="IE1"))
for (s in sets){ base<-sprintf("%s_%s",s$cut,s$spec)
  render_donor(s$cl, sprintf("%s_donor",base), sprintf("et64 v2 merged — %s (%s) — donor", s$spec, s$cut))
  render_continuous(s$cl, sprintf("%s_kNNbias",base), sprintf("et64 v2 merged — %s (%s) — kNN effect", s$spec, s$cut), "knn_effect","diverging")
  for (m in META_LEVELS){ safe<-gsub("[^A-Za-z0-9]","",m)
    render_continuous(s$cl, sprintf("%s_meta%s",base,safe), sprintf("et64 v2 merged — %s (%s) — %% %s", s$spec, s$cut, m),
                      paste0("pct_",m),"sequential",high_col=meta_palette[[m]]) } }
message("\nDone: merged overview plots in ", out_dir)




#
# Constellation v2 — logo-annotated, NODE = unique (TRB,TRA) combo
#   Layout/spacing driven by TCRdist (KK + repel); MST lines NOT drawn.
#   RED line = nodes sharing identical TRB; BLUE line = nodes sharing identical TRA.
#   Node = unique (dominant TRB, dominant TRA); clones merged on exact pair. Unpaired
#   combos kept as disconnected nodes (logo shows "not recovered" for the missing chain).
#   Logo panels have an opaque white backing so dashed leaders terminate under them.
#   2-line title (X cells: + donor/clone detail). Legend exported separately (all 5 classes).
#   Output: exports/TCR_clusters/constellation_v2/logo_annotated/
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr); library(tidyr)
library(cowplot); library(ggseqlogo); library(readr); library(grid)

conga_dir <- "exports/conga/mdata_CMV"; v2_dir <- file.path(conga_dir,"v2")
out_dir   <- "exports/TCR_clusters/constellation_v2/logo_annotated"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
TARGET_CL <- "0"
donor_color_map <- c("D1"="#1b9e77","D2"="#d95f02","D3"="#7570b3","D4"="#e7298a","D5"="#66a61e")

g <- readRDS(file.path(v2_dir,"constellation_v2_graph.rds")); D <- g$D
ord <- read_csv(file.path(v2_dir,"constellation_v2_tcrdist_order.csv"), show_col_types=FALSE)
seq_key_matrix <- paste(ord$cdr3_b_aa, ord$cdr3_a_aa)

# ── settings ─────────────────────────────────────────────────────────────────
DIST_TRANSFORM<-function(d) log1p(d)
NODE_R_MIN<-0.35; NODE_R_MAX<-1.40; GAP_MIN<-0.25
SHARE_TRB_COL<-"#C0392B"; SHARE_TRA_COL<-"#2166AC"; SHARE_LW<-1.5
RING_SCALE<-2.2; LOGO_W<-0.24; LOGO_H<-0.19
META_LEVELS<-c("Tcm","Tem","T-CMC1","NK-like","ISG+")

pull_reads <- function(strings, prefix) unlist(lapply(strings, function(s){
  if (is.na(s)) return(character(0)); parts<-str_split(s,";")[[1]]
  hit<-parts[str_detect(parts,paste0("^",prefix,":"))]; seqs<-str_replace(hit,paste0("^",prefix,":"),"")
  seqs[nchar(seqs)>0] }))
top1 <- function(strings, prefix){ r<-pull_reads(strings,prefix); if(!length(r)) return(NA_character_); names(which.max(table(r))) }

# ── 1. Per-cloneId representative, then MERGE by (trb,tra) into unique-seq nodes ─
md <- mdata_CMV@meta.data %>% dplyr::mutate(cloneId=as.character(cloneId),
        cluster=as.character(TCR_cluster_et64), global_id=as.character(Global_Clone_ID),
        donorS=paste0("D", str_extract(as.character(donor_identity),"\\d")))
cc <- md %>% dplyr::filter(cluster==TARGET_CL, !grepl("notCalled", cloneId))

per_clone <- cc %>% dplyr::group_by(cloneId) %>%
  dplyr::summarise(trb=top1(t_cdr3s_aa,"TRB"), tra=top1(t_cdr3s_aa,"TRA"),
                   is_global=any(!is.na(global_id)),
                   global_id=ifelse(any(!is.na(global_id)), global_id[!is.na(global_id)][1], NA_character_),
                   n_cells=dplyr::n(), .groups="drop") %>%
  dplyr::mutate(nkey = paste(tidyr::replace_na(trb,"NA"), tidyr::replace_na(tra,"NA"), sep="||"))

cell_nkey <- cc %>% dplyr::left_join(per_clone %>% dplyr::select(cloneId,nkey,trb,tra,is_global), by="cloneId")

nodes <- cell_nkey %>% dplyr::group_by(nkey) %>%
  dplyr::summarise(
    trb=dplyr::first(trb), tra=dplyr::first(tra),
    total_cells=dplyr::n(),
    n_major=dplyr::n_distinct(global_id[!is.na(global_id)]),
    n_minor=dplyr::n_distinct(cloneId[is.na(global_id)]),
    n_donor=dplyr::n_distinct(donorS),
    d1=sum(donorS=="D1"), d2=sum(donorS=="D2"), d3=sum(donorS=="D3"),
    d4=sum(donorS=="D4"), d5=sum(donorS=="D5"), .groups="drop") %>%
  dplyr::mutate(paired = !is.na(trb) & !is.na(tra)) %>%
  dplyr::filter(total_cells >= 5)

lc_min <- log(min(g$nodes$total_cells)); lc_max <- log(max(g$nodes$total_cells))
nodes <- nodes %>% dplyr::mutate(
  node_r = NODE_R_MIN + (NODE_R_MAX-NODE_R_MIN)*(log(total_cells)-lc_min)/(lc_max-lc_min),
  mrow = ifelse(paired, match(paste(trb, tra), seq_key_matrix), NA_integer_))
message(sprintf("cluster %s: %d unique-seq nodes (%d paired, %d unpaired)",
                TARGET_CL, nrow(nodes), sum(nodes$paired), sum(!nodes$paired)))

# ── 2. Layout (paired KK on TCRdist; unpaired ring) + repel — SPACING kept ───
paired <- nodes %>% dplyr::filter(paired); unp <- nodes %>% dplyr::filter(!paired)
coords <- tibble(nkey=character(), x=numeric(), y=numeric())
if (nrow(paired)>=2){
  sub<-D[paired$mrow,paired$mrow,drop=FALSE]; sub_t<-DIST_TRANSFORM(sub)
  gg<-graph_from_adjacency_matrix(sub_t,mode="undirected",weighted=TRUE,diag=FALSE)
  set.seed(1); xy<-layout_with_kk(gg,weights=E(gg)$weight); xy<-scale(xy,TRUE,FALSE)
  coords<-tibble(nkey=paired$nkey, x=xy[,1], y=xy[,2])
} else if (nrow(paired)==1){ coords<-tibble(nkey=paired$nkey, x=0, y=0) }
if (nrow(unp)>0){ R0<-if(nrow(coords)) max(sqrt(coords$x^2+coords$y^2))+1 else 0
  a<-seq(0,2*pi,length.out=nrow(unp)+1)[-1]
  coords<-bind_rows(coords, tibble(nkey=unp$nkey, x=R0*cos(a), y=R0*sin(a))) }
rad<-nodes$node_r[match(coords$nkey,nodes$nkey)]
rp<-packcircles::circleRepelLayout(data.frame(x=coords$x,y=coords$y,radius=rad+GAP_MIN/2),
     xysizecols=1:3, sizetype="radius", maxiter=1000)$layout
coords$x<-rp$x; coords$y<-rp$y
cx<-mean(coords$x); cy<-mean(coords$y); coords$x<-coords$x-cx; coords$y<-coords$y-cy
nodes<-nodes %>% dplyr::left_join(coords, by="nkey")

# ── 3. Chain-sharing edges: red = shared TRB, blue = shared TRA (all pairs) ──
p_all <- coords %>% dplyr::select(nkey, x, y)
share_edges <- function(chain_col, colour_tag){
  nn <- nodes %>% dplyr::filter(!is.na(.data[[chain_col]])) %>%
    dplyr::select(nkey, val = dplyr::all_of(chain_col))
  grp <- nn %>% dplyr::group_by(val) %>% dplyr::filter(dplyr::n() >= 2) %>% dplyr::group_split()
  if (length(grp) == 0) return(NULL)
  do.call(rbind, lapply(grp, function(gd){
    cb <- t(combn(gd$nkey, 2)); tibble(from=cb[,1], to=cb[,2], share=colour_tag) }))
}
edges <- dplyr::bind_rows(share_edges("trb","TRB"), share_edges("tra","TRA"))
if (nrow(edges) > 0){
  edges <- edges %>%
    dplyr::left_join(p_all, by=c("from"="nkey")) %>%
    dplyr::left_join(p_all, by=c("to"="nkey"), suffix=c("",".e")) %>%
    dplyr::rename(xend=x.e, yend=y.e)
}

# ── 4. Ring positions (angular order) ────────────────────────────────────────
const_R<-max(sqrt(nodes$x^2+nodes$y^2)+nodes$node_r)
nodes<-nodes %>% dplyr::mutate(ang=atan2(y,x)) %>% dplyr::arrange(ang)
n<-nrow(nodes); ring_ang<-seq(0,2*pi,length.out=n+1)[-1]; ring_ang<-ring_ang+(nodes$ang[1]-ring_ang[1])
ring_R<-const_R*RING_SCALE; nodes$lx<-ring_R*cos(ring_ang); nodes$ly<-ring_R*sin(ring_ang)

# ── 5. Logo panels: strict membership; 2-line title; opaque white backing ────
keep_modal<-function(s){ if(!length(s)) return(character(0))
  l<-nchar(s); m<-as.integer(names(sort(table(l),decreasing=TRUE))[1]); s[l==m] }
make_logo_panel <- function(nd){
  sub <- cell_nkey %>% dplyr::filter(nkey==nd$nkey)
  cell_line   <- sprintf("%d cells:", nd$total_cells)
  detail_line <- sprintf("%d Donor%s, %d major clone%s, %d minor clone%s",
                 nd$n_donor, ifelse(nd$n_donor==1,"","s"),
                 nd$n_major, ifelse(nd$n_major==1,"","s"),
                 nd$n_minor, ifelse(nd$n_minor==1,"","s"))
  ttl_grob <- cowplot::plot_grid(
    cowplot::ggdraw()+cowplot::draw_label(cell_line,   size=16, fontface="bold"),
    cowplot::ggdraw()+cowplot::draw_label(detail_line, size=13),
    ncol=1)
  mk <- function(has, prefix, lab){
    if (!has) return(cowplot::ggdraw()+cowplot::draw_label(paste0(lab,": not recovered"), size=9, colour="grey50"))
    seqs <- keep_modal(pull_reads(sub$t_cdr3s_aa, prefix))
    if (!length(seqs)) return(cowplot::ggdraw()+cowplot::draw_label(paste0(lab,": not recovered"), size=9, colour="grey50"))
    ggseqlogo(seqs, method="bits") + ggtitle(lab) +
      theme(legend.position="none", plot.title=element_text(size=11,face="bold"),
            axis.text=element_text(size=7), axis.title=element_blank(), plot.margin=margin(2,2,2,2))
  }
  panel <- cowplot::plot_grid(ttl_grob, mk(!is.na(nd$trb),"TRB","TRB"), mk(!is.na(nd$tra),"TRA","TRA"),
                     ncol=1, rel_heights=c(0.30, 0.35, 0.35))
  # opaque white backing so leader lines terminate under the panel, not over it
  cowplot::ggdraw() +
    cowplot::draw_grob(grid::rectGrob(gp=grid::gpar(fill="white", col=NA))) +
    cowplot::draw_plot(panel)
}
logo_panels <- lapply(seq_len(n), function(i) make_logo_panel(nodes[i,]))

# ── 6. Base donor-pie constellation (chain-sharing lines) ────────────────────
lead<-nodes %>% dplyr::transmute(x,y,lx,ly)
base<-ggplot()+
  { if(nrow(edges) > 0) list(
      geom_segment(data=edges %>% dplyr::filter(share=="TRB"),
                   aes(x,y,xend=xend,yend=yend), color=SHARE_TRB_COL, linewidth=SHARE_LW, alpha=0.85),
      geom_segment(data=edges %>% dplyr::filter(share=="TRA"),
                   aes(x,y,xend=xend,yend=yend), color=SHARE_TRA_COL, linewidth=SHARE_LW, alpha=0.85)) } +
  geom_segment(data=lead,aes(x=x,y=y,xend=lx,yend=ly),linetype="dashed",color="grey60",linewidth=0.3)+
  ggforce::geom_circle(data=nodes,aes(x0=x,y0=y,r=node_r),fill="white",color="grey30",linewidth=0.2)+
  scatterpie::geom_scatterpie(data=nodes,aes(x=x,y=y,r=node_r),cols=c("d1","d2","d3","d4","d5"),color=NA)+
  scale_fill_manual(values=setNames(unname(donor_color_map),c("d1","d2","d3","d4","d5")),
                    labels=c("D1","D2","D3","D4","D5"),name="Donor")+
  coord_equal(clip="off")+theme_void()+
  ggtitle(sprintf("Cluster %s — logo-annotated constellation", TARGET_CL))+
  theme(plot.title=element_text(size=16,hjust=0.5),legend.position="right",
        legend.title=element_text(size=12),legend.text=element_text(size=10),plot.margin=margin(8,8,8,8))

# ── 7. Compose ───────────────────────────────────────────────────────────────
xr<-range(c(nodes$lx,nodes$x))+c(-1,1)*const_R*0.7; yr<-range(c(nodes$ly,nodes$y))+c(-1,1)*const_R*0.7
to01<-function(v,r)(v-r[1])/(r[2]-r[1])
canvas<-cowplot::ggdraw(base)
for (i in seq_len(n)){ fx<-to01(nodes$lx[i],xr); fy<-to01(nodes$ly[i],yr)
  canvas<-canvas+cowplot::draw_plot(logo_panels[[i]],x=fx-LOGO_W/2,y=fy-LOGO_H/2,width=LOGO_W,height=LOGO_H) }
fn<-file.path(out_dir,sprintf("cluster_%s_logo_annotated.png",TARGET_CL))
ggsave(fn,canvas,width=20,height=20,dpi=200,bg="white")
message(sprintf("cluster %s: %d nodes -> %s", TARGET_CL, n, basename(fn)))

# ── 8. Export chemistry legend ONCE (all 5 chemistry classes) ────────────────
all_aa <- paste(c("A","R","N","D","C","Q","E","G","H","I",
                  "L","K","M","F","P","S","T","W","Y","V"), collapse="")
legend_plot <- ggseqlogo(c(all_aa, all_aa), method="bits") +
  theme(legend.position="right", legend.title=element_text(size=12), legend.text=element_text(size=11))
leg <- cowplot::get_legend(legend_plot)
ggsave(file.path(out_dir,"logo_chemistry_legend.png"),
       cowplot::ggdraw(leg), width=3, height=2.2, dpi=200, bg="white")
message("Wrote logo_chemistry_legend.png (all 5 chemistry classes)")




#
# Constellation v2 — logo-annotated plots, batch (NEW LINE LOGIC, no donor legend)
#   SOLID red  = identical TRB      | SOLID blue  = identical TRA
#   DOTTED red = TRB 1 aa different | DOTTED blue = TRA 1 aa different
#   (1 aa diff = same length, exactly one position differs; indels not counted)
#   Node = unique (TRB,TRA) combo, clones merged, >=5 cells. Unpaired combos kept
#   as nodes (logo shows "not recovered" for the missing chain).
#   Donor legend suppressed (single-cluster panels produce truncated keys).
#   Output: exports/TCR_clusters/constellation_v2/logo_annotated/
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr); library(tidyr)
library(cowplot); library(ggseqlogo); library(readr); library(grid)

conga_dir <- "exports/conga/mdata_CMV"; v2_dir <- file.path(conga_dir,"v2")
out_dir   <- "exports/TCR_clusters/constellation_v2/logo_annotated"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
donor_color_map <- c("D1"="#1b9e77","D2"="#d95f02","D3"="#7570b3","D4"="#e7298a","D5"="#66a61e")

g <- readRDS(file.path(v2_dir,"constellation_v2_graph.rds")); D <- g$D
ord <- read_csv(file.path(v2_dir,"constellation_v2_tcrdist_order.csv"), show_col_types=FALSE)
seq_key_matrix <- paste(ord$cdr3_b_aa, ord$cdr3_a_aa)

# ── settings ─────────────────────────────────────────────────────────────────
DIST_TRANSFORM<-function(d) log1p(d)
NODE_R_MIN<-0.35; NODE_R_MAX<-1.40; GAP_MIN<-0.25
SHARE_TRB_COL<-"#C0392B"; SHARE_TRA_COL<-"#2166AC"
SHARE_LW      <- 1.5       # solid (identical)
SHARE_LW_DASH <- 1.5       # dotted (1 aa difference)
DASH_PATTERN  <- "12"      # dotted
RING_SCALE<-2.2; LOGO_W<-0.24; LOGO_H<-0.19

pull_reads <- function(strings, prefix) unlist(lapply(strings, function(s){
  if (is.na(s)) return(character(0)); parts<-str_split(s,";")[[1]]
  hit<-parts[str_detect(parts,paste0("^",prefix,":"))]; seqs<-str_replace(hit,paste0("^",prefix,":"),"")
  seqs[nchar(seqs)>0] }))
top1 <- function(strings, prefix){ r<-pull_reads(strings,prefix); if(!length(r)) return(NA_character_); names(which.max(table(r))) }
keep_modal<-function(s){ if(!length(s)) return(character(0))
  l<-nchar(s); m<-as.integer(names(sort(table(l),decreasing=TRUE))[1]); s[l==m] }

one_aa_diff <- function(a, b){
  if (is.na(a) || is.na(b)) return(FALSE)
  if (nchar(a) != nchar(b)) return(FALSE)
  sum(strsplit(a,"")[[1]] != strsplit(b,"")[[1]]) == 1L
}

lc_min <- log(min(g$nodes$total_cells)); lc_max <- log(max(g$nodes$total_cells))

md_all <- mdata_CMV@meta.data %>% dplyr::mutate(cloneId=as.character(cloneId),
        cluster=as.character(TCR_cluster_et64), global_id=as.character(Global_Clone_ID),
        donorS=paste0("D", str_extract(as.character(donor_identity),"\\d")),
        spec=as.character(specificity))

# ── the per-cluster builder ──────────────────────────────────────────────────
make_logo_annotated <- function(TARGET_CL){
  cc <- md_all %>% dplyr::filter(cluster==TARGET_CL, !grepl("notCalled", cloneId))
  if (nrow(cc)==0){ message(sprintf("cluster %s: no cells, skipped", TARGET_CL)); return(invisible()) }
  dom_spec <- { w <- table(cc$spec[!is.na(cc$spec)]); if(length(w)) names(which.max(w)) else "NA" }

  per_clone <- cc %>% dplyr::group_by(cloneId) %>%
    dplyr::summarise(trb=top1(t_cdr3s_aa,"TRB"), tra=top1(t_cdr3s_aa,"TRA"),
                     global_id=ifelse(any(!is.na(global_id)), global_id[!is.na(global_id)][1], NA_character_),
                     .groups="drop") %>%
    dplyr::mutate(nkey = paste(tidyr::replace_na(trb,"NA"), tidyr::replace_na(tra,"NA"), sep="||"))
  cell_nkey <- cc %>% dplyr::left_join(per_clone %>% dplyr::select(cloneId,nkey,trb,tra), by="cloneId")

  nodes <- cell_nkey %>% dplyr::group_by(nkey) %>%
    dplyr::summarise(trb=dplyr::first(trb), tra=dplyr::first(tra), total_cells=dplyr::n(),
      n_major=dplyr::n_distinct(global_id[!is.na(global_id)]),
      n_minor=dplyr::n_distinct(cloneId[is.na(global_id)]),
      n_donor=dplyr::n_distinct(donorS),
      d1=sum(donorS=="D1"), d2=sum(donorS=="D2"), d3=sum(donorS=="D3"),
      d4=sum(donorS=="D4"), d5=sum(donorS=="D5"), .groups="drop") %>%
    dplyr::mutate(paired = !is.na(trb) & !is.na(tra)) %>%
    dplyr::filter(total_cells >= 5) %>%
    dplyr::mutate(node_r = NODE_R_MIN + (NODE_R_MAX-NODE_R_MIN)*(log(total_cells)-lc_min)/(lc_max-lc_min),
                  mrow = ifelse(paired, match(paste(trb,tra), seq_key_matrix), NA_integer_))
  if (nrow(nodes) < 1){ message(sprintf("cluster %s: no nodes >=5 cells, skipped", TARGET_CL)); return(invisible()) }

  # layout (spacing driven by TCRdist)
  paired <- nodes %>% dplyr::filter(paired); unp <- nodes %>% dplyr::filter(!paired)
  coords <- tibble(nkey=character(), x=numeric(), y=numeric())
  if (nrow(paired)>=2){
    sub<-D[paired$mrow,paired$mrow,drop=FALSE]; sub_t<-DIST_TRANSFORM(sub)
    gg<-graph_from_adjacency_matrix(sub_t,mode="undirected",weighted=TRUE,diag=FALSE)
    set.seed(1); xy<-layout_with_kk(gg,weights=E(gg)$weight); xy<-scale(xy,TRUE,FALSE)
    coords<-tibble(nkey=paired$nkey, x=xy[,1], y=xy[,2])
  } else if (nrow(paired)==1){ coords<-tibble(nkey=paired$nkey, x=0, y=0) }
  if (nrow(unp)>0){ R0<-if(nrow(coords)) max(sqrt(coords$x^2+coords$y^2))+1 else 0
    a<-seq(0,2*pi,length.out=nrow(unp)+1)[-1]
    coords<-bind_rows(coords, tibble(nkey=unp$nkey, x=R0*cos(a), y=R0*sin(a))) }
  rad<-nodes$node_r[match(coords$nkey,nodes$nkey)]
  rp<-packcircles::circleRepelLayout(data.frame(x=coords$x,y=coords$y,radius=rad+GAP_MIN/2),
       xysizecols=1:3, sizetype="radius", maxiter=1000)$layout
  coords$x<-rp$x; coords$y<-rp$y
  cx<-mean(coords$x); cy<-mean(coords$y); coords$x<-coords$x-cx; coords$y<-coords$y-cy
  nodes<-nodes %>% dplyr::left_join(coords, by="nkey")

  # ── chain-sharing edges: exact (solid) + one-aa (dotted) ──────────────────
  edges <- tibble(x=numeric(),y=numeric(),xend=numeric(),yend=numeric(),
                  share=character(),kind=character())
  if (nrow(nodes) >= 2){
    cb <- t(combn(nrow(nodes), 2))
    edges <- bind_rows(lapply(seq_len(nrow(cb)), function(k){
      i <- cb[k,1]; j <- cb[k,2]; rows <- list()
      tb_i <- nodes$trb[i]; tb_j <- nodes$trb[j]
      if (!is.na(tb_i) && !is.na(tb_j)){
        if (identical(tb_i, tb_j))       rows[[length(rows)+1]] <- tibble(share="TRB", kind="exact")
        else if (one_aa_diff(tb_i,tb_j)) rows[[length(rows)+1]] <- tibble(share="TRB", kind="one_aa") }
      ta_i <- nodes$tra[i]; ta_j <- nodes$tra[j]
      if (!is.na(ta_i) && !is.na(ta_j)){
        if (identical(ta_i, ta_j))       rows[[length(rows)+1]] <- tibble(share="TRA", kind="exact")
        else if (one_aa_diff(ta_i,ta_j)) rows[[length(rows)+1]] <- tibble(share="TRA", kind="one_aa") }
      if (!length(rows)) return(NULL)
      bind_rows(rows) %>% dplyr::mutate(x=nodes$x[i], y=nodes$y[i],
                                        xend=nodes$x[j], yend=nodes$y[j]) }))
  }

  # ring positions (angular order)
  const_R<-max(sqrt(nodes$x^2+nodes$y^2)+nodes$node_r)
  nodes<-nodes %>% dplyr::mutate(ang=atan2(y,x)) %>% dplyr::arrange(ang)
  n<-nrow(nodes); ring_ang<-seq(0,2*pi,length.out=n+1)[-1]; ring_ang<-ring_ang+(nodes$ang[1]-ring_ang[1])
  ring_R<-const_R*RING_SCALE; nodes$lx<-ring_R*cos(ring_ang); nodes$ly<-ring_R*sin(ring_ang)

  # logo panels
  make_panel <- function(nd){
    sub <- cell_nkey %>% dplyr::filter(nkey==nd$nkey)
    cell_line <- sprintf("%d cells:", nd$total_cells)
    detail_line <- sprintf("%d Donor%s, %d major clone%s, %d minor clone%s",
                   nd$n_donor,ifelse(nd$n_donor==1,"","s"), nd$n_major,ifelse(nd$n_major==1,"","s"),
                   nd$n_minor,ifelse(nd$n_minor==1,"","s"))
    ttl_grob <- cowplot::plot_grid(
      cowplot::ggdraw()+cowplot::draw_label(cell_line,size=16,fontface="bold"),
      cowplot::ggdraw()+cowplot::draw_label(detail_line,size=13), ncol=1)
    mk <- function(has,prefix,lab){
      if(!has) return(cowplot::ggdraw()+cowplot::draw_label(paste0(lab,": not recovered"),size=9,colour="grey50"))
      seqs<-keep_modal(pull_reads(sub$t_cdr3s_aa,prefix))
      if(!length(seqs)) return(cowplot::ggdraw()+cowplot::draw_label(paste0(lab,": not recovered"),size=9,colour="grey50"))
      ggseqlogo(seqs,method="bits")+ggtitle(lab)+
        theme(legend.position="none",plot.title=element_text(size=11,face="bold"),
              axis.text=element_text(size=7),axis.title=element_blank(),plot.margin=margin(2,2,2,2)) }
    panel<-cowplot::plot_grid(ttl_grob, mk(!is.na(nd$trb),"TRB","TRB"), mk(!is.na(nd$tra),"TRA","TRA"),
                              ncol=1, rel_heights=c(0.30,0.35,0.35))
    cowplot::ggdraw()+cowplot::draw_grob(grid::rectGrob(gp=grid::gpar(fill="white",col=NA)))+cowplot::draw_plot(panel)
  }
  logo_panels <- lapply(seq_len(n), function(i) make_panel(nodes[i,]))

  # base + compose (NO donor legend)
  lead<-nodes %>% dplyr::transmute(x,y,lx,ly)
  base<-ggplot()+
    { if(nrow(edges)>0) list(
        geom_segment(data=edges %>% dplyr::filter(share=="TRB", kind=="exact"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRB_COL, linewidth=SHARE_LW, alpha=0.85),
        geom_segment(data=edges %>% dplyr::filter(share=="TRA", kind=="exact"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRA_COL, linewidth=SHARE_LW, alpha=0.85),
        geom_segment(data=edges %>% dplyr::filter(share=="TRB", kind=="one_aa"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRB_COL,
                     linewidth=SHARE_LW_DASH, linetype=DASH_PATTERN, alpha=0.8),
        geom_segment(data=edges %>% dplyr::filter(share=="TRA", kind=="one_aa"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRA_COL,
                     linewidth=SHARE_LW_DASH, linetype=DASH_PATTERN, alpha=0.8)) }+
    geom_segment(data=lead,aes(x=x,y=y,xend=lx,yend=ly),linetype="dashed",color="grey60",linewidth=0.3)+
    ggforce::geom_circle(data=nodes,aes(x0=x,y0=y,r=node_r),fill="white",color="grey30",linewidth=0.2)+
    scatterpie::geom_scatterpie(data=nodes,aes(x=x,y=y,r=node_r),cols=c("d1","d2","d3","d4","d5"),color=NA)+
    scale_fill_manual(values=setNames(unname(donor_color_map),c("d1","d2","d3","d4","d5")),
                      labels=c("D1","D2","D3","D4","D5"),name="Donor", guide="none")+
    coord_equal(clip="off")+theme_void()+
    ggtitle(sprintf("Cluster %s — logo-annotated constellation", TARGET_CL))+
    theme(plot.title=element_text(size=16,hjust=0.5), plot.margin=margin(8,8,8,8))

  xr<-range(c(nodes$lx,nodes$x))+c(-1,1)*const_R*0.7; yr<-range(c(nodes$ly,nodes$y))+c(-1,1)*const_R*0.7
  to01<-function(v,r)(v-r[1])/(r[2]-r[1])
  canvas<-cowplot::ggdraw(base)
  for (i in seq_len(n)){ fx<-to01(nodes$lx[i],xr); fy<-to01(nodes$ly[i],yr)
    canvas<-canvas+cowplot::draw_plot(logo_panels[[i]],x=fx-LOGO_W/2,y=fy-LOGO_H/2,width=LOGO_W,height=LOGO_H) }
  fn<-file.path(out_dir,sprintf("%s_cluster_%s_logo_annotated.png", dom_spec, TARGET_CL))
  ggsave(fn,canvas,width=20,height=20,dpi=200,bg="white")
  message(sprintf("[%s] cluster %s: %d nodes, %d edges (%d exact / %d 1aa) -> %s",
                  dom_spec, TARGET_CL, n, nrow(edges),
                  sum(edges$kind=="exact"), sum(edges$kind=="one_aa"), basename(fn)))
}

# ── run for all annotated clusters ───────────────────────────────────────────
for (cl in c("0","2","4","8","18","14","17")) make_logo_annotated(cl)

# cluster 2 only — logos at 50% size (denser cluster)
LOGO_W_saved <- LOGO_W; LOGO_H_saved <- LOGO_H
LOGO_W <- 0.12; LOGO_H <- 0.095
make_logo_annotated("2")
LOGO_W <- LOGO_W_saved; LOGO_H <- LOGO_H_saved

# ── legend once (all 5 chemistry classes) ────────────────────────────────────
all_aa <- paste(c("A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"), collapse="")
leg <- cowplot::get_legend(ggseqlogo(c(all_aa,all_aa),method="bits") +
  theme(legend.position="right",legend.title=element_text(size=12),legend.text=element_text(size=11)))
ggsave(file.path(out_dir,"logo_chemistry_legend.png"), cowplot::ggdraw(leg), width=3, height=2.2, dpi=200, bg="white")
message("Wrote logo_chemistry_legend.png")




#
# Manuscript figure: pp65 constellations (left grid) | dashed divider | IE1 (right column)
#   Each cluster = its own tight panel (merged unique-seq nodes, donor pies).
#   Layout: pp65 row1=[0], row2=[2,4], row3=[8,18]; IE1 column=[17,14]. Big titles.
#   NODE SIZING CONSISTENT ACROSS PANELS: all panels share one coordinate limit
#   (GLOBAL_LIM), so equal on-page size = equal cell count everywhere.
#   Output: exports/TCR_clusters/constellation_v2/manuscript_donor_figure.png
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr); library(cowplot)

conga_dir <- "exports/conga/mdata_CMV"; v2_dir <- file.path(conga_dir,"v2")
out_dir   <- "exports/TCR_clusters/constellation_v2"
g <- readRDS(file.path(v2_dir,"constellation_v2_graph_merged.rds"))
node_tbl <- g$nodes; D <- g$D
donor_color_map <- c("D1"="#1b9e77","D2"="#d95f02","D3"="#7570b3","D4"="#e7298a","D5"="#66a61e")

DIST_TRANSFORM<-function(d) log1p(d)
NODE_R_MIN<-0.35; NODE_R_MAX<-1.40; GAP_MIN<-0.25
EDGE_CUT<-250; EDGE_FAR_COL<-"grey80"; EDGE_FAR_LW<-0.3; EDGE_NEAR_COL<-"grey35"; EDGE_NEAR_LW<-0.8
TITLE_SIZE<-22
lc_min<-log(min(node_tbl$total_cells)); lc_max<-log(max(node_tbl$total_cells))

# ── compute ONE cluster's node layout (shared by extent + render) ────────────
layout_cluster <- function(cl){
  sn <- node_tbl %>% dplyr::filter(as.character(cluster)==cl) %>%
    dplyr::mutate(node_r=NODE_R_MIN+(NODE_R_MAX-NODE_R_MIN)*(log(total_cells)-lc_min)/(lc_max-lc_min),
                  d1=D1,d2=D2,d3=D3,d4=D4,d5=D5, seq_node=uid)
  paired<-sn %>% dplyr::filter(paired); unp<-sn %>% dplyr::filter(!paired)
  coords<-tibble(seq_node=character(),x=numeric(),y=numeric())
  if (nrow(paired)>=2){ sub<-D[paired$mrow,paired$mrow,drop=FALSE]; sub_t<-DIST_TRANSFORM(sub)
    gg<-graph_from_adjacency_matrix(sub_t,mode="undirected",weighted=TRUE,diag=FALSE)
    set.seed(1); xy<-layout_with_kk(gg,weights=E(gg)$weight); xy<-scale(xy,TRUE,FALSE)
    coords<-tibble(seq_node=paired$seq_node,x=xy[,1],y=xy[,2])
  } else if (nrow(paired)==1){ coords<-tibble(seq_node=paired$seq_node,x=0,y=0) }
  if (nrow(unp)>0){ R0<-if(nrow(coords)) max(sqrt(coords$x^2+coords$y^2))+1 else 0
    a<-seq(0,2*pi,length.out=nrow(unp)+1)[-1]
    coords<-bind_rows(coords, tibble(seq_node=unp$seq_node,x=R0*cos(a),y=R0*sin(a))) }
  rad<-sn$node_r[match(coords$seq_node,sn$seq_node)]
  rp<-packcircles::circleRepelLayout(data.frame(x=coords$x,y=coords$y,radius=rad+GAP_MIN/2),
       xysizecols=1:3,sizetype="radius",maxiter=1000)$layout
  coords$x<-rp$x; coords$y<-rp$y
  cx<-mean(coords$x); cy<-mean(coords$y); coords$x<-coords$x-cx; coords$y<-coords$y-cy
  sn<-sn %>% dplyr::left_join(coords,by="seq_node")
  edges<-NULL
  if (nrow(paired)>=2){ sub<-D[paired$mrow,paired$mrow,drop=FALSE]
    gg<-graph_from_adjacency_matrix(sub+1,mode="undirected",weighted=TRUE,diag=FALSE)
    m<-igraph::mst(gg); el<-igraph::as_edgelist(m); i1<-as.integer(el[,1]); i2<-as.integer(el[,2])
    p<-coords %>% dplyr::select(seq_node,x,y)
    edges<-tibble(from=paired$seq_node[i1],to=paired$seq_node[i2],near=sub[cbind(i1,i2)]<=EDGE_CUT) %>%
      dplyr::left_join(p,by=c("from"="seq_node")) %>%
      dplyr::left_join(p,by=c("to"="seq_node"),suffix=c("",".e")) %>% dplyr::rename(xend=x.e,yend=y.e) }
  list(nodes=sn, edges=edges)
}

# ── global coordinate extent (largest cluster) so all panels share scale ─────
all_cl <- c("0","2","4","8","18","17","14")
GLOBAL_LIM <- max(sapply(all_cl, function(cl){
  lay <- layout_cluster(cl)$nodes
  max(sqrt(lay$x^2 + lay$y^2) + lay$node_r) })) * 1.05

# ── render ONE cluster panel with the SHARED coordinate limit ────────────────
render_panel <- function(cl){
  lay <- layout_cluster(cl); sn <- lay$nodes; edges <- lay$edges
  ggplot()+
    { if(!is.null(edges)) list(
        geom_segment(data=edges %>% dplyr::filter(!near),aes(x,y,xend=xend,yend=yend),color=EDGE_FAR_COL,linewidth=EDGE_FAR_LW),
        geom_segment(data=edges %>% dplyr::filter(near), aes(x,y,xend=xend,yend=yend),color=EDGE_NEAR_COL,linewidth=EDGE_NEAR_LW)) }+
    ggforce::geom_circle(data=sn,aes(x0=x,y0=y,r=node_r),fill="white",color="grey30",linewidth=0.2)+
    scatterpie::geom_scatterpie(data=sn,aes(x=x,y=y,r=node_r),cols=c("d1","d2","d3","d4","d5"),color=NA)+
    scale_fill_manual(values=setNames(unname(donor_color_map),c("d1","d2","d3","d4","d5")),
                      labels=c("D1","D2","D3","D4","D5"),name="Donor",guide="none")+
    coord_equal(xlim=c(-GLOBAL_LIM,GLOBAL_LIM), ylim=c(-GLOBAL_LIM,GLOBAL_LIM), clip="off")+
    theme_void()+ ggtitle(sprintf("cluster %s", cl))+
    theme(plot.title=element_text(size=TITLE_SIZE,hjust=0.5,margin=margin(b=2)), plot.margin=margin(2,2,2,2))
}

p0 <-render_panel("0"); p2 <-render_panel("2"); p4 <-render_panel("4")
p8 <-render_panel("8"); p18<-render_panel("18"); p17<-render_panel("17"); p14<-render_panel("14")

# ── compose ──────────────────────────────────────────────────────────────────
pp65_grid <- cowplot::plot_grid(
  p0,
  cowplot::plot_grid(p2, p4, ncol=2),
  cowplot::plot_grid(p8, p18, ncol=2),
  ncol=1, rel_heights=c(1,1,1))
ie1_col <- cowplot::plot_grid(p17, p14, ncol=1, rel_heights=c(1,1))

# donor legend built from a dummy that contains ALL five donors, fixed mapping
leg_df <- data.frame(x = 1:5, y = 1,
                     donor = factor(c("D1","D2","D3","D4","D5"),
                                    levels = c("D1","D2","D3","D4","D5")))
leg_src <- ggplot(leg_df, aes(x, y, fill = donor)) +
  geom_tile() +
  scale_fill_manual(values = donor_color_map, name = "Donor",
                    breaks = c("D1","D2","D3","D4","D5"), drop = FALSE) +
  theme_void() + theme(legend.position = "right",
                       legend.title = element_text(size = 12),
                       legend.text  = element_text(size = 10))
legend <- cowplot::get_legend(leg_src)

divider <- cowplot::ggdraw() +
  cowplot::draw_line(x=c(0.5,0.5), y=c(0.02,0.98), linetype="dashed", color="grey40", linewidth=0.6)

fig <- cowplot::plot_grid(pp65_grid, divider, ie1_col, legend,
                          ncol=4, rel_widths=c(1.0, 0.04, 0.42, 0.14))
fn <- file.path(out_dir, "manuscript_donor_figure.png")
ggsave(fn, fig, width=20, height=16, dpi=300, bg="white")
message(sprintf("Wrote %s (GLOBAL_LIM=%.2f)", basename(fn), GLOBAL_LIM))




#
# Manuscript figure v2: chain-sharing lines instead of MST
#   SOLID red  = nodes with IDENTICAL TRB     | SOLID blue  = identical TRA
#   DASHED red = TRB differing by exactly 1aa | DASHED blue = TRA differing by 1aa
#   (1aa diff = same length, exactly one position differs; indels not counted)
#   Layout: pp65 row1=[0], row2=[2,4], row3=[8,18] | dashed divider | IE1 [17,14]
#   Shared coordinate limit so equal on-page size = equal cell count everywhere.
#   Output: exports/TCR_clusters/constellation_v2/manuscript_donor_figure_v2.png
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr); library(cowplot)

conga_dir <- "exports/conga/mdata_CMV"; v2_dir <- file.path(conga_dir,"v2")
out_dir   <- "exports/TCR_clusters/constellation_v2"
g <- readRDS(file.path(v2_dir,"constellation_v2_graph_merged.rds"))
node_tbl <- g$nodes; D <- g$D
donor_color_map <- c("D1"="#1b9e77","D2"="#d95f02","D3"="#7570b3","D4"="#e7298a","D5"="#66a61e")

DIST_TRANSFORM<-function(d) log1p(d)
NODE_R_MIN<-0.35; NODE_R_MAX<-1.40; GAP_MIN<-0.25
SHARE_TRB_COL<-"#C0392B"; SHARE_TRA_COL<-"#2166AC"
SHARE_LW      <- 1.5     # solid (identical) lines
SHARE_LW_DASH <- 1.5     # dashed (1aa-diff) lines
DASH_PATTERN  <- "12"    # dotted: 1 unit on, 2 off
TITLE_SIZE<-22
lc_min<-log(min(node_tbl$total_cells)); lc_max<-log(max(node_tbl$total_cells))

# exactly one amino-acid substitution (same length, one position differs)
one_aa_diff <- function(a, b){
  if (is.na(a) || is.na(b)) return(FALSE)
  if (nchar(a) != nchar(b)) return(FALSE)
  sum(strsplit(a,"")[[1]] != strsplit(b,"")[[1]]) == 1L
}

# ── layout ONE cluster (nodes + coords) ──────────────────────────────────────
layout_cluster <- function(cl){
  sn <- node_tbl %>% dplyr::filter(as.character(cluster)==cl) %>%
    dplyr::mutate(node_r=NODE_R_MIN+(NODE_R_MAX-NODE_R_MIN)*(log(total_cells)-lc_min)/(lc_max-lc_min),
                  d1=D1,d2=D2,d3=D3,d4=D4,d5=D5, seq_node=uid)
  paired<-sn %>% dplyr::filter(paired); unp<-sn %>% dplyr::filter(!paired)
  coords<-tibble(seq_node=character(),x=numeric(),y=numeric())
  if (nrow(paired)>=2){ sub<-D[paired$mrow,paired$mrow,drop=FALSE]; sub_t<-DIST_TRANSFORM(sub)
    gg<-graph_from_adjacency_matrix(sub_t,mode="undirected",weighted=TRUE,diag=FALSE)
    set.seed(1); xy<-layout_with_kk(gg,weights=E(gg)$weight); xy<-scale(xy,TRUE,FALSE)
    coords<-tibble(seq_node=paired$seq_node,x=xy[,1],y=xy[,2])
  } else if (nrow(paired)==1){ coords<-tibble(seq_node=paired$seq_node,x=0,y=0) }
  if (nrow(unp)>0){ R0<-if(nrow(coords)) max(sqrt(coords$x^2+coords$y^2))+1 else 0
    a<-seq(0,2*pi,length.out=nrow(unp)+1)[-1]
    coords<-bind_rows(coords, tibble(seq_node=unp$seq_node,x=R0*cos(a),y=R0*sin(a))) }
  rad<-sn$node_r[match(coords$seq_node,sn$seq_node)]
  rp<-packcircles::circleRepelLayout(data.frame(x=coords$x,y=coords$y,radius=rad+GAP_MIN/2),
       xysizecols=1:3,sizetype="radius",maxiter=1000)$layout
  coords$x<-rp$x; coords$y<-rp$y
  cx<-mean(coords$x); cy<-mean(coords$y); coords$x<-coords$x-cx; coords$y<-coords$y-cy
  sn %>% dplyr::left_join(coords,by="seq_node")
}

# ── chain-sharing edges for one cluster's laid-out nodes ─────────────────────
chain_edges <- function(sn){
  if (nrow(sn) < 2) return(tibble(x=numeric(),y=numeric(),xend=numeric(),yend=numeric(),
                                  share=character(),kind=character()))
  cb <- t(combn(nrow(sn), 2))
  out <- lapply(seq_len(nrow(cb)), function(k){
    i <- cb[k,1]; j <- cb[k,2]; rows <- list()
    tb_i <- sn$trb[i]; tb_j <- sn$trb[j]
    if (!is.na(tb_i) && !is.na(tb_j)){
      if (identical(tb_i, tb_j))       rows[[length(rows)+1]] <- tibble(share="TRB", kind="exact")
      else if (one_aa_diff(tb_i,tb_j)) rows[[length(rows)+1]] <- tibble(share="TRB", kind="one_aa")
    }
    ta_i <- sn$tra[i]; ta_j <- sn$tra[j]
    if (!is.na(ta_i) && !is.na(ta_j)){
      if (identical(ta_i, ta_j))       rows[[length(rows)+1]] <- tibble(share="TRA", kind="exact")
      else if (one_aa_diff(ta_i,ta_j)) rows[[length(rows)+1]] <- tibble(share="TRA", kind="one_aa")
    }
    if (!length(rows)) return(NULL)
    bind_rows(rows) %>% dplyr::mutate(x=sn$x[i], y=sn$y[i], xend=sn$x[j], yend=sn$y[j])
  })
  bind_rows(out)
}

# ── global coordinate extent so all panels share scale ──────────────────────
all_cl <- c("0","2","4","8","18","17","14")
GLOBAL_LIM <- max(sapply(all_cl, function(cl){
  lay <- layout_cluster(cl); max(sqrt(lay$x^2 + lay$y^2) + lay$node_r) })) * 1.05

# ── render ONE panel ─────────────────────────────────────────────────────────
render_panel <- function(cl){
  sn <- layout_cluster(cl); ed <- chain_edges(sn)
  ggplot()+
    { if(nrow(ed)>0) list(
        geom_segment(data=ed %>% dplyr::filter(share=="TRB", kind=="exact"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRB_COL,
                     linewidth=SHARE_LW, alpha=0.85),
        geom_segment(data=ed %>% dplyr::filter(share=="TRA", kind=="exact"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRA_COL,
                     linewidth=SHARE_LW, alpha=0.85),
        geom_segment(data=ed %>% dplyr::filter(share=="TRB", kind=="one_aa"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRB_COL,
                     linewidth=SHARE_LW_DASH, linetype=DASH_PATTERN, alpha=0.8),
        geom_segment(data=ed %>% dplyr::filter(share=="TRA", kind=="one_aa"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRA_COL,
                     linewidth=SHARE_LW_DASH, linetype=DASH_PATTERN, alpha=0.8)) }+
    ggforce::geom_circle(data=sn,aes(x0=x,y0=y,r=node_r),fill="white",color="grey30",linewidth=0.2)+
    scatterpie::geom_scatterpie(data=sn,aes(x=x,y=y,r=node_r),cols=c("d1","d2","d3","d4","d5"),color=NA)+
    scale_fill_manual(values=setNames(unname(donor_color_map),c("d1","d2","d3","d4","d5")),
                      labels=c("D1","D2","D3","D4","D5"),name="Donor",guide="none")+
    coord_equal(xlim=c(-GLOBAL_LIM,GLOBAL_LIM), ylim=c(-GLOBAL_LIM,GLOBAL_LIM), clip="off")+
    theme_void()+ ggtitle(sprintf("cluster %s", cl))+
    theme(plot.title=element_text(size=TITLE_SIZE,hjust=0.5,margin=margin(b=2)), plot.margin=margin(2,2,2,2))
}

p0 <-render_panel("0"); p2 <-render_panel("2"); p4 <-render_panel("4")
p8 <-render_panel("8"); p18<-render_panel("18"); p17<-render_panel("17"); p14<-render_panel("14")

# ── compose ──────────────────────────────────────────────────────────────────
pp65_grid <- cowplot::plot_grid(
  p0,
  cowplot::plot_grid(p2, p4, ncol=2),
  cowplot::plot_grid(p8, p18, ncol=2),
  ncol=1, rel_heights=c(1,1,1))
ie1_col <- cowplot::plot_grid(p17, p14, ncol=1, rel_heights=c(1,1))

# donor legend from a dummy containing ALL five donors (fixed, correct mapping)
leg_df <- data.frame(x=1:5, y=1,
                     donor=factor(c("D1","D2","D3","D4","D5"), levels=c("D1","D2","D3","D4","D5")))
donor_leg <- cowplot::get_legend(
  ggplot(leg_df, aes(x,y,fill=donor)) + geom_tile() +
    scale_fill_manual(values=donor_color_map, name="Donor",
                      breaks=c("D1","D2","D3","D4","D5"), drop=FALSE) +
    theme_void() + theme(legend.position="right",
                         legend.title=element_text(size=12), legend.text=element_text(size=10)))

# line-type legend (what red/blue and solid/dashed mean)
line_df <- data.frame(
  x=1, xend=2, y=4:1,
  lab=factor(c("Identical TRB","TRB: 1 aa difference","Identical TRA","TRA: 1 aa difference"),
             levels=c("Identical TRB","TRB: 1 aa difference","Identical TRA","TRA: 1 aa difference")))
line_leg <- cowplot::get_legend(
  ggplot(line_df) +
    geom_segment(aes(x=x,xend=xend,y=y,yend=y,color=lab,linetype=lab), linewidth=1.5) +
    scale_color_manual(values=c("Identical TRB"=SHARE_TRB_COL,"TRB: 1 aa difference"=SHARE_TRB_COL,
                                "Identical TRA"=SHARE_TRA_COL,"TRA: 1 aa difference"=SHARE_TRA_COL),
                       name="Shared chain") +
    scale_linetype_manual(values=c("Identical TRB"="solid","TRB: 1 aa difference"=DASH_PATTERN,
                                   "Identical TRA"="solid","TRA: 1 aa difference"=DASH_PATTERN),
                          name="Shared chain") +
    theme_void() + theme(legend.position="right",
                         legend.title=element_text(size=12), legend.text=element_text(size=10)))

legends <- cowplot::plot_grid(donor_leg, line_leg, ncol=1, rel_heights=c(1,1))

divider <- cowplot::ggdraw() +
  cowplot::draw_line(x=c(0.5,0.5), y=c(0.02,0.98), linetype="dashed", color="grey40", linewidth=0.6)

fig <- cowplot::plot_grid(pp65_grid, divider, ie1_col, legends,
                          ncol=4, rel_widths=c(1.0, 0.04, 0.42, 0.18))
fn <- file.path(out_dir, "manuscript_donor_figure_v2.png")
ggsave(fn, fig, width=20, height=16, dpi=300, bg="white")
message(sprintf("Wrote %s", basename(fn)))



#
# Manuscript figure — kNN-bias fill version
#   Same layout/lines as v2, but node fill = the cluster's kNN same-cluster effect
#   (white -> red), with limits computed GLOBALLY across all displayed clusters
#   (both pp65 and IE1) so colors are comparable everywhere.
#   SOLID red/blue = identical TRB/TRA; DOTTED = 1 aa difference.
#   Output: exports/TCR_clusters/constellation_v2/manuscript_donor_figure_kNN.png
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(ggforce); library(stringr); library(cowplot); library(readr)

conga_dir <- "exports/conga/mdata_CMV"; v2_dir <- file.path(conga_dir,"v2")
out_dir   <- "exports/TCR_clusters/constellation_v2"
g <- readRDS(file.path(v2_dir,"constellation_v2_graph_merged.rds"))
node_tbl <- g$nodes; D <- g$D

DIST_TRANSFORM<-function(d) log1p(d)
NODE_R_MIN<-0.35; NODE_R_MAX<-1.40; GAP_MIN<-0.25
SHARE_TRB_COL<-"#C0392B"; SHARE_TRA_COL<-"#2166AC"
SHARE_LW      <- 1.5
SHARE_LW_DASH <- 1.5
DASH_PATTERN  <- "12"          # dotted
KNN_LOW<-"white"; KNN_HIGH<-"#B2182B"; KNN_NA<-"grey85"
TITLE_SIZE<-22
lc_min<-log(min(node_tbl$total_cells)); lc_max<-log(max(node_tbl$total_cells))

# ── attach each cluster's kNN same-cluster effect to its nodes ───────────────
knn <- read_csv(file.path(conga_dir,"..","conga_exports","tcrdist_phenotype_test",
                          "tcrdist_phenotype_specstrat_results.csv"),
                comment="#", show_col_types=FALSE) %>%
  dplyr::mutate(cluster=as.character(cluster)) %>%
  dplyr::select(cluster, knn_effect = same_cluster_effect)
node_tbl <- node_tbl %>% dplyr::mutate(cluster=as.character(cluster)) %>%
  dplyr::left_join(knn, by="cluster")

one_aa_diff <- function(a, b){
  if (is.na(a) || is.na(b)) return(FALSE)
  if (nchar(a) != nchar(b)) return(FALSE)
  sum(strsplit(a,"")[[1]] != strsplit(b,"")[[1]]) == 1L
}

# ── layout ONE cluster ───────────────────────────────────────────────────────
layout_cluster <- function(cl){
  sn <- node_tbl %>% dplyr::filter(cluster==cl) %>%
    dplyr::mutate(node_r=NODE_R_MIN+(NODE_R_MAX-NODE_R_MIN)*(log(total_cells)-lc_min)/(lc_max-lc_min),
                  seq_node=uid)
  paired<-sn %>% dplyr::filter(paired); unp<-sn %>% dplyr::filter(!paired)
  coords<-tibble(seq_node=character(),x=numeric(),y=numeric())
  if (nrow(paired)>=2){ sub<-D[paired$mrow,paired$mrow,drop=FALSE]; sub_t<-DIST_TRANSFORM(sub)
    gg<-graph_from_adjacency_matrix(sub_t,mode="undirected",weighted=TRUE,diag=FALSE)
    set.seed(1); xy<-layout_with_kk(gg,weights=E(gg)$weight); xy<-scale(xy,TRUE,FALSE)
    coords<-tibble(seq_node=paired$seq_node,x=xy[,1],y=xy[,2])
  } else if (nrow(paired)==1){ coords<-tibble(seq_node=paired$seq_node,x=0,y=0) }
  if (nrow(unp)>0){ R0<-if(nrow(coords)) max(sqrt(coords$x^2+coords$y^2))+1 else 0
    a<-seq(0,2*pi,length.out=nrow(unp)+1)[-1]
    coords<-bind_rows(coords, tibble(seq_node=unp$seq_node,x=R0*cos(a),y=R0*sin(a))) }
  rad<-sn$node_r[match(coords$seq_node,sn$seq_node)]
  rp<-packcircles::circleRepelLayout(data.frame(x=coords$x,y=coords$y,radius=rad+GAP_MIN/2),
       xysizecols=1:3,sizetype="radius",maxiter=1000)$layout
  coords$x<-rp$x; coords$y<-rp$y
  cx<-mean(coords$x); cy<-mean(coords$y); coords$x<-coords$x-cx; coords$y<-coords$y-cy
  sn %>% dplyr::left_join(coords,by="seq_node")
}

# ── chain-sharing edges ──────────────────────────────────────────────────────
chain_edges <- function(sn){
  if (nrow(sn) < 2) return(tibble(x=numeric(),y=numeric(),xend=numeric(),yend=numeric(),
                                  share=character(),kind=character()))
  cb <- t(combn(nrow(sn), 2))
  bind_rows(lapply(seq_len(nrow(cb)), function(k){
    i <- cb[k,1]; j <- cb[k,2]; rows <- list()
    tb_i <- sn$trb[i]; tb_j <- sn$trb[j]
    if (!is.na(tb_i) && !is.na(tb_j)){
      if (identical(tb_i, tb_j))       rows[[length(rows)+1]] <- tibble(share="TRB", kind="exact")
      else if (one_aa_diff(tb_i,tb_j)) rows[[length(rows)+1]] <- tibble(share="TRB", kind="one_aa") }
    ta_i <- sn$tra[i]; ta_j <- sn$tra[j]
    if (!is.na(ta_i) && !is.na(ta_j)){
      if (identical(ta_i, ta_j))       rows[[length(rows)+1]] <- tibble(share="TRA", kind="exact")
      else if (one_aa_diff(ta_i,ta_j)) rows[[length(rows)+1]] <- tibble(share="TRA", kind="one_aa") }
    if (!length(rows)) return(NULL)
    bind_rows(rows) %>% dplyr::mutate(x=sn$x[i], y=sn$y[i], xend=sn$x[j], yend=sn$y[j]) }))
}

# ── global scales: coordinate extent + kNN fill limits ──────────────────────
all_cl <- c("0","2","4","8","18","17","14")
GLOBAL_LIM <- max(sapply(all_cl, function(cl){
  lay <- layout_cluster(cl); max(sqrt(lay$x^2 + lay$y^2) + lay$node_r) })) * 1.05

knn_vals <- node_tbl %>% dplyr::filter(cluster %in% all_cl) %>% dplyr::pull(knn_effect)
KNN_LIM  <- range(knn_vals, na.rm=TRUE)
message(sprintf("kNN fill limits (global, both specificities): %.2f to %.2f", KNN_LIM[1], KNN_LIM[2]))

# ── render ONE panel ─────────────────────────────────────────────────────────
render_panel <- function(cl){
  sn <- layout_cluster(cl); ed <- chain_edges(sn)
  ggplot()+
    { if(nrow(ed)>0) list(
        geom_segment(data=ed %>% dplyr::filter(share=="TRB", kind=="exact"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRB_COL, linewidth=SHARE_LW, alpha=0.85),
        geom_segment(data=ed %>% dplyr::filter(share=="TRA", kind=="exact"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRA_COL, linewidth=SHARE_LW, alpha=0.85),
        geom_segment(data=ed %>% dplyr::filter(share=="TRB", kind=="one_aa"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRB_COL,
                     linewidth=SHARE_LW_DASH, linetype=DASH_PATTERN, alpha=0.8),
        geom_segment(data=ed %>% dplyr::filter(share=="TRA", kind=="one_aa"),
                     aes(x,y,xend=xend,yend=yend), color=SHARE_TRA_COL,
                     linewidth=SHARE_LW_DASH, linetype=DASH_PATTERN, alpha=0.8)) }+
    ggforce::geom_circle(data=sn, aes(x0=x,y0=y,r=node_r,fill=knn_effect),
                         color="grey30", linewidth=0.25)+
    scale_fill_gradient(low=KNN_LOW, high=KNN_HIGH, limits=KNN_LIM,
                        na.value=KNN_NA, name="kNN bias", guide="none")+
    coord_equal(xlim=c(-GLOBAL_LIM,GLOBAL_LIM), ylim=c(-GLOBAL_LIM,GLOBAL_LIM), clip="off")+
    theme_void()+ ggtitle(sprintf("cluster %s", cl))+
    theme(plot.title=element_text(size=TITLE_SIZE,hjust=0.5,margin=margin(b=2)), plot.margin=margin(2,2,2,2))
}

p0 <-render_panel("0"); p2 <-render_panel("2"); p4 <-render_panel("4")
p8 <-render_panel("8"); p18<-render_panel("18"); p17<-render_panel("17"); p14<-render_panel("14")

# ── compose ──────────────────────────────────────────────────────────────────
pp65_grid <- cowplot::plot_grid(
  p0,
  cowplot::plot_grid(p2, p4, ncol=2),
  cowplot::plot_grid(p8, p18, ncol=2),
  ncol=1, rel_heights=c(1,1,1))
ie1_col <- cowplot::plot_grid(p17, p14, ncol=1, rel_heights=c(1,1))

# kNN colorbar legend (same fixed limits as the panels)
knn_leg <- cowplot::get_legend(
  ggplot(data.frame(x=1, y=1, v=KNN_LIM), aes(x,y,fill=v)) + geom_tile() +
    scale_fill_gradient(low=KNN_LOW, high=KNN_HIGH, limits=KNN_LIM,
                        na.value=KNN_NA, name="kNN bias\n(same-cluster effect)",
                        guide=guide_colourbar(barheight=grid::unit(3,"cm"))) +
    theme_void() + theme(legend.position="right",
                         legend.title=element_text(size=12), legend.text=element_text(size=10)))

# line-type legend
line_df <- data.frame(x=1, xend=2, y=4:1,
  lab=factor(c("Identical TRB","TRB: 1 aa difference","Identical TRA","TRA: 1 aa difference"),
             levels=c("Identical TRB","TRB: 1 aa difference","Identical TRA","TRA: 1 aa difference")))
line_leg <- cowplot::get_legend(
  ggplot(line_df) +
    geom_segment(aes(x=x,xend=xend,y=y,yend=y,color=lab,linetype=lab), linewidth=1.5) +
    scale_color_manual(values=c("Identical TRB"=SHARE_TRB_COL,"TRB: 1 aa difference"=SHARE_TRB_COL,
                                "Identical TRA"=SHARE_TRA_COL,"TRA: 1 aa difference"=SHARE_TRA_COL),
                       name="Shared chain") +
    scale_linetype_manual(values=c("Identical TRB"="solid","TRB: 1 aa difference"=DASH_PATTERN,
                                   "Identical TRA"="solid","TRA: 1 aa difference"=DASH_PATTERN),
                          name="Shared chain") +
    theme_void() + theme(legend.position="right",
                         legend.title=element_text(size=12), legend.text=element_text(size=10)))

legends <- cowplot::plot_grid(knn_leg, line_leg, ncol=1, rel_heights=c(1,1))

divider <- cowplot::ggdraw() +
  cowplot::draw_line(x=c(0.5,0.5), y=c(0.02,0.98), linetype="dashed", color="grey40", linewidth=0.6)

fig <- cowplot::plot_grid(pp65_grid, divider, ie1_col, legends,
                          ncol=4, rel_widths=c(1.0, 0.04, 0.42, 0.18))
fn <- file.path(out_dir, "manuscript_donor_figure_kNN.png")
ggsave(fn, fig, width=20, height=16, dpi=300, bg="white")
message(sprintf("Wrote %s", basename(fn)))




#
# Jensen-Shannon divergence for donor/specificity clonal dominance over time
# — dots plot version, with reference line at y = 0.85
#

# --- Step 0: subset metadata to valid specificity and clones ---
data <- mdata_CMV@meta.data %>%
  filter(!is.na(Global_Clone_ID),
         specificity %in% c("pp65", "IE1")) %>%
  filter(
    (specificity == "pp65" & grepl("_pp65-clone-", Global_Clone_ID)) |
      (specificity == "IE1"  & grepl("_IE1-clone-",  Global_Clone_ID))
  )

# --- Step 1: count total cells per donor/specificity/timepoint ---
total_cells <- data %>%
  group_by(donor_identity, specificity, timepoint) %>%
  summarise(total_cells = n(), .groups = "drop")

# --- Step 2: count cells per donor/specificity/clone/timepoint ---
clone_counts <- data %>%
  group_by(donor_identity, specificity, Global_Clone_ID, timepoint) %>%
  summarise(n_cells = n(), .groups = "drop") %>%
  left_join(total_cells, by = c("donor_identity", "specificity", "timepoint")) %>%
  mutate(freq = n_cells / total_cells) %>%
  dplyr::select(-total_cells)

# --- Step 3: pivot to wide format ---
clone_wide <- clone_counts %>%
  pivot_wider(
    id_cols = c(donor_identity, specificity, Global_Clone_ID),
    names_from = timepoint,
    values_from = c(n_cells, freq),
    values_fill = 0
  )

# --- Step 4: compute JSD per donor/specificity ---
calc_jsd <- function(early, late) {
  mat <- rbind(early, late)
  JSD(mat, unit = "log2")
}

jsd_donor_spec <- clone_wide %>%
  group_by(donor_identity, specificity) %>%
  summarise(jsd = calc_jsd(freq_Early, freq_Late), .groups = "drop") %>%
  mutate(one_minus_jsd = 1 - jsd)

# --- Step 5: prep plotting frame ---
plot_df <- jsd_donor_spec %>%
  mutate(
    donor = donor_identity,
    specificity = factor(specificity, levels = c("pp65", "IE1"))
  )

# --- Step 6: dots plot ---
JS_grouped_plot <- ggplot(
  plot_df,
  aes(
    x = specificity,
    y = one_minus_jsd,
    color = specificity,
    fill  = specificity,
    shape = donor
  )
) +
  # reference line at y = 0.85 (drawn first so points sit on top)
  geom_hline(
    yintercept = 0.85,
    linetype   = "dashed",
    color      = "grey40",
    linewidth  = 0.6
  ) +
  geom_jitter(
    width  = 0.15,
    alpha  = 0.9,
    stroke = 0.6,
    size   = 4
  ) +
  scale_color_manual(
    name = "Specificity",
    values = c("IE1" = "red", "pp65" = "blue")
  ) +
  scale_fill_manual(
    name = "Specificity",
    values = c("IE1" = "red", "pp65" = "blue")
  ) +
  scale_shape_manual(
    name = "Donor",
    values = c(
      "Donor 1" = 16,
      "Donor 2" = 15,
      "Donor 3" = 17,
      "Donor 4" = 25,
      "Donor 5" = 18
    )
  ) +
  guides(
    shape = guide_legend(
      override.aes = list(
        fill  = rep("black", 5),
        color = rep("black", 5),
        size  = 4
      )
    )
  ) +
  labs(
    x = "Specificity",
    y = "1 - Jensen–Shannon divergence (stability)"
  ) +
  scale_y_continuous(
    limits = c(0, 1.1),
    breaks = seq(0, 1, by = 0.2),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    axis.title.x = element_text(size = 16),
    axis.title.y = element_text(size = 16),
    axis.text.x  = element_text(size = 14),
    axis.text.y  = element_text(size = 12)
  )

# --- Step 7: export ---
ggsave("exports/integrated_v5/timepoint_stability_plots/Donor_clone_spread_Jensen-Shannon_divergence_dots_plot.png",
       JS_grouped_plot, width = 5.9, height = 5.9, dpi = 300)




#
# Per-metacluster clone frequency — Early vs Late scatter plots
#   One panel per metacluster (Tcm, Tem, T-CMC1, NK-like, ISG+).
#   x = clone's frequency in that metacluster at Early
#   y = same clone's frequency in that metacluster at Late
#   Point aesthetics match Clonal_phenotype_Jensen-Shannon_grouped_by_specificity.png:
#     size       = clone size (% of that donor x specificity pool, both timepoints)
#     shape      = donor
#     color/fill = specificity
#   Stats per panel:
#     R^2 = goodness of fit of the OLS best-fit line
#     p   = Pearson correlation test (H0: no association between Early and Late)
#           p < 0.05 => Early and Late frequencies are significantly correlated
#   Dotted grey line = y = x (perfect stability reference).
#

library(dplyr); library(tidyr); library(ggplot2); library(stringr)
library(patchwork); library(readr); library(purrr); library(tibble); library(scales)

# ── 0. Settings ──────────────────────────────────────────────────────────────
MIN_CELLS_PER_TIMEPOINT <- 10   # clone must have >= this many cells at BOTH timepoints
metacluster_order <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")

export_dir <- "exports/integrated_v5/metacluster_freq_scatter"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

spec_colors  <- c("pp65" = "blue", "IE1" = "red")   # matches the referenced PNG
donor_shapes <- c("Donor 1" = 16,   # solid circle
                  "Donor 2" = 15,   # solid square
                  "Donor 3" = 17,   # solid triangle
                  "Donor 4" = 25,   # inverted triangle (fillable)
                  "Donor 5" = 18)   # solid diamond

# ── 1. Parse metadata; donor/specificity from Global_Clone_ID ────────────────
# Global_Clone_ID is the reliable source for specificity (the metadata column
# can be mislabeled for a handful of clones).
meta <- mdata_CMV@meta.data %>%
  dplyr::filter(!is.na(Global_Clone_ID),
                !is.na(metacluster),
                timepoint %in% c("Early", "Late")) %>%
  dplyr::mutate(
    donor       = paste("Donor", str_remove(str_extract(Global_Clone_ID, "D\\d"), "D")),
    specificity = str_extract(Global_Clone_ID, "pp65|IE1"),
    metacluster = factor(metacluster, levels = metacluster_order),
    timepoint   = factor(timepoint,   levels = c("Early", "Late"))
  )

# ── 2. Cells per clone x timepoint; apply the both-timepoints cell floor ─────
clone_tp_cells <- meta %>%
  dplyr::count(donor, specificity, Global_Clone_ID, timepoint, name = "n_tp")

clone_tp_wide <- clone_tp_cells %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = n_tp, values_fill = 0)

keep_clones <- clone_tp_wide %>%
  dplyr::filter(Early >= MIN_CELLS_PER_TIMEPOINT,
                Late  >= MIN_CELLS_PER_TIMEPOINT) %>%
  dplyr::pull(Global_Clone_ID)

message(sprintf("Clones retained: %d of %d (>= %d cells at both timepoints)",
                length(keep_clones), nrow(clone_tp_wide), MIN_CELLS_PER_TIMEPOINT))

# ── 3. Metacluster frequency per clone x timepoint ───────────────────────────
# complete() expands the factor levels so a metacluster with zero cells
# registers as freq = 0 rather than dropping out.
freq_df <- meta %>%
  dplyr::filter(Global_Clone_ID %in% keep_clones) %>%
  dplyr::count(donor, specificity, Global_Clone_ID, timepoint, metacluster, name = "n") %>%
  tidyr::complete(
    tidyr::nesting(donor, specificity, Global_Clone_ID),
    timepoint,
    metacluster,
    fill = list(n = 0)
  ) %>%
  dplyr::group_by(Global_Clone_ID, timepoint) %>%
  dplyr::mutate(freq = n / sum(n)) %>%
  dplyr::ungroup()

scatter_df <- freq_df %>%
  dplyr::select(donor, specificity, Global_Clone_ID, metacluster, timepoint, freq) %>%
  tidyr::pivot_wider(names_from = timepoint, values_from = freq, names_prefix = "freq_") %>%
  dplyr::rename(freq_early = freq_Early, freq_late = freq_Late)

# ── 4. Clone size (percent dominance) — denominator uses ALL clones ──────────
# Matches the referenced plot: (n_early + n_late) / total for that donor x specificity.
clone_size <- clone_tp_cells %>%
  dplyr::group_by(donor, specificity, Global_Clone_ID) %>%
  dplyr::summarise(n_total = sum(n_tp), .groups = "drop") %>%
  dplyr::group_by(donor, specificity) %>%
  dplyr::mutate(clone_size = n_total / sum(n_total) * 100) %>%
  dplyr::ungroup() %>%
  dplyr::select(Global_Clone_ID, clone_size)

scatter_df <- scatter_df %>%
  dplyr::left_join(clone_size, by = "Global_Clone_ID") %>%
  dplyr::mutate(
    specificity = factor(specificity, levels = c("pp65", "IE1")),
    donor       = factor(donor, levels = names(donor_shapes)),
    metacluster = factor(metacluster, levels = metacluster_order)
  )

# ── 5. Correlation tests + best-fit line ─────────────────────────────────────
# R^2         : goodness of fit of the OLS line freq_late ~ freq_early.
#               (For a single predictor this equals pearson_r^2.)
# p_pearson   : cor.test H0 = no linear association. Identical to the lm slope
#               p-value. p < 0.05 => Early and Late are significantly correlated.
# p_spearman  : rank-based equivalent; robust to outliers and to the bounded,
#               zero-heavy shape of proportion data.
# slope/intercept are reported so the DIRECTION of any drift is readable.
format_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p)))
}

# cor.test errors if either variable has zero variance (e.g. a metacluster that
# is 0 for every clone at one timepoint), so fail soft to NA rather than abort.
safe_cor <- function(x, y, method) {
  out <- tryCatch(
    suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(out)) return(list(est = NA_real_, p = NA_real_))
  list(est = unname(out$estimate), p = out$p.value)
}

fit_stats <- scatter_df %>%
  dplyr::group_by(metacluster) %>%
  dplyr::group_modify(~ {
    m  <- stats::lm(freq_late ~ freq_early, data = .x)
    s  <- summary(m)
    pe <- safe_cor(.x$freq_early, .x$freq_late, "pearson")
    sp <- safe_cor(.x$freq_early, .x$freq_late, "spearman")

    tibble::tibble(
      n            = nrow(.x),
      intercept    = unname(coef(m)[1]),
      slope        = unname(coef(m)[2]),
      r_squared    = s$r.squared,
      pearson_r    = pe$est,
      p_pearson    = pe$p,
      spearman_rho = sp$est,
      p_spearman   = sp$p
    )
  }) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    p_pearson_adj  = p.adjust(p_pearson,  method = "BH"),
    p_spearman_adj = p.adjust(p_spearman, method = "BH"),
    label = sprintf("R\u00B2 = %.3f\np = %s", r_squared, format_p(p_pearson))
  )

print(fit_stats %>%
        dplyr::select(metacluster, n, slope, r_squared,
                      pearson_r, p_pearson, spearman_rho, p_spearman))

write_csv(fit_stats %>% dplyr::select(-label),
          file.path(export_dir, "metacluster_freq_early_vs_late_regression_stats.csv"))

# ── 6. Plot builder ──────────────────────────────────────────────────────────
# Fixed size limits so a given clone is the same size in every panel/file.
size_limits <- range(scatter_df$clone_size, na.rm = TRUE)

make_scatter <- function(mc) {
  d  <- scatter_df %>% dplyr::filter(metacluster == mc)
  st <- fit_stats  %>% dplyr::filter(metacluster == mc)

  ggplot(d, aes(x = freq_early, y = freq_late)) +
    # y = x reference: points on this line are unchanged between timepoints
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", color = "grey60", linewidth = 0.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                color = "black", linewidth = 0.7) +
    geom_point(aes(color = specificity, fill = specificity,
                   shape = donor, size = clone_size),
               alpha = 0.9, stroke = 0.6) +
    geom_text(data = st, aes(x = 0.02, y = 0.98, label = label),
              inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4) +
    scale_color_manual(name = "Specificity", values = spec_colors) +
    scale_fill_manual(name  = "Specificity", values = spec_colors) +
    scale_shape_manual(name = "Donor", values = donor_shapes,
                       breaks = names(donor_shapes), limits = names(donor_shapes)) +
    scale_size_continuous(name = "Clone Size",
                          range = c(2, 8), limits = size_limits) +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(accuracy = 1)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25),
                       labels = scales::percent_format(accuracy = 1)) +
    coord_fixed() +
    guides(
      color = guide_legend(order = 1),
      fill  = guide_legend(order = 1),
      shape = guide_legend(order = 2, override.aes = list(
        fill  = rep("black", 5),
        color = rep("black", 5),
        size  = 4)),
      size  = guide_legend(order = 3)
    ) +
    labs(
      x     = sprintf("%s Frequency at Early Timepoint", mc),
      y     = sprintf("%s Frequency at Late Timepoint",  mc),
      title = sprintf("%s Phenotype Stability", mc)
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
}

# ── 7. Build and export ──────────────────────────────────────────────────────
plots <- purrr::map(metacluster_order, make_scatter)
names(plots) <- metacluster_order

for (mc in metacluster_order) {
  safe_mc <- str_replace_all(mc, "[^A-Za-z0-9]", "_")   # T-CMC1 -> T_CMC1, ISG+ -> ISG_
  ggsave(
    filename = file.path(export_dir,
                         sprintf("metacluster_freq_early_vs_late_%s.png", safe_mc)),
    plot   = plots[[mc]],
    width  = 7.5, height = 6, dpi = 300, bg = "white"
  )
  message(sprintf("Saved panel: %s", mc))
}

# Combined 3 x 2 grid with a single shared legend
combined_plot <- patchwork::wrap_plots(plots, ncol = 3, guides = "collect") &
  theme(legend.position = "right")

ggsave(
  filename = file.path(export_dir, "metacluster_freq_early_vs_late_ALL.png"),
  plot   = combined_plot,
  width  = 19, height = 11, dpi = 300, bg = "white"
)

# Underlying values, for reference
write_csv(scatter_df,
          file.path(export_dir, "metacluster_freq_early_vs_late_data.csv"))

message("Done — outputs in ", export_dir)




# ============================================================================
#  Cross-object concordance of pp65 vs IE1 differential expression
#  mdata_CMV (D1-D5)  vs  D0_data_CMV (Visits 1-3 | Visit 4 | all visits)
# ============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(SeuratObject); library(Matrix)
  library(dplyr);  library(tibble);  library(purrr)
  library(ggplot2); library(ggrepel); library(patchwork)
})

# ------------------------------------------------------------------ CONFIG --
SPEC_COL_M   <- "specificity"     # <- adjust to your metadata column names
SPEC_COL_D0  <- "specificity"
TP_COL_D0    <- "timepoint"
SPEC_1       <- "pp65"            # ident.1 -> POSITIVE log2FC
SPEC_2       <- "IE1"             # ident.2 -> NEGATIVE log2FC
VISITS_EARLY <- c("Visit 1", "Visit 2", "Visit 3")
VISITS_LATE  <- c("Visit 4")

RNA_ASSAY <- "RNA"
MIN_PCT   <- 0.10   # detected in >=10% of cells of >=1 specificity, per object
LFC_CUT   <- 0.25   # |log2FC| threshold for concordance calls
FDR_CUT   <- 0.05

GENES_LABEL <- c("GZMB","GZMK","IFNG","PRF1","NKG7","KLRG1","CX3CR1","ZEB2",
                 "TCF7","SELL","CCR7","IL7R","TOX","PDCD1","LAG3","ISG15")
AUTO_LABEL_N <- 0   # bump to auto-label top-N concordant genes by magnitude

PAL <- c("Higher in pp65 (both)" = "#1976D2",
         "Higher in IE1 (both)"  = "#8B0000",
         "Discordant"            = "grey45",
         "Below threshold"       = "grey85")

OUT_DIR <- "exports/cross_object_DE"
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# =============================== SETUP ======================================
# Seurat uses presto::wilcoxauc when available; without it, WilcoxDETest falls
# back to future.apply and exports the whole expression matrix to each worker.
HAVE_PRESTO <- requireNamespace("presto", quietly = TRUE)
if (!HAVE_PRESTO)
  message("NOTE: presto not installed - Wilcoxon DE will be much slower.\n",
          "      install.packages('presto')  (or, if unavailable:\n",
          "      remotes::install_github('immunogenomics/presto'))")

# Force DE to run serially regardless of the session-wide future plan, then
# restore the caller's plan and option. Avoids shipping ~800 MiB to workers.
with_serial_future <- function(expr) {
  if (!requireNamespace("future", quietly = TRUE)) return(force(expr))
  old_plan <- future::plan()
  old_max  <- getOption("future.globals.maxSize")
  future::plan(future::sequential)
  options(future.globals.maxSize = Inf)   # size check fires even when serial
  on.exit({
    future::plan(old_plan)
    options(future.globals.maxSize = old_max)
  }, add = TRUE)
  force(expr)
}


# ============================== HELPERS =====================================

`%or%` <- function(a, b) if (is.null(a)) b else a

## --- object integrity -------------------------------------------------------

# Which layers exist, and do they cover every cell in the object?
layer_report <- function(obj, assay = RNA_ASSAY) {
  a <- obj@assays[[assay]]; cl <- Cells(obj)
  do.call(rbind, lapply(Layers(a), function(l) {
    m <- LayerData(a, layer = l)
    data.frame(layer = l, nfeat = nrow(m), ncell = ncol(m),
               covers_all_cells = all(cl %in% colnames(m)))
  }))
}

# Realign assay feature meta.data so nrow == n(features), rownames == features.
# Writes through the S4 slot directly, bypassing [[<- validation.
repair_feature_meta <- function(obj, assay = RNA_ASSAY, verbose = TRUE) {
  a <- obj@assays[[assay]]
  if (!inherits(a, "Assay5")) return(obj)

  feats <- rownames(slot(a, "features"))
  md    <- slot(a, "meta.data")

  if (is.data.frame(md) && nrow(md) == length(feats) &&
      identical(rownames(md), feats)) {
    if (verbose) message("  feature meta.data OK (", length(feats), " features)")
    return(obj)
  }

  n_stale   <- sum(!rownames(md) %in% feats)
  n_missing <- sum(!feats %in% rownames(md))
  if (verbose)
    message("  repairing feature meta.data: ", nrow(md), " rows -> ",
            length(feats), " features  (", n_stale, " stale, ",
            n_missing, " features absent from meta.data)")

  md_new <- if (!is.data.frame(md) || ncol(md) == 0L) {
    data.frame(row.names = feats)
  } else {
    tmp <- md[match(feats, rownames(md)), , drop = FALSE]
    rownames(tmp) <- feats
    # NA in a logical var-feature flag confuses VariableFeatures(); force FALSE
    for (j in which(vapply(tmp, is.logical, logical(1))))
      tmp[[j]][is.na(tmp[[j]])] <- FALSE
    tmp
  }
  obj@assays[[assay]]@meta.data <- md_new
  obj
}

# Remove redundant split layers from a WORKING COPY (never the session object).
drop_split_layers <- function(obj, assay = RNA_ASSAY, verbose = TRUE) {
  extra <- grep("^(counts|data)\\.", Layers(obj@assays[[assay]]), value = TRUE)
  for (l in extra) {
    ok <- try({ LayerData(obj, assay = assay, layer = l) <- NULL }, silent = TRUE)
    if (inherits(ok, "try-error") && verbose)
      message("  could not drop layer '", l, "' - leaving in place")
  }
  if (verbose && length(extra))
    message("  dropped redundant split layers: ", paste(extra, collapse = ", "))
  obj
}

prep_obj <- function(obj, assay = RNA_ASSAY, verbose = TRUE,
                     need = c("counts", "data"), tidy_layers = TRUE) {
  DefaultAssay(obj) <- assay
  obj <- repair_feature_meta(obj, assay, verbose)

  a <- obj@assays[[assay]]
  if (!inherits(a, "Assay5")) return(obj)

  cl <- Cells(obj)
  complete <- vapply(need, function(l) {
    l %in% Layers(a) && all(cl %in% colnames(LayerData(a, layer = l)))
  }, logical(1))

  if (!all(complete)) {
    if (verbose)
      message("  incomplete: ", paste(need[!complete], collapse = ", "),
              " - attempting JoinLayers")
    joined <- try(JoinLayers(obj, assay = assay), silent = TRUE)
    if (inherits(joined, "try-error"))
      stop("JoinLayers failed after meta.data repair; run layer_report(obj) ",
           "and rebuild the assay manually from the layers that do exist.")
    obj <- joined
  } else if (verbose) {
    message("  layers '", paste(need, collapse = "', '"), "' complete - no join needed")
  }

  if (tidy_layers) obj <- drop_split_layers(obj, assay, verbose)
  obj
}

## --- metadata validation ----------------------------------------------------

check_levels <- function(obj, col, wanted, nm) {
  if (!col %in% colnames(obj@meta.data))
    stop("column '", col, "' not found in ", nm, " metadata")
  have <- unique(as.character(obj@meta.data[[col]]))
  miss <- setdiff(wanted, have)
  if (length(miss))
    stop(nm, " column '", col, "': not found -> ", paste(miss, collapse = ", "),
         "\n  observed values: ", paste(sort(have), collapse = " | "))
  invisible(TRUE)
}

## --- differential expression ------------------------------------------------

# genes detected in >= min_pct of cells in AT LEAST ONE specificity group
detected_genes <- function(obj, spec_col, min_pct = MIN_PCT, assay = RNA_ASSAY) {
  cnts <- GetAssayData(obj, assay = assay, layer = "counts")
  md   <- as.character(obj@meta.data[[spec_col]])
  p1 <- Matrix::rowMeans(cnts[, which(md == SPEC_1), drop = FALSE] > 0)
  p2 <- Matrix::rowMeans(cnts[, which(md == SPEC_2), drop = FALSE] > 0)
  rownames(cnts)[pmax(p1, p2, na.rm = TRUE) >= min_pct]
}

lfc_table <- function(obj, spec_col, features, assay = RNA_ASSAY, label = "") {
  n <- table(factor(as.character(obj@meta.data[[spec_col]]),
                    levels = c(SPEC_1, SPEC_2)))
  if (any(n < 3))
    stop(label, ": too few cells (", SPEC_1, " = ", n[[1]], ", ",
         SPEC_2, " = ", n[[2]], ") for a stable comparison")

  Idents(obj) <- spec_col
  res <- with_serial_future(
    FindMarkers(obj, ident.1 = SPEC_1, ident.2 = SPEC_2,
                features = features, assay = assay,
                logfc.threshold = 0, min.pct = 0,
                min.cells.feature = 0, min.cells.group = 1,
                test.use = "wilcox", verbose = FALSE)
  )
  res |>
    rownames_to_column("gene") |>
    transmute(gene, lfc = avg_log2FC, fdr = p_val_adj,
              pct1 = pct.1, pct2 = pct.2) |>
    as_tibble()
}

# total-least-squares slope (both axes carry measurement error)
orth_fit <- function(x, y) {
  v <- prcomp(cbind(x, y), center = TRUE, scale. = FALSE)$rotation[, 1]
  s <- unname(v[2] / v[1])
  c(slope = s, intercept = mean(y) - s * mean(x))
}


# ============================ BUILD SUBSETS =================================

check_levels(mdata_CMV,   SPEC_COL_M,  c(SPEC_1, SPEC_2), "mdata_CMV")
check_levels(D0_data_CMV, SPEC_COL_D0, c(SPEC_1, SPEC_2), "D0_data_CMV")
check_levels(D0_data_CMV, TP_COL_D0, c(VISITS_EARLY, VISITS_LATE), "D0_data_CMV")

message("Preparing mdata_CMV")
m_prep  <- prep_obj(mdata_CMV)
message("Preparing D0_data_CMV")
d0_prep <- prep_obj(D0_data_CMV)

print(layer_report(d0_prep))

# restrict both objects to cells carrying one of the two specificities
keep_spec <- function(obj, col)
  subset(obj, cells = colnames(obj)[obj@meta.data[[col]] %in% c(SPEC_1, SPEC_2)])

m_prep  <- keep_spec(m_prep,  SPEC_COL_M)
d0_prep <- keep_spec(d0_prep, SPEC_COL_D0)

cat("\nD0 specificity x timepoint:\n")
print(table(d0_prep@meta.data[[SPEC_COL_D0]], d0_prep@meta.data[[TP_COL_D0]]))

cells_by <- function(obj, col, vals) colnames(obj)[obj@meta.data[[col]] %in% vals]

d0_sets <- list(
  "D0 Visits 1-3" = subset(d0_prep, cells = cells_by(d0_prep, TP_COL_D0, VISITS_EARLY)),
  "D0 Visit 4"    = subset(d0_prep, cells = cells_by(d0_prep, TP_COL_D0, VISITS_LATE)),
  "D0 all visits" = d0_prep
)


# =============================== RUN DE =====================================

genes_m <- detected_genes(m_prep, SPEC_COL_M)
message("mdata_CMV: ", length(genes_m), " detected genes - running DE",
        if (!HAVE_PRESTO) " (no presto; this may take several minutes)" else "")
lfc_m <- lfc_table(m_prep, SPEC_COL_M, genes_m, label = "mdata_CMV")

results <- imap(d0_sets, function(obj, nm) {
  g  <- detected_genes(obj, SPEC_COL_D0)
  message(nm, ": ", ncol(obj), " cells | ", length(g), " detected - running DE")
  y  <- lfc_table(obj, SPEC_COL_D0, g, label = nm)
  df <- inner_join(
    rename(lfc_m, lfc_x = lfc, fdr_x = fdr, pct1_x = pct1, pct2_x = pct2),
    rename(y,     lfc_y = lfc, fdr_y = fdr, pct1_y = pct1, pct2_y = pct2),
    by = "gene") |>
    mutate(comparison = nm)
  message("  -> ", nrow(df), " shared with mdata_CMV")
  df
})


# ============================= PLOT MAKER ===================================

concordance_plot <- function(df, y_lab, title,
                             x_lab = sprintf("log2FC %s vs %s \u2014 mdata_CMV (D1-D5)", SPEC_1, SPEC_2),
                             label_genes = GENES_LABEL, auto_n = AUTO_LABEL_N) {

  df <- df |> mutate(
    cat = case_when(
      lfc_x >  LFC_CUT & lfc_y >  LFC_CUT           ~ "Higher in pp65 (both)",
      lfc_x < -LFC_CUT & lfc_y < -LFC_CUT           ~ "Higher in IE1 (both)",
      abs(lfc_x) > LFC_CUT & abs(lfc_y) > LFC_CUT   ~ "Discordant",
      TRUE                                          ~ "Below threshold"),
    cat = factor(cat, levels = names(PAL)))

  resp <- filter(df, fdr_x < FDR_CUT | fdr_y < FDR_CUT)
  conc <- filter(df, abs(lfc_x) > LFC_CUT & abs(lfc_y) > LFC_CUT)

  r_all   <- cor(df$lfc_x, df$lfc_y, method = "pearson")
  rho_all <- cor(df$lfc_x, df$lfc_y, method = "spearman")
  r_sig   <- if (nrow(resp) > 3) cor(resp$lfc_x, resp$lfc_y) else NA_real_
  ols     <- unname(coef(lm(lfc_y ~ lfc_x, data = df))[2])
  os      <- orth_fit(df$lfc_x, df$lfc_y)
  pc      <- if (nrow(conc)) 100 * mean(sign(conc$lfc_x) == sign(conc$lfc_y)) else NA_real_

  stat_lab <- sprintf(
    "Pearson r = %.2f  (r = %.2f, FDR-sig)\nSpearman \u03c1 = %.2f\nTLS slope = %.2f (black) | OLS = %.2f (grey)\nn = %s shared genes\n%.0f%% sign-concordant (n = %s)",
    r_all, r_sig, rho_all, os[["slope"]], ols,
    format(nrow(df), big.mark = ","), pc, format(nrow(conc), big.mark = ","))

  lab_df <- filter(df, gene %in% label_genes)
  if (auto_n > 0) {
    top <- conc |> filter(sign(lfc_x) == sign(lfc_y)) |>
      slice_max(abs(lfc_x) + abs(lfc_y), n = auto_n)
    lab_df <- bind_rows(lab_df, top) |> distinct(gene, .keep_all = TRUE)
  }

  lim <- max(abs(c(df$lfc_x, df$lfc_y)), na.rm = TRUE) * 1.06

  ggplot(df, aes(lfc_x, lfc_y)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_point(data = filter(df, cat == "Below threshold"),
               colour = PAL[["Below threshold"]], size = 0.7, alpha = 0.55) +
    geom_point(data = filter(df, cat != "Below threshold"),
               aes(colour = cat), size = 1.3, alpha = 0.85) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                colour = "grey30", linewidth = 0.45) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                colour = "grey55", linewidth = 0.5) +
    geom_abline(slope = os[["slope"]], intercept = os[["intercept"]],
                colour = "black", linewidth = 0.6) +
    geom_point(data = lab_df, shape = 21, fill = NA, colour = "black",
               size = 2.1, stroke = 0.5) +
    geom_text_repel(data = lab_df, aes(label = gene), size = 3,
                    min.segment.length = 0, segment.size = 0.3,
                    segment.colour = "grey40", box.padding = 0.4,
                    max.overlaps = Inf, seed = 42) +
    annotate("text", x = -lim * 0.96, y = lim * 0.96, label = stat_lab,
             hjust = 0, vjust = 1, size = 2.8, lineheight = 1.2, colour = "grey15") +
    scale_colour_manual(values = PAL, drop = FALSE, name = NULL) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = x_lab, y = y_lab, title = title) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 11, hjust = 0),
      axis.text       = element_text(colour = "black"),
      axis.title      = element_text(size = 9.5),
      legend.position = "bottom",
      legend.text     = element_text(size = 8.5),
      legend.key.size = unit(0.4, "cm"),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line       = element_blank()
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.4), nrow = 2))
}


# ============================== ASSEMBLE ====================================

plots <- imap(results, ~ concordance_plot(
  .x,
  y_lab = sprintf("log2FC %s vs %s \u2014 %s", SPEC_1, SPEC_2, .y),
  title = sprintf("mdata_CMV  vs  %s", .y)))

fig <- wrap_plots(plots, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

ggsave(file.path(OUT_DIR, "pp65_vs_IE1_crossobject_concordance.png"),
       fig, width = 16, height = 6.2, dpi = 400, bg = "white")
ggsave(file.path(OUT_DIR, "pp65_vs_IE1_crossobject_concordance.pdf"),
       fig, width = 16, height = 6.2, device = cairo_pdf, bg = "white")

readr::write_csv(bind_rows(results),
                 file.path(OUT_DIR, "pp65_vs_IE1_lfc_shared_genes.csv"))

fig




#
# Clonal dominance — Early vs Late scatter plot
#   x = clone's % of that donor's pp65- (or IE1-) specific cells at Early
#   y = same clone's % at Late
#   Point aesthetics match the metacluster stability scatters:
#     size       = clone size (% of that donor x specificity pool, both timepoints)
#     shape      = donor
#     color/fill = specificity
#   Stats:
#     R^2 = goodness of fit of the OLS best-fit line
#     p   = Pearson correlation test (H0: no association between Early and Late)
#           p < 0.05 => Early and Late dominance are significantly correlated
#   Dotted grey line = y = x (dominance unchanged between timepoints).
#

library(dplyr); library(tidyr); library(ggplot2); library(stringr)
library(readr); library(tibble); library(scales)

# ── 0. Settings ──────────────────────────────────────────────────────────────
# DENOMINATOR: what counts as "that donor's pp65-/IE1-specific cells"
#   "all_cells"    = every cell of that donor x specificity x timepoint,
#                    including cells not assigned to a major (Global) clone
#   "major_clones" = only cells belonging to Global clones (matches the
#                    denominator used in the JSD / clone-heterogeneity figures)
DENOMINATOR <- "major_clones"

# Clone must have >= this many cells summed across BOTH timepoints.
# Deliberately NOT a per-timepoint floor: a clone at 25% early and 0% late is a
# real result, and a both-timepoints filter would silently delete it.
MIN_CELLS_TOTAL <- 00

export_dir <- "exports/integrated_v5/clonal_dominance_scatter"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

spec_colors  <- c("pp65" = "blue", "IE1" = "red")
donor_shapes <- c("Donor 1" = 16,   # solid circle
                  "Donor 2" = 15,   # solid square
                  "Donor 3" = 17,   # solid triangle
                  "Donor 4" = 25,   # inverted triangle (fillable)
                  "Donor 5" = 18)   # solid diamond

# ── 1. Clone-level cells; donor/specificity from Global_Clone_ID ─────────────
# Global_Clone_ID is the reliable source for a CLONE's specificity (the
# per-cell metadata column can be mislabeled for a handful of clones).
clone_cells <- mdata_CMV@meta.data %>%
  dplyr::filter(!is.na(Global_Clone_ID),
                timepoint %in% c("Early", "Late")) %>%
  dplyr::mutate(
    donor       = paste("Donor", str_remove(str_extract(Global_Clone_ID, "D\\d"), "D")),
    specificity = str_extract(Global_Clone_ID, "pp65|IE1"),
    timepoint   = factor(timepoint, levels = c("Early", "Late"))
  ) %>%
  dplyr::filter(!is.na(specificity)) %>%
  dplyr::count(donor, specificity, Global_Clone_ID, timepoint, name = "n_cells")

# Fill in absent timepoints as 0 so expansions/contractions survive
clone_cells <- clone_cells %>%
  tidyr::complete(
    tidyr::nesting(donor, specificity, Global_Clone_ID),
    timepoint = factor(c("Early", "Late"), levels = c("Early", "Late")),
    fill = list(n_cells = 0)
  )

# ── 2. Denominator: the donor's specificity-defined pool per timepoint ───────
# NOTE: for "all_cells" the denominator necessarily uses the per-cell
# `specificity` metadata column (cells without a Global_Clone_ID have no other
# source), while clone specificity comes from Global_Clone_ID. For the handful
# of clones where those two disagree, a clone's cells are counted against the
# pool of its Global_Clone_ID specificity. Effect is small; noted for the record.
if (DENOMINATOR == "all_cells") {
  denom <- mdata_CMV@meta.data %>%
    dplyr::filter(specificity %in% c("pp65", "IE1"),
                  donor_identity %in% names(donor_shapes),
                  timepoint %in% c("Early", "Late")) %>%
    dplyr::mutate(timepoint = factor(timepoint, levels = c("Early", "Late"))) %>%
    dplyr::count(donor = donor_identity, specificity, timepoint,
                 name = "denom_cells")
} else if (DENOMINATOR == "major_clones") {
  denom <- clone_cells %>%
    dplyr::group_by(donor, specificity, timepoint) %>%
    dplyr::summarise(denom_cells = sum(n_cells), .groups = "drop")
} else {
  stop("DENOMINATOR must be 'all_cells' or 'major_clones'")
}

# ── 3. Percent of response, per clone per timepoint ──────────────────────────
dominance <- clone_cells %>%
  dplyr::left_join(denom, by = c("donor", "specificity", "timepoint")) %>%
  dplyr::mutate(pct_of_response = dplyr::if_else(denom_cells > 0,
                                                 n_cells / denom_cells * 100,
                                                 NA_real_))

scatter_df <- dominance %>%
  dplyr::select(donor, specificity, Global_Clone_ID, timepoint,
                n_cells, pct_of_response) %>%
  tidyr::pivot_wider(names_from  = timepoint,
                     values_from = c(n_cells, pct_of_response)) %>%
  dplyr::rename(n_early   = n_cells_Early,        n_late   = n_cells_Late,
                pct_early = pct_of_response_Early, pct_late = pct_of_response_Late) %>%
  dplyr::mutate(n_total = n_early + n_late) %>%
  dplyr::filter(n_total >= MIN_CELLS_TOTAL,
                !is.na(pct_early), !is.na(pct_late))

message(sprintf("Clones plotted: %d (>= %d cells across both timepoints; denominator = %s)",
                nrow(scatter_df), MIN_CELLS_TOTAL, DENOMINATOR))

# ── 4. Clone size for the point-size aesthetic (both timepoints pooled) ──────
denom_pooled <- denom %>%
  dplyr::group_by(donor, specificity) %>%
  dplyr::summarise(denom_total = sum(denom_cells), .groups = "drop")

scatter_df <- scatter_df %>%
  dplyr::left_join(denom_pooled, by = c("donor", "specificity")) %>%
  dplyr::mutate(
    clone_size  = n_total / denom_total * 100,
    specificity = factor(specificity, levels = c("pp65", "IE1")),
    donor       = factor(donor, levels = names(donor_shapes))
  )

# ── 5. Correlation test + best-fit line ──────────────────────────────────────
format_p <- function(p) {
  ifelse(is.na(p), "NA",
         ifelse(p < 0.001, sprintf("%.2e", p), sprintf("%.3f", p)))
}

safe_cor <- function(x, y, method) {
  out <- tryCatch(
    suppressWarnings(stats::cor.test(x, y, method = method, exact = FALSE)),
    error = function(e) NULL
  )
  if (is.null(out)) return(list(est = NA_real_, p = NA_real_))
  list(est = unname(out$estimate), p = out$p.value)
}

compute_stats <- function(d) {
  m  <- stats::lm(pct_late ~ pct_early, data = d)
  s  <- summary(m)
  pe <- safe_cor(d$pct_early, d$pct_late, "pearson")
  sp <- safe_cor(d$pct_early, d$pct_late, "spearman")
  tibble::tibble(
    n            = nrow(d),
    intercept    = unname(coef(m)[1]),
    slope        = unname(coef(m)[2]),
    r_squared    = s$r.squared,
    pearson_r    = pe$est,
    p_pearson    = pe$p,
    spearman_rho = sp$est,
    p_spearman   = sp$p
  )
}

fit_pooled <- compute_stats(scatter_df) %>%
  dplyr::mutate(label = sprintf("R\u00B2 = %.3f\np = %s",
                                r_squared, format_p(p_pearson)))

fit_by_spec <- scatter_df %>%
  dplyr::group_by(specificity) %>%
  dplyr::group_modify(~ compute_stats(.x)) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(label = sprintf("R\u00B2 = %.3f\np = %s",
                                r_squared, format_p(p_pearson)))

print(fit_pooled  %>% dplyr::select(-label))
print(fit_by_spec %>% dplyr::select(-label))

write_csv(
  dplyr::bind_rows(
    fit_pooled  %>% dplyr::mutate(group = "all clones", .before = 1),
    fit_by_spec %>% dplyr::mutate(group = as.character(specificity), .before = 1) %>%
      dplyr::select(-specificity)
  ) %>% dplyr::select(-label),
  file.path(export_dir, "clonal_dominance_early_vs_late_regression_stats.csv")
)

# ── 6. Shared scale settings ─────────────────────────────────────────────────
axis_max    <- max(c(scatter_df$pct_early, scatter_df$pct_late), na.rm = TRUE)
axis_lim    <- ceiling(axis_max / 10) * 10          # round up to a clean decade
size_limits <- range(scatter_df$clone_size, na.rm = TRUE)

base_layers <- function() {
  list(
    geom_abline(slope = 1, intercept = 0,
                linetype = "dotted", color = "grey60", linewidth = 0.5),
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
                color = "black", linewidth = 0.7),
    geom_point(aes(color = specificity, fill = specificity,
                   shape = donor, size = clone_size),
               alpha = 0.9, stroke = 0.6),
    scale_color_manual(name = "Specificity", values = spec_colors),
    scale_fill_manual(name  = "Specificity", values = spec_colors),
    scale_shape_manual(name = "Donor", values = donor_shapes,
                       breaks = names(donor_shapes), limits = names(donor_shapes)),
    scale_size_continuous(name = "Clone Size", range = c(2, 8), limits = size_limits),
    scale_x_continuous(limits = c(0, axis_lim)),
    scale_y_continuous(limits = c(0, axis_lim)),
    coord_fixed(),
    guides(
      color = guide_legend(order = 1),
      fill  = guide_legend(order = 1),
      shape = guide_legend(order = 2, override.aes = list(
        fill  = rep("black", 5),
        color = rep("black", 5),
        size  = 4)),
      size  = guide_legend(order = 3)
    ),
    labs(x = "% of Response at Early Timepoint",
         y = "% of Response at Late Timepoint"),
    theme_bw(base_size = 14),
    theme(
      plot.title       = element_text(face = "bold", hjust = 0.5),
      legend.position  = "right",
      panel.grid.minor = element_blank()
    )
  )
}

# ── 7. Pooled plot ───────────────────────────────────────────────────────────
dominance_plot <- ggplot(scatter_df, aes(x = pct_early, y = pct_late)) +
  base_layers() +
  geom_text(data = fit_pooled,
            aes(x = axis_lim * 0.02, y = axis_lim * 0.98, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4) +
  ggtitle("Clonal Dominance Stability")

ggsave(file.path(export_dir, "clonal_dominance_early_vs_late.png"),
       dominance_plot, width = 7.5, height = 6, dpi = 300, bg = "white")

# ── 8. Split by specificity ──────────────────────────────────────────────────
dominance_plot_split <- ggplot(scatter_df, aes(x = pct_early, y = pct_late)) +
  base_layers() +
  geom_text(data = fit_by_spec,
            aes(x = axis_lim * 0.02, y = axis_lim * 0.98, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 4) +
  facet_wrap(~ specificity, ncol = 2) +
  ggtitle("Clonal Dominance Stability") +
  theme(strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 14))

ggsave(file.path(export_dir, "clonal_dominance_early_vs_late_by_specificity.png"),
       dominance_plot_split, width = 13, height = 6, dpi = 300, bg = "white")

# Underlying values, for reference
write_csv(scatter_df,
          file.path(export_dir, "clonal_dominance_early_vs_late_data.csv"))

message("Done — outputs in ", export_dir)





# ============================================================================
#  Clone panels v3 — cluster / specificity / timepoint UMAPs + cluster bars
#
#  Four figure families per object:
#    <label>/cluster_colors/                        dots by seurat_clusters
#    <label>/cluster_colors_with_bar_plots/         dots by seurat_clusters + bar
#    <label>/with_bar_plots/                        dots by specificity      + bar
#    <label>/timepoint_colors_with_bar_plots/       dots by timepoint        + bar
#
#  Bar = % of that clone's cells per seurat_cluster, y-axis on the RIGHT.
#
#  Five objects (all under exports/.../clone_panels/<label>/):
#    full                       D0_data_CMV, original UMAP + clusters
#    V1_3                       Visit 1-3 subset, RE-clustered, new UMAP
#    V4                         Visit 4   subset, RE-clustered, new UMAP
#    V1_3_original_clustering   Visit 1-3 cells on the ORIGINAL UMAP/clusters
#    V4_original_clustering     Visit 4   cells on the ORIGINAL UMAP/clusters
#
#  pp65 = 20 clones / ncol 5 ; IE1 = 7 clones / ncol 3 (skinny)
#  Self-contained: needs only D0_data_CMV in the session.
# ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(stringr)
library(purrr)
library(tibble)
library(readr)

# ── 0. Shared settings (reuse from main script if already defined) ──────────
if (!exists("base_export")) base_export <- "Donor_0/exports/timepoint_split"

if (!exists("save_jpeg")) {
  save_jpeg <- function(plot, path, width, height) {
    ggsave(filename = path, plot = plot, width = width, height = height,
           dpi = 300, device = "jpeg", limitsize = FALSE)
  }
}

if (!exists("spec_colors")) {
  spec_colors <- c("pp65" = "#1976D2", "IE1" = "#8B0000")   # blue / red
}

timepoint_colors <- c(
  "Visit 1" = "#4DAF4A",
  "Visit 2" = "#377EB8",
  "Visit 3" = "#FF7F00",
  "Visit 4" = "#E41A1C"
)

cluster_colors <- c(
  "0" = "#C49A00",
  "1" = "#53B400",
  "2" = "#8B0000",
  "3" = "#1976D2",
  "4" = "#F8766D",
  "5" = "#A58AFF",
  "6" = "#FFDB58",
  "7" = "#00C094"
)

# Fallback colors if a re-clustered subset yields more than 8 clusters
cluster_colors_extra <- c("#00B6EB", "#FB61D7", "#7F7F7F", "#B79F00",
                          "#00BA38", "#619CFF", "#F564E3", "#4D4D4D")

CLUSTER_COL   <- "seurat_clusters"
SPEC_COL      <- "specificity"
TIME_COL      <- "timepoint"
CLONE_COL     <- "Global_Clone_ID"
BAR_REL_WIDTH <- 0.26    # thin bar + right-side axis, relative to UMAP panel
PANEL_W       <- 2.6     # inches per UMAP cell (matches existing figures)
PANEL_H       <- 2.8
LEGEND_W      <- 1.5     # inches reserved for the right-hand legend column
UMAP_PAD      <- 0.04    # fraction of range padded around fixed axis limits
WRITE_CSV     <- TRUE    # also export clone x cluster % tables
SKIP_SINGLE_TP <- FALSE  # FALSE -> V4 still gets a timepoint-colored figure

# Subset rebuild controls (re-clustered objects only)
FORCE_RECLUSTER <- FALSE   # TRUE -> ignore cached RDS and recompute
SAVE_SUBSETS    <- TRUE    # write rebuilt objects back to objects/
RECLUSTER_RES   <- 0.8     # matches the main timepoint-split script
RECLUSTER_DIMS  <- 1:30
RECLUSTER_NPCS  <- 30
RECLUSTER_NFEAT <- 2000

# ── 1. Guard ───────────────────────────────────────────────────────────────
if (!exists("D0_data_CMV")) {
  stop("D0_data_CMV not found in the session — load it before running this script.")
}

obj_dir <- file.path(base_export, "objects")
dir.create(obj_dir, recursive = TRUE, showWarnings = FALSE)

message("Timepoint values present: ",
        paste(sort(unique(as.character(D0_data_CMV@meta.data[[TIME_COL]]))),
              collapse = ", "))

V1_3_VISITS <- c("Visit 1", "Visit 2", "Visit 3")
V4_VISITS   <- c("Visit 4")

# ── 2. Re-clustered subsets: session -> cached RDS -> rebuild ──────────────
recluster_obj_v3 <- function(obj, label = "object") {
  message(sprintf("  re-clustering %s (%d cells)...", label, ncol(obj)))
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, normalization.method = "LogNormalize",
                       scale.factor = 1e4, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst",
                              nfeatures = RECLUSTER_NFEAT, verbose = FALSE)
  obj <- ScaleData(obj, features = rownames(obj), verbose = FALSE)
  obj <- RunPCA(obj, npcs = RECLUSTER_NPCS, verbose = FALSE)
  obj <- FindNeighbors(obj, dims = RECLUSTER_DIMS, verbose = FALSE)
  obj <- FindClusters(obj, resolution = RECLUSTER_RES, verbose = FALSE)
  obj <- RunUMAP(obj, dims = RECLUSTER_DIMS, verbose = FALSE)
  message(sprintf("  %s: %d clusters", label,
                  length(levels(obj$seurat_clusters))))
  obj
}

get_or_build_subset <- function(var_name, rds_name, visits, label) {
  if (!FORCE_RECLUSTER && exists(var_name, envir = .GlobalEnv)) {
    message(sprintf("[%s] using '%s' already in session", label, var_name))
    return(get(var_name, envir = .GlobalEnv))
  }

  rds_path <- file.path(obj_dir, rds_name)
  if (!FORCE_RECLUSTER && file.exists(rds_path)) {
    message(sprintf("[%s] loading cached %s", label, rds_name))
    obj <- readRDS(rds_path)
    assign(var_name, obj, envir = .GlobalEnv)
    return(obj)
  }

  message(sprintf("[%s] rebuilding subset from D0_data_CMV", label))
  cells <- colnames(D0_data_CMV)[D0_data_CMV@meta.data[[TIME_COL]] %in% visits]
  if (length(cells) == 0)
    stop(sprintf("No cells matched timepoints: %s", paste(visits, collapse = ", ")))
  message(sprintf("  %s subset: %d cells", label, length(cells)))

  obj <- subset(D0_data_CMV, cells = cells)
  obj <- recluster_obj_v3(obj, label = label)

  if (SAVE_SUBSETS) {
    saveRDS(obj, rds_path)
    message(sprintf("  saved -> %s", rds_path))
  }
  assign(var_name, obj, envir = .GlobalEnv)
  obj
}

obj_V1_3 <- get_or_build_subset("obj_V1_3", "obj_V1_3_reclustered.rds",
                                V1_3_VISITS, "V1_3")
obj_V4   <- get_or_build_subset("obj_V4",   "obj_V4_reclustered.rds",
                                V4_VISITS,   "V4")

# ── 3. Original-clustering subsets (no re-clustering, no new UMAP) ─────────
subset_keep_clustering <- function(visits, label) {
  cells <- colnames(D0_data_CMV)[D0_data_CMV@meta.data[[TIME_COL]] %in% visits]
  if (length(cells) == 0)
    stop(sprintf("No cells matched timepoints: %s", paste(visits, collapse = ", ")))
  message(sprintf("[%s] %d cells on the original UMAP / clusters",
                  label, length(cells)))
  subset(D0_data_CMV, cells = cells)
}

obj_V1_3_orig <- subset_keep_clustering(V1_3_VISITS, "V1_3_original_clustering")
obj_V4_orig   <- subset_keep_clustering(V4_VISITS,   "V4_original_clustering")

# ── 4. Palettes and fixed axis limits ──────────────────────────────────────
build_cluster_pal <- function(obj, level_source = NULL) {
  src <- if (is.null(level_source)) obj else level_source
  lv  <- levels(factor(as.character(src@meta.data[[CLUSTER_COL]])))
  num <- suppressWarnings(as.integer(lv))
  lv  <- if (all(!is.na(num))) lv[order(num)] else sort(lv)

  pal     <- cluster_colors[lv]
  missing <- is.na(pal)
  if (any(missing)) {
    n_miss <- sum(missing)
    pal[missing] <- rep(cluster_colors_extra, length.out = n_miss)
    message(sprintf("  note: %d cluster(s) not in cluster_colors, using fallback: %s",
                    n_miss, paste(lv[missing], collapse = ", ")))
  }
  setNames(as.character(pal), lv)
}

build_timepoint_pal <- function(obj) {
  present <- unique(as.character(obj@meta.data[[TIME_COL]]))
  keep    <- names(timepoint_colors)[names(timepoint_colors) %in% present]
  if (length(keep) == 0) return(timepoint_colors)
  timepoint_colors[keep]
}

umap_limits <- function(obj, pad = UMAP_PAD) {
  emb <- Embeddings(obj, "umap")[, 1:2]
  xr  <- range(emb[, 1]); yr <- range(emb[, 2])
  list(x = xr + c(-1, 1) * diff(xr) * pad,
       y = yr + c(-1, 1) * diff(yr) * pad)
}

# ── 5. Barcode-safe metadata accessors ─────────────────────────────────────
get_cluster_vec <- function(obj, pal) {
  v <- as.character(obj@meta.data[[CLUSTER_COL]])
  setNames(factor(v, levels = names(pal)), rownames(obj@meta.data))
}
get_spec_vec <- function(obj) {
  v <- as.character(obj@meta.data[[SPEC_COL]])
  setNames(factor(v, levels = names(spec_colors)), rownames(obj@meta.data))
}
get_timepoint_vec <- function(obj, pal) {
  v <- as.character(obj@meta.data[[TIME_COL]])
  setNames(factor(v, levels = names(pal)), rownames(obj@meta.data))
}
get_clone_vec <- function(obj) {
  setNames(as.character(obj@meta.data[[CLONE_COL]]), rownames(obj@meta.data))
}

# ── 6. UMAP clone panel ────────────────────────────────────────────────────
make_panel_umap_v3 <- function(obj, clone_id, pal_cluster,
                               color_by      = c("cluster", "specificity", "timepoint"),
                               pal_timepoint = NULL,
                               lims          = NULL,
                               tight_right   = FALSE) {
  color_by <- match.arg(color_by)

  emb <- as.data.frame(Embeddings(obj, "umap"))[, 1:2]
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  emb$clone <- get_clone_vec(obj)[rownames(emb)]

  if (color_by == "cluster") {
    emb$grp <- get_cluster_vec(obj, pal_cluster)[rownames(emb)]
    pal <- pal_cluster
  } else if (color_by == "specificity") {
    emb$grp <- get_spec_vec(obj)[rownames(emb)]
    pal <- spec_colors
  } else {
    if (is.null(pal_timepoint)) pal_timepoint <- build_timepoint_pal(obj)
    emb$grp <- get_timepoint_vec(obj, pal_timepoint)[rownames(emb)]
    pal <- pal_timepoint
  }

  emb$highlight <- !is.na(emb$clone) & emb$clone == clone_id
  n <- sum(emb$highlight)

  p <- ggplot() +
    geom_point(data = emb[!emb$highlight, ], aes(UMAP_1, UMAP_2),
               color = "grey85", size = 0.3) +
    geom_point(data = emb[emb$highlight, ], aes(UMAP_1, UMAP_2, color = grp),
               size = 1.4) +
    scale_color_manual(values = pal, drop = FALSE) +
    ggtitle(paste0(clone_id, "\n(n = ", n, ")"))

  if (!is.null(lims)) p <- p + coord_cartesian(xlim = lims$x, ylim = lims$y)

  p +
    theme_bw() +
    theme(
      plot.title      = element_text(size = 8, face = "bold", hjust = 0.5),
      axis.title      = element_blank(),
      axis.text       = element_blank(),
      axis.ticks      = element_blank(),
      legend.position = "none",
      plot.margin     = if (tight_right) margin(2, 0, 2, 2) else margin(2, 2, 2, 2)
    )
}

# ── 7. Thin stacked bar: % of clone's cells per cluster, axis on the RIGHT ──
clone_cluster_pct <- function(obj, clone_id, pal_cluster) {
  cl    <- get_cluster_vec(obj, pal_cluster)
  clone <- get_clone_vec(obj)
  keep  <- !is.na(clone) & clone == clone_id
  tibble(cluster = cl[keep]) %>%
    count(cluster, .drop = FALSE) %>%
    mutate(pct = if (sum(n) > 0) 100 * n / sum(n) else 0)
}

make_panel_bar_v3 <- function(obj, clone_id, pal_cluster) {
  df <- clone_cluster_pct(obj, clone_id, pal_cluster)

  ggplot(df, aes(x = "", y = pct, fill = cluster)) +
    geom_col(width = 0.85) +
    scale_fill_manual(values = pal_cluster, drop = FALSE) +
    scale_y_continuous(
      limits   = c(0, 100),
      breaks   = c(0, 50, 100),
      expand   = expansion(0, 0),
      position = "right"
    ) +
    labs(y = "% of cells") +
    theme_bw() +
    theme(
      axis.title.x       = element_blank(),
      axis.text.x        = element_blank(),
      axis.ticks.x       = element_blank(),
      axis.title.y.left  = element_blank(),
      axis.text.y.left   = element_blank(),
      axis.ticks.y.left  = element_blank(),
      axis.title.y.right = element_text(size = 6, angle = -90,
                                        margin = margin(l = 1)),
      axis.text.y.right  = element_text(size = 5, margin = margin(l = 1)),
      axis.ticks.y.right = element_line(linewidth = 0.2),
      axis.ticks.length  = unit(1.2, "pt"),
      panel.grid         = element_blank(),
      panel.border       = element_rect(linewidth = 0.3, fill = NA),
      legend.position    = "none",
      plot.margin        = margin(2, 2, 2, 0)
    )
}

make_panel_composite_v3 <- function(obj, clone_id, pal_cluster,
                                    color_by      = "specificity",
                                    pal_timepoint = NULL,
                                    lims          = NULL) {
  wrap_plots(
    make_panel_umap_v3(obj, clone_id, pal_cluster,
                       color_by      = color_by,
                       pal_timepoint = pal_timepoint,
                       lims          = lims,
                       tight_right   = TRUE),
    make_panel_bar_v3(obj, clone_id, pal_cluster),
    widths = c(1, BAR_REL_WIDTH)
  )
}

# ── 8. Standalone legend panel ─────────────────────────────────────────────
make_legend_panel_v3 <- function(pal_cluster, dot_pal = NULL, dot_title = NULL) {
  df_c <- data.frame(x = 1, y = seq_along(pal_cluster),
                     cluster = factor(names(pal_cluster),
                                      levels = names(pal_cluster)))

  p <- ggplot(df_c, aes(x, y, fill = cluster)) +
    geom_point(shape = 21, size = 3, alpha = 0, stroke = 0) +
    scale_fill_manual(values = pal_cluster, drop = FALSE,
                      name = "seurat_cluster") +
    guides(fill = guide_legend(order = 2,
                               override.aes = list(alpha = 1, stroke = 0.3)))

  if (!is.null(dot_pal)) {
    df_d <- data.frame(x = 1, y = seq_along(dot_pal),
                       grp = factor(names(dot_pal), levels = names(dot_pal)))
    p <- p +
      geom_point(data = df_d, aes(x, y, color = grp),
                 shape = 16, size = 3, alpha = 0, inherit.aes = FALSE) +
      scale_color_manual(values = dot_pal, drop = FALSE, name = dot_title) +
      guides(color = guide_legend(order = 1,
                                  override.aes = list(alpha = 1)))
  }

  p +
    theme_void() +
    theme(
      legend.position  = "left",
      legend.title     = element_text(size = 8, face = "bold"),
      legend.text      = element_text(size = 7),
      legend.key.size  = unit(0.35, "cm"),
      legend.spacing.y = unit(0.25, "cm"),
      plot.margin      = margin(2, 2, 2, 2)
    )
}

# ── 9. Clone ID discovery (numeric order, optional cap) ────────────────────
get_clone_ids_v3 <- function(obj, specificity, max_clone = NULL) {
  prefix  <- paste0("D0_", specificity, "_clone_")
  clone   <- get_clone_vec(obj)
  present <- clone[!is.na(clone) & str_starts(clone, prefix)]
  ids     <- unique(present)
  if (length(ids) == 0) return(character(0))
  num <- as.integer(str_remove(ids, prefix))
  ids <- ids[order(num)]
  num <- sort(num)
  if (!is.null(max_clone)) ids <- ids[num <= max_clone]
  ids
}

# ── 10. Grid export ────────────────────────────────────────────────────────
export_clone_grid_v3 <- function(obj, specificity, obj_label,
                                 mode      = c("cluster_colors",
                                               "cluster_colors_with_bar_plots",
                                               "with_bar_plots",
                                               "timepoint_colors_with_bar_plots"),
                                 ncol_grid = 5,
                                 max_clone = NULL,
                                 skinny    = FALSE,
                                 pal_cluster = NULL,
                                 lims        = NULL) {
  mode <- match.arg(mode)

  n_tp <- length(unique(na.omit(as.character(obj@meta.data[[TIME_COL]]))))
  if (SKIP_SINGLE_TP && mode == "timepoint_colors_with_bar_plots" && n_tp < 2) {
    message(sprintf("[%s / %s / %s] only 1 timepoint present — skipping",
                    obj_label, specificity, mode))
    return(invisible(NULL))
  }

  clone_ids <- get_clone_ids_v3(obj, specificity, max_clone)
  if (length(clone_ids) == 0) {
    message(sprintf("[%s / %s / %s] no clones found — skipping",
                    obj_label, specificity, mode))
    return(invisible(NULL))
  }

  message(sprintf("[%s / %s / %s] building %d panels",
                  obj_label, specificity, mode, length(clone_ids)))

  if (is.null(pal_cluster)) pal_cluster <- build_cluster_pal(obj)
  pal_tp <- build_timepoint_pal(obj)

  out_dir <- file.path(base_export, "clone_panels", obj_label, mode)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  if (mode == "cluster_colors") {
    panels   <- map(clone_ids, ~ make_panel_umap_v3(obj, .x, pal_cluster,
                                                    color_by = "cluster",
                                                    lims     = lims))
    legend_p <- make_legend_panel_v3(pal_cluster)
    cell_w   <- PANEL_W
    subtitle <- "dots colored by seurat_cluster"
    tag      <- "cluster"

  } else if (mode == "cluster_colors_with_bar_plots") {
    # dots and bar share the cluster palette -> one legend covers both
    panels   <- map(clone_ids, ~ make_panel_composite_v3(obj, .x, pal_cluster,
                                                         color_by = "cluster",
                                                         lims     = lims))
    legend_p <- make_legend_panel_v3(pal_cluster)
    cell_w   <- PANEL_W * (1 + BAR_REL_WIDTH)
    subtitle <- "dots by seurat_cluster; bar = % of clone's cells per seurat_cluster"
    tag      <- "cluster_bars"

  } else if (mode == "with_bar_plots") {
    panels   <- map(clone_ids, ~ make_panel_composite_v3(obj, .x, pal_cluster,
                                                         color_by = "specificity",
                                                         lims     = lims))
    legend_p <- make_legend_panel_v3(pal_cluster, spec_colors, "Specificity")
    cell_w   <- PANEL_W * (1 + BAR_REL_WIDTH)
    subtitle <- "dots by specificity; bar = % of clone's cells per seurat_cluster"
    tag      <- "bars"

  } else {
    panels   <- map(clone_ids, ~ make_panel_composite_v3(obj, .x, pal_cluster,
                                                         color_by      = "timepoint",
                                                         pal_timepoint = pal_tp,
                                                         lims          = lims))
    legend_p <- make_legend_panel_v3(pal_cluster, pal_tp, "Timepoint")
    cell_w   <- PANEL_W * (1 + BAR_REL_WIDTH)
    subtitle <- "dots by timepoint; bar = % of clone's cells per seurat_cluster"
    tag      <- "timepoint_bars"
  }

  nrow_grid <- ceiling(length(panels) / ncol_grid)

  core <- wrap_plots(panels, ncol = ncol_grid)
  grid <- wrap_plots(core, legend_p,
                     widths = c(ncol_grid * cell_w, LEGEND_W)) +
    plot_annotation(
      title    = sprintf("%s — %s clones (n = %d)",
                         obj_label, specificity, length(clone_ids)),
      subtitle = subtitle,
      theme    = theme(
        plot.title    = element_text(size = 13, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 9, hjust = 0.5)
      )
    )

  suffix   <- paste0("_clone_panels", if (skinny) "_skinny" else "", "_", tag)
  out_path <- file.path(out_dir,
                        sprintf("%s_%s%s.jpeg", obj_label, specificity, suffix))

  save_jpeg(grid, out_path,
            width  = ncol_grid * cell_w + LEGEND_W,
            height = nrow_grid * PANEL_H)

  message(sprintf("  -> %s", out_path))

  if (WRITE_CSV && mode == "with_bar_plots") {
    tbl <- map_dfr(clone_ids, function(cid) {
      clone_cluster_pct(obj, cid, pal_cluster) %>%
        mutate(object = obj_label, specificity = specificity, clone = cid,
               .before = 1)
    })
    write_csv(tbl, file.path(out_dir,
                             sprintf("%s_%s_clone_cluster_percent.csv",
                                     obj_label, specificity)))
  }

  invisible(out_path)
}

# ── 11. Object registry ────────────────────────────────────────────────────
full_pal  <- build_cluster_pal(D0_data_CMV)
full_lims <- umap_limits(D0_data_CMV)

registry_v3 <- list(
  list(label = "full",
       obj   = D0_data_CMV,
       pal   = full_pal,
       lims  = full_lims),

  list(label = "V1_3",
       obj   = obj_V1_3,
       pal   = NULL,              # its own re-clustered levels
       lims  = NULL),             # its own UMAP range

  list(label = "V4",
       obj   = obj_V4,
       pal   = NULL,
       lims  = NULL),

  list(label = "V1_3_original_clustering",
       obj   = obj_V1_3_orig,
       pal   = full_pal,
       lims  = full_lims),

  list(label = "V4_original_clustering",
       obj   = obj_V4_orig,
       pal   = full_pal,
       lims  = full_lims)
)

grid_specs_v3 <- list(
  pp65 = list(ncol_grid = 5, max_clone = 20, skinny = FALSE),
  IE1  = list(ncol_grid = 3, max_clone = 7,  skinny = TRUE)
)

modes_v3 <- c("cluster_colors",
              "cluster_colors_with_bar_plots",
              "with_bar_plots",
              "timepoint_colors_with_bar_plots")

# ── 12. Run: objects x specificity x mode ──────────────────────────────────
for (entry in registry_v3) {
  for (spec in names(grid_specs_v3)) {
    sp <- grid_specs_v3[[spec]]
    for (md in modes_v3) {
      export_clone_grid_v3(entry$obj, spec, entry$label,
                           mode        = md,
                           ncol_grid   = sp$ncol_grid,
                           max_clone   = sp$max_clone,
                           skinny      = sp$skinny,
                           pal_cluster = entry$pal,
                           lims        = entry$lims)
    }
  }
}

message("\nDone — up to 40 figures (5 objects x 2 specificities x 4 modes) under ",
        file.path(base_export, "clone_panels"))




# ============================================================================
#  Donor 0 — manuscript figures with seurat clusters 6 and 7 excluded
#  All percentages / z-scores are recomputed over clusters 0-5 only.
#  DEG heatmap is produced in two versions:
#    (a) original markers  — rankings from the FULL object, 6/7 rows dropped
#    (b) recalculated      — FindAllMarkers re-run on the 6/7-free subset
# ============================================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)

# ── 0. Setup ────────────────────────────────────────────────────────────────
out_dir <- "Donor_0/exports/#_clusters_6_7_excluded_figures"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

KEEP_CLUSTERS <- c("0", "1", "2", "3", "4", "5")

cluster_colors_all <- c(
  "0" = "#C49A00",
  "1" = "#53B400",
  "2" = "#8B0000",
  "3" = "#1976D2",
  "4" = "#F8766D",
  "5" = "#A58AFF",
  "6" = "#FFDB58",
  "7" = "#00C094"
)
cluster_colors <- cluster_colors_all[KEEP_CLUSTERS]   # legend drops 6 / 7

# Subset by cell name (robust to seurat_clusters being a factor), then drop the
# unused factor levels so nothing downstream re-introduces empty 6 / 7 slots.
cells_keep <- colnames(D0_data_CMV)[
  as.character(D0_data_CMV$seurat_clusters) %in% KEEP_CLUSTERS
]
D0_sub <- subset(D0_data_CMV, cells = cells_keep)
D0_sub$seurat_clusters <- factor(as.character(D0_sub$seurat_clusters),
                                 levels = KEEP_CLUSTERS)
Idents(D0_sub) <- "seurat_clusters"

message(sprintf("Cells retained: %d of %d (%.1f%%)",
                ncol(D0_sub), ncol(D0_data_CMV),
                100 * ncol(D0_sub) / ncol(D0_data_CMV)))
print(table(D0_sub$seurat_clusters, D0_sub$timepoint))


# ── 1. All cells — cluster composition per timepoint ────────────────────────
meta <- D0_sub@meta.data

cluster_counts <- meta %>%
  group_by(timepoint, seurat_clusters) %>%
  summarise(n_cells = n(), .groups = "drop")

timepoint_totals <- meta %>%                     # denominator = clusters 0-5 only
  group_by(timepoint) %>%
  summarise(total_cells = n(), .groups = "drop")

plot_df_all <- cluster_counts %>%
  left_join(timepoint_totals, by = "timepoint") %>%
  mutate(percent = 100 * n_cells / total_cells)

p_all <- ggplot(plot_df_all, aes(x = timepoint, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cluster_colors, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ylab("Percent of cells") +
  xlab("Timepoint") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(
  filename = file.path(out_dir, "cluster_composition_over_time.png"),
  plot = p_all, width = 6, height = 5, dpi = 300, bg = "white"
)


# ── 2. Faceted by specificity ───────────────────────────────────────────────
meta_sub <- meta %>% filter(specificity %in% c("pp65", "IE1"))

cluster_counts <- meta_sub %>%
  group_by(specificity, timepoint, seurat_clusters) %>%
  summarise(n_cells = n(), .groups = "drop")

totals <- meta_sub %>%
  group_by(specificity, timepoint) %>%
  summarise(total_cells = n(), .groups = "drop")

plot_df <- cluster_counts %>%
  left_join(totals, by = c("specificity", "timepoint")) %>%
  mutate(percent = 100 * n_cells / total_cells)

p_spec_facet <- ggplot(plot_df,
                       aes(x = timepoint, y = percent, fill = seurat_clusters)) +
  geom_bar(stat = "identity", width = 0.8) +
  scale_fill_manual(values = cluster_colors, drop = FALSE) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  facet_wrap(~ specificity, nrow = 1) +
  ylab("Percent of cells") +
  xlab("Visit") +
  theme_classic() +
  theme(
    axis.text.x      = element_text(angle = 0),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold")
  )

ggsave(
  filename = file.path(out_dir, "cluster_composition_over_time_by_specificity.png"),
  plot = p_spec_facet, width = 8, height = 5, dpi = 300, bg = "white"
)


# ── 3. DEG heatmaps — two versions ──────────────────────────────────────────
DefaultAssay(D0_sub) <- "RNA"

PROTEIN_CODING_ONLY <- TRUE   # matches the filter used for the original figure

# --- protein-coding reference (built once, cached in the session) ----------
if (PROTEIN_CODING_ONLY && !exists("coding_genes")) {
  if (requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
    library(EnsDb.Hsapiens.v86)
    coding_genes <- genes(EnsDb.Hsapiens.v86,
                          filter = GeneBiotypeFilter("protein_coding"))$gene_name
    detach("package:EnsDb.Hsapiens.v86", unload = TRUE)
    detach("package:ensembldb", unload = TRUE)
  } else {
    library(biomaRt)
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    coding_genes <- getBM(
      attributes = c("hgnc_symbol", "gene_biotype"),
      filters    = "hgnc_symbol",
      values     = rownames(D0_sub[["RNA"]]),
      mart       = mart
    ) %>%
      filter(gene_biotype == "protein_coding") %>%
      pull(hgnc_symbol)
  }
  coding_genes <- unique(coding_genes[!is.na(coding_genes) & coding_genes != ""])
  message(sprintf("Protein-coding reference: %d symbols", length(coding_genes)))
}

# --- (a) markers from the full object -------------------------------------
if (!exists("d0_cluster_markers")) {
  tmp <- D0_data_CMV
  DefaultAssay(tmp) <- "RNA"
  Idents(tmp) <- "seurat_clusters"
  d0_cluster_markers <- FindAllMarkers(tmp, only.pos = TRUE, logfc.threshold = 0.25)
  rm(tmp)
}

# --- (b) markers recomputed on clusters 0-5 only ---------------------------
Idents(D0_sub) <- "seurat_clusters"
d0_cluster_markers_sub <- FindAllMarkers(D0_sub, only.pos = TRUE,
                                         logfc.threshold = 0.25)

# Shared selection: pct.1 filter -> protein-coding -> top 10 per cluster.
# Filters on the `gene` column, NOT rownames — rownames are made unique, so a
# gene that marks two clusters becomes "GENE.1" and would be dropped silently.
top10_by_cluster <- function(markers, label = "") {
  m <- markers %>%
    filter(pct.1 > 0.1, as.character(cluster) %in% KEEP_CLUSTERS)

  if (PROTEIN_CODING_ONLY) {
    n_before <- nrow(m)
    m <- m %>% filter(gene %in% coding_genes)
    message(sprintf("[%s] protein-coding filter: %d -> %d marker rows",
                    label, n_before, nrow(m)))
  }

  m %>%
    mutate(cluster = factor(as.character(cluster), levels = KEEP_CLUSTERS)) %>%
    group_by(cluster) %>%
    top_n(n = 10, wt = avg_log2FC) %>%
    ungroup() %>%
    arrange(cluster, desc(avg_log2FC))
}

top_orig   <- top10_by_cluster(d0_cluster_markers,     "original markers")
top_recalc <- top10_by_cluster(d0_cluster_markers_sub, "recalculated markers")

# Scale once over the union — per-gene z-scoring is independent of which other
# genes are in the set, so both heatmaps get identical scaling for shared genes.
genes_union <- unique(c(top_orig$gene, top_recalc$gene))
D0_sub <- ScaleData(D0_sub, features = genes_union, verbose = FALSE)

save_deg_heatmap <- function(top_tbl, filename, label) {
  genes   <- top_tbl$gene                 # DoHeatmap de-duplicates internally
  n_genes <- length(unique(genes))
  message(sprintf("[%s] %d unique genes from %d marker rows (clusters 0-5)",
                  label, n_genes, nrow(top_tbl)))

  p <- DoHeatmap(
    D0_sub,
    features     = genes,
    group.colors = unname(cluster_colors)
  ) +
    scale_fill_gradientn(colors = c("blue", "white", "red")) +
    theme(axis.text.y = element_text(size = 18))

  ggsave(
    filename = file.path(out_dir, filename),
    plot = p,
    width = 10, height = max(8, 20 * n_genes / 80), dpi = 300, bg = "white"
  )

  write.csv(top_tbl,
            file.path(out_dir, sub("\\.png$", "_genes.csv", filename)),
            row.names = FALSE)
}

save_deg_heatmap(top_orig,
                 "D0_CMV_DEGs_per_cluster_original_markers.png",
                 "original markers")
save_deg_heatmap(top_recalc,
                 "D0_CMV_DEGs_per_cluster_recalculated.png",
                 "recalculated markers")

# Quick comparison: which genes changed between the two versions
gene_diff <- tibble(
  gene = union(top_orig$gene, top_recalc$gene)
) %>%
  mutate(
    in_original     = gene %in% top_orig$gene,
    in_recalculated = gene %in% top_recalc$gene,
    status = case_when(
      in_original & in_recalculated ~ "shared",
      in_original                   ~ "original only",
      TRUE                          ~ "recalculated only"
    )
  ) %>%
  arrange(status, gene)

print(count(gene_diff, status))
write.csv(gene_diff,
          file.path(out_dir, "D0_CMV_DEGs_per_cluster_version_comparison.csv"),
          row.names = FALSE)
          

# ── 4. CITEseq bubble plot — z-scored across clusters 0-5 only ──────────────
d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_sub[["FB"]])
)

p_d0      <- DotPlot(D0_sub, features = d0_markers,
                     group.by = "seurat_clusters", assay = "FB")
plot_data <- p_d0$data

plot_data$id            <- factor(as.character(plot_data$id), levels = KEEP_CLUSTERS)
plot_data$features.plot <- factor(plot_data$features.plot, levels = rev(d0_markers))

message(sprintf("avg.exp.scaled range (clusters 0-5): %.2f to %.2f",
                min(plot_data$avg.exp.scaled), max(plot_data$avg.exp.scaled)))

p_bubble <- ggplot(plot_data, aes(x = id, y = features.plot)) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6)) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red")) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 12)
  ) +
  labs(x = "Seurat cluster", y = "CITEseq marker")

ggsave(
  filename = file.path(out_dir, "D0_CMV_FB_bubbleplot_by_seurat_cluster_short.png"),
  plot = p_bubble, width = 6.5, height = 3.5, dpi = 300, bg = "white"
) # height was 4.5 in non-short version

message("Done — five figures written to ", out_dir)



# ---- D0: FB markers x specificity x visit block (V1-3 vs V4) ----

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

message("Timepoint values:")
print(table(D0_data_CMV$timepoint, useNA = "ifany"))

tp        <- as.character(D0_data_CMV$timepoint)
spec_chr  <- as.character(D0_data_CMV$specificity)

visit_blk <- ifelse(tp %in% c("Visit 1", "Visit 2", "Visit 3"), "Visits 1-3",
                    ifelse(tp == "Visit 4", "Visit 4", NA))

D0_data_CMV$spec_visit <- ifelse(
  !is.na(spec_chr) & spec_chr %in% c("pp65", "IE1") & !is.na(visit_blk),
  paste0(spec_chr, " (", visit_blk, ")"), NA
)

grp_levels <- c("pp65 (Visits 1-3)", "pp65 (Visit 4)",
                "IE1 (Visits 1-3)",  "IE1 (Visit 4)")

keep_cells <- colnames(D0_data_CMV)[D0_data_CMV$spec_visit %in% grp_levels]
d0_sv <- subset(D0_data_CMV, cells = keep_cells)
d0_sv$spec_visit <- factor(d0_sv$spec_visit, levels = grp_levels)

grp_n <- table(d0_sv$spec_visit)
message("Cells per group:"); print(grp_n)

# scale = FALSE to stay comparable with the 2-row specificity figure
p_sv <- DotPlot(d0_sv, features = d0_markers,
                group.by = "spec_visit", assay = "FB", scale = FALSE)
plot_data <- p_sv$data

plot_data$features.plot <- factor(plot_data$features.plot, levels = d0_markers)
plot_data$id <- factor(as.character(plot_data$id), levels = rev(grp_levels))

y_levels  <- levels(plot_data$id)
y_lab_col <- ifelse(grepl("^pp65", y_levels), "#1976D2", "#8B0000")
y_lab_txt <- setNames(
  paste0(y_levels, "\nn=", as.integer(grp_n[y_levels])), y_levels
)

split_at <- sum(grepl("^pp65", y_levels)) + 0.5

p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_hline(yintercept = split_at, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 8), name = "% cells\nexpressing") +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red"),
                        name = "Mean CLR\nexpression\n(log1p)") +
  scale_y_discrete(labels = y_lab_txt) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y  = element_text(size = 9.5, face = "bold", color = y_lab_col),
        axis.title   = element_text(size = 12),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()) +
  labs(x = "CITEseq marker", y = NULL,
       title = "D0 CMV-specific CD8 T cells: CITEseq markers by specificity and visit")

dir.create("corrected_FB_exports/D0", showWarnings = FALSE, recursive = TRUE)
ggsave("corrected_FB_exports/D0/D0_CMV_FB_bubbleplot_by_specificity_visit_short.png",
       p, width = 8, height = 3.8, dpi = 300)

# ---- D0: FB markers x specificity x visit block (V1-3 vs V4), z-scored ----

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

message("Timepoint values:")
print(table(D0_data_CMV$timepoint, useNA = "ifany"))

tp       <- as.character(D0_data_CMV$timepoint)
spec_chr <- as.character(D0_data_CMV$specificity)

visit_blk <- ifelse(tp %in% c("Visit 1", "Visit 2", "Visit 3"), "Visits 1-3",
                    ifelse(tp == "Visit 4", "Visit 4", NA))

D0_data_CMV$spec_visit <- ifelse(
  !is.na(spec_chr) & spec_chr %in% c("pp65", "IE1") & !is.na(visit_blk),
  paste0(spec_chr, " (", visit_blk, ")"), NA
)

grp_levels <- c("pp65 (Visits 1-3)", "pp65 (Visit 4)",
                "IE1 (Visits 1-3)",  "IE1 (Visit 4)")

keep_cells <- colnames(D0_data_CMV)[D0_data_CMV$spec_visit %in% grp_levels]
d0_sv <- subset(D0_data_CMV, cells = keep_cells)
d0_sv$spec_visit <- factor(d0_sv$spec_visit, levels = grp_levels)

grp_n <- table(d0_sv$spec_visit)
message("Cells per group:"); print(grp_n)

# scale = TRUE -> avg.exp.scaled is a z-score of each marker across the 4 rows
p_sv <- DotPlot(d0_sv, features = d0_markers,
                group.by = "spec_visit", assay = "FB", scale = TRUE)
plot_data <- p_sv$data

plot_data$features.plot <- factor(plot_data$features.plot, levels = d0_markers)
plot_data$id <- factor(as.character(plot_data$id), levels = rev(grp_levels))

message(sprintf("z-score range: %.2f to %.2f",
                min(plot_data$avg.exp.scaled), max(plot_data$avg.exp.scaled)))

# symmetric limits so grey sits exactly at z = 0 (the 4-row mean)
z_lim <- max(abs(plot_data$avg.exp.scaled), na.rm = TRUE)

y_levels  <- levels(plot_data$id)
y_lab_col <- ifelse(grepl("^pp65", y_levels), "#1976D2", "#8B0000")
y_lab_txt <- setNames(
  paste0(y_levels, "\nn=", as.integer(grp_n[y_levels])), y_levels
)

split_at <- sum(grepl("^pp65", y_levels)) + 0.5

p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_hline(yintercept = split_at, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 8), name = "% cells\nexpressing") +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red"),
                        limits = c(-z_lim, z_lim),
                        name = "Scaled mean\nexpression\n(z across rows)") +
  scale_y_discrete(labels = y_lab_txt) +
  theme_bw() +
  theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y  = element_text(size = 9.5, face = "bold", color = y_lab_col),
        axis.title   = element_text(size = 12),
        plot.title   = element_text(face = "bold", hjust = 0.5),
        panel.grid.minor = element_blank()) +
  labs(x = "CITEseq marker", y = NULL,
       title = "D0 CMV-specific CD8 T cells: CITEseq markers by specificity and visit")

dir.create("corrected_FB_exports/D0", showWarnings = FALSE, recursive = TRUE)
ggsave("corrected_FB_exports/D0/D0_CMV_FB_bubbleplot_by_specificity_visit_short_zscore.png",
       p, width = 8, height = 3.8, dpi = 300)




# ============================================================================
#  D0 — CITEseq bubble plot by major clone (IE1 1-7, pp65 1-20)
#  Labels: "pp65 Clone 1" style, no per-clone n
# ============================================================================

library(Seurat)
library(dplyr)
library(ggplot2)
library(stringr)

# ── Settings ────────────────────────────────────────────────────────────────
# FALSE = faithful remake of the original figure (all cells, clusters 0-7)
# TRUE  = restrict to clusters 0-5, so the scaled expression matches the other
#         6/7-excluded panels
EXCLUDE_67 <- FALSE

if (EXCLUDE_67) {
  obj_src  <- D0_sub                                # built earlier in the session
  out_path <- "Donor_0/exports/#_clusters_6_7_excluded_figures/D0_CMV_FB_bubbleplot_by_clone.png"
} else {
  obj_src  <- D0_data_CMV
  out_path <- "corrected_FB_exports/D0/D0_CMV_FB_bubbleplot_by_clone.png"
}
dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(obj_src[["FB"]])
)

# ── Clone ordering: numeric, not lexicographic (clone_2 before clone_10) ────
ie1_clones  <- paste0("D0_IE1_clone_",  1:7)
pp65_clones <- paste0("D0_pp65_clone_", 1:20)
clone_levels_all <- c(ie1_clones, pp65_clones)

obs_clones <- unique(as.character(obj_src$Global_Clone_ID))
missing    <- setdiff(clone_levels_all, obs_clones)
if (length(missing)) message("Not found in object: ", paste(missing, collapse = ", "))
clone_levels <- clone_levels_all[clone_levels_all %in% obs_clones]

keep_cells <- colnames(obj_src)[
  !is.na(obj_src$Global_Clone_ID) &
    as.character(obj_src$Global_Clone_ID) %in% clone_levels
]
d0_clones <- subset(obj_src, cells = keep_cells)
d0_clones$Global_Clone_ID <- factor(as.character(d0_clones$Global_Clone_ID),
                                    levels = clone_levels)

clone_n <- table(d0_clones$Global_Clone_ID)
message("Cells per clone:"); print(clone_n)
message("Total cells retained: ", ncol(d0_clones))

# ── Dot plot ────────────────────────────────────────────────────────────────
p_clone <- DotPlot(d0_clones, features = d0_markers,
                   group.by = "Global_Clone_ID", assay = "FB", scale = TRUE)
plot_data <- p_clone$data

# Markers along x; clones as rows, IE1 block on top, clone_1 at top of each block
plot_data$features.plot <- factor(plot_data$features.plot, levels = d0_markers)
plot_data$id <- factor(as.character(plot_data$id), levels = rev(clone_levels))

# Display labels only — factor LEVELS stay as the raw Global_Clone_IDs so the
# specificity color logic and the block-split line keep working off them.
fmt_clone_label <- function(ids) {
  m <- str_match(ids, "^D0_(pp65|IE1)_clone_(\\d+)$")
  ifelse(is.na(m[, 2]), ids, sprintf("%s Clone %s", m[, 2], m[, 3]))
}

y_levels  <- levels(plot_data$id)
y_lab_col <- ifelse(grepl("_pp65_", y_levels), "#1976D2", "#8B0000")
y_lab_txt <- setNames(fmt_clone_label(y_levels), y_levels)

# Divider between the pp65 block and the IE1 block
split_at <- sum(grepl("_pp65_", y_levels)) + 0.5

p <- ggplot(plot_data, aes(x = features.plot, y = id)) +
  geom_hline(yintercept = split_at, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) +
  scale_size(range = c(0, 6), name = "% cells\nexpressing") +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red"),
                        name = "Scaled mean\nexpression") +
  scale_y_discrete(labels = y_lab_txt) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 9, face = "bold", color = y_lab_col),
    axis.title       = element_text(size = 12),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "CITEseq marker", y = "Clone",
       title = "D0 CMV-specific CD8 T cells: CITEseq markers by major clone")

ggsave(out_path, p, width = 8, height = 10.5, dpi = 300, bg = "white")
message("Wrote ", out_path)



#
# Clonotype highlight plots — colored by metacluster
#

# Output directory
output_dir <- "exports/integrated_v5/Clonotype_plots/Clonotype_plots_colored_by_cluster/"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Standard metacluster palette / order
metacluster_order <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")
metacluster_colors <- c(
  "Tcm"     = "#F8766D",
  "Tem"     = "#A3A500",
  "T-CMC1"  = "#00BF7D",
  "NK-like" = "#00B0F6",
  "ISG+"    = "#E76BF3"
)

BG_COLOR <- "grey92"   # background (non-highlighted) point color
BG_SIZE  <- 0.5        # background point size
HL_SIZE  <- 1.3       # matches sizes.highlight in the original plots

# Pull UMAP coords + metadata once, rather than per plot
umap_df <- as.data.frame(Embeddings(mdata_CMV, reduction = "umap"))
colnames(umap_df)[1:2] <- c("UMAP_1", "UMAP_2")
umap_df$Global_Clone_ID <- as.character(mdata_CMV$Global_Clone_ID)
umap_df$metacluster <- factor(as.character(mdata_CMV$metacluster),
                              levels = metacluster_order)

# Single clone panel: light grey backdrop + highlighted cells by metacluster
make_clone_plot <- function(clonotype_id, plot_title) {

  hl <- !is.na(umap_df$Global_Clone_ID) & umap_df$Global_Clone_ID == clonotype_id
  if (!any(hl)) return(NULL)

  bg <- umap_df[!hl, ]
  fg <- umap_df[ hl, ]

  ggplot() +
    geom_point(data = bg, aes(UMAP_1, UMAP_2),
               color = BG_COLOR, size = BG_SIZE) +
    geom_point(data = fg, aes(UMAP_1, UMAP_2, color = metacluster),
               size = HL_SIZE) +
    # drop = FALSE keeps colors fixed per metacluster across every plot;
    # na.value catches highlighted cells with no metacluster assignment
    scale_color_manual(values = metacluster_colors,
                       drop = FALSE, na.value = "grey40") +
    ggtitle(plot_title) +
    theme_cowplot() +
    theme(
      legend.position  = "none",
      # explicit white fills so the PNG is opaque even on older ggplot2
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )
}

# Loop through donors and clones
for (x in 1:5) {
  for (spec in c("IE1", "pp65")) {
    for (i in 1:20) {

      clonotype_id <- paste0("D", x, "_", spec, "-clone-", i)
      plot_title   <- paste0("Donor ", x, " - ", spec, " clone ", i)

      p <- make_clone_plot(clonotype_id, plot_title)

      if (!is.null(p)) {
        ggsave(
          filename = file.path(output_dir, paste0(clonotype_id, ".png")),
          plot = p, width = 6, height = 6, dpi = 300
        )
      }
    }
  }
}

# ---- Standalone metacluster legend (plots themselves have none) ----
legend_src <- ggplot(
  data.frame(x = 1, y = seq_along(metacluster_order),
             metacluster = factor(metacluster_order, levels = metacluster_order)),
  aes(x, y, color = metacluster)
) +
  geom_point(size = 4) +
  scale_color_manual(values = metacluster_colors, name = "Metacluster") +
  theme_void() +
  theme(legend.position = "right")

ggsave(
  filename = file.path(output_dir, "metacluster_color_legend.png"),
  plot = as.ggplot(cowplot::get_legend(legend_src)),
  width = 2.5, height = 2, dpi = 300, bg = "white"
)




# ---- D0: CITEseq marker violins x specificity x visit block ----

library(tidyr)

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

tp       <- as.character(D0_data_CMV$timepoint)
spec_chr <- as.character(D0_data_CMV$specificity)

visit_blk <- ifelse(tp %in% c("Visit 1", "Visit 2", "Visit 3"), "Visits 1-3",
                    ifelse(tp == "Visit 4", "Visit 4", NA))

D0_data_CMV$spec_visit <- ifelse(
  !is.na(spec_chr) & spec_chr %in% c("pp65", "IE1") & !is.na(visit_blk),
  paste0(spec_chr, " (", visit_blk, ")"), NA
)

grp_levels <- c("pp65 (Visits 1-3)", "pp65 (Visit 4)",
                "IE1 (Visits 1-3)",  "IE1 (Visit 4)")

keep_cells <- colnames(D0_data_CMV)[D0_data_CMV$spec_visit %in% grp_levels]
d0_sv <- subset(D0_data_CMV, cells = keep_cells)
d0_sv$spec_visit <- factor(d0_sv$spec_visit, levels = grp_levels)

grp_n <- table(d0_sv$spec_visit)
message("Cells per group:"); print(grp_n)

# CLR values straight from the FB assay (avoids the RNA/FB name collisions)
expr <- GetAssayData(d0_sv, assay = "FB", layer = "data")[d0_markers, , drop = FALSE]
grp_vec <- setNames(as.character(d0_sv$spec_visit), colnames(d0_sv))

vln_df <- as.data.frame(t(as.matrix(expr))) %>%
  tibble::rownames_to_column("cell") %>%
  mutate(group = grp_vec[cell]) %>%
  pivot_longer(all_of(d0_markers), names_to = "marker", values_to = "expr") %>%
  mutate(
    marker = factor(marker, levels = d0_markers),
    group  = factor(group,  levels = grp_levels),
    spec   = factor(ifelse(grepl("^pp65", group), "pp65", "IE1"),
                    levels = c("pp65", "IE1"))
  )

spec_cols <- c(pp65 = "#1976D2", IE1 = "#8B0000")

# x labels carry the group n so the reader can weight the V4 rows appropriately
x_labs <- setNames(
  sprintf("%s\nn = %d", grp_levels, as.integer(grp_n[grp_levels])),
  grp_levels
)

p <- ggplot(vln_df, aes(x = group, y = expr, fill = spec)) +
  geom_violin(scale = "width", trim = TRUE, color = "black",
              linewidth = 0.3, alpha = 0.8) +
  geom_boxplot(width = 0.12, outlier.shape = NA,
               fill = "white", alpha = 0.7, linewidth = 0.3) +
  geom_vline(xintercept = 2.5, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  facet_grid(marker ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = spec_cols, name = "Specificity") +
  scale_x_discrete(labels = x_labs) +
  labs(x = NULL, y = "CLR expression",
       title = "D0 CMV-specific CD8 T cells: CITEseq markers by specificity and visit") +
  theme_bw(base_size = 12) +
  theme(
    strip.placement  = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
    axis.text.x      = element_text(size = 9.5),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    legend.position  = "right"
  )

dir.create("corrected_FB_exports/D0", showWarnings = FALSE, recursive = TRUE)
ggsave("corrected_FB_exports/D0/D0_CMV_FB_violin_by_specificity_visit.png",
       p, width = 6.5, height = 1.4 * length(d0_markers) + 1.2, dpi = 300)



# ---- D0: CITEseq marker violins, 8 markers x 4 groups grid ----

library(tidyr)

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

tp       <- as.character(D0_data_CMV$timepoint)
spec_chr <- as.character(D0_data_CMV$specificity)

visit_blk <- ifelse(tp %in% c("Visit 1", "Visit 2", "Visit 3"), "Visits 1-3",
                    ifelse(tp == "Visit 4", "Visit 4", NA))

D0_data_CMV$spec_visit <- ifelse(
  !is.na(spec_chr) & spec_chr %in% c("pp65", "IE1") & !is.na(visit_blk),
  paste0(spec_chr, " (", visit_blk, ")"), NA
)

grp_levels <- c("pp65 (Visits 1-3)", "pp65 (Visit 4)",
                "IE1 (Visits 1-3)",  "IE1 (Visit 4)")

keep_cells <- colnames(D0_data_CMV)[D0_data_CMV$spec_visit %in% grp_levels]
d0_sv <- subset(D0_data_CMV, cells = keep_cells)
d0_sv$spec_visit <- factor(d0_sv$spec_visit, levels = grp_levels)

grp_n <- table(d0_sv$spec_visit)
message("Cells per group:"); print(grp_n)

# CLR values straight from the FB assay (avoids the RNA/FB name collisions)
expr <- GetAssayData(d0_sv, assay = "FB", layer = "data")[d0_markers, , drop = FALSE]
grp_vec <- setNames(as.character(d0_sv$spec_visit), colnames(d0_sv))

# column strip labels carry the group n
grp_labs <- setNames(
  sprintf("%s\nn = %d", grp_levels, as.integer(grp_n[grp_levels])),
  grp_levels
)

vln_df <- as.data.frame(t(as.matrix(expr))) %>%
  tibble::rownames_to_column("cell") %>%
  mutate(group = grp_vec[cell]) %>%
  pivot_longer(all_of(d0_markers), names_to = "marker", values_to = "expr") %>%
  mutate(
    marker = factor(marker, levels = d0_markers),
    spec   = factor(ifelse(grepl("^pp65", group), "pp65", "IE1"),
                    levels = c("pp65", "IE1")),
    # facet on the labelled version so strips show the n
    group  = factor(grp_labs[group], levels = unname(grp_labs[grp_levels]))
  )

spec_cols <- c(pp65 = "#1976D2", IE1 = "#8B0000")

p <- ggplot(vln_df, aes(x = spec, y = expr, fill = spec)) +
  geom_violin(scale = "width", trim = TRUE, color = "black",
              linewidth = 0.3, alpha = 0.8, width = 0.85) +
  geom_boxplot(width = 0.14, outlier.shape = NA,
               fill = "white", alpha = 0.7, linewidth = 0.3) +
  facet_grid(marker ~ group, scales = "free_y", switch = "y") +
  scale_fill_manual(values = spec_cols, name = "Specificity") +
  labs(x = NULL, y = "CLR expression",
       title = "D0 CMV-specific CD8 T cells: CITEseq markers by specificity and visit") +
  theme_bw(base_size = 11) +
  theme(
    strip.placement   = "outside",
    strip.background  = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold", size = 10),
    strip.text.x      = element_text(face = "bold", size = 9.5,
                                     margin = margin(b = 3)),
    axis.text.x       = element_blank(),
    axis.ticks.x      = element_blank(),
    axis.text.y       = element_text(size = 8),
    panel.grid.minor  = element_blank(),
    panel.spacing.x   = unit(0.4, "lines"),
    panel.spacing.y   = unit(0.3, "lines"),
    plot.title        = element_text(face = "bold", hjust = 0.5),
    legend.position   = "right"
  )

dir.create("corrected_FB_exports/D0", showWarnings = FALSE, recursive = TRUE)
ggsave("corrected_FB_exports/D0/D0_CMV_FB_violin_by_specificity_visit_grid.png",
       p, width = 9, height = 1.15 * length(d0_markers) + 1.4, dpi = 300)



# ---- D0: CITEseq marker violins, 2 blocks of 4 markers side by side ----

library(tidyr)
library(patchwork)

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

tp       <- as.character(D0_data_CMV$timepoint)
spec_chr <- as.character(D0_data_CMV$specificity)

visit_blk <- ifelse(tp %in% c("Visit 1", "Visit 2", "Visit 3"), "Visits 1-3",
                    ifelse(tp == "Visit 4", "Visit 4", NA))

D0_data_CMV$spec_visit <- ifelse(
  !is.na(spec_chr) & spec_chr %in% c("pp65", "IE1") & !is.na(visit_blk),
  paste0(spec_chr, " (", visit_blk, ")"), NA
)

grp_levels <- c("pp65 (Visits 1-3)", "pp65 (Visit 4)",
                "IE1 (Visits 1-3)",  "IE1 (Visit 4)")

keep_cells <- colnames(D0_data_CMV)[D0_data_CMV$spec_visit %in% grp_levels]
d0_sv <- subset(D0_data_CMV, cells = keep_cells)
d0_sv$spec_visit <- factor(d0_sv$spec_visit, levels = grp_levels)

grp_n <- table(d0_sv$spec_visit)
message("Cells per group:"); print(grp_n)

# CLR values straight from the FB assay (avoids the RNA/FB name collisions)
expr <- GetAssayData(d0_sv, assay = "FB", layer = "data")[d0_markers, , drop = FALSE]
grp_vec <- setNames(as.character(d0_sv$spec_visit), colnames(d0_sv))

vln_df <- as.data.frame(t(as.matrix(expr))) %>%
  tibble::rownames_to_column("cell") %>%
  mutate(group = grp_vec[cell]) %>%
  pivot_longer(all_of(d0_markers), names_to = "marker", values_to = "expr") %>%
  mutate(
    group = factor(group, levels = grp_levels),
    spec  = factor(ifelse(grepl("^pp65", group), "pp65", "IE1"),
                   levels = c("pp65", "IE1"))
  )

spec_cols <- c(pp65 = "#1976D2", IE1 = "#8B0000")

# x labels carry the group n so the reader can weight the V4 columns
x_labs <- setNames(
  sprintf("%s\nn = %d", grp_levels, as.integer(grp_n[grp_levels])),
  grp_levels
)

# split into two blocks; ceiling keeps the left block fuller if length is odd
half   <- ceiling(length(d0_markers) / 2)
blocks <- list(d0_markers[seq_len(half)],
               d0_markers[seq(half + 1, length(d0_markers))])

make_block <- function(markers) {
  d <- vln_df %>%
    filter(marker %in% markers) %>%
    mutate(marker = factor(marker, levels = markers))

  ggplot(d, aes(x = group, y = expr, fill = spec)) +
    geom_violin(scale = "width", trim = TRUE, color = "black",
                linewidth = 0.3, alpha = 0.8) +
    geom_boxplot(width = 0.12, outlier.shape = NA,
                 fill = "white", alpha = 0.7, linewidth = 0.3) +
    geom_vline(xintercept = 2.5, linetype = "dashed",
               color = "grey40", linewidth = 0.4) +
    facet_grid(marker ~ ., scales = "free_y", switch = "y") +
    scale_fill_manual(values = spec_cols, name = "Specificity") +
    scale_x_discrete(labels = x_labs) +
    labs(x = NULL, y = "CLR expression") +
    theme_bw(base_size = 12) +
    theme(
      strip.placement   = "outside",
      strip.background  = element_blank(),
      strip.text.y.left = element_text(angle = 0, face = "bold", size = 11),
      axis.text.x       = element_text(size = 9),
      panel.grid.minor  = element_blank(),
      legend.position   = "right"
    )
}

p_left  <- make_block(blocks[[1]])
p_right <- make_block(blocks[[2]])

p <- (p_left | p_right) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "D0 CMV-specific CD8 T cells: CITEseq markers by specificity and visit",
    theme = theme(plot.title = element_text(face = "bold", hjust = 0.5, size = 13))
  ) &
  theme(legend.position = "right")

dir.create("corrected_FB_exports/D0", showWarnings = FALSE, recursive = TRUE)
ggsave("corrected_FB_exports/D0/D0_CMV_FB_violin_by_specificity_visit_short.png",
       p, width = 13, height = 1.4 * half + 1.6, dpi = 300)



#
#
# Dummy clonal pie charts — same donut style as the clone_piecharts_* figures
#   Two sets of three panels, meant to read as a timecourse:
#     Set 1 (polyclonal start):  12 clones -> 8 clones -> jittered repeat
#     Set 2 (oligoclonal start):  7 clones -> 5 clones -> jittered repeat
#
#   Clone identity is FIXED across the three panels of a set: a given clone keeps
#   its color in every panel, mirroring the ET-vs-LT pies where the clone mapping
#   is built once per donor and then applied to each timepoint.
#
#   Slice ORDER is set per panel via `order_by_size`:
#     TRUE  -> largest to smallest, clockwise from 12 o'clock (used for set 1)
#     FALSE -> fixed clone-number order (used for set 2)
#   Because the order differs between panels, each panel is built as its own
#   ggplot and the combined view is assembled with patchwork rather than facets.
#
#   Exports 6 individual panel PNGs plus 2 combined 1x3 figures.
#
#

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

export_dir <- "exports/dummy_pie_charts"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)


# ── 1. Fixed clone palette (identical to the real clone pie charts) ──────────

clone_colors_fixed <- c(
  "Clone 1"  = "#1b9e77",
  "Clone 2"  = "#d95f02",
  "Clone 3"  = "#7570b3",
  "Clone 4"  = "#e7298a",
  "Clone 5"  = "#66a61e",
  "Clone 6"  = "#e6ab02",
  "Clone 7"  = "#a6761d",
  "Clone 8"  = "#666666",
  "Clone 9"  = "#1f78b4",
  "Clone 10" = "#b2df8a",
  "Clone 11" = "#33a02c",
  "Clone 12" = "#fb9a99",
  "Clone 13" = "#e31a1c",
  "Clone 14" = "#fdbf6f",
  "Clone 15" = "#ff7f00",
  "Clone 16" = "#cab2d6",
  "Clone 17" = "#6a3d9a",
  "Clone 18" = "#ffff99",
  "Clone 19" = "#b15928",
  "Clone 20" = "#8dd3c7",
  "Other (<1%)" = "grey95"
)

clone_levels_fixed <- c(paste0("Clone ", 1:20), "Other (<1%)")


# ── 2. Dummy percentages ────────────────────────────────────────────────────
# Each named vector is one panel: clone -> percent. Values are hard-coded (not
# drawn at runtime) so the figure is byte-identical every time it is rebuilt.
# Percentages within a panel sum to 100; the plotting code renormalizes anyway.

set1_panels <- list(
  "Timepoint 1" = c(
    "Clone 1"  = 30,   "Clone 2"  = 15,
    "Clone 3"  = 9,    "Clone 4"  = 8,    "Clone 5"  = 7,
    "Clone 6"  = 6,    "Clone 7"  = 5,
    "Clone 8"  = 5,    "Clone 9"  = 4.5,  "Clone 10" = 4,
    "Clone 11" = 3.5,  "Clone 12" = 3
  ),
  # Clones 6, 9, 11, 12 drop out. Clone 1 dominates; clones 3, 5 and 10 expand
  # into the 10-20% range (10 is a dramatic expansion off a 4% baseline).
  "Timepoint 2" = c(
    "Clone 1"  = 45,
    "Clone 2"  = 7,
    "Clone 3"  = 17,
    "Clone 4"  = 4,
    "Clone 5"  = 12,
    "Clone 7"  = 2.5,
    "Clone 8"  = 1.5,
    "Clone 10" = 11
  ),
  # Same 8 clones, values nudged a few percent. Deliberately kept small enough
  # that no two clones swap rank, so the size-ordered panel 3 is a near-copy of
  # panel 2 rather than a visible reshuffle.
  "Timepoint 3" = c(
    "Clone 1"  = 44.3,
    "Clone 2"  = 6.6,
    "Clone 3"  = 18.2,
    "Clone 4"  = 4.5,
    "Clone 5"  = 11.6,
    "Clone 7"  = 2.5,
    "Clone 8"  = 1.7,
    "Clone 10" = 10.6
  )
)

set2_panels <- list(
  "Timepoint 1" = c(
    "Clone 1" = 60,
    "Clone 2" = 12,
    "Clone 3" = 9,
    "Clone 4" = 7,
    "Clone 5" = 5,
    "Clone 6" = 4,
    "Clone 7" = 3
  ),
  # Clones 5 and 7 drop out; Clone 1 tightens its grip on the repertoire.
  "Timepoint 2" = c(
    "Clone 1" = 75,
    "Clone 2" = 7,
    "Clone 3" = 10,
    "Clone 4" = 3,
    "Clone 6" = 5
  ),
  "Timepoint 3" = c(
    "Clone 1" = 73.8,
    "Clone 2" = 6.5,
    "Clone 3" = 10.9,
    "Clone 4" = 3.4,
    "Clone 6" = 5.4
  )
)


# ── 3. Reshape to long format ───────────────────────────────────────────────

panels_to_df <- function(panels) {
  bind_rows(lapply(names(panels), function(p) {
    v <- panels[[p]]
    tibble(panel = p, Clone = names(v), pct = as.numeric(v))
  })) %>%
    group_by(panel) %>%
    mutate(prop = pct / sum(pct)) %>%
    ungroup() %>%
    mutate(panel = factor(panel, levels = names(panels)))
}

set1_df <- panels_to_df(set1_panels)
set2_df <- panels_to_df(set2_panels)


# ── 4. Single-panel donut ───────────────────────────────────────────────────
# Slices run clockwise from 12 o'clock. geom_col stacks the FIRST factor level
# on top of the stack, so the level order is reversed before plotting: whatever
# ends up last is drawn first, i.e. immediately clockwise of 12 o'clock.

make_pie_panel <- function(d, order_by_size = TRUE, title = NULL) {

  slice_order <- if (order_by_size) {
    as.character(d$Clone[order(-d$prop)])            # largest first, clockwise
  } else {
    clone_levels_fixed[clone_levels_fixed %in% as.character(d$Clone)]
  }

  d <- d %>% mutate(Clone = factor(Clone, levels = rev(slice_order)))

  p <- ggplot(d, aes(x = 2, y = prop, fill = Clone)) +
    geom_col(color = "black", width = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    scale_fill_manual(values = clone_colors_fixed, drop = FALSE) +
    theme_void(base_size = 14) +
    theme(
      legend.position  = "none",
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  if (!is.null(title)) p <- p + ggtitle(title)
  p
}


# ── 5. Build a set: 3 individual PNGs + 1 combined 1x3 PNG ──────────────────

make_dummy_pie_set <- function(df, file_tag, order_by_size = TRUE,
                               panel_size = 4, combined_height = 4.5) {

  panel_names <- levels(df$panel)

  plots <- lapply(panel_names, function(p) {
    make_pie_panel(filter(df, panel == p),
                   order_by_size = order_by_size,
                   title         = p)
  })
  names(plots) <- panel_names

  # Individual exports — one file per panel
  for (p in panel_names) {
    fn <- file.path(export_dir,
                    sprintf("%s_%s.png", file_tag,
                            gsub("[^A-Za-z0-9]+", "_", p)))
    ggsave(fn, plots[[p]],
           width = panel_size, height = panel_size, dpi = 300, bg = "white")
    message(sprintf("  wrote %s", fn))
  }

  # Combined 1 x 3
  combined <- wrap_plots(plots, nrow = 1)
  fn_comb  <- file.path(export_dir, sprintf("%s_combined.png", file_tag))
  ggsave(fn_comb, combined,
         width = panel_size * length(plots), height = combined_height,
         dpi = 300, bg = "white")
  message(sprintf("  wrote %s", fn_comb))

  invisible(plots)
}


# ── 6. Build both sets ──────────────────────────────────────────────────────
# Set 1 is ordered largest-to-smallest within each panel; set 2 keeps the fixed
# clone-number order. Flip set 2's flag to TRUE to size-order it as well.

message("Set 1 (polyclonal, size-ordered):")
set1_plots <- make_dummy_pie_set(set1_df, "dummy_set1_polyclonal",
                                 order_by_size = TRUE)

message("Set 2 (oligoclonal, clone-number order):")
set2_plots <- make_dummy_pie_set(set2_df, "dummy_set2_oligoclonal",
                                 order_by_size = FALSE)


# ── 7. Optional: legend export (the pies themselves carry none) ──────────────

export_dummy_legend <- function(df, file_tag) {
  present <- clone_levels_fixed[clone_levels_fixed %in% as.character(df$Clone)]

  legend_src <- ggplot(
    tibble(x = seq_along(present), y = 1,
           Clone = factor(present, levels = present)),
    aes(x, y, fill = Clone)
  ) +
    geom_tile() +
    scale_fill_manual(values = clone_colors_fixed, name = "Clone") +
    guides(fill = guide_legend(ncol = 1)) +
    theme_void() +
    theme(legend.position = "right",
          legend.title    = element_text(size = 12, face = "bold"),
          legend.text     = element_text(size = 10))

  ggsave(
    filename = file.path(export_dir, paste0(file_tag, "_legend.png")),
    plot   = cowplot::ggdraw(cowplot::get_legend(legend_src)),
    width  = 2.2, height = 0.32 * length(present) + 0.6,
    dpi    = 300, bg = "white"
  )
}

# export_dummy_legend(set1_df, "dummy_set1_polyclonal")
# export_dummy_legend(set2_df, "dummy_set2_oligoclonal")


# ── 8. Optional: regenerate panel 3 as a live jitter of panel 2 ─────────────
# Replaces the hard-coded "Timepoint 3" values with a random +/- 10% wobble.
# Note this can flip the rank of two similar clones, which is visible once the
# panel is size-ordered; re-run with a different seed if that happens.

jitter_panel <- function(v, sd_frac = 0.10, seed = 42) {
  set.seed(seed)
  out <- v * runif(length(v), 1 - sd_frac, 1 + sd_frac)
  out / sum(out) * 100
}

# set1_panels[["Timepoint 3"]] <- jitter_panel(set1_panels[["Timepoint 2"]])
# set2_panels[["Timepoint 3"]] <- jitter_panel(set2_panels[["Timepoint 2"]])


#
#
# Dummy metacluster stacked bars — same style as the clone composition bar plots
#   Two sets of three panels, meant to read as a timecourse:
#     Panel 1 -> one composition
#     Panel 2 -> a different composition
#     Panel 3 -> panel 2 with the percentages nudged a few points
#
#   Single vertical bar per panel, percentages summing to 100, no axis titles.
#   Stacking follows the project convention (metacluster_order sets the levels,
#   so Tcm sits at the top of the stack and ISG+ at the bottom).
#
#   Exports 6 individual panel PNGs plus 2 combined 1x3 figures.
#
#

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

export_dir <- "exports/dummy_stacked_bars"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# TRUE  -> y axis shows 0/25/50/75/100 tick numbers
# FALSE -> bare axis lines, no tick numbers at all
SHOW_AXIS_NUMBERS <- TRUE


# ── 1. Fixed metacluster palette (same as the real stacked bar figures) ─────

metacluster_order <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")

metacluster_colors <- c(
  "Tcm"     = "#F8766D",
  "Tem"     = "#A3A500",
  "T-CMC1"  = "#00BF7D",
  "NK-like" = "#00B0F6",
  "ISG+"    = "#E76BF3"
)


# ── 2. Dummy percentages ────────────────────────────────────────────────────
# Each named vector is one panel: metacluster -> percent, summing to 100.
# Hard-coded rather than drawn at runtime so the figure rebuilds identically.

set1_panels <- list(
  "Timepoint 1" = c("Tcm" = 32,   "Tem" = 28,   "T-CMC1" = 18,
                    "NK-like" = 14,  "ISG+" = 8),
  "Timepoint 2" = c("Tcm" = 12,   "Tem" = 41,   "T-CMC1" = 9,
                    "NK-like" = 31,  "ISG+" = 7),
  # panel 2 nudged a couple of points per metacluster
  "Timepoint 3" = c("Tcm" = 13.4, "Tem" = 39.6, "T-CMC1" = 8.2,
                    "NK-like" = 32.5, "ISG+" = 6.3)
)

set2_panels <- list(
  "Timepoint 1" = c("Tcm" = 9,    "Tem" = 16,   "T-CMC1" = 35,
                    "NK-like" = 11,  "ISG+" = 29),
  "Timepoint 2" = c("Tcm" = 24,   "Tem" = 11,   "T-CMC1" = 20,
                    "NK-like" = 38,  "ISG+" = 7),
  "Timepoint 3" = c("Tcm" = 22.7, "Tem" = 12.1, "T-CMC1" = 21.3,
                    "NK-like" = 36.4, "ISG+" = 7.5)
)


# ── 3. Reshape to long format ───────────────────────────────────────────────

panels_to_df <- function(panels) {
  bind_rows(lapply(names(panels), function(p) {
    v <- panels[[p]]
    tibble(panel = p, metacluster = names(v), pct = as.numeric(v))
  })) %>%
    group_by(panel) %>%
    mutate(percent = pct / sum(pct) * 100) %>%     # renormalize, just in case
    ungroup() %>%
    mutate(
      panel       = factor(panel, levels = names(panels)),
      metacluster = factor(metacluster, levels = metacluster_order)
    )
}

set1_df <- panels_to_df(set1_panels)
set2_df <- panels_to_df(set2_panels)


# ── 4. Single-panel stacked bar ─────────────────────────────────────────────

make_bar_panel <- function(d, title = NULL) {

  p <- ggplot(d, aes(x = "", y = percent, fill = metacluster)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = metacluster_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      legend.position  = "none",
      axis.title       = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  if (!SHOW_AXIS_NUMBERS) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  if (!is.null(title)) p <- p + ggtitle(title)
  p
}


# ── 5. Build a set: 3 individual PNGs + 1 combined 1x3 PNG ──────────────────

make_dummy_bar_set <- function(df, file_tag,
                               panel_width = 2, panel_height = 4.5) {

  panel_names <- levels(df$panel)

  plots <- lapply(panel_names, function(p) {
    make_bar_panel(filter(df, panel == p), title = p)
  })
  names(plots) <- panel_names

  # Individual exports — one file per panel
  for (p in panel_names) {
    fn <- file.path(export_dir,
                    sprintf("%s_%s.png", file_tag,
                            gsub("[^A-Za-z0-9]+", "_", p)))
    ggsave(fn, plots[[p]],
           width = panel_width, height = panel_height, dpi = 300, bg = "white")
    message(sprintf("  wrote %s", fn))
  }

  # Combined 1 x 3
  combined <- wrap_plots(plots, nrow = 1)
  fn_comb  <- file.path(export_dir, sprintf("%s_combined.png", file_tag))
  ggsave(fn_comb, combined,
         width = panel_width * length(plots) + 1, height = panel_height,
         dpi = 300, bg = "white")
  message(sprintf("  wrote %s", fn_comb))

  invisible(plots)
}


# ── 6. Build both sets ──────────────────────────────────────────────────────

message("Set 1:")
set1_plots <- make_dummy_bar_set(set1_df, "dummy_bars_set1")

message("Set 2:")
set2_plots <- make_dummy_bar_set(set2_df, "dummy_bars_set2")


# ── 7. Optional: legend export (the bars themselves carry none) ─────────────

export_metacluster_legend <- function(file_tag = "dummy_bars") {
  legend_src <- ggplot(
    tibble(x = seq_along(metacluster_order), y = 1,
           metacluster = factor(metacluster_order, levels = metacluster_order)),
    aes(x, y, fill = metacluster)
  ) +
    geom_tile() +
    scale_fill_manual(values = metacluster_colors, name = "Metacluster") +
    guides(fill = guide_legend(ncol = 1)) +
    theme_void() +
    theme(legend.position = "right",
          legend.title    = element_text(size = 12, face = "bold"),
          legend.text     = element_text(size = 10))

  ggsave(
    filename = file.path(export_dir, paste0(file_tag, "_legend.png")),
    plot   = cowplot::ggdraw(cowplot::get_legend(legend_src)),
    width  = 2.2, height = 2.2, dpi = 300, bg = "white"
  )
}

# export_metacluster_legend()


# ── 8. Optional: regenerate panel 3 as a live jitter of panel 2 ────────────

jitter_panel <- function(v, sd_frac = 0.10, seed = 42) {
  set.seed(seed)
  out <- v * runif(length(v), 1 - sd_frac, 1 + sd_frac)
  out / sum(out) * 100
}

# set1_panels[["Timepoint 3"]] <- jitter_panel(set1_panels[["Timepoint 2"]])
# set2_panels[["Timepoint 3"]] <- jitter_panel(set2_panels[["Timepoint 2"]])



#
#
# Dummy metacluster stacked bars — second batch (re-randomized percentages)
#   Identical structure to dummy_metacluster_bars.R, new compositions:
#     Panel 1 -> one composition
#     Panel 2 -> a different composition
#     Panel 3 -> panel 2 with the percentages nudged a few points
#
#   File tags are set3 / set4 so these sit alongside the first batch rather
#   than overwriting it.
#
#

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

export_dir <- "exports/dummy_stacked_bars"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# TRUE  -> y axis shows 0/25/50/75/100 tick numbers
# FALSE -> bare axis lines, no tick numbers at all
SHOW_AXIS_NUMBERS <- TRUE


# ── 1. Fixed metacluster palette (same as the real stacked bar figures) ─────

metacluster_order <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")

metacluster_colors <- c(
  "Tcm"     = "#F8766D",
  "Tem"     = "#A3A500",
  "T-CMC1"  = "#00BF7D",
  "NK-like" = "#00B0F6",
  "ISG+"    = "#E76BF3"
)


# ── 2. Dummy percentages ────────────────────────────────────────────────────
# Each named vector is one panel: metacluster -> percent, summing to 100.
# Hard-coded rather than drawn at runtime so the figure rebuilds identically.

set3_panels <- list(
  "Timepoint 1" = c("Tcm" = 21,   "Tem" = 12,   "T-CMC1" = 26,
                    "NK-like" = 33,  "ISG+" = 8),
  "Timepoint 2" = c("Tcm" = 38,   "Tem" = 24,   "T-CMC1" = 11,
                    "NK-like" = 9,   "ISG+" = 18),
  # panel 2 nudged a couple of points per metacluster
  "Timepoint 3" = c("Tcm" = 36.5, "Tem" = 25.4, "T-CMC1" = 11.8,
                    "NK-like" = 8.4, "ISG+" = 17.9)
)

set4_panels <- list(
  "Timepoint 1" = c("Tcm" = 15,   "Tem" = 37,   "T-CMC1" = 14,
                    "NK-like" = 6,   "ISG+" = 28),
  "Timepoint 2" = c("Tcm" = 7,    "Tem" = 19,   "T-CMC1" = 42,
                    "NK-like" = 27,  "ISG+" = 5),
  "Timepoint 3" = c("Tcm" = 7.6,  "Tem" = 18.1, "T-CMC1" = 43.4,
                    "NK-like" = 25.6, "ISG+" = 5.3)
)


# ── 3. Reshape to long format ───────────────────────────────────────────────

panels_to_df <- function(panels) {
  bind_rows(lapply(names(panels), function(p) {
    v <- panels[[p]]
    tibble(panel = p, metacluster = names(v), pct = as.numeric(v))
  })) %>%
    group_by(panel) %>%
    mutate(percent = pct / sum(pct) * 100) %>%     # renormalize, just in case
    ungroup() %>%
    mutate(
      panel       = factor(panel, levels = names(panels)),
      metacluster = factor(metacluster, levels = metacluster_order)
    )
}

set3_df <- panels_to_df(set3_panels)
set4_df <- panels_to_df(set4_panels)


# ── 4. Single-panel stacked bar ─────────────────────────────────────────────

make_bar_panel <- function(d, title = NULL) {

  p <- ggplot(d, aes(x = "", y = percent, fill = metacluster)) +
    geom_bar(stat = "identity", position = "stack", width = 0.6) +
    scale_fill_manual(values = metacluster_colors, drop = FALSE) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25),
                       expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme_classic(base_size = 14) +
    theme(
      legend.position  = "none",
      axis.title       = element_blank(),
      axis.text.x      = element_blank(),
      axis.ticks.x     = element_blank(),
      plot.title       = element_text(size = 14, face = "bold", hjust = 0.5),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  if (!SHOW_AXIS_NUMBERS) {
    p <- p + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank())
  }

  if (!is.null(title)) p <- p + ggtitle(title)
  p
}


# ── 5. Build a set: 3 individual PNGs + 1 combined 1x3 PNG ──────────────────

make_dummy_bar_set <- function(df, file_tag,
                               panel_width = 2, panel_height = 4.5) {

  panel_names <- levels(df$panel)

  plots <- lapply(panel_names, function(p) {
    make_bar_panel(filter(df, panel == p), title = p)
  })
  names(plots) <- panel_names

  # Individual exports — one file per panel
  for (p in panel_names) {
    fn <- file.path(export_dir,
                    sprintf("%s_%s.png", file_tag,
                            gsub("[^A-Za-z0-9]+", "_", p)))
    ggsave(fn, plots[[p]],
           width = panel_width, height = panel_height, dpi = 300, bg = "white")
    message(sprintf("  wrote %s", fn))
  }

  # Combined 1 x 3
  combined <- wrap_plots(plots, nrow = 1)
  fn_comb  <- file.path(export_dir, sprintf("%s_combined.png", file_tag))
  ggsave(fn_comb, combined,
         width = panel_width * length(plots) + 1, height = panel_height,
         dpi = 300, bg = "white")
  message(sprintf("  wrote %s", fn_comb))

  invisible(plots)
}


# ── 6. Build both sets ──────────────────────────────────────────────────────

message("Set 3:")
set3_plots <- make_dummy_bar_set(set3_df, "dummy_bars_set3")

message("Set 4:")
set4_plots <- make_dummy_bar_set(set4_df, "dummy_bars_set4")


# ── 7. Optional: regenerate panel 3 as a live jitter of panel 2 ────────────

jitter_panel <- function(v, sd_frac = 0.10, seed = 7) {
  set.seed(seed)
  out <- v * runif(length(v), 1 - sd_frac, 1 + sd_frac)
  out / sum(out) * 100
}

# set3_panels[["Timepoint 3"]] <- jitter_panel(set3_panels[["Timepoint 2"]])
# set4_panels[["Timepoint 3"]] <- jitter_panel(set4_panels[["Timepoint 2"]])



#
#
# Recreating clone pies including only global clones
#
#

library(dplyr); library(tidyr); library(stringr); library(ggplot2)

export_dir <- "exports/integrated_v5/updated_clone_frequency_plots"
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# ── Settings ────────────────────────────────────────────────────────────────

# TRUE  -> cells whose cloneId is "notCalled" are counted in the grey slice
# FALSE -> dropped entirely, matching the denominator of your existing 1% pies
INCLUDE_NOTCALLED <- FALSE

# TRUE  -> a global clone's cells are pooled by the specificity in its
#          Global_Clone_ID, so cross-specificity spillover follows the clone
# FALSE -> pooled by the per-cell `specificity` metadata column (old behaviour)
USE_ID_SPECIFICITY <- TRUE

OTHER_LABEL <- "No global clone"

clone_colors_fixed <- setNames(
  c("#1b9e77","#d95f02","#7570b3","#e7298a","#66a61e",
    "#e6ab02","#a6761d","#666666","#1f78b4","#b2df8a",
    "#33a02c","#fb9a99","#e31a1c","#fdbf6f","#ff7f00",
    "#cab2d6","#6a3d9a","#ffff99","#b15928","#8dd3c7",
    "grey95"),
  c(paste0("Clone ", 1:20), OTHER_LABEL)
)
slice_levels <- names(clone_colors_fixed)

# ── 1. Per-cell table with slice assignment ─────────────────────────────────

cell_tbl <- mdata_CMV@meta.data %>%
  transmute(
    donor_identity = as.character(donor_identity),
    spec_meta      = as.character(specificity),
    cloneId        = as.character(cloneId),
    timepoint      = as.character(timepoint),
    global_id      = as.character(Global_Clone_ID)
  ) %>%
  mutate(
    is_notcalled = is.na(cloneId) | grepl("notCalled", cloneId),
    spec_id      = str_extract(global_id, "pp65|IE1"),
    specificity  = if (USE_ID_SPECIFICITY) coalesce(spec_id, spec_meta) else spec_meta,
    clone_num    = as.integer(str_extract(global_id, "\\d+$")),
    slice        = if_else(!is.na(global_id),
                           paste0("Clone ", clone_num),
                           OTHER_LABEL)
  ) %>%
  filter(
    donor_identity %in% paste("Donor", 1:5),
    specificity %in% c("pp65", "IE1"),
    INCLUDE_NOTCALLED | !is_notcalled
  ) %>%
  mutate(
    donor_identity = factor(donor_identity, levels = paste("Donor", 1:5)),
    specificity    = factor(specificity, levels = c("pp65", "IE1")),
    timepoint      = factor(timepoint, levels = c("Early", "Late"))
  )

stopifnot(max(cell_tbl$clone_num, na.rm = TRUE) <= 20)   # palette has 20 slots

# ── 2. Pie builder ──────────────────────────────────────────────────────────
# Slices run clockwise from 12 o'clock in clone-number order, grey last.
# geom_col stacks the first factor level on top, so levels are reversed here.

make_pie <- function(df, facet_formula, filename,
                     width = 12, height = 6, show_legend = FALSE) {

  df <- df %>% mutate(slice = factor(slice, levels = rev(slice_levels)))

  p <- ggplot(df, aes(x = 2, y = prop, fill = slice)) +
    geom_col(color = "black", width = 1) +
    coord_polar(theta = "y") +
    xlim(0.5, 2.5) +
    facet_grid(facet_formula) +
    scale_fill_manual(values = clone_colors_fixed, drop = FALSE,
                      name = "Global clone") +
    theme_void(base_size = 14) +
    theme(
      legend.position  = if (show_legend) "right" else "none",
      strip.text       = element_text(size = 13),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background  = element_rect(fill = "white", color = NA)
    )

  ggsave(file.path(export_dir, filename), p,
         width = width, height = height, dpi = 300, bg = "white")
  message("  wrote ", filename)
  invisible(p)
}

# ── 3. Combined-timepoint pies (specificity ~ donor) ────────────────────────

pie_combined <- cell_tbl %>%
  count(donor_identity, specificity, slice, name = "cell_count") %>%
  group_by(donor_identity, specificity) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

make_pie(pie_combined, specificity ~ donor_identity,
         "clone_piecharts_global_clones.png")

# ── 4. Early vs Late pies, one figure per specificity (timepoint ~ donor) ────
# Denominator is that donor x specificity x timepoint pool, so each panel sums
# to 100%. A clone absent at a timepoint simply has no slice in that panel.

pie_by_tp <- cell_tbl %>%
  filter(!is.na(timepoint)) %>%
  count(donor_identity, specificity, timepoint, slice, name = "cell_count") %>%
  group_by(donor_identity, specificity, timepoint) %>%
  mutate(prop = cell_count / sum(cell_count)) %>%
  ungroup()

for (sp in c("pp65", "IE1")) {
  make_pie(filter(pie_by_tp, specificity == sp),
           timepoint ~ donor_identity,
           sprintf("%s_clone_piecharts_global_clones.png", sp))
}

# ── 5. Legend (the pies themselves carry none) ──────────────────────────────

present <- slice_levels[slice_levels %in% unique(pie_combined$slice)]
legend_src <- ggplot(
  tibble(x = seq_along(present), y = 1,
         slice = factor(present, levels = present)),
  aes(x, y, fill = slice)
) +
  geom_tile() +
  scale_fill_manual(values = clone_colors_fixed, name = "Global clone") +
  guides(fill = guide_legend(ncol = 1)) +
  theme_void() +
  theme(legend.position = "right",
        legend.title    = element_text(size = 12, face = "bold"),
        legend.text     = element_text(size = 10))

ggsave(file.path(export_dir, "clone_piecharts_global_clones_legend.png"),
       cowplot::ggdraw(cowplot::get_legend(legend_src)),
       width = 2.4, height = 0.32 * length(present) + 0.6,
       dpi = 300, bg = "white")

# ── 6. Sanity check ─────────────────────────────────────────────────────────
# named_slices should reproduce your C3 counts: D1 8/1, D2 5/1, D3 8/2,
# D4 10/5, D5 3/5, summing to 48.

cat("\n=== Combined-timepoint pies ===\n")
pie_combined %>%
  group_by(donor_identity, specificity) %>%
  summarise(
    named_slices = sum(slice != OTHER_LABEL),
    pct_named    = round(100 * sum(prop[slice != OTHER_LABEL]), 1),
    pct_grey     = round(100 * sum(prop[slice == OTHER_LABEL]), 1),
    pool_cells   = sum(cell_count),
    .groups = "drop"
  ) %>% as.data.frame() %>% print()

cat(sprintf("\nTotal named slices across all pies: %d\n",
            sum(pie_combined$slice != OTHER_LABEL)))

cat("\n=== Early vs Late pies ===\n")
pie_by_tp %>%
  group_by(donor_identity, specificity, timepoint) %>%
  summarise(
    named_slices = sum(slice != OTHER_LABEL),
    pct_grey     = round(100 * sum(prop[slice == OTHER_LABEL]), 1),
    .groups = "drop"
  ) %>% as.data.frame() %>% print()




#
#
# RE-EXPORT — Donor 0 pp65 clone frequency over time (log scale)
#   Clones 1, 2, 3, 4, 19 highlighted: drawn above all other clones, with
#   lines and points at 1.5x the base weight. x-axis relabeled Visit -> Week.
#   Output: re-exports_for_figures/pp65_clone_freq_over_time_log.png
#
#

library(dplyr)
library(ggplot2)

# ── 0. Output folder (shared by all re-exported figures) ────────────────────
reexport_dir <- "re-exports_for_figures"
dir.create(reexport_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Settings ─────────────────────────────────────────────────────────────
clone_colors <- c(
  "D0_pp65_clone_1"  = "#1b9e77","D0_pp65_clone_2"  = "#d95f02","D0_pp65_clone_3"  = "#7570b3",
  "D0_pp65_clone_4"  = "#e7298a","D0_pp65_clone_5"  = "#66a61e","D0_pp65_clone_6"  = "#e6ab02",
  "D0_pp65_clone_7"  = "#a6761d","D0_pp65_clone_8"  = "#666666","D0_pp65_clone_9"  = "#1f78b4",
  "D0_pp65_clone_10" = "#b2df8a","D0_pp65_clone_11" = "#33a02c","D0_pp65_clone_12" = "#fb9a99",
  "D0_pp65_clone_13" = "#e31a1c","D0_pp65_clone_14" = "#fdbf6f","D0_pp65_clone_15" = "#ff7f00",
  "D0_pp65_clone_16" = "#cab2d6","D0_pp65_clone_17" = "#6a3d9a","D0_pp65_clone_18" = "#ffff99",
  "D0_pp65_clone_19" = "#b15928","D0_pp65_clone_20" = "#8dd3c7"
)

# Highlighted clones, listed BOTTOM -> TOP: the last entry is drawn last and
# therefore sits on top of everything else (clone 1 here)
highlight_ids <- paste0("D0_pp65_clone_", c(19, 4, 3, 2, 1))

BASE_LINEWIDTH <- 1      # original geom_line size
BASE_POINTSIZE <- 3      # original geom_point size
HL_SCALE       <- 1.3    # weight multiplier for highlighted clones

visit_to_week <- c(
  "Visit 1" = "Week 20",
  "Visit 2" = "Week 21",
  "Visit 3" = "Week 32",
  "Visit 4" = "Week 69"
)

# ── 2. Data (same calculation as the original figure) ───────────────────────
meta_df <- D0_data_CMV@meta.data %>%
  dplyr::select(specificity, timepoint, Global_Clone_ID) %>%
  dplyr::filter(specificity == "pp65",
                Global_Clone_ID %in% names(clone_colors))

percent_df <- meta_df %>%
  dplyr::group_by(timepoint, Global_Clone_ID) %>%
  dplyr::summarise(count = dplyr::n(), .groups = "drop") %>%
  dplyr::group_by(timepoint) %>%
  dplyr::mutate(percent = count / sum(count) * 100) %>%
  dplyr::ungroup()

# Relabel visits as weeks (fails loudly if an unexpected timepoint appears)
stopifnot(all(as.character(percent_df$timepoint) %in% names(visit_to_week)))
percent_df <- percent_df %>%
  dplyr::mutate(timepoint = factor(unname(visit_to_week[as.character(timepoint)]),
                                   levels = unname(visit_to_week)))

missing_hl <- setdiff(highlight_ids, unique(percent_df$Global_Clone_ID))
if (length(missing_hl))
  warning("Highlighted clones not found: ", paste(missing_hl, collapse = ", "))

# Background = every plotted clone that is NOT highlighted
bg_df <- percent_df %>% dplyr::filter(!Global_Clone_ID %in% highlight_ids)

# One line + point layer pair per highlighted clone, so each clone's dots and
# line sit above every clone drawn before it (flattened into a single list)
hl_layers <- unlist(lapply(highlight_ids, function(cid) {
  d <- dplyr::filter(percent_df, Global_Clone_ID == cid)
  list(
    geom_line(data = d,  linewidth = BASE_LINEWIDTH * HL_SCALE),
    geom_point(data = d, size      = BASE_POINTSIZE * HL_SCALE)
  )
}), recursive = FALSE)

# ── 3. Plot ─────────────────────────────────────────────────────────────────
pp65_plot <- ggplot(percent_df,
                    aes(x = timepoint, y = percent,
                        group = Global_Clone_ID, color = Global_Clone_ID)) +
  geom_line(data = bg_df,  linewidth = BASE_LINEWIDTH) +
  geom_point(data = bg_df, size      = BASE_POINTSIZE) +
  hl_layers +
  scale_color_manual(values = clone_colors) +
  scale_y_log10() +
  labs(x = NULL, y = "Clone Frequency (log scale)") +
  ggtitle("pp65 Major Clones") +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x  = element_text(size = 14),
    axis.title.y = element_text(size = 16, margin = margin(r = 15)),
    plot.title   = element_text(size = 18, face = "bold", hjust = 0.5)
  )

# ── 4. Export ───────────────────────────────────────────────────────────────
ggsave(file.path(reexport_dir, "pp65_clone_freq_over_time_log.png"),
       plot = pp65_plot, width = 8, height = 6, dpi = 300, bg = "white")


#
#
# RE-EXPORT — D0 CITEseq bubble plot by seurat_cluster (short)
#   All bubbles one fixed size; color = unscaled CLR expression, computed the
#   same way as the specificity bubble plot (DotPlot with scale = FALSE).
#   Output: re-exports_for_figures/D0_CMV_FB_bubbleplot_by_seurat_cluster_short.png
#
#

library(Seurat)
library(dplyr)
library(ggplot2)

# ── 0. Output folder (shared by all re-exported figures) ────────────────────
reexport_dir <- "re-exports_for_figures"
dir.create(reexport_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Settings ─────────────────────────────────────────────────────────────
EXCLUDE_67   <- TRUE        # TRUE = manuscript version (clusters 0-5)
                            # FALSE = all clusters (original corrected_FB_exports version)
COLOR_METRIC <- "dotplot"   # "dotplot": log1p(mean(expm1(CLR))) — identical to what
                            #            DotPlot(scale = FALSE) returns, i.e. what the
                            #            specificity bubble plot shows
                            # "mean":    plain arithmetic mean of the CLR values
BUBBLE_SIZE  <- 6           # fixed size for every bubble

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

# ── 2. Cells and clusters to plot ───────────────────────────────────────────
cluster_vec  <- setNames(as.character(D0_data_CMV$seurat_clusters), colnames(D0_data_CMV))
all_clusters <- as.character(sort(as.numeric(unique(cluster_vec))))
keep_clusters <- if (EXCLUDE_67) setdiff(all_clusters, c("6", "7")) else all_clusters

cells_keep <- names(cluster_vec)[cluster_vec %in% keep_clusters]
message(sprintf("Clusters plotted: %s | cells: %d of %d",
                paste(keep_clusters, collapse = ", "),
                length(cells_keep), length(cluster_vec)))

# ── 3. Per-cluster CLR summary straight from the FB assay ───────────────────
expr <- as.matrix(
  GetAssayData(D0_data_CMV, assay = "FB", layer = "data")[d0_markers, cells_keep, drop = FALSE]
)
grp <- factor(cluster_vec[colnames(expr)], levels = keep_clusters)

plot_data <- bind_rows(lapply(levels(grp), function(g) {
  x <- expr[, grp == g, drop = FALSE]
  tibble(
    marker      = rownames(x),
    cluster     = g,
    mean_clr    = rowMeans(x),
    dotplot_val = log1p(rowMeans(expm1(x))),
    n_cells     = ncol(x)
  )
})) %>%
  mutate(
    value   = if (COLOR_METRIC == "dotplot") dotplot_val else mean_clr,
    cluster = factor(cluster, levels = keep_clusters),
    marker  = factor(marker,  levels = rev(d0_markers))
  )

message(sprintf("Color range (%s): %.2f to %.2f",
                COLOR_METRIC, min(plot_data$value), max(plot_data$value)))

legend_title <- if (COLOR_METRIC == "dotplot") {
  "Mean CLR\nexpression\n(log1p)"
} else {
  "Mean CLR\nexpression"
}

# ── 4. Plot ─────────────────────────────────────────────────────────────────
p_bubble <- ggplot(plot_data, aes(x = cluster, y = marker)) +
  geom_point(aes(color = value), size = BUBBLE_SIZE) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red"),
                        name = legend_title) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title  = element_text(size = 12)
  ) +
  labs(x = "Seurat cluster", y = "CITEseq marker")

# ── 5. Export ───────────────────────────────────────────────────────────────
out_name <- if (EXCLUDE_67) {
  "D0_CMV_FB_bubbleplot_by_seurat_cluster_short.png"
} else {
  "D0_CMV_FB_bubbleplot_by_seurat_cluster_short_all_clusters.png"
}
plot_w <- if (EXCLUDE_67) 6.5 else 7
plot_h <- if (EXCLUDE_67) 3.5 else 4

ggsave(file.path(reexport_dir, out_name),
       plot = p_bubble, width = plot_w, height = plot_h, dpi = 300, bg = "white")


#
#
# RE-EXPORT — D0 CITEseq bubble plot by specificity x visit block (short)
#   All bubbles one fixed size; color = unscaled CLR (same metric as the
#   cluster bubble plot). y labels: (Visits 1-3) -> (Weeks 20-32),
#   (Visit 4) -> (Week 69).
#   Output: re-exports_for_figures/D0_CMV_FB_bubbleplot_by_specificity_visit_short.png
#
#

library(Seurat)
library(dplyr)
library(ggplot2)

# ── 0. Output folder (shared by all re-exported figures) ────────────────────
reexport_dir <- "re-exports_for_figures"
dir.create(reexport_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Settings ─────────────────────────────────────────────────────────────
EXCLUDE_67   <- FALSE       # FALSE = all D0 cells, as in the original figure
                            # TRUE  = drop seurat clusters 6 and 7 first
COLOR_METRIC <- "dotplot"   # "dotplot": log1p(mean(expm1(CLR))) — DotPlot(scale = FALSE)
                            # "mean":    plain arithmetic mean of the CLR values
BUBBLE_SIZE  <- 8           # fixed size for every bubble

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

visit_block_map <- c(
  "Visit 1" = "Weeks 20-32",
  "Visit 2" = "Weeks 20-32",
  "Visit 3" = "Weeks 20-32",
  "Visit 4" = "Week 69"
)

grp_levels <- c("pp65 (Weeks 20-32)", "pp65 (Week 69)",
                "IE1 (Weeks 20-32)",  "IE1 (Week 69)")

spec_cols <- c(pp65 = "#1976D2", IE1 = "#8B0000")

# ── 2. Assign cells to specificity x week-block groups ──────────────────────
md        <- D0_data_CMV@meta.data
spec_chr  <- as.character(md$specificity)
block_chr <- unname(visit_block_map[as.character(md$timepoint)])
clus_chr  <- as.character(md$seurat_clusters)

grp_vec <- setNames(
  ifelse(spec_chr %in% c("pp65", "IE1") & !is.na(block_chr),
         paste0(spec_chr, " (", block_chr, ")"),
         NA_character_),
  rownames(md)
)
if (EXCLUDE_67) grp_vec[clus_chr %in% c("6", "7")] <- NA_character_

cells_keep <- names(grp_vec)[!is.na(grp_vec)]
grp        <- factor(grp_vec[cells_keep], levels = grp_levels)

grp_n <- table(grp)
message("Cells per group:"); print(grp_n)
stopifnot(all(grp_n > 0))

# ── 3. Per-group CLR summary straight from the FB assay ─────────────────────
expr <- as.matrix(
  GetAssayData(D0_data_CMV, assay = "FB", layer = "data")[d0_markers, cells_keep, drop = FALSE]
)

plot_data <- bind_rows(lapply(grp_levels, function(g) {
  x <- expr[, grp == g, drop = FALSE]
  tibble(
    marker      = rownames(x),
    group       = g,
    mean_clr    = rowMeans(x),
    dotplot_val = log1p(rowMeans(expm1(x)))
  )
})) %>%
  mutate(
    value  = if (COLOR_METRIC == "dotplot") dotplot_val else mean_clr,
    marker = factor(marker, levels = d0_markers),
    group  = factor(group,  levels = rev(grp_levels))   # pp65 block on top
  )

message(sprintf("Color range (%s): %.2f to %.2f",
                COLOR_METRIC, min(plot_data$value), max(plot_data$value)))

legend_title <- if (COLOR_METRIC == "dotplot") {
  "Mean CLR\nexpression\n(log1p)"
} else {
  "Mean CLR\nexpression"
}

# ── 4. y-axis labels: colored by specificity, with group n ──────────────────
y_levels  <- levels(plot_data$group)
y_lab_col <- unname(ifelse(grepl("^pp65", y_levels), spec_cols["pp65"], spec_cols["IE1"]))
y_lab_txt <- setNames(
  paste0(y_levels, "\nn=", as.integer(grp_n[y_levels])),
  y_levels
)

# divider between the pp65 block and the IE1 block
split_at <- sum(grepl("^IE1", y_levels)) + 0.5

# ── 5. Plot ─────────────────────────────────────────────────────────────────
p <- ggplot(plot_data, aes(x = marker, y = group)) +
  geom_hline(yintercept = split_at, linetype = "dashed",
             color = "grey40", linewidth = 0.4) +
  geom_point(aes(color = value), size = BUBBLE_SIZE) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red"),
                        name = legend_title) +
  scale_y_discrete(labels = y_lab_txt) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 9.5, face = "bold", color = y_lab_col),
    axis.title       = element_text(size = 12),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "CITEseq marker", y = NULL,
       title = "D0 CMV-specific CD8 T cells: CITEseq markers by specificity and visit")

# ── 6. Export ───────────────────────────────────────────────────────────────
ggsave(file.path(reexport_dir, "D0_CMV_FB_bubbleplot_by_specificity_visit_short.png"),
       plot = p, width = 8, height = 3.8, dpi = 300, bg = "white")



#
#
# RE-EXPORT — D0 CITEseq bubble plot by major clone (IE1 1-7, pp65 1-20)
#   All bubbles one fixed size; color = unscaled CLR (same metric as the
#   cluster and specificity bubble plots). Labels "pp65 Clone 1" style, no n.
#   Output: re-exports_for_figures/D0_CMV_FB_bubbleplot_by_clone.png
#
#

library(Seurat)
library(dplyr)
library(stringr)
library(ggplot2)

# ── 0. Output folder (shared by all re-exported figures) ────────────────────
reexport_dir <- "re-exports_for_figures"
dir.create(reexport_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Settings ─────────────────────────────────────────────────────────────
EXCLUDE_67   <- FALSE         # FALSE = all D0 cells, as in the original figure
                              # TRUE  = drop seurat clusters 6 and 7 first
COLOR_METRIC <- "dotplot"     # "dotplot": log1p(mean(expm1(CLR))) — DotPlot(scale = FALSE)
                              # "mean":    plain arithmetic mean of the CLR values
BUBBLE_SIZE  <- 6             # fixed size for every bubble
BLOCK_ORDER  <- c("IE1", "pp65")   # specificity blocks, top -> bottom

d0_markers <- intersect(
  c("CD45RA", "CD95", "CD57", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(D0_data_CMV[["FB"]])
)

# Numeric ordering, not lexicographic (clone_2 before clone_10)
clone_sets <- list(
  IE1  = paste0("D0_IE1_clone_",  1:7),
  pp65 = paste0("D0_pp65_clone_", 1:20)
)

spec_cols <- c(pp65 = "#1976D2", IE1 = "#8B0000")

# ── 2. Assign cells to clones ───────────────────────────────────────────────
md        <- D0_data_CMV@meta.data
clone_vec <- setNames(as.character(md$Global_Clone_ID), rownames(md))
clus_chr  <- as.character(md$seurat_clusters)
if (EXCLUDE_67) clone_vec[clus_chr %in% c("6", "7")] <- NA_character_

clone_levels_all <- unlist(clone_sets[BLOCK_ORDER], use.names = FALSE)   # top -> bottom
present <- clone_levels_all[clone_levels_all %in% clone_vec]
missing <- setdiff(clone_levels_all, present)
if (length(missing))
  message("Not found (or no cells left): ", paste(missing, collapse = ", "))

cells_keep <- names(clone_vec)[!is.na(clone_vec) & clone_vec %in% present]
grp        <- factor(clone_vec[cells_keep], levels = present)

clone_n <- table(grp)
message("Cells per clone:"); print(clone_n)
message("Total cells retained: ", length(cells_keep))

# ── 3. Per-clone CLR summary straight from the FB assay ─────────────────────
expr <- as.matrix(
  GetAssayData(D0_data_CMV, assay = "FB", layer = "data")[d0_markers, cells_keep, drop = FALSE]
)

plot_data <- bind_rows(lapply(present, function(g) {
  x <- expr[, grp == g, drop = FALSE]
  tibble(
    marker      = rownames(x),
    clone       = g,
    mean_clr    = rowMeans(x),
    dotplot_val = log1p(rowMeans(expm1(x)))
  )
})) %>%
  mutate(
    value  = if (COLOR_METRIC == "dotplot") dotplot_val else mean_clr,
    marker = factor(marker, levels = d0_markers),
    clone  = factor(clone,  levels = rev(present))   # first entry drawn at the top
  )

message(sprintf("Color range (%s): %.2f to %.2f",
                COLOR_METRIC, min(plot_data$value), max(plot_data$value)))

legend_title <- if (COLOR_METRIC == "dotplot") {
  "Mean CLR\nexpression\n(log1p)"
} else {
  "Mean CLR\nexpression"
}

# ── 4. y-axis labels: "pp65 Clone 1" style, colored by specificity ──────────
fmt_clone_label <- function(ids) {
  m <- str_match(ids, "^D0_(pp65|IE1)_clone_(\\d+)$")
  ifelse(is.na(m[, 2]), ids, sprintf("%s Clone %s", m[, 2], m[, 3]))
}

y_levels  <- levels(plot_data$clone)
y_spec    <- str_extract(y_levels, "pp65|IE1")
y_lab_col <- unname(spec_cols[y_spec])
y_lab_txt <- setNames(fmt_clone_label(y_levels), y_levels)

# divider between the two specificity blocks (skipped if only one is present)
bottom_spec <- BLOCK_ORDER[length(BLOCK_ORDER)]
split_at    <- sum(y_spec == bottom_spec) + 0.5
divider     <- if (length(unique(y_spec)) == 2) {
  geom_hline(yintercept = split_at, linetype = "dashed",
             color = "grey40", linewidth = 0.4)
} else NULL

# ── 5. Plot ─────────────────────────────────────────────────────────────────
p <- ggplot(plot_data, aes(x = marker, y = clone)) +
  divider +
  geom_point(aes(color = value), size = BUBBLE_SIZE) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red"),
                        name = legend_title) +
  scale_y_discrete(labels = y_lab_txt) +
  theme_bw() +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 9, face = "bold", color = y_lab_col),
    axis.title       = element_text(size = 12),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank()
  ) +
  labs(x = "CITEseq marker", y = "Clone",
       title = "D0 CMV-specific CD8 T cells: CITEseq markers by major clone")

# ── 6. Export ───────────────────────────────────────────────────────────────
ggsave(file.path(reexport_dir, "D0_CMV_FB_bubbleplot_by_clone.png"),
       plot = p, width = 8, height = 10.5, dpi = 300, bg = "white")


#
#
# RE-EXPORT — D1-5 CITEseq bubble plot by metacluster
#   All bubbles one fixed size; color = unscaled CLR (same metric as the D0
#   bubble plots). Layout matches the original: markers on x, metaclusters on y.
#   Output: re-exports_for_figures/D1-5_CMV_FB_bubbleplot_by_metacluster.png
#
#

library(Seurat)
library(dplyr)
library(ggplot2)

# ── 0. Output folder (shared by all re-exported figures) ────────────────────
reexport_dir <- "re-exports_for_figures"
dir.create(reexport_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Settings ─────────────────────────────────────────────────────────────
COLOR_METRIC <- "dotplot"   # "dotplot": log1p(mean(expm1(CLR))) — DotPlot(scale = FALSE)
                            # "mean":    plain arithmetic mean of the CLR values
BUBBLE_SIZE  <- 6           # fixed size for every bubble

mdata_markers <- intersect(
  c("CD45RA", "CCR7", "CD95", "CD57", "KLRG1", "CXCR3", "CD226", "PD-1", "TIGIT", "Tim-3"),
  rownames(mdata_CMV[["FB"]])
)

metacluster_order <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")

# ── 2. Cells with a metacluster assignment ──────────────────────────────────
md       <- mdata_CMV@meta.data
meta_vec <- setNames(as.character(md$metacluster), rownames(md))

n_drop <- sum(!meta_vec %in% metacluster_order)
if (n_drop) message(sprintf("Dropping %d cells with no / unexpected metacluster", n_drop))

cells_keep <- names(meta_vec)[meta_vec %in% metacluster_order]
grp        <- factor(meta_vec[cells_keep], levels = metacluster_order)

grp_n <- table(grp)
message("Cells per metacluster:"); print(grp_n)
stopifnot(all(grp_n > 0))

# ── 3. Per-metacluster CLR summary straight from the FB assay ───────────────
expr <- as.matrix(
  GetAssayData(mdata_CMV, assay = "FB", layer = "data")[mdata_markers, cells_keep, drop = FALSE]
)

plot_data <- bind_rows(lapply(metacluster_order, function(g) {
  x <- expr[, grp == g, drop = FALSE]
  tibble(
    marker      = rownames(x),
    metacluster = g,
    mean_clr    = rowMeans(x),
    dotplot_val = log1p(rowMeans(expm1(x)))
  )
})) %>%
  mutate(
    value       = if (COLOR_METRIC == "dotplot") dotplot_val else mean_clr,
    marker      = factor(marker,      levels = mdata_markers),
    metacluster = factor(metacluster, levels = rev(metacluster_order))   # Tcm on top
  )

message(sprintf("Color range (%s): %.2f to %.2f",
                COLOR_METRIC, min(plot_data$value), max(plot_data$value)))

legend_title <- if (COLOR_METRIC == "dotplot") {
  "Mean CLR\nexpression\n(log1p)"
} else {
  "Mean CLR\nexpression"
}

# ── 4. Plot ─────────────────────────────────────────────────────────────────
p <- ggplot(plot_data, aes(x = marker, y = metacluster)) +
  geom_point(aes(color = value), size = BUBBLE_SIZE) +
  scale_color_gradientn(colors = c("blue", "lightgrey", "red"),
                        name = legend_title) +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title  = element_blank()
  )

# ── 5. Export ───────────────────────────────────────────────────────────────
ggsave(file.path(reexport_dir, "D1-5_CMV_FB_bubbleplot_by_metacluster.png"),
       plot = p, width = 7, height = 2.5, dpi = 300, bg = "white")



#
#
# RE-EXPORT — pp65 vs IE1 cross-object concordance (Latent cohort vs Seroconverter)
#   Panels (left -> right): all visits | weeks 20-32 | week 69
#   One fitted line (TLS). Annotation: n shared genes, % sign-concordant,
#   TLS slope + r². Reads the saved DE table, so DE is not re-run.
#   Output: re-exports_for_figures/pp65_vs_IE1_crossobject_concordance.{png,pdf}
#
#

suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(purrr); library(readr)
  library(ggplot2); library(ggrepel); library(patchwork)
})

# ── 0. Output folder (shared by all re-exported figures) ────────────────────
reexport_dir <- "re-exports_for_figures"
dir.create(reexport_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Settings (thresholds unchanged from the original figure) ─────────────
SPEC_1  <- "pp65"   # positive log2FC
SPEC_2  <- "IE1"    # negative log2FC
LFC_CUT <- 0.25     # |log2FC| threshold for concordance calls

SHOW_IDENTITY <- TRUE   # dashed y = x reference line (not a fitted line)

GENES_LABEL  <- c("GZMB","GZMK","IFNG","PRF1","NKG7","KLRG1","CX3CR1","ZEB2",
                  "TCF7","SELL","CCR7","IL7R","TOX","PDCD1","LAG3","ISG15")
AUTO_LABEL_N <- 0       # bump to auto-label top-N concordant genes by magnitude

PAL <- c("Higher in pp65 (both)" = "#1976D2",
         "Higher in IE1 (both)"  = "#8B0000",
         "Discordant"            = "grey45",
         "Below threshold"       = "grey85")

X_LAB <- sprintf("log2FC %s vs %s - Latent cohort",  SPEC_1, SPEC_2)
Y_LAB <- sprintf("log2FC %s vs %s - Seroconverter", SPEC_1, SPEC_2)

# Panel order (left -> right) and titles
panel_spec <- tribble(
  ~comparison,      ~title,
  "D0 all visits",  "Latent cohort vs Seroconverter (all visits)",
  "D0 Visits 1-3",  "Latent cohort vs Seroconverter (weeks 20-32)",
  "D0 Visit 4",     "Latent cohort vs Seroconverter (week 69)"
)

# ── 2. Load the saved DE results ────────────────────────────────────────────
de_csv <- "exports/cross_object_DE/pp65_vs_IE1_lfc_shared_genes.csv"

de_all <- if (file.exists(de_csv)) {
  read_csv(de_csv, show_col_types = FALSE)
} else if (exists("results") && is.list(results) &&
           all(panel_spec$comparison %in% names(results))) {
  bind_rows(results)
} else {
  stop("DE table not found — run the cross-object DE block first, ",
       "or point de_csv at the saved CSV.")
}

stopifnot(all(c("gene", "lfc_x", "lfc_y", "comparison") %in% names(de_all)))
missing_cmp <- setdiff(panel_spec$comparison, unique(de_all$comparison))
if (length(missing_cmp))
  stop("Comparisons missing from DE table: ", paste(missing_cmp, collapse = ", "))

# ── 3. Helpers ──────────────────────────────────────────────────────────────
# total-least-squares slope (both axes carry measurement error)
orth_fit <- function(x, y) {
  v <- prcomp(cbind(x, y), center = TRUE, scale. = FALSE)$rotation[, 1]
  s <- unname(v[2] / v[1])
  c(slope = s, intercept = mean(y) - s * mean(x))
}

panel_stats <- function(df) {
  conc <- filter(df, abs(lfc_x) > LFC_CUT & abs(lfc_y) > LFC_CUT)
  os   <- orth_fit(df$lfc_x, df$lfc_y)
  r    <- cor(df$lfc_x, df$lfc_y, method = "pearson")
  tibble(
    n_shared            = nrow(df),
    n_concordance_calls = nrow(conc),
    pct_sign_concordant = if (nrow(conc)) 100 * mean(sign(conc$lfc_x) == sign(conc$lfc_y)) else NA_real_,
    tls_slope           = os[["slope"]],
    tls_intercept       = os[["intercept"]],
    pearson_r           = r,
    r_squared           = r^2
  )
}

# ── 4. Single-panel plot ────────────────────────────────────────────────────
concordance_plot_v2 <- function(df, title) {

  df <- df %>%
    filter(!is.na(lfc_x), !is.na(lfc_y)) %>%
    mutate(
      cat = case_when(
        lfc_x >  LFC_CUT & lfc_y >  LFC_CUT         ~ "Higher in pp65 (both)",
        lfc_x < -LFC_CUT & lfc_y < -LFC_CUT         ~ "Higher in IE1 (both)",
        abs(lfc_x) > LFC_CUT & abs(lfc_y) > LFC_CUT ~ "Discordant",
        TRUE                                        ~ "Below threshold"),
      cat = factor(cat, levels = names(PAL))
    )

  st <- panel_stats(df)

  stat_lab <- sprintf(
    "n = %s shared genes\n%.0f%% sign-concordant (n = %s)\nTLS slope = %.2f, r\u00B2 = %.2f",
    format(st$n_shared, big.mark = ","),
    st$pct_sign_concordant,
    format(st$n_concordance_calls, big.mark = ","),
    st$tls_slope, st$r_squared
  )

  lab_df <- filter(df, gene %in% GENES_LABEL)
  if (AUTO_LABEL_N > 0) {
    top <- df %>%
      filter(abs(lfc_x) > LFC_CUT & abs(lfc_y) > LFC_CUT,
             sign(lfc_x) == sign(lfc_y)) %>%
      slice_max(abs(lfc_x) + abs(lfc_y), n = AUTO_LABEL_N)
    lab_df <- bind_rows(lab_df, top) %>% distinct(gene, .keep_all = TRUE)
  }

  lim <- max(abs(c(df$lfc_x, df$lfc_y)), na.rm = TRUE) * 1.06

  # Legend-only layer: one invisible point per category. This gives every key
  # (including "Below threshold") a dot, and keeps the legend identical in all
  # three panels so patchwork merges them into a single legend.
  legend_df <- tibble(lfc_x = 0, lfc_y = 0,
                      cat = factor(names(PAL), levels = names(PAL)))

  ggplot(df, aes(lfc_x, lfc_y)) +
    geom_hline(yintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_vline(xintercept = 0, colour = "grey80", linewidth = 0.3) +
    geom_point(data = filter(df, cat == "Below threshold"),
               colour = PAL[["Below threshold"]], size = 0.7, alpha = 0.55) +
    geom_point(data = filter(df, cat != "Below threshold"),
               aes(colour = cat), size = 1.3, alpha = 0.85, show.legend = FALSE) +
    geom_point(data = legend_df, aes(colour = cat), size = 2.4, alpha = 0) +
    { if (SHOW_IDENTITY)
        geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                    colour = "grey30", linewidth = 0.45) } +
    geom_abline(slope = st$tls_slope, intercept = st$tls_intercept,
                colour = "black", linewidth = 0.6) +
    geom_point(data = lab_df, shape = 21, fill = NA, colour = "black",
               size = 2.1, stroke = 0.5) +
    geom_text_repel(data = lab_df, aes(label = gene), size = 3,
                    min.segment.length = 0, segment.size = 0.3,
                    segment.colour = "grey40", box.padding = 0.4,
                    max.overlaps = Inf, seed = 42) +
    annotate("text", x = -lim * 0.96, y = lim * 0.96, label = stat_lab,
             hjust = 0, vjust = 1, size = 2.8, lineheight = 1.2, colour = "grey15") +
    scale_colour_manual(values = PAL, drop = FALSE, name = NULL) +
    coord_fixed(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
    labs(x = X_LAB, y = Y_LAB, title = title) +
    theme_classic(base_size = 11) +
    theme(
      plot.title      = element_text(face = "bold", size = 11, hjust = 0),
      axis.text       = element_text(colour = "black"),
      axis.title      = element_text(size = 9.5),
      legend.position = "bottom",
      legend.text     = element_text(size = 8.5),
      legend.key.size = unit(0.4, "cm"),
      panel.border    = element_rect(colour = "black", fill = NA, linewidth = 0.5),
      axis.line       = element_blank()
    ) +
    guides(colour = guide_legend(override.aes = list(size = 2.4, alpha = 1), nrow = 2))
}

# ── 5. Assemble ─────────────────────────────────────────────────────────────
plots <- map2(panel_spec$comparison, panel_spec$title, function(cmp, ttl) {
  concordance_plot_v2(filter(de_all, comparison == cmp), ttl)
})

fig <- wrap_plots(plots, nrow = 1) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "bottom")

# ── 6. Export ───────────────────────────────────────────────────────────────
ggsave(file.path(reexport_dir, "pp65_vs_IE1_crossobject_concordance.png"),
       fig, width = 16, height = 6.2, dpi = 400, bg = "white")
ggsave(file.path(reexport_dir, "pp65_vs_IE1_crossobject_concordance.pdf"),
       fig, width = 16, height = 6.2, device = cairo_pdf, bg = "white")

# Per-panel stats, for the figure legend / methods
stats_tbl <- map2_dfr(panel_spec$comparison, panel_spec$title, function(cmp, ttl) {
  d <- filter(de_all, comparison == cmp, !is.na(lfc_x), !is.na(lfc_y))
  panel_stats(d) %>% mutate(panel = ttl, comparison = cmp, .before = 1)
})
write_csv(stats_tbl, file.path(reexport_dir, "pp65_vs_IE1_crossobject_concordance_stats.csv"))
print(stats_tbl)




#
# RE-EXPORT — Constellation v2: constellation and combined logo as separate images
#
#   Lines between nodes (colour = chain, pattern = degree of similarity):
#     SOLID  red  = identical TRB            | SOLID  blue = identical TRA
#     DOTTED red  = TRB 1 aa difference      | DOTTED blue = TRA 1 aa difference
#     DASHED red  = TRB 3 aa differences  (extra version, THREE_AA_CLUSTERS only)
#   aa differences = same length, number of mismatched positions (indels not counted)
#
#   Node = unique (TRB, TRA) combo, clones merged, >= 5 cells (unchanged).
#   Per cluster:
#     <spec>_cluster_<id>_constellation.png           — the constellation
#     <spec>_cluster_<id>_constellation_with_3aa.png  — same, plus TRB 3-aa lines
#                                                       (THREE_AA_CLUSTERS only)
#     <spec>_cluster_<id>_logo.png                    — summary text above one TRB
#                                                       and one TRA logo, one
#                                                       sequence per node
#
#   Output: re-exports_for_figures/constellation_logo_annotated/
#
library(dplyr); library(igraph); library(ggplot2); library(packcircles)
library(scatterpie); library(ggforce); library(stringr); library(tidyr)
library(cowplot); library(ggseqlogo); library(readr); library(grid)

# ── 0. Paths ─────────────────────────────────────────────────────────────────
conga_dir <- "exports/conga/mdata_CMV"
v2_dir    <- file.path(conga_dir, "v2")

reexport_dir <- "re-exports_for_figures"
out_dir      <- file.path(reexport_dir, "constellation_logo_annotated")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Settings ──────────────────────────────────────────────────────────────
CLUSTERS_TO_PLOT <- c("0", "2", "4", "8", "18", "14", "17")
MIN_NODE_CELLS   <- 5

donor_color_map <- c("D1" = "#1b9e77", "D2" = "#d95f02", "D3" = "#7570b3",
                     "D4" = "#e7298a", "D5" = "#66a61e")

DIST_TRANSFORM <- function(d) log1p(d)
NODE_R_MIN <- 0.35; NODE_R_MAX <- 1.40; GAP_MIN <- 0.25

TRB_COL <- "#C0392B"
TRA_COL <- "#2166AC"

LW_EXACT <- 1.5
LW_1AA   <- 1.5
LW_3AA   <- 1.5

DASH_1AA <- "12"   # dotted, as in the original figures
DASH_3AA <- "44"   # long, even dashes — clearly distinct from the dotted 1-aa lines

# Clusters that also get a second constellation with TRB 3-aa lines
THREE_AA_CLUSTERS <- c("17")

# Logo height mode:
#   "bits" = column height is conservation (information content), letters scaled
#            by frequency within the column; variable positions are shorter
#   "prob" = every column full height; letters show composition only
LOGO_METHOD <- "prob"

# Logo titles: TRUE = always "x of y nodes" (e.g. "8 of 8 nodes");
# FALSE = "8 nodes" when every node is included, "x of y nodes" otherwise
ALWAYS_X_OF_Y <- TRUE

# Coordinate scaling of the constellation.
# FALSE = each cluster scaled to fill its own panel, with a floor so very small
#         clusters don't balloon. TRUE = one shared scale across all clusters.
SHARED_SCALE    <- FALSE
MIN_HALF_EXTENT <- 3
PAD_FRAC        <- 0.12

# Image sizes (inches)
CONST_W <- 8;   CONST_H <- 8
LOGO_W  <- 5.5             # logo image width
LOGOS_H <- 3.9             # height of the TRB + TRA logos together (unchanged)
FIG_DPI <- 300

CONST_TITLE     <- "Cluster %s"
LOGO_TITLE_SIZE <- 12

# Summary text above the logos.
# SUMMARY_TEXT_SIZE is the target; if the widest line across ALL clusters would
# not fit across the image at that size, every image uses the largest size that
# does fit, so the text is the same size in every logo image.
SUMMARY_TEXT_SIZE   <- 28      # pt
SUMMARY_LINEHEIGHT  <- 1.15
SUMMARY_BLANK_LINES <- 1       # blank line(s) left between the text and the TRB logo
SUMMARY_SIDE_PAD_IN <- 0.2     # clear space at each side of the widest line
SUMMARY_VPAD_PT     <- 16      # clear space above + below the text block (excl. blank lines)

# ── 2. Inputs ────────────────────────────────────────────────────────────────
g   <- readRDS(file.path(v2_dir, "constellation_v2_graph.rds")); D <- g$D
ord <- read_csv(file.path(v2_dir, "constellation_v2_tcrdist_order.csv"), show_col_types = FALSE)
seq_key_matrix <- paste(ord$cdr3_b_aa, ord$cdr3_a_aa)

# Node-size scale shared with the other constellation figures
lc_min <- log(min(g$nodes$total_cells)); lc_max <- log(max(g$nodes$total_cells))

md_all <- mdata_CMV@meta.data %>%
  dplyr::mutate(cloneId   = as.character(cloneId),
                cluster   = as.character(TCR_cluster_et64),
                global_id = as.character(Global_Clone_ID),
                donorS    = paste0("D", str_extract(as.character(donor_identity), "\\d")),
                spec      = as.character(specificity))

# ── 3. Helpers ───────────────────────────────────────────────────────────────
pull_reads <- function(strings, prefix) unlist(lapply(strings, function(s) {
  if (is.na(s)) return(character(0))
  parts <- str_split(s, ";")[[1]]
  hit   <- parts[str_detect(parts, paste0("^", prefix, ":"))]
  seqs  <- str_replace(hit, paste0("^", prefix, ":"), "")
  seqs[nchar(seqs) > 0]
}))

top1 <- function(strings, prefix) {
  r <- pull_reads(strings, prefix)
  if (!length(r)) return(NA_character_)
  names(which.max(table(r)))
}

keep_modal <- function(s) {
  if (!length(s)) return(character(0))
  l <- nchar(s); m <- as.integer(names(sort(table(l), decreasing = TRUE))[1])
  s[l == m]
}

# Number of mismatched positions between two equal-length sequences;
# NA if either is missing or the lengths differ
aa_mismatches <- function(a, b) {
  if (is.na(a) || is.na(b) || nchar(a) != nchar(b)) return(NA_integer_)
  sum(strsplit(a, "")[[1]] != strsplit(b, "")[[1]])
}

empty_edges <- function() tibble(x = numeric(), y = numeric(),
                                 xend = numeric(), yend = numeric(),
                                 share = character(), kind = character())

# Rendered width of one line of text, in inches
text_width_in <- function(txt, size) {
  grid::convertWidth(
    grid::grobWidth(grid::textGrob(txt, gp = grid::gpar(fontsize = size))),
    "inches", valueOnly = TRUE
  )
}

# ── 4. Build one cluster: nodes, layout, edges, summary counts ──────────────
build_cluster <- function(cl) {

  cc <- md_all %>% dplyr::filter(cluster == cl, !grepl("notCalled", cloneId))
  if (nrow(cc) == 0) { message(sprintf("cluster %s: no cells, skipped", cl)); return(NULL) }

  dom_spec <- { w <- table(cc$spec[!is.na(cc$spec)]); if (length(w)) names(which.max(w)) else "NA" }

  # Per-cloneId dominant TRB / TRA, then merge clones sharing the exact pair
  per_clone <- cc %>%
    dplyr::group_by(cloneId) %>%
    dplyr::summarise(trb = top1(t_cdr3s_aa, "TRB"),
                     tra = top1(t_cdr3s_aa, "TRA"), .groups = "drop") %>%
    dplyr::mutate(nkey = paste(tidyr::replace_na(trb, "NA"),
                               tidyr::replace_na(tra, "NA"), sep = "||"))

  cell_nkey <- cc %>%
    dplyr::left_join(per_clone %>% dplyr::select(cloneId, nkey, trb, tra), by = "cloneId")

  nodes <- cell_nkey %>%
    dplyr::group_by(nkey) %>%
    dplyr::summarise(trb = dplyr::first(trb), tra = dplyr::first(tra),
                     total_cells = dplyr::n(),
                     d1 = sum(donorS == "D1"), d2 = sum(donorS == "D2"),
                     d3 = sum(donorS == "D3"), d4 = sum(donorS == "D4"),
                     d5 = sum(donorS == "D5"), .groups = "drop") %>%
    dplyr::mutate(paired = !is.na(trb) & !is.na(tra)) %>%
    dplyr::filter(total_cells >= MIN_NODE_CELLS) %>%
    dplyr::mutate(
      node_r = NODE_R_MIN + (NODE_R_MAX - NODE_R_MIN) *
               (log(total_cells) - lc_min) / (lc_max - lc_min),
      mrow     = ifelse(paired, match(paste(trb, tra), seq_key_matrix), NA_integer_),
      has_dist = paired & !is.na(mrow)
    )

  if (nrow(nodes) < 1) {
    message(sprintf("cluster %s: no nodes >= %d cells, skipped", cl, MIN_NODE_CELLS))
    return(NULL)
  }
  n_unmapped <- sum(nodes$paired & !nodes$has_dist)
  if (n_unmapped)
    message(sprintf("cluster %s: %d paired node(s) not in the TCRdist matrix, placed on the outer ring",
                    cl, n_unmapped))

  # Layout: paired nodes by Kamada-Kawai on log1p(TCRdist), others on a ring,
  # then repelled so circles don't overlap (unchanged). The layout does not
  # depend on the lines, so both cluster-17 versions share identical positions.
  lay <- nodes %>% dplyr::filter(has_dist)
  unp <- nodes %>% dplyr::filter(!has_dist)
  coords <- tibble(nkey = character(), x = numeric(), y = numeric())
  if (nrow(lay) >= 2) {
    sub_t <- DIST_TRANSFORM(D[lay$mrow, lay$mrow, drop = FALSE])
    gg <- graph_from_adjacency_matrix(sub_t, mode = "undirected", weighted = TRUE, diag = FALSE)
    set.seed(1); xy <- layout_with_kk(gg, weights = E(gg)$weight); xy <- scale(xy, TRUE, FALSE)
    coords <- tibble(nkey = lay$nkey, x = xy[, 1], y = xy[, 2])
  } else if (nrow(lay) == 1) {
    coords <- tibble(nkey = lay$nkey, x = 0, y = 0)
  }
  if (nrow(unp) > 0) {
    R0 <- if (nrow(coords)) max(sqrt(coords$x^2 + coords$y^2)) + 1 else 0
    a  <- seq(0, 2 * pi, length.out = nrow(unp) + 1)[-1]
    coords <- bind_rows(coords, tibble(nkey = unp$nkey, x = R0 * cos(a), y = R0 * sin(a)))
  }
  rad <- nodes$node_r[match(coords$nkey, nodes$nkey)]
  rp  <- packcircles::circleRepelLayout(
           data.frame(x = coords$x, y = coords$y, radius = rad + GAP_MIN / 2),
           xysizecols = 1:3, sizetype = "radius", maxiter = 1000)$layout
  coords$x <- rp$x - mean(rp$x); coords$y <- rp$y - mean(rp$y)
  nodes <- nodes %>% dplyr::left_join(coords, by = "nkey")

  # Chain-sharing edges. TRB: identical, 1 aa, and (for listed clusters) 3 aa.
  # TRA: identical and 1 aa. 3-aa edges are only drawn in the "_with_3aa" version.
  has_3aa_version <- cl %in% THREE_AA_CLUSTERS

  edges <- empty_edges()
  if (nrow(nodes) >= 2) {
    cb <- t(combn(nrow(nodes), 2))
    edges <- bind_rows(lapply(seq_len(nrow(cb)), function(k) {
      i <- cb[k, 1]; j <- cb[k, 2]; rows <- list()

      mb <- aa_mismatches(nodes$trb[i], nodes$trb[j])
      kind_b <- if (is.na(mb))                       NA_character_
                else if (mb == 0)                    "exact"
                else if (mb == 1)                    "one_aa"
                else if (mb == 3 && has_3aa_version) "three_aa"
                else                                 NA_character_
      if (!is.na(kind_b)) rows[[length(rows) + 1]] <- tibble(share = "TRB", kind = kind_b)

      ma <- aa_mismatches(nodes$tra[i], nodes$tra[j])
      if (!is.na(ma) && ma <= 1)
        rows[[length(rows) + 1]] <- tibble(share = "TRA",
                                           kind = c("exact", "one_aa")[ma + 1])

      if (!length(rows)) return(NULL)
      bind_rows(rows) %>% dplyr::mutate(x = nodes$x[i], y = nodes$y[i],
                                        xend = nodes$x[j], yend = nodes$y[j])
    }))
    if (nrow(edges) == 0) edges <- empty_edges()
  }

  if (has_3aa_version && !any(edges$kind == "three_aa"))
    message(sprintf("cluster %s: no TRB pairs with exactly 3 aa differences", cl))

  # Summary counts over the cells in the plotted nodes
  kept_cells <- cell_nkey %>% dplyr::filter(nkey %in% nodes$nkey)
  n_major <- dplyr::n_distinct(kept_cells$global_id[!is.na(kept_cells$global_id)])
  n_minor <- dplyr::n_distinct(kept_cells$cloneId[is.na(kept_cells$global_id)])
  n_donor <- dplyr::n_distinct(kept_cells$donorS)

  list(
    cl              = cl,
    dom_spec        = dom_spec,
    nodes           = nodes,
    edges           = edges,
    has_3aa_version = has_3aa_version,
    n_major         = n_major,
    n_minor         = n_minor,
    n_donor         = n_donor,
    summary_lab     = sprintf("%d major and %d minor clones\nfrom %d donor%s represented",
                              n_major, n_minor, n_donor, ifelse(n_donor == 1, "", "s")),
    extent          = max(sqrt(nodes$x^2 + nodes$y^2) + nodes$node_r)
  )
}

# ── 5. Combined per-cluster logo (one sequence per node = equal weights) ────
make_chain_logo <- function(seqs, chain_lab) {
  seqs <- seqs[!is.na(seqs) & nchar(seqs) > 0]
  if (!length(seqs))
    return(cowplot::ggdraw() +
             cowplot::draw_label(paste0(chain_lab, ": not recovered"),
                                 size = 11, colour = "grey50"))

  keep    <- keep_modal(seqs)
  n_keep  <- length(keep)
  n_total <- length(seqs)

  # Node count only: "x of y nodes" (or "y nodes" when ALWAYS_X_OF_Y is FALSE
  # and every node is included)
  n_lab <- if (ALWAYS_X_OF_Y || n_keep < n_total) {
    sprintf("%d of %d node%s", n_keep, n_total, ifelse(n_total == 1, "", "s"))
  } else {
    sprintf("%d node%s", n_keep, ifelse(n_keep == 1, "", "s"))
  }

  ggseqlogo(keep, method = LOGO_METHOD) +
    ggtitle(sprintf("%s  (%s)", chain_lab, n_lab)) +
    theme(legend.position = "none",
          plot.title   = element_text(size = LOGO_TITLE_SIZE, face = "bold"),
          axis.text    = element_text(size = 8),
          axis.title.y = element_text(size = 9),
          axis.title.x = element_blank(),
          plot.margin  = margin(4, 8, 4, 8))
}

# ── 6. Constellation plot (with or without the 3-aa lines) ──────────────────
make_const_plot <- function(cd, half_lim, show_3aa = FALSE) {

  nodes <- cd$nodes; edges <- cd$edges

  edge_layer <- function(sh, kd, col, lw, lty, a) {
    d <- dplyr::filter(edges, share == sh, kind == kd)
    if (!nrow(d)) return(NULL)
    geom_segment(data = d, aes(x = x, y = y, xend = xend, yend = yend),
                 colour = col, linewidth = lw, linetype = lty, alpha = a)
  }

  # Solid first, then dashed/dotted on top so lines sharing a segment stay visible
  ggplot() +
    edge_layer("TRB", "exact", TRB_COL, LW_EXACT, "solid", 0.85) +
    edge_layer("TRA", "exact", TRA_COL, LW_EXACT, "solid", 0.85) +
    { if (show_3aa) edge_layer("TRB", "three_aa", TRB_COL, LW_3AA, DASH_3AA, 0.80) } +
    edge_layer("TRB", "one_aa", TRB_COL, LW_1AA, DASH_1AA, 0.80) +
    edge_layer("TRA", "one_aa", TRA_COL, LW_1AA, DASH_1AA, 0.80) +
    ggforce::geom_circle(data = nodes, aes(x0 = x, y0 = y, r = node_r),
                         fill = "white", color = "grey30", linewidth = 0.2) +
    scatterpie::geom_scatterpie(data = nodes, aes(x = x, y = y, r = node_r),
                                cols = c("d1", "d2", "d3", "d4", "d5"), color = NA) +
    scale_fill_manual(values = setNames(unname(donor_color_map), c("d1", "d2", "d3", "d4", "d5")),
                      labels = c("D1", "D2", "D3", "D4", "D5"), name = "Donor", guide = "none") +
    coord_equal(xlim = c(-half_lim, half_lim), ylim = c(-half_lim, half_lim), clip = "off") +
    theme_void() +
    ggtitle(sprintf(CONST_TITLE, cd$cl)) +
    theme(plot.title  = element_text(size = 16, hjust = 0.5),
          plot.margin = margin(8, 8, 8, 8))
}

# ── 7. Render one cluster ────────────────────────────────────────────────────
render_cluster <- function(cd, half_lim, summary_size, text_h_in) {

  nodes <- cd$nodes

  # Constellation(s): the standard version for every cluster, plus a
  # "_with_3aa" version for clusters in THREE_AA_CLUSTERS
  const_base <- sprintf("%s_cluster_%s_constellation", cd$dom_spec, cd$cl)

  ggsave(file.path(out_dir, paste0(const_base, ".png")),
         make_const_plot(cd, half_lim, show_3aa = FALSE),
         width = CONST_W, height = CONST_H, dpi = FIG_DPI, bg = "white")

  if (cd$has_3aa_version) {
    ggsave(file.path(out_dir, paste0(const_base, "_with_3aa.png")),
           make_const_plot(cd, half_lim, show_3aa = TRUE),
           width = CONST_W, height = CONST_H, dpi = FIG_DPI, bg = "white")
  }

  # Logo image: summary text above the TRB and TRA logos. The text is anchored
  # to the top of its block; the space below it (SUMMARY_BLANK_LINES lines)
  # separates it from the TRB logo. The logos stay LOGOS_H tall.
  y_top <- 1 - (SUMMARY_VPAD_PT / 2) / (text_h_in * 72)

  text_panel <- cowplot::ggdraw() +
    cowplot::draw_label(cd$summary_lab, x = 0.5, y = y_top,
                        hjust = 0.5, vjust = 1,
                        size = summary_size, lineheight = SUMMARY_LINEHEIGHT)

  logo_img <- cowplot::plot_grid(
    text_panel,
    make_chain_logo(nodes$trb, "TRB"),
    make_chain_logo(nodes$tra, "TRA"),
    ncol = 1, rel_heights = c(text_h_in, LOGOS_H / 2, LOGOS_H / 2)
  )

  fn_logo <- file.path(out_dir, sprintf("%s_cluster_%s_logo.png", cd$dom_spec, cd$cl))
  ggsave(fn_logo, logo_img, width = LOGO_W, height = text_h_in + LOGOS_H,
         dpi = FIG_DPI, bg = "white")

  e <- cd$edges
  message(sprintf(
    "[%s] cluster %s: %d nodes | %d major, %d minor, %d donors | edges: %d exact, %d 1aa%s",
    cd$dom_spec, cd$cl, nrow(nodes), cd$n_major, cd$n_minor, cd$n_donor,
    sum(e$kind == "exact"), sum(e$kind == "one_aa"),
    if (cd$has_3aa_version)
      sprintf(", %d TRB 3aa (in _with_3aa version only)", sum(e$kind == "three_aa")) else ""))
}

# ── 8. Build all clusters, size the summary text once, then render ──────────
cluster_data <- Filter(Negate(is.null), lapply(CLUSTERS_TO_PLOT, build_cluster))
if (!length(cluster_data)) stop("No clusters could be built")

global_lim <- max(vapply(cluster_data, function(cd) cd$extent, numeric(1))) * (1 + PAD_FRAC)

# One summary font size for every logo image: the target, or the largest size
# at which the widest line of any cluster still fits across the image
all_lines <- unlist(strsplit(vapply(cluster_data, function(cd) cd$summary_lab, character(1)), "\n"))
widest_in <- max(vapply(all_lines, text_width_in, numeric(1), size = SUMMARY_TEXT_SIZE))
avail_in  <- LOGO_W - 2 * SUMMARY_SIDE_PAD_IN

summary_size <- if (widest_in > avail_in) SUMMARY_TEXT_SIZE * avail_in / widest_in else SUMMARY_TEXT_SIZE

# Text block height: two lines of text, the blank line(s) below, and padding
text_h_in <- (SUMMARY_VPAD_PT +
              (2 + SUMMARY_BLANK_LINES) * summary_size * SUMMARY_LINEHEIGHT) / 72

message(sprintf("Summary text: %.1f pt%s | logo image %.2f x %.2f in",
                summary_size,
                if (summary_size < SUMMARY_TEXT_SIZE)
                  sprintf(" (capped from %d pt to fit the widest line)", SUMMARY_TEXT_SIZE) else "",
                LOGO_W, text_h_in + LOGOS_H))

for (cd in cluster_data) {
  half_lim <- if (SHARED_SCALE) global_lim else max(cd$extent * (1 + PAD_FRAC), MIN_HALF_EXTENT)
  render_cluster(cd, half_lim, summary_size, text_h_in)
}

# ── 9. Standalone legends ────────────────────────────────────────────────────
# Amino-acid chemistry key (all 5 chemistry classes)
all_aa <- paste(c("A","R","N","D","C","Q","E","G","H","I",
                  "L","K","M","F","P","S","T","W","Y","V"), collapse = "")
chem_leg <- cowplot::get_legend(
  ggseqlogo(c(all_aa, all_aa), method = "bits") +
    theme(legend.position = "right",
          legend.title = element_text(size = 12), legend.text = element_text(size = 11)))
ggsave(file.path(out_dir, "logo_chemistry_legend.png"),
       cowplot::ggdraw(chem_leg), width = 3, height = 2.2, dpi = FIG_DPI, bg = "white")

# Line-type keys: one standard, one for the "_with_3aa" constellations
line_styles <- tibble::tribble(
  ~lab,                     ~col,     ~lty,     ~lw,
  "Identical TRB",          TRB_COL,  "solid",  LW_EXACT,
  "TRB: 1 aa difference",   TRB_COL,  DASH_1AA, LW_1AA,
  "TRB: 3 aa differences",  TRB_COL,  DASH_3AA, LW_3AA,
  "Identical TRA",          TRA_COL,  "solid",  LW_EXACT,
  "TRA: 1 aa difference",   TRA_COL,  DASH_1AA, LW_1AA
)

save_line_legend <- function(labels, filename) {
  st  <- line_styles %>% dplyr::filter(lab %in% labels)
  lvl <- st$lab
  df  <- tibble(lab = factor(lvl, levels = lvl), y = rev(seq_along(lvl)))

  leg <- cowplot::get_legend(
    ggplot(df) +
      geom_segment(aes(x = 0, xend = 1, y = y, yend = y,
                       colour = lab, linetype = lab, linewidth = lab)) +
      scale_colour_manual(name = "Shared chain",    values = setNames(st$col, lvl)) +
      scale_linetype_manual(name = "Shared chain",  values = setNames(st$lty, lvl)) +
      scale_linewidth_manual(name = "Shared chain", values = setNames(st$lw,  lvl)) +
      theme_void() +
      theme(legend.position  = "right",
            legend.key.width = unit(2, "cm"),
            legend.title     = element_text(size = 12),
            legend.text      = element_text(size = 11)))

  ggsave(file.path(out_dir, filename), cowplot::ggdraw(leg),
         width = 4.5, height = 0.4 * length(lvl) + 0.6, dpi = FIG_DPI, bg = "white")
}

save_line_legend(c("Identical TRB", "TRB: 1 aa difference",
                   "Identical TRA", "TRA: 1 aa difference"),
                 "line_type_legend.png")

save_line_legend(c("Identical TRB", "TRB: 1 aa difference", "TRB: 3 aa differences",
                   "Identical TRA", "TRA: 1 aa difference"),
                 sprintf("line_type_legend_cluster%s.png",
                         paste(THREE_AA_CLUSTERS, collapse = "_")))

message("\nDone: outputs in ", out_dir)




#
#
# RE-EXPORT — Cluster cards v2
#   Header : one line of large text,
#            "Cluster x - pp65 - kNN bias: Effect size x, p = y",
#            followed by a blank line before the plots
#   Left   : per-clone KDE contours over the full-object UMAP (unchanged)
#   Right  : per-clone metacluster composition bars (unchanged)
#   Clone labels in the UMAP legend and on the bar x axis:
#            "Donor 3 - pp65 clone 1"
#   The bar legend always shows all five metaclusters.
#   Cards for clusters with >= 2 unique TCR sequences and >= 1 major clone with
#   >= MIN_CELLS cells in the cluster (same gate as before).
#   Output: re-exports_for_figures/cluster_cards_v2/
#
#
library(Seurat); library(dplyr); library(stringr); library(tidyr)
library(ggplot2); library(cowplot); library(readr)

# ── 0. Paths ─────────────────────────────────────────────────────────────────
conga_dir <- "exports/conga/mdata_CMV"
knn_file  <- "exports/conga/conga_exports/tcrdist_phenotype_test/tcrdist_phenotype_specstrat_results.csv"

reexport_dir <- "re-exports_for_figures"
out_dir      <- file.path(reexport_dir, "cluster_cards_v2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── 1. Settings ──────────────────────────────────────────────────────────────
MIN_CELLS   <- 5
META_LEVELS <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")

# Project-standard metacluster palette (matches the metacluster UMAP and the
# corrected_colors stacked bars)
meta_palette <- c(
  "Tcm"     = "#F8766D",
  "Tem"     = "#A3A500",
  "T-CMC1"  = "#00BF7D",
  "NK-like" = "#00B0F6",
  "ISG+"    = "#E76BF3"
)

# kNN bias values shown in the header (specificity-stratified test).
# Same-cluster metric by default; swap to "cross_donor_effect" /
# "cross_donor_p_adj" for the cross-donor version, or "_p_emp" for unadjusted p.
KNN_EFFECT_COL <- "same_cluster_effect"
KNN_P_COL      <- "same_cluster_p_adj"    # BH-adjusted across clusters

# Header text and spacing.
# HEADER_LINE_H is the row height the text sits in (as before); the blank
# line(s) are added below it, so the text keeps its distance from the top edge.
HEADER_TEXT_SIZE   <- 24                    # pt
HEADER_LINE_H      <- 1.8                   # text row height, in multiples of the font size
HEADER_BLANK_LINES <- 1                     # blank line(s) between the text and the plots
HEADER_LINEHEIGHT  <- 1.2                   # height of one blank line, in multiples of the font size

header_h_pt <- HEADER_TEXT_SIZE * (HEADER_LINE_H + HEADER_BLANK_LINES * HEADER_LINEHEIGHT)
HEADER_H_IN <- header_h_pt / 72
HEADER_Y    <- 1 - (HEADER_TEXT_SIZE * HEADER_LINE_H / 2) / header_h_pt   # text centre (npc)

# Bottom row and card size (inches)
BOTTOM_H_IN <- 6.0
CARD_W      <- 13
CARD_DPI    <- 200

# Bar x-axis labels
BAR_LABEL_ANGLE <- 45      # steeper than the original 30 so longer labels stay clear of the y axis

# ── 2. Inputs ────────────────────────────────────────────────────────────────
md <- mdata_CMV@meta.data

cell_tbl <- tibble(
  cluster     = as.character(md$TCR_cluster_et64),
  cloneId     = as.character(md$cloneId),
  spec        = as.character(md$specificity),
  metacluster = as.character(md$metacluster),
  global_id   = as.character(md$Global_Clone_ID)) %>%
  dplyr::filter(!is.na(cluster), !grepl("notCalled", cloneId)) %>%
  dplyr::mutate(is_global = !is.na(global_id))

# Clusters with >= 2 unique TCR sequences (the constellation min2seq set)
g_data <- readRDS(file.path(conga_dir, "constellation_graph.rds"))
keep_clusters <- g_data$nodes %>%
  dplyr::distinct(cluster, mrow) %>%
  dplyr::count(cluster, name = "n_seq") %>%
  dplyr::filter(n_seq >= 2) %>%
  dplyr::mutate(cluster = as.character(cluster)) %>%
  dplyr::pull(cluster)
message(sprintf("Clusters with >=2 unique sequences: %d", length(keep_clusters)))

# Dominant specificity per cluster (cell-weighted), used in the header and filename
cluster_spec <- cell_tbl %>%
  dplyr::filter(!is.na(spec)) %>%
  dplyr::count(cluster, spec, name = "n") %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(dom_spec = spec[which.max(n)], .groups = "drop")

# kNN bias results
knn_res <- read_csv(knn_file, comment = "#", show_col_types = FALSE) %>%
  dplyr::mutate(cluster = as.character(cluster))
stopifnot(all(c(KNN_EFFECT_COL, KNN_P_COL) %in% names(knn_res)))

knn_lookup <- knn_res %>%
  dplyr::transmute(cluster,
                   effect = .data[[KNN_EFFECT_COL]],
                   p      = .data[[KNN_P_COL]])

# UMAP coordinates for every cell (backdrop + clone contours)
umap_all <- as.data.frame(Embeddings(mdata_CMV, "umap"))
colnames(umap_all)[1:2] <- c("UMAP_1", "UMAP_2")
umap_all$global_id <- as.character(md$Global_Clone_ID)

# ── 3. Helpers ───────────────────────────────────────────────────────────────
# "D3_pp65-clone-1" -> "Donor 3 - pp65 clone 1", ordered by donor, then
# specificity, then clone number
clone_labels <- function(major_ids) {
  info <- tibble(
    global_id = major_ids,
    don  = str_match(major_ids, "^D(\\d+)")[, 2],
    spec = str_match(major_ids, "_(pp65|IE1)-")[, 2],
    num  = suppressWarnings(as.integer(str_match(major_ids, "clone-(\\d+)")[, 2]))
  ) %>%
    dplyr::mutate(
      label = ifelse(is.na(don) | is.na(spec) | is.na(num), global_id,
                     sprintf("Donor %s - %s clone %d", don, spec, num))
    ) %>%
    dplyr::arrange(as.integer(don), factor(spec, levels = c("pp65", "IE1")), num)
  setNames(info$label, info$global_id)
}

fmt_p <- function(p) {
  if (is.na(p)) "p = NA" else if (p < 0.001) "p < 0.001" else sprintf("p = %.3f", p)
}

make_header <- function(cl, dom) {
  row <- knn_lookup[knn_lookup$cluster == cl, ]
  knn_txt <- if (!nrow(row) || is.na(row$effect[1])) {
    "kNN bias: not tested"
  } else {
    sprintf("kNN bias: Effect size %.2f, %s", row$effect[1], fmt_p(row$p[1]))
  }
  sprintf("Cluster %s - %s - %s", cl, dom, knn_txt)
}

# ── 4. Contour UMAP (unchanged apart from the clone labels) ─────────────────
make_umap <- function(major_ids) {
  lab <- clone_labels(major_ids)
  clone_cells <- umap_all %>%
    dplyr::filter(global_id %in% major_ids) %>%
    dplyr::mutate(clone = factor(lab[global_id], levels = unname(lab)))

  ggplot() +
    geom_point(data = umap_all, aes(UMAP_1, UMAP_2),
               color = "grey70", size = 0.35, stroke = 0) +
    geom_density_2d(data = clone_cells, aes(UMAP_1, UMAP_2, color = clone),
                    linewidth = 0.4, bins = 10) +
    scale_color_brewer(palette = "Set1", name = "Major clone", drop = FALSE) +
    theme_classic(base_size = 11) +
    theme(legend.position = "right",
          legend.key.size = unit(0.8, "lines"),
          aspect.ratio    = 1) +
    labs(title = NULL, x = "UMAP 1", y = "UMAP 2") +
    guides(color = guide_legend(override.aes = list(linewidth = 1.2)))
}

# ── 5. Metacluster proportion bars ───────────────────────────────────────────
make_bars <- function(cl, major_ids) {
  lab <- clone_labels(major_ids)

  bd <- cell_tbl %>%
    dplyr::filter(cluster == cl, global_id %in% major_ids,
                  metacluster %in% META_LEVELS) %>%
    dplyr::mutate(clone       = factor(lab[global_id], levels = unname(lab)),
                  metacluster = factor(metacluster, levels = META_LEVELS)) %>%
    dplyr::count(clone, metacluster, name = "n") %>%
    dplyr::group_by(clone) %>%
    dplyr::mutate(prop = 100 * n / sum(n)) %>%
    dplyr::ungroup()

  # Legend-only layer: one zero-height, invisible bar per metacluster, so all
  # five keys are filled even when a cluster has no cells in some of them
  legend_df <- tibble(
    clone       = factor(as.character(bd$clone[1]), levels = levels(bd$clone)),
    prop        = 0,
    metacluster = factor(META_LEVELS, levels = META_LEVELS)
  )

  ggplot(bd, aes(clone, prop, fill = metacluster)) +
    geom_col(width = 0.75, color = "white", linewidth = 0.2) +
    geom_col(data = legend_df, width = 0.75, alpha = 0) +
    scale_fill_manual(values = meta_palette, breaks = META_LEVELS,
                      limits = META_LEVELS, drop = FALSE, name = "Metacluster") +
    guides(fill = guide_legend(override.aes = list(alpha = 1, colour = "white",
                                                   linewidth = 0.2))) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(title = NULL, x = NULL, y = "% of cells") +
    theme_classic(base_size = 11) +
    theme(axis.text.x     = element_text(angle = BAR_LABEL_ANGLE, hjust = 1, vjust = 1),
          legend.key.size = unit(0.8, "lines"),
          plot.margin     = margin(5.5, 5.5, 5.5, 20))
}

# ── 6. Compose one card ──────────────────────────────────────────────────────
render_card <- function(cl) {

  major_ids <- cell_tbl %>%
    dplyr::filter(cluster == cl, is_global) %>%
    dplyr::count(global_id, name = "n") %>%
    dplyr::filter(n >= MIN_CELLS) %>%
    dplyr::pull(global_id)

  if (length(major_ids) == 0) {
    message(sprintf("cluster %s: no major clones >= %d cells, skipped", cl, MIN_CELLS))
    return(invisible())
  }

  dom <- cluster_spec$dom_spec[match(cl, cluster_spec$cluster)]
  dom <- ifelse(is.na(dom), "NA", dom)

  # Header text sits in the top HEADER_LINE_H of the row; the blank line(s)
  # below it separate it from the plots
  header_txt   <- make_header(cl, dom)
  header_panel <- cowplot::ggdraw() +
    cowplot::draw_label(header_txt, x = 0.5, y = HEADER_Y,
                        hjust = 0.5, vjust = 0.5,
                        size = HEADER_TEXT_SIZE, fontface = "bold")

  bottom <- cowplot::plot_grid(make_umap(major_ids), make_bars(cl, major_ids),
                               ncol = 2, rel_widths = c(0.65, 0.35))

  card <- cowplot::plot_grid(header_panel, bottom, ncol = 1,
                             rel_heights = c(HEADER_H_IN, BOTTOM_H_IN))

  fn <- file.path(out_dir, sprintf("%s_cluster_%s_card.png", dom, cl))
  ggsave(fn, card, width = CARD_W, height = HEADER_H_IN + BOTTOM_H_IN,
         dpi = CARD_DPI, bg = "white")

  message(sprintf("[%s] cluster %s: %d major clones | %s -> %s",
                  dom, cl, length(major_ids), header_txt, basename(fn)))
}

# ── 7. Render all qualifying clusters ────────────────────────────────────────
kc <- keep_clusters[order(suppressWarnings(as.integer(keep_clusters)), keep_clusters)]
invisible(lapply(kc, render_card))
message(sprintf("\nDone: cards in %s", out_dir))





library(Seurat)
library(ggplot2)

# ---- Settings ----
obj             <- mdata_CMV        # swap to D0_data_CMV if needed
gene            <- "ZEB2"           # human symbol, uppercase
metacluster_col <- "metacluster"    # meta.data column holding the 5 metaclusters
reduction       <- "umap"
out_dir         <- "ZEB2_plots"

mc_levels <- c("Tcm", "Tem", "T-CMC1", "NK-like", "ISG+")
mc_colors <- c(                     # replace with your project metacluster hex codes
  "Tcm"     = "#4DAF4A",
  "Tem"     = "#377EB8",
  "T-CMC1"  = "#E41A1C",
  "NK-like" = "#984EA3",
  "ISG+"    = "#FF7F00"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Checks ----
DefaultAssay(obj) <- "RNA"          # gene expression, not FB
stopifnot(gene %in% rownames(obj))
stopifnot(metacluster_col %in% colnames(obj@meta.data))

mc_vals    <- as.character(obj@meta.data[[metacluster_col]])
unexpected <- setdiff(na.omit(unique(mc_vals)), mc_levels)
if (length(unexpected)) {
  stop("Labels not in mc_levels: ", paste(unexpected, collapse = ", "))
}
obj@meta.data[[metacluster_col]] <- factor(mc_vals, levels = mc_levels)

# ---- UMAP: ZEB2 highlighted ----
p_umap <- FeaturePlot(
  obj,
  features   = gene,
  reduction  = reduction,
  order      = TRUE,                # expressing cells drawn on top
  pt.size    = 0.4,
  max.cutoff = "q99",               # stop a few outliers compressing the scale
  cols       = c("grey85", "#B2182B")
) +
  ggtitle(gene) +
  theme(aspect.ratio = 1)

# ---- Violins: ZEB2 by metacluster ----
obj_vln <- if (anyNA(obj@meta.data[[metacluster_col]])) {
  subset(obj, cells = colnames(obj)[!is.na(obj@meta.data[[metacluster_col]])])
} else {
  obj
}

p_vln <- VlnPlot(
  obj_vln,
  features = gene,
  group.by = metacluster_col,
  cols     = mc_colors[mc_levels],
  pt.size  = 0
) +
  stat_summary(fun = median, geom = "point", size = 2, colour = "black") +
  labs(x = NULL, y = paste(gene, "expression (log-normalized)")) +
  NoLegend()

# ---- Show and save ----
print(p_umap)
print(p_vln)

ggsave(file.path(out_dir, paste0(gene, "_UMAP.png")),
       p_umap, width = 5, height = 5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0(gene, "_violin_by_metacluster.png")),
       p_vln, width = 5, height = 4, dpi = 300, bg = "white")


library(Seurat)
library(ggplot2)
library(dplyr)

# ---- Settings ----
obj       <- D0_data_CMV
gene      <- "ZEB2"
visit_col <- "timepoint"          # meta.data column holding the visit
spec_col  <- "specificity"    # meta.data column holding pp65 / IE1
clone_col <- "Global_Clone_ID"        # meta.data column holding clone IDs
cmv_specs <- c("pp65", "IE1")
n_top     <- 20
out_dir   <- "ZEB2_plots/D0"

visit_levels  <- c("Visit 1", "Visit 2", "Visit 3", "Visit 4")   # values as stored in visit_col
visit_display <- c("Visit 1\n(Wk 20)", "Visit 2\n(Wk 21)", "Visit 3\n(Wk 32)", "Visit 4\n(Wk 69)")
visit_colors  <- setNames(                                        # placeholder: swap in your timepoint colours
  c("#C6DBEF", "#6BAED6", "#2171B5", "#08306B"), visit_display
)

pp65_clone_palette <- c(
  "#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
  "#e6ab02", "#a6761d", "#666666", "#1f78b4", "#b2df8a",
  "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
  "#cab2d6", "#6a3d9a", "#ffff99", "#b15928", "#8dd3c7"
)

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Checks ----
DefaultAssay(obj) <- "RNA"
stopifnot(gene %in% rownames(obj))

missing_cols <- setdiff(c(visit_col, spec_col, clone_col), colnames(obj@meta.data))
if (length(missing_cols)) {
  stop("Missing meta.data columns: ", paste(missing_cols, collapse = ", "))
}

fetched <- FetchData(obj, vars = c(gene, visit_col, spec_col, clone_col))
df <- data.frame(
  expr  = fetched[[gene]],
  visit = as.character(fetched[[visit_col]]),
  spec  = as.character(fetched[[spec_col]]),
  clone = as.character(fetched[[clone_col]])
)

unexpected_visits <- setdiff(na.omit(unique(df$visit)), visit_levels)
if (length(unexpected_visits)) {
  stop("Visit values not in visit_levels: ", paste(unexpected_visits, collapse = ", "))
}
if (!all(cmv_specs %in% df$spec)) {
  stop("Specificity values found: ", paste(unique(df$spec), collapse = ", "))
}

df <- df |>
  filter(spec %in% cmv_specs, !is.na(visit)) |>
  mutate(visit = factor(visit, levels = visit_levels, labels = visit_display))

# ---- Top pp65 clones (all visits combined) ----
top_clones <- df |>
  filter(spec == "pp65", !is.na(clone), tolower(clone) != "na") |>
  count(clone, name = "n_cells") |>
  arrange(desc(n_cells), clone) |>
  slice_head(n = n_top) |>
  mutate(clone_label = paste("Clone", row_number()))

print(top_clones, n = n_top)
write.csv(top_clones, file.path(out_dir, "top20_pp65_clones.csv"), row.names = FALSE)

clone_levels <- top_clones$clone_label
clone_colors <- setNames(pp65_clone_palette[seq_along(clone_levels)], clone_levels)

df_clone <- df |>
  filter(spec == "pp65") |>
  inner_join(select(top_clones, clone, clone_label), by = "clone") |>
  mutate(clone_label = factor(clone_label, levels = clone_levels))

# ---- Violin helper ----
violin_plot <- function(d, group, fills, title, angle = 45) {
  y_lab <- max(d$expr, na.rm = TRUE) * 1.05

  labs_df <- d |>
    group_by(across(all_of(group)), .drop = FALSE) |>
    summarise(n = n(), pct = 100 * mean(expr > 0), .groups = "drop") |>
    mutate(
      label = ifelse(n > 0, sprintf("n=%d\n%.0f%%", n, pct), "n=0"),
      y     = y_lab
    )

  ggplot(d, aes(x = .data[[group]], y = expr, fill = .data[[group]])) +
    geom_violin(scale = "width", colour = "grey25", linewidth = 0.3) +
    stat_summary(fun = median, geom = "point", size = 1.5, colour = "black") +
    geom_text(
      data = labs_df,
      aes(x = .data[[group]], y = y, label = label),
      inherit.aes = FALSE, size = 2.5, vjust = 0, lineheight = 0.9
    ) +
    scale_x_discrete(drop = FALSE) +
    scale_fill_manual(values = fills, drop = FALSE) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.2))) +
    labs(title = title, x = NULL, y = paste(gene, "expression (log-normalized)")) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.x     = element_text(angle = angle, hjust = if (angle == 0) 0.5 else 1)
    )
}

# ---- Plots ----
p_visit <- violin_plot(
  df, "visit", visit_colors,
  paste(gene, "- all CMV-specific cells by visit"),
  angle = 0
)

p_clone_all <- violin_plot(
  df_clone, "clone_label", clone_colors,
  paste(gene, "- top", n_top, "pp65 clones, all visits")
)

p_clone_v4 <- violin_plot(
  filter(df_clone, visit == visit_display[4]), "clone_label", clone_colors,
  paste(gene, "- top", n_top, "pp65 clones, Visit 4 only")
)

# ---- Show and save ----
print(p_visit)
print(p_clone_all)
print(p_clone_v4)

ggsave(file.path(out_dir, paste0(gene, "_CMV_by_visit.png")),
       p_visit, width = 5, height = 4.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0(gene, "_top20_pp65_clones_all_visits.png")),
       p_clone_all, width = 10, height = 4.5, dpi = 300, bg = "white")
ggsave(file.path(out_dir, paste0(gene, "_top20_pp65_clones_visit4.png")),
       p_clone_v4, width = 10, height = 4.5, dpi = 300, bg = "white")

