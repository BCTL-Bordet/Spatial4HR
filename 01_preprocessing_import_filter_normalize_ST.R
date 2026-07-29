# 01_import_filter_normalize_ST.R

# reading the annotation-fraction files
# merging annotations into metadata
# histopathology-based spot filtering


####################################
## ST analysis using semla
####################################

# input : spaceranger output
# output: object_merged

library(Seurat)
library(SeuratObject)
library(semla)
library(hdf5r)
library(stringr)
library(data.table)
library(Polychrome)
library(dplyr)
library(survival)
library(magrittr)
library(tibble)
library(patchwork)
library(parallel)


root <- "/mnt/bctl/.../ductals_st_data/spaceranger/ST"
annot_classes <- c(
  "Tumor", "Necrosis", "Fat_tissue", "High_TILs_stroma",
  "Cellular_stroma", "Acellular_stroma", "Vessels",
  "Artefact", "Canal_galactophore", "Nodule_lymphoid",
  "In_situ", "Nerve", "Lymphocyte", "Hole",
  "Microcalcification", "Out", "Apocrine metaplasia"
)
info_colnames <- c(
  "barcode", "in_tissue", "array_row", "array_col",
  "col_pxl", "row_pxl"
)
coords_colnames <- c(
  "barcode", "in_tissue", "coord1", "coord2",
  "pxl_row_in_fullres", "pxl_col_in_fullres"
)


read_samples <- function(st_id) {
  # ======================== Read Data from Filesystem ========================
  # load sample
  samples <- paste0(root, st_id, "/outs/filtered_feature_bc_matrix.h5")
  spotfiles <- paste0(root, st_id, "/outs/spatial/tissue_positions_list.csv")
  imgs <- paste0(root, st_id, "/outs/spatial/tissue_hires_image.png")
  json <- paste0(root, st_id, "/outs/spatial/scalefactors_json.json")
  info_table <- as.data.frame(cbind(samples, spotfiles, imgs, json))
  
  # append into list of samples
  st_sample <- ReadVisiumData(
    info_table,
    platform = "Visium",
    min.cells = 5,
    min.features = 200,
    type = "2"
  ) 
  
  # set orig ident
  st_sample@meta.data$orig.ident <- st_id
  
  # read annot data
  annotations <- fread(
    paste0(root, st_id, "/outs/spatial/tissue_positions_list_annotation.csv")
  )
  
  # read coordinates
  coordinates <- read.csv(
    paste0(root, st_id, "/outs/spatial/tissue_positions_list.csv"),
    header = FALSE,
  )
  
  # ======================== Merge Metadata ========================
  
  # filter annotation spots by st spots
  st_sample@meta.data %<>%
    mutate(barcode = rownames(.))
  
  st_spots <- rownames(st_sample@meta.data)
  
  annotations %<>%
    set_colnames(c(info_colnames, annot_classes)) %>%
    filter(barcode %in% st_spots) %>%
    select(all_of(c(annot_classes, "barcode"))) %>%
    arrange(match(barcode, st_spots))
  
  coordinates %<>%
    set_colnames(coords_colnames) %>%
    filter(barcode %in% st_spots) %>%
    arrange(match(barcode, st_spots))
  
  st_sample@meta.data <- Reduce(
    function(x, y) {
      merge(x, y, by = "barcode", all = TRUE)
    },
    list(st_sample@meta.data, annotations, coordinates)
  )
  rownames(st_sample@meta.data) <- st_sample@meta.data$barcode
  
  # filtering spots
  st_sample %<>%
    SubsetSTData(Hole < 0.3) %>%
    SubsetSTData(Artefact < 0.3) %>%
    SubsetSTData(Out < 0.3)
  
  # filtering genes and loading images
  genes <- rownames(st_sample)
  non_meta_genes <- genes[!(grepl("RPL", genes) | grepl("RPS", genes) | grepl("MT-", genes) | grepl("MTRNR", genes))]
  
  st_sample %<>%
    SubsetSTData(features = non_meta_genes)
  
  st_sample <- LoadImages(st_sample, verbose = TRUE, time.resolve = FALSE)
  st_sample <- MaskImages(st_sample)
  
  return(st_sample)
}

res <- lapply(c(73:95, 97:159), read_samples)


list_objects_ductals <- res

for (st_id in 1:length(list_objects_ductals)) {
  # Ensure the correct assay is set
  DefaultAssay(list_objects_ductals[[st_id]]) <- "Spatial"
  
  # Perform SCTransform normalization
  list_objects_ductals[[st_id]] <- SCTransform(list_objects_ductals[[st_id]], assay = "Spatial")
}

# merging data
object_merged <- MergeSTData(
  x = list_objects_ductals[[1]], y = list_objects_ductals[2:length(list_objects_ductals)], # nolint: line_length_linter.
  add.spot.ids = paste0("sample", c(73:95, 97:159))
)

object_merged <- RunPCA(object_merged, assay = "SCT", verbose = FALSE, features = rownames(object_merged@assays$SCT@scale.data))
ElbowPlot(object_merged)

object_merged <- FindNeighbors(object_merged, reduction = "pca", dims = 1:10)
object_merged <- FindClusters(object_merged, resolution = 0.4, verbose = FALSE)
object_merged <- RunUMAP(object_merged, reduction = "pca", dims = 1:10)

DimPlot(object_merged, reduction = "umap", group.by = "seurat_clusters", label = F)
DimPlot(object_merged, reduction = "umap", group.by = "orig.ident", label = F)
table(object_merged@meta.data$orig.ident, object_merged@meta.data$seurat_clusters)

# PrepSCTFindMarkers {Seurat} R Documentation
object_merged <- PrepSCTFindMarkers(object_merged, assay = "SCT", verbose = FALSE)

# find markers of all clusters
all_markers <- FindAllMarkers(object_merged)

# saveRDS(object_merged, 'object_merged.RDS')