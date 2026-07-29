# 05_CARD.R


####################################
## Cell type deconvolution using CARD
####################################

## using wu-tang reference
## (adipocytes from tang and celltype minor from wu, epithelials merged together)

library(CARD)
library(Seurat)
library(SeuratObject)
library(semla)
library(hdf5r)
library(stringr)
library(data.table)
library(Polychrome)
library(dplyr)
library(pgirmess)
library(survival)
library(patchwork)
library(tidyverse)
library(magrittr)
library(tibble)
library(patchwork)
library(parallel)
library(TOAST)
library(MuSiC)


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
    # minGenesPerSpot = 200,
    # minSpotsPerGene = 5,
    min.cells = 5,
    min.features = 200
  ) # filters for UMIs, genes, …from ST utility # nolint
  
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

list_objects_ductals <- lapply(c(73:95, 97:159), read_samples)


## WUTANG

wutang <- readRDS("~/Desktop/final_scripts/00_data/Wu_Tang_reannotated_2.rds")
Idents(wutang) <- wutang$celltype_subset

table(wutang@meta.data[wutang@meta.data$celltype_major == "Normal Epithelial", ]$celltype_subset)
table(wutang@meta.data[wutang@meta.data$celltype_major == "Cancer Epithelial", ]$celltype_subset)

# Luminal Progenitors      Mature Luminal        Myoepithelial
# Cancer Basal SC  Cancer Cycling  Cancer Her2 SC  Cancer LumA SC  Cancer LumB SC

wutang <- RenameIdents(wutang, "Luminal Progenitors" = "Epithelial")
wutang <- RenameIdents(wutang, "Mature Luminal" = "Epithelial")
wutang <- RenameIdents(wutang, "Myoepithelial" = "Epithelial")

wutang <- RenameIdents(wutang, "Cancer Basal SC" = "Epithelial")
wutang <- RenameIdents(wutang, "Cancer Cycling" = "Epithelial")
wutang <- RenameIdents(wutang, "Cancer Her2 SC" = "Epithelial")
wutang <- RenameIdents(wutang, "Cancer LumA SC" = "Epithelial")
wutang <- RenameIdents(wutang, "Cancer LumB SC" = "Epithelial")

table(Idents(wutang))

wutang$CARD_wutang <- Idents(wutang)


#  1 - A scRNA-seq gene expression file: - scRNA_gene_exp_filt

wutang@meta.data <- wutang@meta.data %>%
  mutate(subtype = ifelse(is.na(subtype), "adipo", subtype))

wutang <- subset(wutang, subset = subtype != "HER2+" & subtype != "TNBC")

wutang_down <- subset(wutang, downsample = 300)

scRNA_gene_exp <- as.matrix(wutang_down[["RNA"]]$counts)



#  2 - A cell type label file: - cell_type_label

cellID <- rownames(wutang_down@meta.data)
cellType <- as.vector(wutang_down@meta.data[["CARD_wutang"]])
sampleInfo <- wutang_down@meta.data[["Patient"]]

sampleInfo[is.na(sampleInfo)] <- "Tang"

cell_type_label <- cbind(cellID, cellType, sampleInfo)
rownames(cell_type_label) <- cellID
head(cell_type_label)
cell_type_label <- as.data.frame(cell_type_label)


## run CARD

for (i in 1:length(list_objects_ductals)) {
  ST_ptx <- list_objects_ductals[[i]]
  
  # location
  
  coord_loc <- paste0(
    "/mnt/bctl/bkar0016/ductals_st_data/spaceranger/ST",
    unique(list_objects_ductals[[i]]@meta.data[["orig.ident"]]),
    "/outs/spatial/tissue_positions_list.csv"
  )
  coord <- read.csv(coord_loc, header = F)
  rownames(coord) <- coord$V1
  coord_filt <- coord[match(colnames(ST_ptx), rownames(coord)), ]
  coord_x_y <- coord_filt[, c("V3", "V4")]
  coord_pxl <- coord_filt[, c("V5", "V6")]
  colnames(coord_x_y) <- c("x", "y")
  
  colnames(ST_ptx@assays$Spatial@layers$counts) <- rownames(ST_ptx@meta.data)
  rownames(ST_ptx@assays$Spatial@layers$counts) <- row.names(ST_ptx)
  
  
  #### deconvolution ####
  CARD_obj <- createCARDObject(
    sc_count = scRNA_gene_exp,
    sc_meta = cell_type_label,
    spatial_count = ST_ptx@assays$Spatial@layers$counts,
    spatial_location = coord_x_y,
    ct.varname = "cellType",
    ct.select = unique(cell_type_label$cellType),
    sample.varname = "sampleInfo",
    minCountGene = 100,
    minCountSpot = 5
  )
  # ## QC on scRNASeq dataset! ...
  # ## QC on spatially-resolved dataset! ...
  
  
  # length(intersect(row.names(ST_ptx), row.names(erpos))) # Common names
  
  CARD_obj <- CARD_deconvolution(CARD_object = CARD_obj)
  
  file_save <- paste0(
    "/mnt/bctl/myloc/CARD_obj_",
    unique(list_objects_ductals[[i]]@meta.data[["orig.ident"]]), ".RDS"
  )
  saveRDS(CARD_obj, file_save)
  
  
  prop_card <- CARD_obj@Proportion_CARD
  file_save_prop <- paste0(
    "/mnt/bctl/myloc/CARD_prop_",
    unique(list_objects_ductals[[i]]@meta.data[["orig.ident"]]), ".txt"
  )
  write.delim(prop_card, file_save_prop, row.names = T)
}






####################################
## Extended Data Fig. 5 – Correlation Heatmap of Cell-Type Proportions Across ST Spots
####################################

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

mat <- cor(ST_ubermeta[,celltype_cols], use = "pairwise.complete.obs", method = "spearman")
diag(mat) <- NA  # Use NA if you want it to be blank (transparent)

# Define a color scale
col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

# Generate the heatmap
Heatmap(
  mat,
  name = "Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_title = "Correlation Heatmap of Cell-Type Proportions Across ST Spots",
  column_title_gp = gpar(fontsize = 14, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    title_position = "topcenter"
  )
)















