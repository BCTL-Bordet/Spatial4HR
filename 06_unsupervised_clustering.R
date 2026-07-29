# 06_unsupervised_clustering.R

# use the output from 01_preprocessing_import_filter_normalize_ST.R (all_markers)


###############################################################
## Fig 3 - Cluster analysis
###############################################################

####################################
## Fig. 3a - Spatial distribution of cluster identities across the tissue
####################################

library(Seurat)
library(SeuratObject)
library(semla)
library(hdf5r)
library(stringr)
library(data.table)
library(Polychrome)
library(dplyr)
library(magrittr)
library(tibble)
library(patchwork)
library(parallel)

# object_merged <- readRDS("~/Desktop/final_scripts/00_data/object_merged_final.RDS")
ST_ubermeta <- object_merged@meta.data

ST_ubermeta$seurat_clusters <- factor(ST_ubermeta$seurat_clusters + 1L,
                                      levels = 1:24)
## updated cluster numbering 1-24 instead of 0-23



root <- "~/Desktop/ductals_st_data/spaceranger/ST"
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
  # Load sample data paths
  samples <- paste0(root, st_id, "/outs/filtered_feature_bc_matrix.h5")
  spotfiles <- paste0(root, st_id, "/outs/spatial/tissue_positions_list.csv")
  imgs <- paste0(root, st_id, "/outs/spatial/tissue_hires_image.png")
  json <- paste0(root, st_id, "/outs/spatial/scalefactors_json.json")
  info_table <- as.data.frame(cbind(samples, spotfiles, imgs, json))
  
  # Verify that the image file exists
  if (!file.exists(imgs)) {
    stop(paste("Image file does not exist at:", imgs))
  }
  
  # Load the Visium data
  st_sample <- ReadVisiumData(
    info_table,
    assay = "Spatial",
    min.cells = 5,
    min.features = 200
  )
  st_sample@meta.data$orig.ident <- st_id
  
  # Load annotations and coordinates
  annotations <- fread(
    paste0(root, st_id, "/outs/spatial/tissue_positions_list_annotation.csv")
  )
  coordinates <- read.csv(
    paste0(root, st_id, "/outs/spatial/tissue_positions_list.csv"),
    header = FALSE
  )
  
  # Process metadata
  st_sample@meta.data$barcode <- rownames(st_sample@meta.data)
  st_spots <- rownames(st_sample@meta.data)
  
  annotations <- set_colnames(annotations, c(info_colnames, annot_classes)) %>%
    filter(barcode %in% st_spots) %>%
    select(all_of(c(annot_classes, "barcode"))) %>%
    arrange(match(barcode, st_spots))
  
  coordinates <- set_colnames(coordinates, coords_colnames) %>%
    filter(barcode %in% st_spots) %>%
    arrange(match(barcode, st_spots))
  
  st_sample@meta.data <- Reduce(
    function(x, y) merge(x, y, by = "barcode", all = TRUE),
    list(st_sample@meta.data, annotations, coordinates)
  )
  rownames(st_sample@meta.data) <- st_sample@meta.data$barcode
  
  # Filter spots
  st_sample <- SubsetSTData(st_sample, Hole < 0.3)
  st_sample <- SubsetSTData(st_sample, Artefact < 0.3)
  st_sample <- SubsetSTData(st_sample, Out < 0.3)
  
  # Filter genes
  genes <- rownames(st_sample)
  non_meta_genes <- genes[!(grepl("RPL", genes) | grepl("RPS", genes) | grepl("MT-", genes) | grepl("MTRNR", genes))]
  st_sample <- SubsetSTData(st_sample, features = non_meta_genes)
  
  # Load images with error handling
  
  st_sample <- LoadImages(st_sample, verbose = TRUE, time.resolve = FALSE)
  
  if (is.null(st_sample)) {
    stop("st_sample is NULL after attempting to load images.")
  }
  return(st_sample)
}


# Spatial distribution of annotations/genes/clusters.. across the tissue

st_id <- "131" # Replace with your actual sample ID
st_object <- read_samples(st_id)

st_object <- LoadImages(st_object, verbose = FALSE)
cols <- RColorBrewer::brewer.pal(11, "Spectral") |> rev()

MapFeatures(st_object,
            features = c("Tumor", "Cellular_stroma", "Acellular_stroma", "Necrosis"),
            image_use = 'raw',
            colors = cols)

st131_meta <- ST_ubermeta[ST_ubermeta$orig.ident == '131',]

table(st131_meta$barcode == st_object@meta.data$barcode)

st_object@meta.data <- cbind(st_object@meta.data, st131_meta[,c(29,31:69)])


# Plot Seurat clusters with distinct colors
MapLabels(
  object = st_object,
  column_name = "seurat_clusters",
  pt_size = 2, # Adjust point size
  image_use = "raw",
  label_size = 5
)







####################################
## Fig. 3d - Total spots per cluster and # of patients who have ≥5 spots in that cluster,
####################################

library(dplyr)
library(ggplot2)
library(tidyr)


# ---- fixed cluster palette (0–23) ----
cluster_colors <- c(
  "0"  = "#0072B2",  # blue
  "1"  = "#E69F00",  # orange
  "2"  = "#009E73",  # green
  "3"  = "#D55E00",  # vermillion
  "4"  = "#56B4E9",  # sky blue
  "5"  = "#CC79A7",  # reddish purple
  "6"  = "#F0E442",  # yellow
  "7"  = "#000000",  # black
  "8"  = "#9AD0F3",  # light sky
  "9"  = "#0099CC",  # cyan-blue
  "10" = "#7CAE00",  # lime green
  "11" = "#C77CFF",  # violet
  "12" = "#E41A1C",  # red
  "13" = "#4DAF4A",  # green2
  "14" = "#984EA3",  # purple2
  "15" = "#FF7F00",  # orange2
  "16" = "#A65628",  # brown
  "17" = "#377EB8",  # blue2
  "18" = "#F781BF",  # pink
  "19" = "#999999",  # gray
  "20" = "#1B9E77",  # teal
  "21" = "#E6AB02",  # mustard
  "22" = "#A6CEE3",  # pastel blue
  "23" = "#B2DF8A"   # pastel green
)

## 1) Per–patient-per–cluster spot counts
patient_cluster_counts <- ST_ubermeta %>%
  dplyr::count(seurat_clusters, orig.ident, name = "spots")

## 2) Summarise per cluster
summary_df <- patient_cluster_counts %>%
  group_by(seurat_clusters) %>%
  summarise(
    Patients_ge5 = sum(spots >= 5),
    Total_Spots  = sum(spots),
    .groups = "drop"
  ) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = as.character(0:23)))

## 3) Dual-axis bar + points with fixed colors
scale_factor <- max(summary_df$Total_Spots, na.rm = TRUE) /
  max(summary_df$Patients_ge5, na.rm = TRUE)

p <- ggplot(summary_df, aes(x = seurat_clusters)) +
  geom_col(aes(y = Total_Spots, fill = seurat_clusters), width = 0.8) +
  geom_point(aes(y = Patients_ge5 * scale_factor), size = 2, color = "black") +
  geom_text(aes(y = Patients_ge5 * scale_factor, label = Patients_ge5),
            vjust = -0.8, size = 3) +
  scale_fill_manual(values = cluster_colors) +   # <- use your palette here
  scale_y_continuous(
    name = "Total spots",
    sec.axis = sec_axis(~ . / scale_factor, name = "Patients with ≥5 spots")
  ) +
  labs(x = "Cluster",
       title = "Spots per cluster and patients contributing ≥5 spots") +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.position = "none"
  )

p



## rotated 

## 1) Per–patient-per–cluster spot counts
patient_cluster_counts <- ST_ubermeta %>%
  dplyr::count(seurat_clusters, orig.ident, name = "spots")

## 2) Summarise per cluster
summary_df <- patient_cluster_counts %>%
  group_by(seurat_clusters) %>%
  summarise(
    Patients_ge5 = sum(spots >= 5),
    Total_Spots  = sum(spots),
    .groups = "drop"
  ) %>%
  mutate(seurat_clusters = factor(seurat_clusters, levels = as.character(0:23)))

## 3) Scale factor for dual axis
scale_factor <- max(summary_df$Total_Spots, na.rm = TRUE) /
  max(summary_df$Patients_ge5, na.rm = TRUE)

## Reverse cluster order: 0 at top, 23 at bottom
summary_df <- summary_df %>%
  mutate(seurat_clusters = factor(seurat_clusters,
                                  levels = rev(as.character(0:23))))
## 4) Plot rotated
p <- ggplot(summary_df, aes(x = seurat_clusters)) +
  geom_col(aes(y = Total_Spots, fill = seurat_clusters), width = 0.7) +
  geom_point(aes(y = Patients_ge5 * scale_factor), size = 2, color = "black") +
  geom_text(aes(y = Patients_ge5 * scale_factor, label = Patients_ge5),
            hjust = -0.3, size = 4) +
  scale_fill_manual(values = cluster_colors) +
  scale_y_continuous(
    name = "Total spots",
    sec.axis = sec_axis(~ . / scale_factor, name = "Patients with ≥5 spots")
  ) +
  labs(x = "Cluster",
       title = "Spots per cluster and patients contributing ≥5 spots") +
  coord_flip() +  # <-- rotate the whole plot
  theme_minimal(base_size = 12) +
  theme(
    legend.position = "none"
  )

p




####################################
## Fig. 3g - Association between spatial cluster abundance and RFS
####################################

library(dplyr)
library(ggplot2)
library(paletteer)
library(randomcoloR)
library(tidyr)
library(survival)
library(survminer)
library(readxl)

# object_merged <- readRDS("~/Desktop/final_scripts/00_data/object_merged_final.RDS")
ST_ubermeta <- object_merged@meta.data
head(ST_ubermeta)

sum(is.na(ST_ubermeta$seurat_clusters))
sum(is.na(ST_ubermeta$orig.ident))

ductal_meta <- read.csv("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt", sep="")


# Calculate the percentage of each cluster within each 'orig.ident'
percentage_df <- ST_ubermeta %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

print(percentage_df)

# Create a numeric factor for the seurat_clusters to ensure correct ordering in the plot
percentage_df$seurat_clusters <- factor(percentage_df$seurat_clusters, levels = unique(percentage_df$seurat_clusters))

# Convert orig.ident to factor with correct levels
percentage_df$orig.ident <- factor(percentage_df$orig.ident, levels = unique(percentage_df$orig.ident),
                                   labels = paste0("ST", unique(percentage_df$orig.ident)))

# cluster_colors is defined above

# Create bar plot where each patient has a single bar with segmented percentages by cluster
ggplot(percentage_df, aes(x = orig.ident, y = percentage, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Percentage of Clusters by Patient", x = "Patient", y = "Percentage") +
  scale_fill_manual(name = "Seurat Clusters", values = cluster_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

percentage_df$orig.ident <- as.factor(percentage_df$orig.ident)
ductal_meta$name <- as.factor(ductal_meta$name)

percentage_df <- percentage_df[,-3]
percentage_wide <- percentage_df %>%
  spread(key = seurat_clusters, value = percentage, fill = 0) %>%
  rename_with(~paste0("cluster", .), -orig.ident)  # Rename clusters as cluster0, cluster1, ...

# Merge the reshaped percentage data with pseudobulk_meta
merged_data <- ductal_meta %>%
  left_join(percentage_wide, by = c("name" = "orig.ident"))

cluster_columns <- paste0("cluster", 0:23)

merged_data_ordered <- merged_data %>%
  dplyr::select(-starts_with("cluster"), all_of(cluster_columns), everything())

# forest plots
allForest(merged_data_ordered[,c(51:74)], 
          y = Surv(time = merged_data_ordered$time, event = merged_data_ordered$status), fdr = T, new_page = T)




####################################
## Cluster enrichments - Hallmarks
####################################

library(tidyverse)
library(fgsea)
library(tibble)
library(data.table)
library(stringr)
library(ggplot2)

all_markers_res04 <- readRDS("~/Desktop/final_scripts/00_data/all_markers_res04.RDS")

filtered_markers <- all_markers_res04 %>% filter(abs(avg_log2FC) > 0, p_val_adj < 0.05)

cluster_marker_genes <- filtered_markers %>%
  group_by(cluster) %>%
  summarise(marker_genes = list(gene)) %>%
  deframe()

cluster_marker_genes


pathways.hallmark_0 <- gmtPathways("~/Desktop/final_scripts/00_data/h.all.v2023.2.Hs.symbols.gmt.txt") # hallmarks
pathways.hallmark = c(pathways.hallmark_0)

pathways.hallmark %>%
  head() %>%
  lapply(head)

# Initialize empty lists to store results and plots
all_fgseaResTidy <- list()
all_plots <- list()

# Create a loop for clusters 0 through 23
for (i in 0:23) {
  cluster_data <- all_markers_res04[all_markers_res04$cluster == as.character(i), c(7,2)]
  ranks <- deframe(cluster_data)
  
  set.seed(2)
  fgseaRes <- fgsea(pathways = pathways.hallmark, stats = ranks)
  
  fgseaResTidy = fgseaRes[order(fgseaRes$NES), ]
  
  # Store fgseaResTidy in the list
  all_fgseaResTidy[[paste0("cluster_", i)]] <- fgseaResTidy
  
  # Prepare data for plotting
  a = fgseaResTidy[abs(fgseaResTidy$padj) < 0.25, ]
  a$`pval. adjust.` =  with(a, ifelse(padj < 0.05, "< 0.05", "> 0.05"))
  group.colors <- c(`< 0.05` = "#009E73", `> 0.05` = "#56B4E9")
  
  a$pathway = gsub("_", " ", a$pathway)
  
  pathwaysss = a$pathway
  for (j in 1:length(a$pathway)) {
    if (str_length(a$pathway[j]) > 40) {
      pathwaysss[j] = paste0(substr(a$pathway[j], 1, 39), ".")
    }
  }
  
  a$pathwaysss = pathwaysss
  
  # Create ggplot for the current cluster
  p <- ggplot(a, aes(reorder(pathwaysss, NES), NES)) +
    geom_col(aes(fill = `pval. adjust.`)) + 
    scale_fill_manual(values = group.colors) +
    coord_flip() +
    labs(x = "Pathway", y = "Normalized Enrichment Score",
         title = paste("Hallmark Pathways (GSEA) for cluster", i)) +
    theme_minimal()
  
  # Store the plot in the list
  all_plots[[paste0("cluster_", i)]] <- p
  
  # Show in a nice table:
  fgseaResTidy %>%
    dplyr::select(-leadingEdge, -ES) %>%
    arrange(padj) %>%
    DT::datatable()
}

# Save the list of fgseaResTidy and plots as RDS files
# saveRDS(all_fgseaResTidy, "/Users/bengisukarakose/Desktop/final_scripts/00_data/all_fgseaResTidy_clusters.RDS")
# saveRDS(all_plots, "/Users/bengisukarakose/Desktop/final_scripts/00_data/all_plots_clusters.RDS")


# all_fgseaResTidy <- readRDS("/Users/bengisukarakose/Desktop/final_scripts/00_data/all_fgseaResTidy_clusters.RDS")


####################################
## Fig. 3f - Hallmark gene set enrichment analysis of spatial clusters
####################################

library(dplyr)
library(ggplot2)
library(stringr)
library(data.table)

## 1) Bind all clusters into one data frame ---------------------------------
# Convert the named list -> long data frame and tag cluster id
fgsea_long <- rbindlist(
  lapply(names(all_fgseaResTidy), function(nm) {
    df <- as.data.frame(all_fgseaResTidy[[nm]])
    df$cluster <- gsub("^cluster_", "", nm)
    df
  }),
  fill = TRUE
)

## 2) Clean & order factors --------------------------------------------------
fgsea_long <- fgsea_long %>%
  mutate(
    cluster = factor(cluster, levels = as.character(0:23)),
    pathway_clean = pathway %>%
      str_replace("^HALLMARK_", "") %>%
      str_replace_all("_", " "),
    # (optional) truncate very long pathway labels for plotting only
    pathway_plot = ifelse(nchar(pathway_clean) > 45,
                          paste0(substr(pathway_clean, 1, 44), "…"),
                          pathway_clean),
    neglog10_fdr = -log10(padj),
    sig = padj < 0.05
  )  %>%
  filter(pathway_clean != "PANCREAS BETA CELLS")   # remove this hallmark

# Order pathways on y-axis by overall significance across clusters
pathway_order <- fgsea_long %>%
  group_by(pathway_plot) %>%
  summarise(best_fdr = min(padj, na.rm = TRUE), .groups = "drop") %>%
  arrange(best_fdr) %>%
  pull(pathway_plot)

fgsea_long <- fgsea_long %>%
  mutate(pathway_plot = factor(pathway_plot, levels = rev(pathway_order)))

## 3A) Variant A: show Top N pathways per cluster (by FDR) -------------------
topN <- 10   # change to what you prefer
fgsea_topN <- fgsea_long %>%
  filter(padj < 0.05) %>%               # keep somewhat liberal threshold to rank within cluster
  group_by(cluster) %>%
  slice_min(order_by = padj, n = topN, with_ties = FALSE) %>%
  ungroup()

p_topN <- ggplot(fgsea_topN, aes(x = cluster, y = pathway_plot)) +
  geom_point(aes(size = neglog10_fdr, color = NES), alpha = 0.85) +
  scale_size_continuous(name = expression(-log[10]("FDR")), range = c(2, 8)) +
  scale_color_gradient2(name = "NES", low = "#2c7bb6", mid = "grey90",
                        high = "#d7191c", midpoint = 0) +
  labs(
    x = "Cluster",
    y = "Hallmark pathway",
    # title = paste0("Hallmark enrichment by cluster (Top ", topN, " per cluster)"),
    title = "Hallmark enrichment by cluster",
    subtitle = "Dot size: −log10(FDR). Color: NES."
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    panel.grid.minor = element_blank()
  )

p_topN



## cluster stats

# object_merged <- readRDS("~/Desktop/final_scripts/00_data/object_merged_final.RDS")
ST_ubermeta <- object_merged@meta.data

library(dplyr)
library(ggplot2)
library(paletteer)
library(randomcoloR)

# Calculate the percentage of each cluster within each 'orig.ident'
percentage_df <- ST_ubermeta %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# View the calculated percentages to ensure correctness
print(percentage_df)

core_clusters <- c("0", "1", "2", "3")

# Calculate sum of core cluster percentages per patient
core_cluster_summary <- percentage_df %>%
  filter(seurat_clusters %in% core_clusters) %>%
  group_by(orig.ident) %>%
  summarise(core_cluster_percentage = sum(percentage)) %>%
  arrange(core_cluster_percentage)

# View the range
range(core_cluster_summary$core_cluster_percentage)


# Count number of unique patients and total spots per cluster
cluster_distribution_summary <- ST_ubermeta %>%
  group_by(seurat_clusters) %>%
  summarise(
    n_spots = n(),
    n_patients = n_distinct(orig.ident)
  ) %>%
  arrange(n_patients)

print(cluster_distribution_summary, n=100)

#########

library(Seurat)
library(SeuratObject)
library(semla)

DimPlot(object_merged, group.by = "seurat_clusters", label = TRUE)


####################################
## Fig. 3e - Histological Composition of ST Clusters
####################################

cluster_annotation <- c()

for (i in 0:(length(unique(ST_ubermeta$seurat_clusters))-1)) {
  ST_ubermeta_cluster <- ST_ubermeta[ST_ubermeta$seurat_clusters == i,5:21]
  annotation <- colMeans(ST_ubermeta_cluster)
  cluster_annotation <- rbind(cluster_annotation,annotation)}

rownames(cluster_annotation) <- paste0('cluster', 0:(length(unique(ST_ubermeta$seurat_clusters))-1))
cluster_annotation <- as.data.frame(cluster_annotation)

# Load necessary libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(tibble)

# Convert row names to a column for easy plotting
cluster_annotation <- cluster_annotation %>%
  rownames_to_column(var = "cluster")

# Reshape data from wide to long format for ggplot
cluster_annotation_long <- cluster_annotation %>%
  pivot_longer(cols = -cluster, names_to = "feature", values_to = "percentage")


col_vector = c("#017801", "Black", "Navy", "#e9e900", "#ff9980", "#e8d1bb", "#dc0000", "#6e2400", "#9980e6", "#80801a" , "#ccffcc", "#4d8080", "#c4417f", "#40e5f6", "#41344c", "Grey", "Magenta")
binded <- as.data.frame(cbind(colnames(ST_ubermeta_cluster),col_vector))
binded <- binded[order(binded$V1), ]
binded$col_vector


# # Create palette of colors
col_vector = binded$col_vector

cluster_annotation_long$cluster <- factor(cluster_annotation_long$cluster, 
                                          levels = paste0("cluster", 0:23))


# Plot-specific palette: names lock each color to its feature
plot_cols <- setNames(binded$col_vector, binded$V1)

# Remove unwanted features from the palette
plot_cols <- plot_cols[
  !names(plot_cols) %in% c("Nodule_lymphoid", "Apocrine.metaplasia")
]

# Rename the corresponding palette entry
names(plot_cols)[names(plot_cols) == "Canal_galactophore"] <- "Normal breast"

ggplot(
  cluster_annotation_long %>%
    filter(!feature %in% c("Nodule_lymphoid", "Apocrine.metaplasia")) %>%
    mutate(
      feature = recode(
        feature,
        "Canal_galactophore" = "Normal breast"
      )
    ),
  aes(x = cluster, y = percentage, fill = feature)
) +
  geom_col() +
  labs(
    title = "Histological Composition of ST Clusters",
    x = "Clusters",
    y = "Percentage"
  ) +
  scale_fill_manual(
    values = plot_cols,
    labels = function(x) gsub("_", " ", x)
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      vjust = 0.5,
      hjust = 1
    )
  )



####################################
## Fig. 3c - Percentage of Clusters by Patient
####################################

library(dplyr)
library(paletteer)
library(randomcoloR)
library(Seurat)
library(SeuratObject)
library(semla)
library(ggplot2)

str(ST_ubermeta)
ST_ubermeta$orig.ident <- as.factor(ST_ubermeta$orig.ident )

# Calculate the percentage of each cluster within each 'orig.ident'
percentage_df <- ST_ubermeta %>%
  group_by(orig.ident, seurat_clusters) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(orig.ident) %>%
  mutate(percentage = count / sum(count) * 100) %>%
  ungroup()

# View the calculated percentages to ensure correctness
print(percentage_df)
percentage_df$seurat_clusters <- as.factor(percentage_df$seurat_clusters)

DimPlot(object_merged, group.by = "seurat_clusters", label = TRUE)




# Convert seurat_clusters to factor with consistent levels
ordered_clusters <- sort(as.numeric(levels(percentage_df$seurat_clusters)))

percentage_df$seurat_clusters <- factor(
  percentage_df$seurat_clusters,
  levels = as.character(ordered_clusters)  # ensure character match
)

# Plot
ggplot(percentage_df, aes(x = orig.ident, y = percentage, fill = seurat_clusters)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(title = "Percentage of Clusters by Patient", x = "Patient", y = "Percentage") +
  scale_fill_manual(name = "Seurat Clusters", values = cluster_colors) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))





####################################
## Fig. 3h - METABRIC Validation
####################################

library(GSVA)
library(dplyr)
library(survival)
library(survminer)


filtered_markers_pos <- all_markers_res04[all_markers_res04$avg_log2FC > 1, ]
positive_cluster_genes <- split(filtered_markers_pos$gene, filtered_markers_pos$cluster)
names(positive_cluster_genes) <- paste0('cluster_', unique(filtered_markers_pos$cluster))

counts_metabric <- read.delim("~/Desktop/bc_cox_model/2_metabric/counts_metabric.gct")
counts_metabric <- counts_metabric[,-2]
rownames(counts_metabric) <- counts_metabric$NAME
counts_metabric <- counts_metabric[,-1]

data_clinical_patient <- read.delim("~/Desktop/bc_cox_model/2_metabric/data_clinical_patient.txt", comment.char="#")
data_clinical_sample <- read.delim("~/Desktop/bc_cox_model/2_metabric/data_clinical_sample.txt", comment.char="#")

table(data_clinical_patient$PATIENT_ID == data_clinical_sample$PATIENT_ID)

metabric_meta <- cbind(data_clinical_patient, data_clinical_sample)
colnames(counts_metabric) = sub("\\.", "-", colnames(counts_metabric))
metabric_meta = metabric_meta[metabric_meta$PATIENT_ID %in% colnames(counts_metabric),]

genesets_sig <- positive_cluster_genes[c(7,12,13,18)]



# Create the GSVAParams object for the "gsva" method
params <- gsvaParam(as.matrix(counts_metabric),
                    kcdf="Gaussian", geneSets = genesets_sig)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)

table(rownames(module_sig) == metabric_meta$PATIENT_ID)

metabric_meta <- cbind(metabric_meta, module_sig)


metabric_meta$OS_STATUS = as.numeric(gsub(":.*", "", metabric_meta$OS_STATUS))
metabric_meta$RFS_STATUS = as.numeric(gsub(":.*", "", metabric_meta$RFS_STATUS))

# forest plots
allForest(metabric_meta[,c(38:41)], 
          y = Surv(time = metabric_meta$RFS_MONTHS, event = metabric_meta$RFS_STATUS), fdr = T, new_page = T)

allForest(metabric_meta[,c(38:41)], 
          y = Surv(time = metabric_meta$OS_MONTHS, event = metabric_meta$OS_STATUS), fdr = T, new_page = T)


# figure
colnames(metabric_meta)[38:41] <- gsub('_', '', colnames(metabric_meta)[38:41])

allForest(metabric_meta[,c(38:41)], 
          y = Surv(time = metabric_meta$RFS_MONTHS, event = metabric_meta$RFS_STATUS), fdr = T, new_page = T)




####################################
## Fig. 3i - SCAN-B Validation
####################################

filtered_markers_pos <- all_markers_res04[all_markers_res04$avg_log2FC > 1, ]
positive_cluster_genes <- split(filtered_markers_pos$gene, filtered_markers_pos$cluster)
names(positive_cluster_genes) <- paste0('cluster_', unique(filtered_markers_pos$cluster))

counts_scanb <- read.delim("~/Desktop/bc_cox_model/3_scanb/counts_scanb.gct")
counts_scanb <- counts_scanb[,-2]
rownames(counts_scanb) <- counts_scanb$NAME
counts_scanb <- counts_scanb[,-1]

clin_rev_scanb <- read.delim("~/Desktop/bc_cox_model/3_scanb/clin_rev_scanb.txt")
clin_all_scanb <- read_excel("~/Desktop/bc_cox_model/3_scanb/Supplementary Data Table 1 - 2023-01-13.xlsx")

# Filter clin_all_scanb to keep only rows with GEX.assay in colnames(counts_scanb)
clin_filtered <- clin_all_scanb %>% 
  filter(GEX.assay %in% colnames(counts_scanb))

genesets_sig <- positive_cluster_genes[c(7,12,13,18)]

# Create the GSVAParams object for the "gsva" method
params <- gsvaParam(as.matrix(counts_scanb),  
                    kcdf="Gaussian", geneSets = genesets_sig)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)

table(rownames(module_sig) == clin_rev_scanb$name)
scanb_meta <- cbind(clin_rev_scanb, module_sig)

table(rownames(module_sig) == clin_filtered$GEX.assay)
scanb_big_meta <- cbind(clin_filtered, module_sig)

scanb_meta$status = as.numeric(scanb_meta$status)
scanb_meta$time = as.numeric(scanb_meta$time)

scanb_big_meta$DRFi_days = as.numeric(scanb_big_meta$DRFi_days)
scanb_big_meta$DRFi_event = as.numeric(scanb_big_meta$DRFi_event)

scanb_big_meta$OS_days = as.numeric(scanb_big_meta$OS_days)
scanb_big_meta$OS_event = as.numeric(scanb_big_meta$OS_event)

scanb_big_meta$RFi_days = as.numeric(scanb_big_meta$RFi_days)
scanb_big_meta$RFi_event = as.numeric(scanb_big_meta$RFi_event)

scanb_big_meta$BCFi_days = as.numeric(scanb_big_meta$BCFi_days)
scanb_big_meta$BCFi_event = as.numeric(scanb_big_meta$BCFi_event)



# forest plots

allForest(scanb_big_meta[,c(88:91)],
          y = Surv(time = scanb_big_meta$DRFi_days, event = scanb_big_meta$DRFi_event), fdr = T, new_page = T)

allForest(scanb_big_meta[,c(88:91)],
          y = Surv(time = scanb_big_meta$OS_days, event = scanb_big_meta$OS_event), fdr = T, new_page = T)

allForest(scanb_big_meta[,c(88:91)],
          y = Surv(time = scanb_big_meta$RFi_days, event = scanb_big_meta$RFi_event), fdr = T, new_page = T)

allForest(scanb_big_meta[,c(88:91)],
          y = Surv(time = scanb_big_meta$BCFi_days, event = scanb_big_meta$BCFi_event), fdr = T, new_page = T)


## fig 
colnames(scanb_big_meta)[88:91] <- gsub('_', '', colnames(scanb_big_meta)[88:91])

allForest(scanb_big_meta[,c(88:91)],
          y = Surv(time = scanb_big_meta$RFi_days, event = scanb_big_meta$RFi_event), fdr = T, new_page = T)



####################################
## Extended Data Fig.7 - Cell-type associations with ST clusters (bubbleplot)
####################################

# cluster column
ST_ubermeta$cluster <- as.factor(ST_ubermeta$seurat_clusters)
subtype <- levels(ST_ubermeta$cluster)

# deconvolution columns
signatures <- colnames(ST_ubermeta)[31:69]

# scale (same spirit as genefu::rescale)
ST_ubermeta[, signatures] <- apply(
  ST_ubermeta[, signatures],
  2,
  function(x) scale(x)[,1]
)

effect <- matrix(
  0,
  nrow = length(signatures),
  ncol = length(subtype),
  dimnames = list(signatures, subtype)
)

pvalue <- effect
ci <- effect

for (j in seq_along(signatures)) {
  
  cat("Signature", j, "/", length(signatures), "\n")
  
  for (i in seq_along(subtype)) {
    
    id <- subtype[i]
    ST_ubermeta$tmp <- ifelse(ST_ubermeta$cluster == id, 1, 0)
    
    f <- as.formula(paste("tmp ~", signatures[j]))
    res.logist <- glm(f, data = ST_ubermeta, family = binomial)
    
    p <- wilcox.test(
      ST_ubermeta[[signatures[j]]] ~ ST_ubermeta$tmp
    )$p.value
    
    effect[j, i] <- round(exp(coef(res.logist))[2], 2)
    pvalue[j, i] <- p
  }
}

padjust <- matrix(
  p.adjust(pvalue, method = "fdr"),
  nrow = nrow(pvalue),
  dimnames = dimnames(pvalue)
)

# filtering 
effect[padjust > 0.25] <- NA
effect[effect <= .25] <- .25
effect[effect >= 4] <- 4

# transpose for plotting
e_ductal <- t(effect)

# clean labels (same as old script)
rownames(e_ductal) <- gsub("_", " ", rownames(e_ductal))
colnames(e_ductal) <- gsub("\\.", " ", colnames(e_ductal))

rownames(e_ductal) <- paste0("cluster ", rownames(e_ductal))

plot_legend_1 <- matrix(c("0.25","4"), nrow = 1)
class(plot_legend_1) <- "numeric"
colnames(plot_legend_1) <- c("depletion","enrichment")
rownames(plot_legend_1) <- "Direction"

plot_legend_2 <- matrix(c("0.4","0.2"), nrow = 1)
class(plot_legend_2) <- "numeric"
colnames(plot_legend_2) <- c("large","small")
rownames(plot_legend_2) <- "Effect size"

library(Cairo)
par(mar=c(4.2, 1, 6, 1))

circle <- function(x,y,r,nsteps=100,...){  
  rs <- seq(0,2*pi,len=180)
  xc <- x + r * cos(rs+pi/2)
  yc <- y + r * sin(rs+pi/2)
  polygon(xc,yc,...)
}

plot(
  row(e_ductal)*1,
  col(e_ductal)*1,
  type="n",
  xlab="", ylab="",
  xlim=c(0.5,nrow(e_ductal)*1+0.5),
  ylim=c(0.5,ncol(e_ductal)*1+1),
  axes=FALSE,
  ann=FALSE,
  asp=1
)

for(i in 1:nrow(e_ductal)){
  for(j in 1:ncol(e_ductal)){
    circle(
      i*1,
      j*1,
      abs(log2(e_ductal[i,j]))*0.2,
      col=adjustcolor(
        ifelse(e_ductal[i,j]>1,"#009E73","#D55E00"),
        alpha.f = 0.8
      ),
      border = NA
    )
  }
}

axis(1, at=1:nrow(e_ductal)*1, labels=FALSE, pos=0)
text(
  x=rep(.2,ncol(e_ductal)),
  y=1:ncol(e_ductal),
  labels=colnames(e_ductal),
  cex.axis=1.5,
  las=2,
  pos=2,
  xpd=TRUE
)

axis(2, at=1:ncol(e_ductal)*1, pos=0.5, labels=FALSE)
text(
  x=1:nrow(e_ductal),
  # y=rep(-1,nrow(e_ductal)),
  y=rep(-1.6, nrow(e_ductal)),
  labels=rownames(e_ductal),
  cex.axis=1.5,
  srt=90,
  adj=1,
  xpd=TRUE
)



