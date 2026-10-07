####################################
## 03_tumor_growth_patterns.R
####################################

## Tumor Growth Patterns

## 5 bins are defined based on percentages of tumor annotations of spots
## (0%", "0-25%", "25-50%", "50-75%", "75-100%")

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyr)
library(Seurat)
library(SeuratObject)
library(semla)
library(hdf5r)
library(stringr)
library(data.table)
library(Polychrome)
library(magrittr)
library(tibble)
library(patchwork)
library(parallel)


####################################
## Fig 2a - Tumor content binning
####################################

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

st_id <- "159"
st_object <- read_samples(st_id)

st_object <- LoadImages(st_object, verbose = FALSE)
cols <- c("0" = "white", "1" = "red")

ST_ubermeta <- st_object@meta.data

ST_ubermeta$tumor_group <- cut(
  ST_ubermeta$Tumor,
  breaks = c(-Inf, 0, 0.25, 0.50, 0.75, 1.00),
  labels = c("0%", "0-25%", "25-50%", "50-75%", "75-100%"),
  right = TRUE,  # Include the right endpoint in intervals
  include.lowest = TRUE
)

ST_ubermeta$tumor_group[ST_ubermeta$Tumor == 0] <- "0%"
head(ST_ubermeta)

# Add new columns based on tumor_group
ST_ubermeta <- ST_ubermeta %>%
  mutate(
    q0 = ifelse(tumor_group == "0%", 1, 0),
    q1 = ifelse(tumor_group == "0-25%", 1, 0),
    q2 = ifelse(tumor_group == "25-50%", 1, 0),
    q3 = ifelse(tumor_group == "50-75%", 1, 0),
    q4 = ifelse(tumor_group == "75-100%", 1, 0)
  )


st_object@meta.data <- ST_ubermeta


MapFeatures(st_object,
            features = c("q0", 'q1', 'q2', 'q3', 'q4'),
            image_use = "raw",
            colors = cols)



# ---- Show and save each bin map ----
feat_list <- c("q0", "q1", "q2", "q3", "q4")
labels    <- c("0", "0-25", "25-50", "50-75", "75-100")

plots <- list()

for (i in seq_along(feat_list)) {
  f   <- feat_list[i]
  lab <- labels[i]
  
  p <- MapFeatures(
    st_object,
    features  = f,
    image_use = "raw",
    colors    = cols,
    pt_size = 1.3
  ) +
    ggtitle(lab) +
    theme(
      plot.title = element_text(size = 16, face = "bold"),
      legend.position = "none" # hides redundant legend
    )
  
  # Store in list so you can view them together if needed
  plots[[lab]] <- p
  
  # Show in RStudio / plotting window
  print(p)
  
  # Save
  # ggsave(paste0("Fig2A_", lab, ".pdf"), p,
  #        width = 4, height = 5, units = "in", useDingbats = FALSE, bg = "white")
  # ggsave(paste0("Fig2A_", lab, ".png"), p,
  #        width = 4, height = 5, units = "in", dpi = 600, bg = "white")
}


####################################
## Extended Data Fig.6 - Visualizing the amount of spots in a dot plot
####################################

# object_merged <- readRDS("~/Desktop/final_scripts/00_data/object_merged_final.RDS")
ST_ubermeta <- object_merged@meta.data
ST_ubermeta$id <- paste(ST_ubermeta$barcode, ST_ubermeta$orig.ident, sep = "_")
colnames(ST_ubermeta)

# create tumor_group based on tumor annotations
ST_ubermeta$tumor_group <- cut(
  ST_ubermeta$Tumor,
  breaks = c(-Inf, 0, 0.25, 0.50, 0.75, 1.00),
  labels = c("0%", "0-25%", "25-50%", "50-75%", "75-100%"),
  right = TRUE,  # Include the right endpoint in intervals
  include.lowest = TRUE
)

ST_ubermeta$tumor_group[ST_ubermeta$Tumor == 0] <- "0%"


spot_count_per_group <- as.data.frame(table(ST_ubermeta[,c(2,71)]))

spot_count_per_group <- spot_count_per_group %>%
  group_by(orig.ident) %>%  # Group by orig.ident
  mutate(Freq_percentage = Freq / sum(Freq) * 100 )  # Divide by total Freq for that orig.ident

spot_count_per_group$orig.ident <- paste0('ST', spot_count_per_group$orig.ident)

# Create the dot plot ( with freq percentages )
ggplot(spot_count_per_group, aes(x = tumor_group, y = factor(orig.ident), size = Freq_percentage)) +
  geom_point(alpha = 0.7, color = "grey") +  # Use semi-transparent dots
  scale_size_continuous(name = "Spot Count", range = c(1,4)) +  # Adjust the size range
  labs(
    title = "Dot Plot of Spot Counts by Tumor Group",
    x = "Tumor Group",
    y = "Patient ID"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "right"
  )




############################################################
## Tumor growth pattern analysis
##
## This script:
## 1. Classifies Visium spots according to tumor content
## 2. Computes per-sample tumor composition profiles
## 3. Identifies recurrent tumor growth patterns
##    using hierarchical clustering
############################################################

library(dplyr)
library(tidyr)
library(reshape2)
library(ggplot2)
library(NbClust)
library(dendextend)

## Build tumor composition matrix
spot_table <- table(
  ST_ubermeta$orig.ident,
  ST_ubermeta$tumor_group
)

spot_table <- as.data.frame(spot_table)

spot_matrix <- dcast(
  spot_table,
  Var1 ~ Var2,
  value.var = "Freq"
)

rownames(spot_matrix) <- paste0(
  "ST",
  spot_matrix$Var1
)

spot_matrix <- spot_matrix[,-1]

spot_matrix <- sweep(
  as.matrix(spot_matrix),
  1,
  rowSums(spot_matrix),
  "/"
) * 100



## Hierarchical clustering
distance_matrix <- dist(
  spot_matrix,
  method = "euclidean")

hc <- hclust(
  distance_matrix,
  method = "ward.D2")


## Determine optimal number of clusters
NbClust(
  spot_matrix,
  distance = "euclidean",
  method = "ward.D2",
  min.nc = 2,
  max.nc = 10)

## Majority rule identified 3 clusters
cluster_assignments <- cutree(hc, k = 3)

## Assign descriptive names
cluster_names <- c(
  "1" = "Nests",
  "2" = "Scattered",
  "3" = "Cell-dense")

tumor_patterns <- data.frame(
  Sample = rownames(spot_matrix),
  Cluster = cluster_assignments,
  TumorPattern = cluster_names[
    as.character(cluster_assignments)])

## Visualize clustering
plot(as.dendrogram(hc) %>%
       color_branches(k = 3),
     main = "Tumor pattern clustering")

## Save results


####################################
## Fig 2b - Representative annotated tissue sections illustrating the three tumor growth patterns
####################################


####################################
## Fig 2c - Histology - tumor patterns
####################################

library(ComplexHeatmap)
library(viridis)
library(grid)
library(circlize)
library(textshape)

ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")
annotations <- colnames(ductal_meta)[grepl("_annotation$", colnames(ductal_meta))]

### Heatmap with the significance brackets

####################################
## Statistical significance
##
## For each histopathological feature:
## - Kruskal-Wallis test across the
##   three tumor architectural patterns
## - Benjamini-Hochberg correction
## - FDR < 0.05 considered significant
##
## For each significant feature,
## the tumor pattern with the most
## extreme absolute z-score is outlined.
####################################

# --- aggregate by tumor pattern ---

heatmap_data <- ductal_meta %>%
  group_by(tumor_pattern) %>%
  summarise(
    across(
      all_of(annotations),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  column_to_rownames("tumor_pattern")


row_order <- c(
  "Cell-dense",
  "Scattered",
  "Nests"
)


heatmap_matrix <- scale(heatmap_data)

heatmap_matrix <- heatmap_matrix[
  row_order,
  ,
  drop = FALSE
]

histology_stats <- lapply(
  annotations,
  function(feature_i) {
    
    dat <- ductal_meta %>%
      select(
        tumor_pattern,
        all_of(feature_i)
      )
    
    colnames(dat)[2] <- "value"
    
    dat <- dat %>%
      filter(
        !is.na(tumor_pattern),
        is.finite(value)
      )
    
    kw <- kruskal.test(
      value ~ tumor_pattern,
      data = dat
    )
    
    data.frame(
      feature = feature_i,
      p_value = kw$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  )


# significant features

significant_features <- histology_stats %>%
  filter(FDR < 0.05)


# identify the tumor pattern with the
# most extreme z-score for each significant feature

box_results <- lapply(
  significant_features$feature,
  function(feature_i) {
    
    z_values <- heatmap_matrix[
      row_order,
      feature_i
    ]
    
    selected_pattern <- row_order[
      which.max(abs(z_values))
    ]
    
    data.frame(
      feature = feature_i,
      tumor_pattern = selected_pattern,
      z_score = z_values[selected_pattern]
    )
  }
) %>%
  bind_rows()


# build TRUE/FALSE matrix for black outlines

sig_matrix <- matrix(
  FALSE,
  nrow = nrow(heatmap_matrix),
  ncol = ncol(heatmap_matrix),
  dimnames = dimnames(heatmap_matrix)
)

for (k in seq_len(nrow(box_results))) {
  
  sig_matrix[
    box_results$tumor_pattern[k],
    box_results$feature[k]
  ] <- TRUE
}


####################################
## Clean names
####################################

clean_names <- colnames(heatmap_matrix)

clean_names <- gsub(
  "_annotation$",
  "",
  clean_names
)

clean_names <- gsub(
  "[._]",
  " ",
  clean_names
)

# Canal galactophore corresponds to
# the normal breast compartment

clean_names[
  clean_names == "Canal galactophore"
] <- "Normal breast"

colnames(heatmap_matrix) <- clean_names
colnames(sig_matrix) <- clean_names


####################################
## Colors
####################################

rng <- max(
  abs(
    range(
      heatmap_matrix,
      na.rm = TRUE
    )
  )
)

brks <- seq(
  -rng,
  rng,
  length.out = 101
)

color_fun <- colorRamp2(
  c(-1, 0, 1),
  c("#4575b4", "white", "#d73027")
)

red_white_blue <- color_fun(
  seq(-1, 1, length.out = 100)
)



####################################
## Final heatmap + significance boxes
####################################

ComplexHeatmap::Heatmap(
  heatmap_matrix,
  
  name = "Z score",
  
  col = color_fun,
  
  ####################################
  ## Keep tumor-pattern order fixed
  ####################################
  
  cluster_rows = FALSE,
  
  ####################################
  ## Cluster histological features
  ## using same settings as original
  ## pheatmap
  ####################################
  
  cluster_columns = TRUE,
  
  clustering_distance_columns = "euclidean",
  
  clustering_method_columns = "ward.D2",
  
  ####################################
  ## Labels
  ####################################
  
  row_names_side = "right",
  
  column_names_rot = 90,
  
  ####################################
  ## No default cell borders
  ####################################
  
  rect_gp = gpar(
    col = NA
  ),
  
  ####################################
  ## Add black significance boxes
  ####################################
  
  cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    if (
      isTRUE(
        sig_matrix[i, j]
      )
    ) {
      
      grid.rect(
        x = x,
        y = y,
        width = width,
        height = height,
        gp = gpar(
          fill = NA,
          col = "black",
          lwd = 1.5
        )
      )
    }
  }
)




####################################
## Fig 2d - Deconvolution - tumor patterns
####################################

library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyr)

CARD_subset <- colnames(ST_ubermeta)[31:69]
columns_to_plot <- CARD_subset

col_means_by_group <- ST_ubermeta %>%
  group_by(orig.ident) %>%
  summarise(across(all_of(CARD_subset), mean, na.rm = TRUE), .groups = "drop")

head(col_means_by_group)

col_means_by_group$name <- paste0('ST', col_means_by_group$orig.ident)

tumor_patterns <- as.data.frame(cbind(ductal_meta$name, ductal_meta$tumor_pattern))
colnames(tumor_patterns) <- c('name', 'tumor_pattern')

CARD_data <- left_join(col_means_by_group, tumor_patterns, by='name')

library(dplyr)
library(ComplexHeatmap)
library(viridis)
library(grid)
library(tibble)
library(circlize)

### Heatmap with the significance brackets
####################################
## Aggregate by tumor pattern
## for heatmap
####################################

agg_df <- CARD_data %>%
  group_by(tumor_pattern) %>%
  summarise(
    across(
      all_of(columns_to_plot),
      mean,
      na.rm = TRUE
    ),
    .groups = "drop"
  )


row_order <- c(
  "Cell-dense",
  "Scattered",
  "Nests"
)


####################################
## Build heatmap matrix
####################################

heatmap_matrix <- agg_df %>%
  column_to_rownames(
    "tumor_pattern"
  ) %>%
  as.matrix()


heatmap_matrix <- heatmap_matrix[
  row_order,
  ,
  drop = FALSE
]


####################################
## Column-wise z-score
####################################

heatmap_matrix <- scale(
  heatmap_matrix
)


####################################
## Global Kruskal-Wallis tests
## one test per CARD cell type
####################################

CARD_global_stats <- lapply(
  columns_to_plot,
  function(feature_i) {
    
    dat <- CARD_data %>%
      select(
        tumor_pattern,
        all_of(feature_i)
      )
    
    colnames(dat)[2] <- "value"
    
    dat <- dat %>%
      filter(
        !is.na(tumor_pattern),
        is.finite(value)
      )
    
    kw <- kruskal.test(
      value ~ tumor_pattern,
      data = dat
    )
    
    data.frame(
      feature = feature_i,
      p_value = kw$p.value
    )
  }
) %>%
  bind_rows() %>%
  mutate(
    FDR = p.adjust(
      p_value,
      method = "BH"
    )
  )


####################################
## Print all statistics
####################################

CARD_global_stats %>%
  arrange(FDR) %>%
  print


####################################
## Significant features
####################################

significant_features <- CARD_global_stats %>%
  filter(
    FDR < 0.05
  )


cat(
  "\n\nSIGNIFICANT CARD FEATURES FDR < 0.05:\n"
)

significant_features %>%
  arrange(FDR) %>%
  print


####################################
## For each significant feature,
## find the tumor pattern with
## the largest absolute z-score
####################################

box_results <- lapply(
  significant_features$feature,
  function(feature_i) {
    
    z_values <- heatmap_matrix[
      row_order,
      feature_i
    ]
    
    selected_pattern <- row_order[
      which.max(
        abs(z_values)
      )
    ]
    
    data.frame(
      feature = feature_i,
      tumor_pattern = selected_pattern,
      z_score = z_values[
        selected_pattern
      ]
    )
  }
) %>%
  bind_rows()


####################################
## Clean feature names
####################################

box_results <- box_results %>%
  mutate(
    
    feature_clean = gsub(
      "[._]",
      " ",
      feature
    ),
    
    feature_clean = gsub(
      "  ",
      " ",
      feature_clean
    )
  )


####################################
## Print predicted boxes
####################################

cat(
  "\n\nPREDICTED BOXES:\n"
)

box_results %>%
  select(
    tumor_pattern,
    feature_clean,
    z_score
  ) %>%
  arrange(
    tumor_pattern,
    feature_clean
  ) %>%
  print


####################################
## Build significance matrix
####################################

sig_matrix <- matrix(
  FALSE,
  nrow = nrow(heatmap_matrix),
  ncol = ncol(heatmap_matrix),
  dimnames = dimnames(heatmap_matrix)
)


for (
  k in seq_len(
    nrow(box_results)
  )
) {
  
  sig_matrix[
    box_results$tumor_pattern[k],
    box_results$feature[k]
  ] <- TRUE
}


####################################
## Clean heatmap names
####################################

clean_names <- gsub(
  "[._]",
  " ",
  columns_to_plot
)

clean_names <- gsub(
  "  ",
  " ",
  clean_names
)

colnames(
  heatmap_matrix
) <- clean_names

colnames(
  sig_matrix
) <- clean_names


####################################
## Heatmap colors
####################################

rng <- max(
  abs(
    range(
      heatmap_matrix,
      na.rm = TRUE
    )
  )
)

brks <- seq(
  -rng,
  rng,
  length.out = 101
)

color_fun <- colorRamp2(
  c(
    -1,
    0,
    1
  ),
  c(
    "#4575b4",
    "white",
    "#d73027"
  )
)

red_white_blue <- color_fun(
  seq(
    -1,
    1,
    length.out = 100
  )
)


####################################
## Final heatmap + significance boxes
####################################

ComplexHeatmap::Heatmap(
  heatmap_matrix,
  
  name = "Z score",
  
  col = color_fun,
  
  ####################################
  ## Keep tumor-pattern order fixed
  ####################################
  
  cluster_rows = FALSE,
  
  ####################################
  ## Cluster CARD cell types
  ## using same settings as original
  ## pheatmap defaults
  ####################################
  
  cluster_columns = TRUE,
  
  clustering_distance_columns = "euclidean",
  
  clustering_method_columns = "complete",
  
  ####################################
  ## Labels
  ####################################
  
  row_names_side = "right",
  
  column_names_rot = 90,
  
  row_names_gp = gpar(
    fontsize = 10
  ),
  
  column_names_gp = gpar(
    fontsize = 7
  ),
  
  ####################################
  ## No default cell borders
  ####################################
  
  rect_gp = gpar(
    col = NA
  ),
  
  ####################################
  ## Add black significance boxes
  ####################################
  
  cell_fun = function(
    j,
    i,
    x,
    y,
    width,
    height,
    fill
  ) {
    
    if (
      isTRUE(
        sig_matrix[i, j]
      )
    ) {
      
      grid.rect(
        x = x,
        y = y,
        width = width,
        height = height,
        gp = gpar(
          fill = NA,
          col = "black",
          lwd = 1.5
        )
      )
    }
  }
)




####################################
## Fig 2e,f,g,h - Distribution of clinicopathologic features across tumor growth patterns
####################################

library(dplyr)
library(ggplot2)
library(scales)
library(vcd) 
library(ggplot2)
library(scales)

ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")

####################################
## Plotting dataframe
####################################

df <- ductal_meta %>%
  filter(!is.na(tumor_pattern), !is.na(pam50))

df$tumor_pattern <- factor(
  df$tumor_pattern,
  levels = c("Cell-dense", "Scattered", "Nests")
)

df$GRADE <- as.factor(df$GRADE)

df$KI67_CATEGORIES <- factor(
  df$KI67_CATEGORIES,
  levels = c("≤10", "10-20", ">20")
)

####################################
## Palettes
####################################

pam50_colors <- c(
  "LumA"        = "#F2D86D",
  "LumB"        = "#82AEEB",
  "Her2"        = "#96D4AC",
  "Basal"       = "#E48A87",
  "Normal"      = "#BDBDBD",
  "claudin-low" = "#A569BD",
  "NC"          = "#999999"
)

nodal_cols <- c(
  "N0" = "#C7C7C7",
  "N+" = "#009E73"
)

grade_pal <- c(
  "1" = "#1B9E77",
  "2" = "#D95F02",
  "3" = "#7570B3"
)

ki67_cols <- c(
  "≤10"   = "#56B4E9",
  "10-20" = "#E69F00",
  ">20"   = "#009E73"
)

####################################
## Helper: add global p value
####################################

add_global_pvalue <- function(plot, data, xvar, fillvar, y_pos = 1.08) {
  
  tmp <- data %>%
    filter(!is.na(.data[[xvar]]), !is.na(.data[[fillvar]]))
  
  tab <- table(tmp[[xvar]], tmp[[fillvar]])
  
  chi <- suppressWarnings(chisq.test(tab))
  
  if (any(chi$expected < 5)) {
    test <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
    test_name <- "Fisher"
  } else {
    test <- chi
    test_name <- "Chi-square"
  }
  
  p_lab <- paste0(" p = ", signif(test$p.value, 3))
  
  plot +
    annotate(
      "text",
      x = 2,
      y = y_pos,
      label = p_lab,
      size = 4.5
    ) +
    coord_cartesian(
      ylim = c(0, y_pos + 0.05),
      clip = "off"
    )
}

####################################
## Helper: common stacked bar plot
####################################

make_stacked_bar <- function(data, fill_var, fill_label, fill_colors) {
  
  ggplot(data, aes(x = tumor_pattern, fill = .data[[fill_var]])) +
    geom_bar(position = "fill", color = "white", linewidth = 0.2) +
    scale_y_continuous(labels = percent_format()) +
    scale_fill_manual(values = fill_colors, na.translate = FALSE) +
    labs(
      x = NULL,
      y = "Percentage",
      fill = fill_label
    ) +
    theme_bw() +
    theme(
      text = element_text(size = 14),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
      axis.text.y = element_text(size = 12),
      axis.title = element_text(size = 13),
      legend.title = element_text(size = 12),
      legend.text = element_text(size = 11),
      plot.margin = margin(10, 10, 10, 10)
    )
}

####################################
## PAM50
####################################

df_pam50 <- df %>%
  filter(!is.na(pam50))

p_pam50 <- make_stacked_bar(
  data = df_pam50,
  fill_var = "pam50",
  fill_label = "PAM50",
  fill_colors = pam50_colors
)

p_pam50 <- add_global_pvalue(
  p_pam50,
  data = df_pam50,
  xvar = "tumor_pattern",
  fillvar = "pam50"
)

####################################
## Nodal status
####################################

df_nodal <- df %>%
  filter(!is.na(N_STATUS))

p_nodal <- make_stacked_bar(
  data = df_nodal,
  fill_var = "N_STATUS",
  fill_label = "Nodal status",
  fill_colors = nodal_cols
)

p_nodal <- add_global_pvalue(
  p_nodal,
  data = df_nodal,
  xvar = "tumor_pattern",
  fillvar = "N_STATUS"
)

####################################
## Grade
####################################

df_grade <- df %>%
  filter(!is.na(GRADE))

p_grade <- make_stacked_bar(
  data = df_grade,
  fill_var = "GRADE",
  fill_label = "Grade",
  fill_colors = grade_pal
)

p_grade <- add_global_pvalue(
  p_grade,
  data = df_grade,
  xvar = "tumor_pattern",
  fillvar = "GRADE"
)

####################################
## Ki67
####################################

df_ki67 <- df %>%
  filter(!is.na(KI67_CATEGORIES))

p_ki67 <- make_stacked_bar(
  data = df_ki67,
  fill_var = "KI67_CATEGORIES",
  fill_label = "Ki67",
  fill_colors = ki67_cols
)

p_ki67 <- add_global_pvalue(
  p_ki67,
  data = df_ki67,
  xvar = "tumor_pattern",
  fillvar = "KI67_CATEGORIES"
)

####################################
## Print individual plots
####################################

print(p_pam50)
print(p_nodal)
print(p_grade)
print(p_ki67)

####################################
## Save individual plots
####################################

ggsave(
  "Fig2f_PAM50_by_tumor_pattern_with_pvalue.pdf",
  p_pam50,
  width = 5,
  height = 5
)

ggsave(
  "Fig2f_Nodal_status_by_tumor_pattern_with_pvalue.pdf",
  p_nodal,
  width = 5,
  height = 5
)

ggsave(
  "Fig2f_Grade_by_tumor_pattern_with_pvalue.pdf",
  p_grade,
  width = 5,
  height = 5
)

ggsave(
  "Fig2f_Ki67_by_tumor_pattern_with_pvalue.pdf",
  p_ki67,
  width = 5,
  height = 5
)

####################################
## Combined figure: 1 row, 4 columns
####################################

fig2f_combined <- p_pam50 + p_nodal + p_grade + p_ki67 +
  plot_layout(nrow = 1)

print(fig2f_combined)

ggsave(
  "/Users/bengisukarakose/Desktop/NATCOM REVISION/0_last revision files 7.10/new_figs/Fig2efgh_clinicopathologic_features_by_tumor_pattern_with_pvalues.pdf",
  fig2f_combined,
  width = 18,
  height = 4.8
)




####################################
## Print statistical results
####################################

print_global_test <- function(data, xvar, fillvar, label) {
  
  tmp <- data %>%
    filter(!is.na(.data[[xvar]]), !is.na(.data[[fillvar]]))
  
  tab <- table(tmp[[xvar]], tmp[[fillvar]])
  
  chi <- suppressWarnings(chisq.test(tab))
  
  if (any(chi$expected < 5)) {
    test <- fisher.test(tab, simulate.p.value = TRUE, B = 10000)
    test_name <- "Fisher exact test simulated"
  } else {
    test <- chi
    test_name <- "Chi-square test"
  }
  
  cat("\n\n===== ", label, " =====\n", sep = "")
  cat(test_name, "\n")
  print(test)
  
  if (all(dim(tab) > 1)) {
    cat("Cramer's V:", assocstats(tab)$cramer, "\n")
  }
}

print_global_test(df_pam50, "tumor_pattern", "pam50", "PAM50")
print_global_test(df_nodal, "tumor_pattern", "N_STATUS", "Nodal status")
print_global_test(df_grade, "tumor_pattern", "GRADE", "Grade")
print_global_test(df_ki67, "tumor_pattern", "KI67_CATEGORIES", "Ki67")

