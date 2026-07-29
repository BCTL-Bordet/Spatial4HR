# 02_fig1_histopathology.R

####################################
## Fig 1a is BioRender image 
####################################

####################################
## Fig 1b – Cohort characteristics
####################################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

## Read data
ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")

df <- ductal_meta

## Define factors
df_plot <- df %>%
  mutate(
    GRADE_F   = ifelse(is.na(GRADE), NA, paste0("G", as.integer(GRADE))),
    GRADE_F   = factor(GRADE_F, levels = c("G1", "G2", "G3")),
    T_STAGE_F = factor(T_STAGE, levels = c("T1","T2","T3","T4")),
    NODAL_F   = factor(N_STATUS, levels = c("N0","N+")),
    RELAPSE_F = factor(
      case_when(
        is.na(Relapse) ~ NA_character_,
        Relapse == 1   ~ "Yes",
        TRUE           ~ "No"
      ),
      levels = c("No","Yes")
    )
  )


grade_cols   <- c("G1"="#1B9E77","G2"="#D95F02","G3"="#7570B3")
stage_cols   <- c("T1"="#56B4E9","T2"="#009E73","T3"="#F0E442","T4"="#D55E00")
nodal_cols   <- c("N0"="#E69F00","N+"="#0072B2")
relapse_cols <- c("No"="#7F7F7F","Yes"="#CC0000")

theme_nc_bar <- theme_classic(base_size = 11) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "plain", size = 13),
    axis.text.x  = element_text(size = 10, color = "black"),
    axis.text.y  = element_text(size = 11, face = "bold", color = "black"),
    axis.title   = element_blank(),
    axis.ticks   = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    legend.position   = "none",
    plot.margin       = margin(6, 14, 6, 6)
  )

plots_info <- list(
  list(var="GRADE_F",   levels=c("G1","G2","G3"),      title="Tumor grade",   cols=grade_cols),
  list(var="T_STAGE_F", levels=c("T1","T2","T3","T4"), title="Tumor size",   cols=stage_cols),
  list(var="NODAL_F",   levels=c("N0","N+"),           title="Nodal status",  cols=nodal_cols),
  list(var="RELAPSE_F", levels=c("No","Yes"),          title="Relapse",       cols=relapse_cols)
)

plots_bar <- list()



## Build plots

for (info in plots_info) {
  
  tmp <- df_plot %>%
    mutate(.x = .data[[info$var]]) %>%
    filter(!is.na(.x)) %>%
    dplyr::count(.x, name = "n") %>%
    complete(.x = factor(info$levels, levels = info$levels), fill = list(n = 0)) %>%
    mutate(
      N   = sum(n),
      pct = 100 * n / N,
      lab = paste0("n=", n),
      .x  = factor(.x, levels = info$levels)
    )
  
  
  p <- ggplot(tmp, aes(x = .x, y = pct, fill = .x)) +
    geom_col(width = 0.55) +
    geom_text(
      aes(x = .x, y = pct + 2.5, label = lab),
      hjust = 0.5, vjust = 0,
      size = 3.8
    ) +
    scale_fill_manual(values = info$cols) +
    scale_y_continuous(
      limits = c(0, max(tmp$pct) + 8),
      breaks = pretty_breaks(4),
      labels = function(x) paste0(x, "%")
    ) +
    labs(title = sprintf("%s (n=%d)", info$title, unique(tmp$N))) +
    theme_nc_bar
  
  plots_bar[[info$var]] <- p
}


panel_bar <- (plots_bar$GRADE_F | plots_bar$T_STAGE_F) /
  (plots_bar$NODAL_F | plots_bar$RELAPSE_F)

panel_bar

# ggsave(
#   "Fig1b_cohort_barplot.pdf",
#   panel_bar,
#   width  = 5.5,
#   height = 5,
#   dpi    = 600
# )



####################################
## Fig 1c – H&E and respective annotations
####################################


####################################
## Extended Data Fig.2  – Prevalence of histopathological compartments across tumors.
####################################

annotations <- ductal_meta %>% dplyr::select(contains("_annotation"))
head(annotations)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

normalize_key <- function(x) {
  x |>
    stringr::str_replace_all("_annotation$", "") |>  # drop suffix if present
    stringr::str_to_lower() |>
    stringr::str_replace_all("[^a-z0-9]+", "")       # keep letters/digits only
}

# build a named color vector from 'custer_cols' table
custer_cols <- readRDS("~/Desktop/bc_cox_model/103_paper_updates/FIGURES AND SCRIPTS/FIG1_cohort/fig1d//histo_colors.RDS")
col_key_norm <- normalize_key(custer_cols$V1)
ann_colors   <- setNames(custer_cols$col_vector, col_key_norm)

# long-format  annotations and count presence
presence_threshold <- 0   

ann_long <- annotations %>%
  select(where(is.numeric)) %>%                      # drop SEQ_ID or other non-numeric columns
  pivot_longer(
    cols = everything(),
    names_to   = "annotation_raw",
    values_to  = "value"
  ) %>%
  mutate(
    key_norm   = normalize_key(annotation_raw),      # normalized key to join colors
    ann_label  = annotation_raw |>
      str_replace("_annotation$", "") |> # pretty label for plotting
      str_replace_all("_", " ") |>
      str_replace_all("\\.", " ")
  ) %>%
  group_by(key_norm, ann_label) %>%
  summarise(n_patients = sum(value > presence_threshold, na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(n_patients))

# keep only annotations that are actually present at least once
ann_long <- ann_long %>% filter(n_patients > 0)

# order factors by count (descending)
ann_long <- ann_long %>%
  mutate(ann_label = factor(ann_label, levels = ann_label))

# match colors (fallback to grey if a color is missing)
plot_colors <- ann_colors[normalize_key(levels(ann_long$ann_label))]
names(plot_colors) <- levels(ann_long$ann_label)
plot_colors[is.na(plot_colors)] <- "#BDBDBD" 

p <- ggplot(ann_long, aes(x = ann_label, y = n_patients, fill = ann_label)) +
  geom_col(width = 0.6, color = "black", linewidth = 0.3) +       
  scale_fill_manual(values = plot_colors, guide = "none") +
  labs(x = NULL, y = "Number of samples") +
  geom_text(aes(label = n_patients),
            vjust = -0.4, size = 5.5) +         # bigger labels
  theme_classic(base_size = 14, base_family = "Helvetica") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 16),
    axis.text.y = element_text(size = 15),
    axis.title.y = element_text(size = 17, margin = margin(r = 10)),
    plot.title  = element_text(size = 18, hjust = 0.5),
    plot.margin = margin(8, 8, 8, 8)
  ) +
  expand_limits(y = max(ann_long$n_patients) * 1.15)

p
# ggsave("annotation_counts_bar_largefont.png", p, width = 7, height = 5, dpi = 600)
# ggsave("annotation_counts_bar_largefont.pdf", p, width = 7, height = 5, dpi = 600)





####################################
## Extended Data Fig. 1 - Distribution of annotated histopathological compartments.
####################################

# object_merged <- readRDS("~/Desktop/final_scripts/00_data/object_merged_final.RDS")
annotation_cols_prev <- c(
  "Tumor",
  "Necrosis",
  "Fat_tissue",
  "High_TILs_stroma",
  "Cellular_stroma",
  "Acellular_stroma",
  "Vessels",
  "Artefact",
  "Canal_galactophore",
  "Nodule_lymphoid",
  "In_situ",
  "Nerve",
  "Lymphocyte",
  "Hole",
  "Microcalcification",
  "Out",
  "Apocrine.metaplasia"
)

ann_mat <- object_merged@meta.data %>%
  dplyr::select(dplyr::all_of(annotation_cols_prev))

normalize_key <- function(x) {
  x |>
    stringr::str_replace_all("[._ ]+", "_") |>
    stringr::str_replace_all("_+$", "") |>
    stringr::str_to_lower()
}

col_key <- normalize_key(custer_cols$V1)
ann_colors <- setNames(custer_cols$col_vector, col_key)

totals <- colSums(ann_mat, na.rm = TRUE)
grand_total <- sum(totals)

plot_df <- tibble(
  annotation_raw = names(totals),
  total = as.numeric(totals)
) %>%
  mutate(
    pct = 100 * total / grand_total,
    
    # Keep original annotation name for colour matching
    key_norm = normalize_key(annotation_raw),
    
    # Change only the displayed label
    ann_label = dplyr::recode(
      annotation_raw,
      "Canal_galactophore" = "Normal_breast"
    ) |>
      stringr::str_replace_all("[._]", " ")
  ) %>%
  arrange(desc(pct))

plot_cols <- ann_colors[plot_df$key_norm]
plot_cols[is.na(plot_cols)] <- "#BDBDBD"
names(plot_cols) <- plot_df$ann_label

plot_df$ann_label <- factor(
  plot_df$ann_label,
  levels = plot_df$ann_label
)

p_horz <- ggplot(
  plot_df,
  aes(y = ann_label, x = pct, fill = ann_label)
) +
  geom_col(
    width = 0.55,
    color = "black",
    linewidth = 0.3
  ) +
  scale_fill_manual(
    values = plot_cols,
    guide = "none"
  ) +
  scale_x_continuous(
    labels = function(x) paste0(x, "%")
  ) +
  labs(
    x = "Total annotation percentage",
    y = NULL
  ) +
  geom_text(
    aes(label = sprintf("%.1f%%", pct)),
    hjust = -0.2,
    size = 5.2
  ) +
  coord_cartesian(
    xlim = c(0, max(plot_df$pct) * 1.15)
  ) +
  theme_classic(
    base_size = 14,
    base_family = "Helvetica"
  ) +
  theme(
    axis.text = element_text(size = 16),
    axis.title.x = element_text(
      size = 17,
      margin = margin(t = 10)
    ),
    plot.margin = margin(8, 18, 8, 8)
  )

print(p_horz)


####################################
## Fig 1d – Correlation heatmap (histology & deconvolution)
####################################

# object_merged <- readRDS("~/Desktop/final_scripts/00_data/object_merged_final.RDS")
ST_ubermeta <- object_merged@meta.data

library(ComplexHeatmap)
library(circlize)
library(grid)

# Cell-type proportion columns
celltype_cols <- c(
  "adipocyte",
  "Endothelial.ACKR1",
  "Endothelial.RGS5",
  "Endothelial.CXCL12",
  "CAFs.MSC.iCAF.like.s1",
  "CAFs.MSC.iCAF.like.s2",
  "CAFs.myCAF.like.s4",
  "PVL.Immature.s1",
  "Endothelial.Lymphatic.LYVE1",
  "B.cells.Memory",
  "T_cells_c4_CD8._ZFP36",
  "T_cells_c6_IFIT1",
  "T_cells_c7_CD8._IFNG",
  "T_cells_c8_CD8._LAG3",
  "T_cells_c0_CD4._CCR7",
  "T_cells_c1_CD4._IL7R",
  "T_cells_c2_CD4._T.regs_FOXP3",
  "T_cells_c3_CD4._Tfh_CXCL13",
  "T_cells_c9_NK_cells_AREG",
  "T_cells_c10_NKT_cells_FCGR3A",
  "Myeloid_c10_Macrophage_1_EGR1",
  "Myeloid_c12_Monocyte_1_IL1B",
  "Myeloid_c1_LAM1_FABP5",
  "Myeloid_c8_Monocyte_2_S100A9",
  "Epithelial",
  "CAFs.Transitioning.s3",
  "CAFs.myCAF.like.s5",
  "PVL.Differentiated.s3",
  "PVL_Immature.s2",
  "B.cells.Naive",
  "Plasmablasts",
  "Myeloid_c2_LAM2_APOE",
  "Myeloid_c9_Macrophage_2_CXCL10",
  "Myeloid_c11_cDC2_CD1C",
  "Myeloid_c4_DCs_pDC_IRF7",
  "Myeloid_c3_cDC1_CLEC9A",
  "Myeloid_c0_DC_LAMP3",
  "T_cells_c5_CD8._GZMK",
  "Myeloid_c7_Monocyte_3_FCGR3A"
)

# Histology annotation columns
annotation_cols <- c(
  "Tumor",
  "Necrosis",
  "Fat_tissue",
  "High_TILs_stroma",
  "Cellular_stroma",
  "Acellular_stroma",
  "Vessels",
  "Canal_galactophore",
  "In_situ",
  "Nerve",
  "Lymphocyte"
)

mat <- cor(
  object_merged@meta.data[, celltype_cols],
  object_merged@meta.data[, annotation_cols],
  use = "pairwise.complete.obs",
  method = "spearman"
)


mat_t <- t(mat)

# Clean labels: replace underscores and dots with spaces
rownames(mat_t) <- gsub("_", " ", rownames(mat_t))
colnames(mat_t) <- gsub("_", " ", colnames(mat_t))

rownames(mat_t) <- gsub("\\.", " ", rownames(mat_t))
colnames(mat_t) <- gsub("\\.", " ", colnames(mat_t))

rownames(mat_t) <- gsub("\\s+", " ", rownames(mat_t))
colnames(mat_t) <- gsub("\\s+", " ", colnames(mat_t))

rownames(mat_t) <- gsub("Canal galactophore", "Normal breast", rownames(mat_t))

col_fun <- colorRamp2(c(-0.2, 0, 0.2), c("#4575b4", "white", "#d73027"))

# Heatmap
ht <- Heatmap(
  mat_t,
  name = "Correlation",
  col = col_fun,
  cluster_rows = TRUE,
  cluster_columns = TRUE,
  show_row_names = TRUE,
  show_column_names = TRUE,
  row_names_gp = gpar(fontsize = 10),
  column_names_gp = gpar(fontsize = 10),
  column_names_rot = 45,
  row_title = "Morphological compartments",
  column_title = "Deconvoluted cell types",
  row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  column_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(
    color_bar = "continuous",
    legend_direction = "vertical",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 11, fontface = "bold"),
    labels_gp = gpar(fontsize = 10)
  )
)


draw(ht)




####################################
## Extended Data Fig. 3 – Associations between histology & clinical variables
####################################

library(reshape2)
library(dplyr)
library(ggplot2)

annotation_cols_clin <- c(
  "Tumor_annotation",
  "Necrosis_annotation",
  "Fat_tissue_annotation",
  "High_TILs_stroma_annotation",
  "Cellular.stroma_annotation",
  "Acellular.stroma_annotation",
  "Vessels_annotation",
  "Canal_galactophore_annotation",
  "In_situ_annotation",
  "Nerve_annotation"
)

scaled_data <- as.data.frame(apply(ductal_meta[,annotation_cols_clin], 2, genefu::rescale))

scaled_data <- scaled_data %>%
  dplyr::rename(
    Acellular_stroma_annotation = Acellular.stroma_annotation,
    Cellular_stroma_annotation = Cellular.stroma_annotation,
    Normal_breast_annotation = Canal_galactophore_annotation
  )


# ----------------
# Relapse 
# ----------------

scaled_relapse <- cbind(scaled_data, Relapse = ductal_meta$status)

relapse_df <- melt(scaled_relapse, id.vars = 'Relapse')
names(relapse_df) <- c("Relapse", "Variable", "Value")
relapse_df$Variable <- sub('_annotation', '', relapse_df$Variable)

# Perform Wilcoxon tests for each continuous column
test_results <- lapply(unique(relapse_df$Variable), function(var_name) {
  test_result <- wilcox.test(Value ~ Relapse, data = relapse_df[relapse_df$Variable == var_name, ])
  return(data.frame(Variable = var_name, p_value = test_result$p.value))
})

# Combine the results into a single dataframe
results_df <- do.call(rbind, test_results)

# Join p-values with the original data
relapse_df <- merge(relapse_df, results_df, by = "Variable", all.x = TRUE)
relapse_df$Variable <- gsub('_', ' ', relapse_df$Variable)

## ----- significance labels (BH-adjusted p) -----
p_values <- relapse_df %>%
  group_by(Variable) %>%
  summarise(p_value   = unique(p_value), .groups = "drop") %>%
  mutate(p_corrected = p.adjust(p_value, method = "BH"),
         p.signif = case_when(
           p_corrected <= 1e-4 ~ "****",
           p_corrected <= 1e-3 ~ "***",
           p_corrected <= 1e-2 ~ "**",
           p_corrected <= 5e-2 ~ "*",
           TRUE                ~ "ns"
         ))

## y-position for labels = a bit above the max within each Variable
y_pos <- relapse_df %>%
  group_by(Variable) %>%
  summarise(y = max(Value, na.rm = TRUE) + 0.03, .groups = "drop")
p_values <- left_join(p_values, y_pos, by = "Variable")

## ----- palette similar to reference (blue / orange) -----
pal <- c("0" = "#4C78A8",  # No relapse
         "1" = "#F58518")  # Relapse


## ----- the plot -----
p <- ggplot(relapse_df, aes(Variable, Value, fill = factor(Relapse))) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,                # don’t draw outlier dots from boxplot
    position = position_dodge2(width = 0.75, preserve = "single")
  ) +
  geom_point(
    aes(fill = factor(Relapse)),
    shape = 21, stroke = 0.3, alpha = 0.65, size = 1.4,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)
  ) +
  # significance stars
  geom_text(
    data = p_values,
    aes(x = Variable, y = y, label = p.signif),
    inherit.aes = FALSE, vjust = -0.2, size = 4
  ) +
  scale_fill_manual(values = pal, labels = c("No", "Relapse")) +
  labs(
    # title = "Histology Annotation Differences Between Relapse (–) vs (+) Samples",
    x = "",
    y = "Fraction",
    fill = "Relapse status"
  ) +
  coord_cartesian(ylim = c(0, NA)) +         # keep points; crop only display
  theme_classic(base_size = 12) +
  theme(
    panel.border   = element_blank(),        # ← removes the black border
    axis.text.x    = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title   = element_text(size = 11),
    legend.text    = element_text(size = 10),
    plot.title     = element_text(hjust = 0.5, face = "bold")
  )

p




# ----------------
# Grade
# ----------------

scaled_grade <- cbind(scaled_data, Grade = ductal_meta$GRADE)

# Melt the data for easy plotting
grade_df <- melt(scaled_grade, id.vars = 'Grade')
names(grade_df) <- c("Grade", "Variable", "Value")
grade_df$Variable <- sub('_annotation', '', grade_df$Variable)

# Perform Wilcoxon tests for each continuous column
test_results <- lapply(unique(grade_df$Variable), function(var_name) {
  test_result <- kruskal.test(Value ~ Grade, data = grade_df[grade_df$Variable == var_name, ])
  return(data.frame(Variable = var_name, p_value = test_result$p.value))
})

# Combine the results into a single dataframe
results_df <- do.call(rbind, test_results)

# Join p-values with the original data
grade_df <- merge(grade_df, results_df, by = "Variable", all.x = TRUE)
grade_df$Variable <- gsub('_', ' ', grade_df$Variable)

## ----- significance labels (BH-adjusted p) -----
p_values <- grade_df %>%
  group_by(Variable) %>%
  summarise(p_value   = unique(p_value), .groups = "drop") %>%
  mutate(p_corrected = p.adjust(p_value, method = "BH"),
         p.signif = case_when(
           p_corrected <= 1e-4 ~ "****",
           p_corrected <= 1e-3 ~ "***",
           p_corrected <= 1e-2 ~ "**",
           p_corrected <= 5e-2 ~ "*",
           TRUE                ~ "ns"
         ))

## y-position for labels = a bit above the max within each Variable
y_pos <- grade_df %>%
  group_by(Variable) %>%
  summarise(y = max(Value, na.rm = TRUE) + 0.03, .groups = "drop")
p_values <- left_join(p_values, y_pos, by = "Variable")

## ----- palette  -----
pal <- c("1" = "#1B9E77", 
         "2" = "#D95F02", 
         "3" = "#7570B3")

## ----- the plot -----
p <- ggplot(grade_df, aes(Variable, Value, fill = factor(Grade))) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,                # don’t draw outlier dots from boxplot
    position = position_dodge2(width = 0.75, preserve = "single")
  ) +
  geom_point(
    aes(fill = factor(Grade)),
    shape = 21, stroke = 0.3, alpha = 0.65, size = 1.4,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)
  ) +
  # significance stars
  geom_text(
    data = p_values,
    aes(x = Variable, y = y, label = p.signif),
    inherit.aes = FALSE, vjust = -0.2, size = 4
  ) +
  scale_fill_manual(values = pal, labels = c("1", "2", "3")) +
  labs(
    # title = "Histology Annotation Differences Between Grade 1, 2 vs 3 Samples",
    x = "",
    y = "Fraction",
    fill = "Grade"
  ) +
  coord_cartesian(ylim = c(0, NA)) +         # keep points; crop only display
  theme_classic(base_size = 12) +
  theme(
    panel.border   = element_blank(),        # ← removes the black border
    axis.text.x    = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title   = element_text(size = 11),
    legend.text    = element_text(size = 10),
    plot.title     = element_text(hjust = 0.5, face = "bold")
  )

p


# ----------------
# Age 
# ----------------

ductal_meta$age_class <- ifelse(ductal_meta$AGE < 50, "<50", "≥50")
scaled_age <- cbind(scaled_data, Age = ductal_meta$age_class)

# Melt the data for easy plotting
age_df <- melt(scaled_age, id.vars = 'Age')
names(age_df) <- c("Age", "Variable", "Value")
age_df$Variable <- sub('_annotation', '', age_df$Variable)

# Perform tests for each continuous column
test_results <- lapply(unique(age_df$Variable), function(var_name) {
  test_result <- wilcox.test(Value ~ Age, data = age_df[age_df$Variable == var_name, ])
  return(data.frame(Variable = var_name, p_value = test_result$p.value))
})

# Combine the results into a single dataframe
results_df <- do.call(rbind, test_results)

# Join p-values with the original data
age_df <- merge(age_df, results_df, by = "Variable", all.x = TRUE)
age_df$Variable <- gsub('_', ' ', age_df$Variable)

# Convert to ordered factor
age_df$Age <- factor(age_df$Age, levels = c("<50", "≥50"), ordered = TRUE)

## ----- significance labels (BH-adjusted p) -----
p_values <- age_df %>%
  group_by(Variable) %>%
  summarise(p_value   = unique(p_value), .groups = "drop") %>%
  mutate(p_corrected = p.adjust(p_value, method = "BH"),
         p.signif = case_when(
           p_corrected <= 1e-4 ~ "****",
           p_corrected <= 1e-3 ~ "***",
           p_corrected <= 1e-2 ~ "**",
           p_corrected <= 5e-2 ~ "*",
           TRUE                ~ "ns"
         ))

## y-position for labels = a bit above the max within each Variable
y_pos <- age_df %>%
  group_by(Variable) %>%
  summarise(y = max(Value, na.rm = TRUE) + 0.03, .groups = "drop")
p_values <- left_join(p_values, y_pos, by = "Variable")

## ----- palette  -----
pal <- c("<50" = "#0072B2", "≥50" = "#E69F00")

## ----- the plot -----
p <- ggplot(age_df, aes(Variable, Value, fill = factor(Age))) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,                # don’t draw outlier dots from boxplot
    position = position_dodge2(width = 0.75, preserve = "single")
  ) +
  geom_point(
    aes(fill = factor(Age)),
    shape = 21, stroke = 0.3, alpha = 0.65, size = 1.4,
    position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)
  ) +
  # significance stars
  geom_text(
    data = p_values,
    aes(x = Variable, y = y, label = p.signif),
    inherit.aes = FALSE, vjust = -0.2, size = 4
  ) +
  scale_fill_manual(values = pal, labels = c("<50", "≥50")) +
  labs(
    # title = "Histology Annotation Differences Between Ages",
    x = "",
    y = "Fraction",
    fill = "Age"
  ) +
  coord_cartesian(ylim = c(0, NA)) +         # keep points; crop only display
  theme_classic(base_size = 12) +
  theme(
    panel.border   = element_blank(),        # ← removes the black border
    axis.text.x    = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title   = element_text(size = 11),
    legend.text    = element_text(size = 10),
    plot.title     = element_text(hjust = 0.5, face = "bold")
  )

p





####################################
## Extended Data Fig. 4 – Spot-level correlations between gene expression and histopathological annotations
####################################

library(Seurat)
library(dplyr)
library(ggplot2)
library(ggrepel)


## 1. Calculate gene–annotation correlations

# object_merged <- readRDS("~/Desktop/final_scripts/00_data/object_merged_final.RDS")

expr_mat <- GetAssayData(
  object_merged,
  assay = "SCT",
  slot = "data")

## Select the 4,000 most variable genes
gene_variances <- apply(expr_mat, 1, var)
top_genes <- names(sort(gene_variances, decreasing = TRUE))[1:4000]
expr_mat_top <- expr_mat[top_genes, ]
gene_names <- rownames(expr_mat_top)

## Transpose to spots × genes
expr_mat_t <- t(as.matrix(expr_mat_top))

## Extract spot-level histopathological annotation fractions
annotations <- object_merged@meta.data[,c(5:11, 13, 15:17)]

## Confirm identical spot ordering
stopifnot(
  identical(
    rownames(expr_mat_t),
    rownames(annotations)))

## Calculate Spearman correlations between gene expression and annotations
cor_results <- sapply(
  colnames(annotations),
  function(annotation_name) {
    cor(
      expr_mat_t,
      annotations[[annotation_name]],
      method = "spearman")})

rownames(cor_results) <- gene_names



## 2. Select the strongest positive and negative correlates

## Top 30 positively correlated genes per annotation
top30_genes_per_annotation <- apply(
  cor_results,
  2,
  function(x) {
    names(
      sort(x, decreasing = TRUE)
    )[1:30]
  }
)

top30_gene_list <- lapply(
  colnames(top30_genes_per_annotation),
  function(annotation_name) {
    top30_genes_per_annotation[, annotation_name]
  }
)

names(top30_gene_list) <- colnames(top30_genes_per_annotation)

top30_cor_values <- apply(
  cor_results,
  2,
  function(x) {
    sort(x, decreasing = TRUE)[1:30]
  }
)

## Top 30 negatively correlated genes per annotation
anticor30_genes_per_annotation <- apply(
  cor_results,
  2,
  function(x) {
    names(
      sort(x, decreasing = FALSE)
    )[1:30]
  }
)

anticor30_gene_list <- lapply(
  colnames(anticor30_genes_per_annotation),
  function(annotation_name) {
    anticor30_genes_per_annotation[, annotation_name]
  }
)

names(anticor30_gene_list) <- colnames(
  anticor30_genes_per_annotation
)

anticor30_cor_values <- apply(
  cor_results,
  2,
  function(x) {
    sort(x, decreasing = FALSE)[1:30]
  }
)

top30_cor_values <- as.data.frame(top30_cor_values)
anticor30_cor_values <- as.data.frame(anticor30_cor_values)



## 3. Prepare plotting table


create_data_table <- function(
    gene_list,
    cor_values,
    annotation_name
) {
  direction <- ifelse(
    cor_values > 0,
    "positive_cor",
    "negative_cor"
  )
  
  data.frame(
    gene = gene_list,
    cor_value = cor_values,
    cluster = annotation_name,
    direction = direction
  )
}

df_list <- list()

## Positive correlates
for (annotation_name in names(top30_gene_list)) {
  df_list[[annotation_name]] <- create_data_table(
    gene_list = top30_gene_list[[annotation_name]],
    cor_values = top30_cor_values[[annotation_name]],
    annotation_name = annotation_name
  )
}

## Negative correlates
for (annotation_name in names(anticor30_gene_list)) {
  anticor_name <- paste0(
    annotation_name,
    "_anticor"
  )
  
  df_list[[anticor_name]] <- create_data_table(
    gene_list = anticor30_gene_list[[annotation_name]],
    cor_values = anticor30_cor_values[[annotation_name]],
    annotation_name = anticor_name
  )
}

df_plot <- do.call(
  rbind,
  df_list
)

df_plot$cluster <- gsub("_anticor", "",  df_plot$cluster)

df_plot$cluster <- factor(df_plot$cluster)



## 4. Plot strongest positive and negative correlates


col_map <- c("positive_cor" = "#D73027", "negative_cor" = "#4575B4")

p <- ggplot(
  df_plot,
  aes(
    x = cluster,
    y = cor_value,
    color = direction
  )
) +
  geom_jitter(
    width = 0.25,
    size = 2,
    alpha = 0.8
  ) +
  scale_color_manual(
    values = col_map
  ) +
  theme_minimal(
    base_size = 14
  ) +
  labs(
    x = "Histopathological annotations",
    y = "Spearman correlation coefficient"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

## Label the five strongest positive and negative genes per annotation
top_genes_to_label <- df_plot %>%
  group_by(cluster) %>%
  arrange(
    cluster,
    desc(cor_value)
  ) %>%
  mutate(
    rank = row_number()
  ) %>%
  filter(
    rank <= 5 |
      rank > n() - 5
  ) %>%
  ungroup()

p_labeled <- p +
  geom_text_repel(
    data = top_genes_to_label,
    aes(label = gene),
    max.overlaps = 100,
    size = 3.5,
    box.padding = 0.25,
    show.legend = FALSE
  )

p_labeled



## 5. Export positively correlated gene sets in GMT format

gmt_path <- "top30_gene_sets_from_cor.gmt"

con <- file(gmt_path, open = "wt")

for (set_name in names(top30_gene_list)) {
  gene_vector <- top30_gene_list[[set_name]]
  line <- paste( c(set_name, "na", gene_vector), collapse = "\t")
  writeLines(line, con)}

close(con)


