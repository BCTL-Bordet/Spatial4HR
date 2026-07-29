# 13_CNV_mutation.R

##################################################
## Global CNV burden across Spatial4HR+ subtypes in METABRIC 
##################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2)
library(pheatmap)
library(broom)
library(RColorBrewer)
library(ggpubr)

data_cna <- read.delim("~/Desktop/bc_cox_model/brca_metabric/data_cna.txt")[,-2]
metabric_meta <- readRDS("~/Desktop/final_scripts/00_data/metabric_meta.RDS")[,-25]

colnames(data_cna)[-1] <- gsub("\\.", "-", colnames(data_cna)[-1])

# Subset for matched samples
common_samples <- intersect(colnames(data_cna)[-1], metabric_meta$PATIENT_ID)
data_cna_filtered <- data_cna[, c("Hugo_Symbol", common_samples)]
metabric_meta_filtered <- metabric_meta %>% filter(PATIENT_ID %in% common_samples)

# Melt to long format and join subtypes
data_cna_long <- melt(data_cna_filtered, id.vars = "Hugo_Symbol",
                      variable.name = "PATIENT_ID", value.name = "CNV_Status") %>%
  left_join(metabric_meta_filtered[, c("PATIENT_ID", "predicted_subtype")], by = "PATIENT_ID") %>%
  filter(!is.na(predicted_subtype)) %>%
  mutate(Altered = ifelse(CNV_Status != 0, 1, 0))

##################################################
## Fig. 7f - CNV Burden per Subtype
##################################################

cna_long <- data_cna_long
colnames(cna_long) <- c("Hugo_Symbol", "name", "CNV",  "predicted_subtype" ,"Altered" )

cnv_burden <- data_cna_long %>%
  mutate(altered = CNV_Status != 0) %>%
  group_by(PATIENT_ID, predicted_subtype) %>%
  summarise(
    cnv_burden = sum(altered, na.rm = TRUE),
    .groups = "drop"
  )

colnames(cnv_burden) <- c('name', 'Subtype', 'cnv_burden')

pairwise.wilcox.test(cnv_burden$cnv_burden, cnv_burden$Subtype,
                     p.adjust.method = "fdr")  # Adjust for multiple testing

subtype_colors <- c(
  "SF" = "#E46B71",   # red
  "EnR" = "#8CBF40",  # green
  "PI" = "#38A7B0",   # teal
  "HL" = "#8C6BB1"    # purple
)

cnv_burden$Subtype <- factor(
  cnv_burden$Subtype,
  levels = c("EnR", "SF", "PI", "HL"))

ggplot(cnv_burden, aes(x = Subtype, y = cnv_burden, fill = Subtype)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_jitter(aes(color = Subtype), width = 0.2, size = 1.5, alpha = 0.6) +
  scale_fill_manual(values = subtype_colors) +
  scale_color_manual(values = subtype_colors) +
  stat_compare_means(
    method = "wilcox.test", 
    comparisons = list(
      c("HL", "SF"),
      c("HL", "PI"),
      c("HL", "EnR"),
      c("SF", "PI"),
      c("SF", "EnR"),
      c("PI", "EnR")
    ),
    label = "p.signif"
  ) +
  theme_minimal(base_size = 14) +
  labs(
    title = "CNV Burden (Gain+Loss) Across Subtypes",
    y = "CNV Burden (Number of Altered Genes)",
    x = "Spatial4HR+ Subtype"
  )




##################################################
## Mutations - Spatial4HR+ subtypes in METABRIC 
##################################################

data_mutations <- read.delim("~/Desktop/bc_cox_model/brca_metabric/data_mutations.txt", comment.char="#")
data_mutations$PATIENT_ID <- gsub("MTS-T", "MB-", data_mutations$Tumor_Sample_Barcode)

metabric_meta <- readRDS("~/Desktop/final_scripts/00_data/metabric_meta.RDS")[,-25]

length(unique(data_mutations$PATIENT_ID))
length(unique(metabric_meta$PATIENT_ID))

mutations_matched <- data_mutations %>%
  filter(PATIENT_ID %in% metabric_meta$PATIENT_ID)

# Keep only functional mutations
functional_classes <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Ins", "Frame_Shift_Del",
                        "Splice_Site", "In_Frame_Ins", "In_Frame_Del")

mut_filtered <- mutations_matched %>%
  filter(Variant_Classification %in% functional_classes)

# Create binary mutation matrix
mut_matrix <- mut_filtered %>%
  dplyr::select(PATIENT_ID, Hugo_Symbol) %>%
  distinct() %>%
  mutate(Mutated = 1) %>%
  pivot_wider(names_from = Hugo_Symbol, values_from = Mutated, values_fill = 0)

# Convert to data frame with rownames
mut_matrix_df <- as.data.frame(mut_matrix)
rownames(mut_matrix_df) <- mut_matrix_df$PATIENT_ID
mut_matrix_df$PATIENT_ID <- NULL

table(metabric_meta$predicted_subtype)

merged_data <- metabric_meta %>%
  dplyr::select(PATIENT_ID, predicted_subtype = predicted_subtype) %>% 
  inner_join(mutate(mut_matrix_df, PATIENT_ID = rownames(mut_matrix_df)), by = "PATIENT_ID")

rownames(merged_data) <- merged_data$PATIENT_ID
merged_data$PATIENT_ID <- NULL

(merged_data)[1:10,1:10]


####################################
# Fig. 7g - Differential mutation frequencies across Spatial4HR+ subtypes in METABRIC
####################################

library(purrr)
library(broom)

# Fisher's test for one gene
fisher_by_gene <- function(gene_name) {
  tbl <- table(merged_data[[gene_name]], merged_data$predicted_subtype)
  if (nrow(tbl) == 2 && all(colSums(tbl) > 0)) {
    test <- fisher.test(tbl)
    tibble(
      Gene = gene_name,
      p.value = test$p.value,
      odds_ratio = test$estimate[[1]]
    )
  } else {
    # Skip genes with 0 counts in any group
    tibble(Gene = gene_name, p.value = NA, odds_ratio = NA)
  }
}

# Apply to all gene columns
gene_cols <- setdiff(colnames(merged_data), "predicted_subtype")
gene_stats <- map_dfr(gene_cols, fisher_by_gene)

# Adjust p-values for multiple testing
gene_stats <- gene_stats %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  arrange(p_adj)

# Filter for significant genes (e.g., FDR < 0.05)
sig_mutations <- gene_stats %>%
  filter(p_adj < 0.05)

head(sig_mutations)
# Gene    p.value    p_adj
# 1 PIK3CA 1.53e-12 2.63e-10
# 2 GATA3  1.09e- 7 9.35e- 6
# 3 CBFB   2.19e- 4 1.26e- 2


top_genes <- sig_mutations$Gene

plot_data <- merged_data %>%
  dplyr::select(predicted_subtype, all_of(top_genes)) %>%
  pivot_longer(cols = all_of(top_genes), names_to = "Gene", values_to = "Mutation_Status") %>%
  group_by(predicted_subtype, Gene) %>%
  summarise(Freq = mean(Mutation_Status), .groups = "drop")

subtype_colors <- c(
  "SF" = "#E46B71",   # red
  "EnR" = "#8CBF40",  # green
  "PI" = "#38A7B0",   # teal
  "HL" = "#8C6BB1"    # purple
)


plot_data$predicted_subtype <- factor(plot_data$predicted_subtype, levels = c('EnR', 'SF', 'PI', 'HL'))

ggplot(plot_data, aes(x = predicted_subtype, y = Freq, fill = predicted_subtype)) +
  geom_col(position = "dodge", width = 0.7) +
  facet_wrap(~ Gene, scales = "free_y") +
  ylab("Mutation Frequency") +
  xlab("Spatial Subtype") +
  ggtitle("Subtype-specific Mutation Frequencies") +
  scale_fill_manual(values = subtype_colors) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1)
  )



