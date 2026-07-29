# 11_genomic_risk.R


##################################################
## Fig. 6a -  RFS stratified by Spatial4HR+ subtype
##################################################

ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")

# 1) Order
ductal_meta$ductal_subtype <- factor(
  ductal_meta$ductal_subtype,
  levels = c("EnR","SF","PI","HL")           # display names
)

# 2) Fit KM using years on the x-axis
km_fit <- survfit(Surv(time/365, status) ~ ductal_subtype, data = ductal_meta)

# 3) Plot with unnamed palette in the same order as factor levels
pal <- c("#8CBF40", "#E46B71", "#38A7B0", "#8C6BB1")  # EnR, SF, PI, HL

ggsurvplot(
  km_fit,
  data = ductal_meta,
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = TRUE,          # show row labels
  risk.table.y.text.col = TRUE,      # color row labels by strata (optional)
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (years)",
  ylab = "Survival Probability",
  title = "Kaplan–Meier Survival Curves by Spatial Subtypes",
  legend.title = "Spatial Subtype",
  legend.labs = levels(ductal_meta$ductal_subtype),  # EnR, SF, PI, HL
  palette = pal,
  break.time.by = 3                  # ticks every 3 years (adjust as you like)
)


################################################
## Molecular subtyping and genomic risk scoring
################################################

library(AIMS)
library(genefu)
library(org.Hs.eg.db)
library(AnnotationDbi)

######################
## 1. PAM50 intrinsic molecular subtyping using AIMS
######################

## Load whole-section pseudobulk expression matrix
expression_matrix <- readRDS(
  "/Users/bengisukarakose/Desktop/final_scripts/00_data/pb_ductals_not_norm.RDS")

## Map gene symbols to Entrez Gene identifiers
hs <- org.Hs.eg.db
my.symbols <- rownames(expression_matrix)

annotation <- AnnotationDbi::select(
  hs,
  keys = my.symbols,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL")

## Remove duplicated gene-symbol mappings
annotation_not_dup <- annotation[
  !duplicated(annotation$SYMBOL),]

## Run AIMS
matrice2 <- as.matrix(expression_matrix)

pam50_aims <- applyAIMS(
  matrice2,
  annotation_not_dup$ENTREZID
)

pam50class <- pam50_aims$cl

table(pam50class)




######################
## 2. Genomic risk scoring using genefu
######################

## Load log2-transformed CPM whole-section pseudobulk expression matrix
counts_ductal <- readRDS("~/Desktop/final_scripts/00_data/ductal_pb_rpm_log.RDS")

## Map gene symbols to Entrez Gene identifiers
my.symbols <- rownames(counts_ductal)

annotation <- AnnotationDbi::select(
  hs,
  keys = my.symbols,
  columns = c("ENTREZID", "SYMBOL"),
  keytype = "SYMBOL"
)

colnames(annotation) <- c(
  "name",
  "EntrezGene.ID"
)

annotation <- na.omit(annotation)
annotation <- annotation[
  !duplicated(annotation$name),
]

rownames(annotation) <- annotation$name

######################
## 3. Genomic Grade Index
######################

data("sig.ggi")

ggi_calc <- ggi(
  data = t(as.matrix(counts_ductal)),
  annot = annotation,
  do.mapping = TRUE
)

score_ggi_genefu <- ggi_calc$score

######################
## 4. MammaPrint 70-gene score and risk classification
######################

data("sig.gene70")

gene70_res <- gene70(
  data = t(as.matrix(counts_ductal)),
  annot = annotation,
  do.mapping = TRUE
)

score_gene70 <- gene70_res$score
risk_gene70 <- gene70_res$risk

######################
## 5. Oncotype DX recurrence score and risk classification
######################

data("sig.oncotypedx")

oncotype_res <- oncotypedx(
  data = t(as.matrix(counts_ductal)),
  annot = annotation,
  do.mapping = TRUE
)

score_oncotypedx <- oncotype_res$score
risk_oncotypedx <- oncotype_res$risk

## inspect outputs
score_ggi_genefu
score_gene70
risk_gene70
score_oncotypedx
risk_oncotypedx




##################################################
## Fig. 6b - Patient-level genomic risk scores derived from 
## Oncotype DX and MammaPrint across Spatial4HR+ subtypes
##################################################

ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")

score_gene70
score_ggi_genefu
score_oncotypedx

library(dplyr)
library(ggplot2)
library(ggpubr)
library(FSA)

subtype_colors <- c(
  "SF"  = "#E46B71",  # red/pink
  "EnR" = "#8CBF40",  # green
  "PI"  = "#38A7B0",  # teal
  "HL"  = "#8C6BB1"   # purple
)

# ---- Define score list ----
score_list <- list(
  Gene70 = score_gene70,
  GGI_Genefu = score_ggi_genefu,
  OncotypeDX = score_oncotypedx
)

# ---- Iterate over each score ----
for (score_name in names(score_list)) {
  cat("\nProcessing:", score_name, "\n")

  # Prepare score dataframe
  score_df <- data.frame(
    name = names(score_list[[score_name]]),
    score = as.numeric(score_list[[score_name]]),
    stringsAsFactors = FALSE
  )

  # Merge with subtype data
  ductal_meta_scored <- left_join(score_df, ductal_meta, by = "name")
  ductal_meta_scored$ductal_subtype <- factor(ductal_meta_scored$ductal_subtype,
                                              levels = c('EnR', 'SF', 'PI', 'HL'))


  # Kruskal-Wallis test
  kruskal_res <- kruskal.test(score ~ ductal_subtype, data = ductal_meta_scored)
  print(kruskal_res)

  # Dunn's test
  dunn_res <- dunnTest(score ~ ductal_subtype, data = ductal_meta_scored, method = "bh")
  print(dunn_res)

  # Set pairwise comparisons
  subtype_levels <- levels(ductal_meta_scored$ductal_subtype)
  comparisons <- combn(subtype_levels, 2, simplify = FALSE)

  # Max y for annotation positioning
  max_y <- max(ductal_meta_scored$score, na.rm = TRUE)

  # Create plot
  p <- ggplot(ductal_meta_scored, aes(x = ductal_subtype, y = score, fill = ductal_subtype)) +
    geom_boxplot() +
    geom_jitter(width = 0.2, alpha = 0.6) +
    stat_compare_means(comparisons = comparisons, method = "wilcox.test", label = "p.signif") +
    scale_fill_manual(values = subtype_colors) +
    labs(
      title = paste(score_name, "Score"),
      x = "Spatial Subtype",
      y = paste(score_name, "Score")
    ) +
    theme_minimal(base_size = 16) +   # increase base font size
    theme(
      plot.title = element_text(size = 14),
      axis.title = element_text(size = 14),
      axis.text = element_text(size = 12),
      legend.position = "none"   # removes the legend
    )

  # Show plot
  print(p)
}



## =========================================================
## Fig. 6c - Distribution of genomic risk categories assigned by 
## MammaPrint and Oncotype DX across Spatial4HR+ subtypes
## =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(binom)
})

ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")

## -------- 0) Output directory ----------
out_dir <- "~/Desktop/bc_cox_model/denememememe/risk_scores"
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

## -------- 1) Build long data frame ----------
risk_list <- list(
  OncotypeDX  = risk_oncotypedx,
  Gene70      = risk_gene70
)

risk_long <- bind_rows(
  tibble(
    sample = names(risk_list$OncotypeDX),
    assay  = "Oncotype DX",
    call   = as.numeric(risk_list$OncotypeDX) # may include 0.5 for Intermediate
  ),
  tibble(
    sample = names(risk_list$Gene70),
    assay  = "MammaPrint",
    call   = as.numeric(risk_list$Gene70)
  )
) %>% filter(!is.na(sample))

## -------- 2) Attach subtype info ----------
risk_long <- risk_long %>%
  left_join(ductal_meta, by = c("sample" = "name")) %>%
  mutate(
    ductal_subtype = factor(ductal_subtype, levels = c("EnR","SF","PI","HL"))
  ) %>%
  filter(!is.na(ductal_subtype), !is.na(call))

## -------- 3) Define risk categories ----------
risk_cats <- risk_long %>%
  mutate(
    risk_cat = case_when(
      call >= 1    ~ "High",
      call == 0.5  ~ "Intermediate",  # only expected for Oncotype
      call <= 0    ~ "Low",
      TRUE         ~ NA_character_
    ),
    high_flag = as.integer(risk_cat == "High"),
    low_flag  = as.integer(risk_cat == "Low")
  ) %>%
  filter(!is.na(risk_cat))

## -------- 4) % High-risk summary with 95% Wilson CI ----------
disc_high <- risk_cats %>%
  group_by(assay, ductal_subtype) %>%
  summarise(
    n = n(),
    k_high = sum(high_flag, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    prop_high = k_high / n,
    lower = binom::binom.confint(k_high, n, method = "wilson")[, "lower"],
    upper = binom::binom.confint(k_high, n, method = "wilson")[, "upper"],
    label = paste0(k_high, "/", n)
  )

## -------- 5) Risk composition (100% stacked) ----------
comp <- risk_cats %>%
  count(assay, ductal_subtype, risk_cat, name = "k") %>%
  group_by(assay, ductal_subtype) %>%
  mutate(pct = k / sum(k)) %>%
  ungroup()

## -------- 6) Colors ----------
subtype_colors <- c("EnR"="#8CBF40","SF"="#E46B71","PI"="#38A7B0","HL"="#8C6BB1")

risk_colors <- c("Low"='#FFDCDB', "Intermediate"='#F49795',  "High"='#800f2f')

## -------- 7) Plot: Risk composition ----------
comp <- comp %>%
  mutate(risk_cat = factor(risk_cat, levels = c("Low", "Intermediate", "High")))

p_stack <- ggplot(comp, aes(x = ductal_subtype, y = pct, fill = risk_cat)) +
  geom_col(width = 0.72, color = "black") +
  facet_wrap(~ assay, nrow = 1) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  scale_fill_manual(
    values = risk_colors,
    breaks = c("Low", "Intermediate", "High"),  # enforce order in legend
    name = "Risk call"
  ) +
  labs(
    title = "Risk-call composition within each subtype",
    x = "Spatial Subtype",
    y = "Composition (%)"
  ) +
  theme_minimal(base_size = 18) +   # global font size
  theme(
    plot.title   = element_text(size = 18, hjust = 0.5),
    axis.title   = element_text(size = 16),
    axis.text    = element_text(size = 14),
    strip.text   = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  )


print(p_stack)





#############################################################
## Fig. 6d - Kaplan–Meier survival curve
## EnR vs SF within MammaPrint (gene70) low-risk tumors
#############################################################

library(dplyr)
library(survival)
library(survminer)

 
## 1) Select MammaPrint low-risk samples (gene70)
mp_lowrisk <- risk_long %>%
  filter(
    assay == "MammaPrint",
    call == 0,                          # gene70 low risk
    ductal_subtype %in% c("EnR", "SF")
  ) %>%
  dplyr::select(
    sample
  ) %>%
  distinct()

 
## 2) Merge with survival metadata
km_df <- ductal_meta %>%
  inner_join(
    mp_lowrisk,
    by = c("name" = "sample")
  )

 
## 3) Set factor levels 
km_df <- km_df %>%
  mutate(
    ductal_subtype = factor(
      ductal_subtype,
      levels = c("EnR", "SF")))

 
## 4) Fit Kaplan–Meier model
fit <- survfit(
  Surv(time / 365.25, status) ~ ductal_subtype,
  data = km_df
)

 
## 5) Colors
pal <- c(
  "EnR" = "#8CBF40",   # green
  "SF"  = "#E46B71"    # red
)
 
## 6) Plot 
p <- ggsurvplot(
  fit,
  data = km_df,
  
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = TRUE,
  risk.table.y.text.col = TRUE,
  
  pval = TRUE,
  conf.int = FALSE,
  
  xlab = "Time (years)",
  ylab = "Relapse-free survival probability",
  
  legend.title = "Spatial program\n(MammaPrint low-risk)",
  legend.labs = c(
    "EnR",
    "SF"
  ),
  palette = pal,
  break.time.by = 3,
  ggtheme = theme_classic(base_size = 12)
)

print(p)



#############################################################
## Fig. 6e - EnR vs SF Stroma Pseudobulk Comparison 
## on MammaPrint low risk patients
#############################################################

library(dplyr)
library(stringr)
library(GSEABase)
library(GSVA)
library(limma)
library(pheatmap)
library(tibble)
library(ggplot2)

tumor_pb  <- readRDS("~/Desktop/final_scripts/00_data/50_percent_tumor_pb.RDS")
stroma_pb <- readRDS("~/Desktop/final_scripts/00_data/100_percent_stroma_pb.RDS")

ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")

## Spatial subtype annotation 
meta <- ductal_meta %>%
  mutate(
    sample = paste0("ST", orig.ident),
    ductal_subtype = factor(
      ductal_subtype,
      levels = c("EnR", "SF", "PI", "HL")
    )
  ) %>%
  dplyr::select(sample, ductal_subtype)

## Keep EnR vs SF only
meta_ES <- meta %>% filter(ductal_subtype %in% c("EnR", "SF"))

library(qusage)

genesets <- read.gmt(
  file = "~/Desktop/final_scripts/00_data/h.all.v2023.2.Hs.symbols.gmt")

## Prepare stroma matrix 
stroma_mat <- as.data.frame(stroma_pb)
colnames(stroma_mat) <- paste0("ST", colnames(stroma_mat))

# remove samples with all-NA stroma (ST80)
stroma_mat <- stroma_mat[, colSums(is.na(stroma_mat)) < nrow(stroma_mat)]

stroma_samples <- intersect(colnames(stroma_mat), meta_ES$sample)
stroma_mat <- stroma_mat[, stroma_samples]

meta_stroma <- meta_ES %>%
  filter(sample %in% stroma_samples) %>%
  arrange(match(sample, stroma_samples))

## GSVA (stroma) 
params_stroma <- gsvaParam(
  as.matrix(stroma_mat),
  kcdf = "Gaussian",
  geneSets = genesets
)

gsva_res_stroma <- gsva(params_stroma)

module_sig_stroma <- t(gsva_res_stroma)
module_sig_stroma_m <- as.matrix(module_sig_stroma)
colnames(module_sig_stroma_m) <- gsub("HALLMARK_", "", colnames(module_sig_stroma_m))

## limma (stroma) 
design_stroma <- model.matrix(~ ductal_subtype, data = meta_stroma)
fit_stroma <- lmFit(t(module_sig_stroma_m), design_stroma)
fit_stroma <- eBayes(fit_stroma)

stats_stroma <- topTable(
  fit_stroma,
  coef = "ductal_subtypeSF",
  number = Inf,
  adjust.method = "fdr"
) %>%
  rownames_to_column("Signature") %>%
  arrange(adj.P.Val)

head(stats_stroma)
stats_stroma[stats_stroma$adj.P.Val<0.05,]


## Stroma pseudobulk plot
## (Top 10 stromal programs, EnR vs SF)

subtype_cols <- c(
  EnR = "#8CBF40",
  SF  = "#E46B71"
)

## Select top 10 stromal signatures (FDR-ranked) ------------
stromal_sigs <- stats_stroma %>%
  filter(adj.P.Val < 0.05) %>%
  arrange(adj.P.Val) %>%
  slice_head(n = 10) %>%
  pull(Signature)

# Explicitly remove UV_RESPONSE_DN if present
stromal_sigs <- setdiff(stromal_sigs, "UV_RESPONSE_DN")

## Subset + scale GSVA matrix -------------------------------
heat_stroma <- module_sig_stroma_m[, stromal_sigs, drop = FALSE]
heat_stroma <- scale(heat_stroma)

## Pretty pathway names -------------------------------------
pretty_names <- c(
  ANGIOGENESIS                      = "Angiogenesis",
  EPITHELIAL_MESENCHYMAL_TRANSITION = "Epithelial–Mesenchymal Transition",
  ESTROGEN_RESPONSE_LATE            = "Estrogen Response (Late)",
  ESTROGEN_RESPONSE_EARLY           = "Estrogen Response (Early)",
  TGF_BETA_SIGNALING                = "TGF-B Signaling",
  HEDGEHOG_SIGNALING                = "Hedgehog Signaling",
  APICAL_JUNCTION                   = "Apical Junction",
  COAGULATION                       = "Coagulation",
  KRAS_SIGNALING_UP                 = "KRAS Signaling (Up)",
  KRAS_SIGNALING_DN                 = "KRAS Signaling (Down)",
  WNT_BETA_CATENIN_SIGNALING        = "WNT/β-Catenin Signaling",
  E2F_TARGETS                       = "E2F Targets",
  MYOGENESIS                        = "Myogenesis",
  MYC_TARGETS_V1                    = "MYC Targets (V1)",
  MYC_TARGETS_V2                    = "MYC Targets (V2)",
  OXIDATIVE_PHOSPHORYLATION         = "Oxidative Phosphorylation",
  NOTCH_SIGNALING                   = "Notch Signaling",
  IL2_STAT5_SIGNALING               = "IL2–STAT5 Signaling",
  HYPOXIA                           = "Hypoxia"
)

# Apply readable labels
colnames(heat_stroma) <- pretty_names[colnames(heat_stroma)]

## Row annotation 
annotation_stroma <- data.frame(
  `Spatial4HR+ subtype` = meta_stroma$ductal_subtype,
  row.names = meta_stroma$sample,
  check.names = FALSE
)

colnames(annotation_stroma) <- "Spatial4HR+ subtype"

ann_colors <- list(
  `Spatial4HR+ subtype` = subtype_cols)

## Heatmap
pheatmap(
  heat_stroma,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  annotation_row = annotation_stroma,
  annotation_names_row = FALSE,
  annotation_colors = ann_colors,
  color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
  main = "Stroma pseudobulk: EnR vs SF",
  fontsize_row = 8,
  fontsize_col = 10,
  angle_col = '45',         # <<< rotate x-axis labels
  border_color = NA,
  legend = F
)


