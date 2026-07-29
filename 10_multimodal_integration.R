# 10_multimodal_integration.R

# ===============================
# MOVICS MultiModal Integration Pipeline
# ===============================

library(MOVICS)
library(dplyr)
library(ggplot2)
library(patchwork)
library(ggpubr)
library(tidyr)
library(survival)
library(survminer)
library(GSVA)
library(cluster)
library(pheatmap)
library(msigdbr)
library(limma)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(grid)

Spatial4HR_sample_metadata <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")

module_scores_meta <- Spatial4HR_sample_metadata[,c('M4', 'M14')]
rownames(module_scores_meta) <- Spatial4HR_sample_metadata$name

cluster_props_matrix <- read.delim("~/Desktop/final_scripts/00_data/cluster_percentages.txt", row.names=1)

annotation_cols <- c(
    "Tumor_annotation",
    "Necrosis_annotation",
    "Fat_tissue_annotation",
    "High_TILs_stroma_annotation",
    "Cellular.stroma_annotation",
    "Acellular.stroma_annotation",
    "Vessels_annotation",
    "Canal_galactophore_annotation",
    "In_situ_annotation",
    "Nerve_annotation",
    "Lymphocyte_annotation",
    "Microcalcification_annotation")

histo_annot <- Spatial4HR_sample_metadata[, annotation_cols, drop = FALSE]
rownames(histo_annot) <- Spatial4HR_sample_metadata$name

# ductal_meta <- readRDS("~/Desktop/PhD/latest data 11jul/ductals_clinics_meta.RDS")
# histo_annot <- ductal_meta[, c(4:11, 13:16)]
# rownames(histo_annot) <- ductal_meta$name

# clinical <- ductal_meta[, c(2,3,20,25,28,29,65,66)]
# rownames(clinical) <- ductal_meta$name

# object_merged <- readRDS("~/Desktop/final_scripts/00_data/object_merged_final.RDS")
ST_ubermeta <- object_merged@meta.data

ST_ubermeta$id <- paste(ST_ubermeta$barcode, ST_ubermeta$orig.ident, sep = "_")
CARD_subset <- colnames(ST_ubermeta)[31:69]

col_means_by_patient <- ST_ubermeta %>%
  group_by(orig.ident) %>%
  summarise(across(all_of(CARD_subset), mean, na.rm = TRUE), .groups = "drop")
col_means_by_patient <- as.data.frame(col_means_by_patient)
rownames(col_means_by_patient) <- paste0('ST', col_means_by_patient$orig.ident)
col_means_by_patient <- col_means_by_patient[, -1]

# --- Build Multi-Omics List ---
omicsList <- list(
  histo              = t(histo_annot),
  module             = t(module_scores_meta),
  clusters           = t(cluster_props_matrix),
  cell_infiltration  = t(col_means_by_patient)
)

# --- Find Optimal Clusters ---
optk.breast <- getClustNum(data = omicsList, try.N.clust = 2:8, scale = TRUE, fig.name = "CLUSTER NUMBER")
# we choose 4 clusters
# Extended Data Fig. 9a

# # --- Run Top 3 Best Methods ---
# best.methods <- c("CIMLR", "MoCluster", "iClusterBayes")
# moic.res <- getMOIC(data = omicsList, methodslist = best.methods, N.clust = 4)
# 
# # --- Cluster Sizes ---
# for (method in names(moic.res)) {
#   cat("\n=== Cluster Sizes -", method, "===\n")
#   print(table(moic.res[[method]]$clust.res$clust))
# }
# 
# # --- Kaplan-Meier Survival Plots ---
# for (method in names(moic.res)) {
#   clust <- moic.res[[method]]$clust.res$clust
#   clinical$Cluster <- factor(clust)
#   surv_obj <- Surv(clinical$time, clinical$status)
#   fit <- survfit(surv_obj ~ Cluster, data = clinical)
#   print(ggsurvplot(fit, data = clinical, pval = TRUE, title = paste("Survival -", method)))
#   logrank <- survdiff(surv_obj ~ Cluster, data = clinical)
#   pval <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)
#   cat("Log-rank p-value:", pval, "\n")
# }
# 
# # we decided on using CIMLR

# --- Run CIMLR clustering ---
moic.cimlr <- getMOIC(
  data = omicsList,
  methodslist = list("CIMLR"),
  N.clust = 4
)

# --- View Cluster Assignments ---
table(moic.cimlr$clust.res$clust)

# --- Kaplan-Meier Survival Analysis --- 
clinical$Cluster <- factor(moic.cimlr$clust.res$clust)
surv_obj <- Surv(clinical$time, clinical$status)
fit <- survfit(surv_obj ~ Cluster, data = clinical)

# Plot KM curve
ggsurvplot(fit, data = clinical, pval = TRUE, title = "CIMLR Clusters")

# Log-rank test
logrank <- survdiff(surv_obj ~ Cluster, data = clinical)
pval <- 1 - pchisq(logrank$chisq, df = length(logrank$n) - 1)
cat("Log-rank p-value:", pval, "\n")
# Log-rank p-value: 0.01378699 



#########################################
## now we have 4 subtypes
#########################################

# Recode Cluster numbers to subtype names
annotation_col$Subtype <- factor(
  annotation_col$Cluster,
  levels = c(1, 2, 3, 4),
  labels = c("SF", "EnR", "PI", "HL")
)

# Define subtype colors
ann_colors <- list(
  Subtype = c(
    "SF"  = "#E41A1C",  # red
    "EnR" = "#4DAF4A",  # green
    "PI"  = "#377EB8",  # blue
    "HL"  = "#984EA3"   # purple
  )
)





####################################
## Fig 5a - multimodal integration heatmap
####################################

# --- Standardize Data for Visualization ---
plotdata <- getStdiz(
  data       = omicsList,
  halfwidth  = c(2, 2, 2, 2),
  centerFlag = c(TRUE, TRUE, TRUE, TRUE),
  scaleFlag  = c(TRUE, TRUE, TRUE, TRUE)
)

# --- Extract Important Features ---
feat   <- moic.cimlr$feat.res
feat1  <- feat[feat$dataset == "histo", ][, "feature"]
feat2  <- feat[feat$dataset == "module", ][, "feature"]
feat3  <- feat[feat$dataset == "clusters", ][, "feature"]
# feat4  <- feat[feat$dataset == "cell_infiltration", ][c(1,3,4,5,6,7,8,10,11,12,14,21,22,23,28), "feature"]
feat4  <- feat[feat$dataset == "cell_infiltration", ][, "feature"]
annRow <- list(feat1, feat2, feat3, feat4)

# --- Define Color Palettes ---
color_histo        <- colorRampPalette(brewer.pal(9, "YlGnBu"))(100)
color_module       <- colorRampPalette(c("blue", "white", "red"))(100)
color_clusters     <- colorRampPalette(c("grey90", "black"))(100)
color_infiltration <- colorRampPalette(brewer.pal(9, "RdYlBu"))(100)
col.list <- list(color_histo, color_module, color_clusters, color_infiltration)

# === Rename Clusters to Custom Subtypes ===
cluster_assignment <- moic.cimlr$clust.res$clust
sample_ids <- rownames(moic.cimlr$clust.res)
names(cluster_assignment) <- sample_ids

# Define desired labels
# correct mapping
subtype_labels <- c("1" = "EnR", "2" = "SF", "3" = "PI", "4" = "HL")
subtypes_named <- setNames(subtype_labels[as.character(cluster_assignment)], sample_ids)


# --- Create Column Annotation with Correct Levels ---
annCol <- data.frame(
  Subtype = factor(subtypes_named, levels = c("EnR","SF","PI","HL")),
  row.names = sample_ids
)

# --- Define Subtype Annotation Colors ---
annotation_colors <- list(
  Subtype = c(
    "SF"  = "#E46B71",  # red/pink
    "EnR" = "#8CBF40",  # green
    "PI"  = "#38A7B0",  # teal
    "HL"  = "#8C6BB1"   # purple
  )
)

# --- Sort Samples by Subtype ---
ordered_samples <- rownames(annCol)[order(annCol$Subtype)]
ordered_plotdata <- lapply(plotdata, function(mat) mat[, ordered_samples])

# --- Update Cluster Labels in `clust.res` to Match Subtypes ---
ordered_clust_res <- moic.cimlr$clust.res
ordered_clust_res$clust <- as.character(annCol[ordered_samples, "Subtype"])
annCol <- annCol[ordered_samples, , drop = FALSE]


#### heatmap

# -----------------------------
# 1) Fix subtype order + colors
# -----------------------------
subtype_colors <- c(
  "SF"  = "#E46B71",  # red/pink
  "EnR" = "#8CBF40",  # green
  "PI"  = "#38A7B0",  # teal
  "HL"  = "#8C6BB1"   # purple
)

annCol$Subtype <- factor(
  annCol$Subtype,
  levels = c("EnR", "SF", "PI", "HL")
)

# reorder samples according to subtype
ordered_samples <- rownames(annCol)[order(annCol$Subtype)]
annCol <- annCol[ordered_samples, , drop = FALSE]
ordered_plotdata <- lapply(ordered_plotdata, function(m) m[, ordered_samples])

# -----------------------------
# 2) Row-wise z-scoring
# -----------------------------
row_z <- function(m) t(scale(t(m)))

histo_z    <- row_z(ordered_plotdata$histo)
module_z   <- row_z(ordered_plotdata$module)
clusters_z <- row_z(ordered_plotdata$clusters)
infil_z    <- row_z(ordered_plotdata$cell_infiltration)


# -----------------------------
# 3) Restrict rows (your feature sets)
# -----------------------------
if (!is.null(annRow[[1]])) histo_z    <- histo_z[rownames(histo_z)    %in% annRow[[1]], , drop = FALSE]
if (!is.null(annRow[[2]])) module_z   <- module_z[rownames(module_z)  %in% annRow[[2]], , drop = FALSE]
if (!is.null(annRow[[3]])) clusters_z <- clusters_z[rownames(clusters_z)%in% annRow[[3]], , drop = FALSE]
if (!is.null(annRow[[4]])) infil_z    <- infil_z[rownames(infil_z)    %in% annRow[[4]], , drop = FALSE]


# --- Relabel selected modules with biological annotation ---
rownames(module_z) <- dplyr::recode(
  rownames(module_z),
  "M14 " = "M14: proliferation",
  "M4 "  = "M4: interferon signalling"
)

# -----------------------------
# 4) Color maps
# -----------------------------
col_histo    <- colorRamp2(c(-2, 0, 2), c("#CDECE7", "white", "#1F9D8B"))
col_module   <- colorRamp2(c(-2, 0, 2), c("#2C7BE5", "white", "#D64545"))
col_clusters <- colorRamp2(c(-2, 0, 2), c("#FFFFFF", "#9E9E9E", "#111111"))
col_infil    <- colorRamp2(c(-2, 0, 2), c("#FDE725", "white", "#355F8D"))

# -----------------------------
# 5) Top annotation
# -----------------------------
top_anno <- HeatmapAnnotation(
  "spatial subtype" = annCol$Subtype,              # <- the name shown above the bar
  col = list("spatial subtype" = subtype_colors),  # <- key must match
  annotation_legend_param = list(
    "spatial subtype" = list(title = "spatial subtype")  # legend title
  ),
  simple_anno_size = unit(4, "mm"),
  annotation_name_gp = gpar(fontsize = 10, fontface = "bold"),
  annotation_name_rot = 0
)

col_split <- annCol$Subtype  # enforce subtype split

# -----------------------------
# 6) Base arguments
# -----------------------------
base_args <- list(
  cluster_columns = FALSE,
  column_split = col_split,
  show_column_names = FALSE,
  column_title = NULL,
  use_raster = TRUE,
  raster_quality = 2,
  column_gap = unit(1, "mm")  # thin white gaps between subtype blocks
)

# -----------------------------
# 7) Build heatmap blocks
# -----------------------------
hm_histo <- do.call(Heatmap, c(list(
  histo_z, name = "histo", col = col_histo,
  top_annotation = top_anno,
  cluster_rows = TRUE, clustering_method_rows = "ward.D2",
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 9),
  row_title = "histology", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "histology", at = c(-2,-1,0,1,2)),
  gap = unit(2, "mm")
), base_args))

hm_module <- do.call(Heatmap, c(list(
  module_z, name = "modules", col = col_module,
  cluster_rows = TRUE, clustering_method_rows = "ward.D2",
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 9),
  row_title = "modules", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "modules", at = c(-2,-1,0,1,2)),
  gap = unit(2, "mm")
), base_args))

hm_clusters <- do.call(Heatmap, c(list(
  clusters_z, name = "clusters", col = col_clusters,
  cluster_rows = TRUE, clustering_method_rows = "ward.D2",
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 9),
  row_title = "spatial clusters", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "spatial clusters", at = c(-2,-1,0,1,2)),
  gap = unit(2, "mm")
), base_args))

hm_infil <- do.call(Heatmap, c(list(
  infil_z, name = "cell infiltration", col = col_infil,
  cluster_rows = TRUE, clustering_method_rows = "ward.D2",
  show_row_names = TRUE, row_names_gp = gpar(fontsize = 9),
  row_title = "cell type deconvolution", row_title_gp = gpar(fontsize = 12, fontface = "bold"),
  heatmap_legend_param = list(title = "cell type deconvolution", at = c(-2,-1,0,1,2))
), base_args))

# -----------------------------
# 8) Combine + save
# -----------------------------

ht <- hm_histo %v% hm_module %v% hm_clusters %v% hm_infil
ht
# pdf("HEATMAP_CIMLR_FINAL_SUBTYPES_complexheatmap.pdf", width = 12, height = 13, useDingbats = FALSE)
# draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right",
#      padding = unit(c(4, 8, 4, 4), "mm"))
# dev.off()





####################################
## Fig 5b - Hallmark pathway activity across Spatial4HR+ subtypes
####################################

library(GSVA)
library(dplyr)
library(qusage)
library(survival)
library(survminer)

ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")
ductal_pb_norm <- readRDS("~/Desktop/final_scripts/00_data/ductal_pb_rpm_log.RDS")

ductal_meta$ductal_subtype <- factor(ductal_meta$ductal_subtype,
                                     levels = c("EnR", "SF", "PI", "HL"))

table(ductal_meta$ductal_subtype)
# EnR  SF  PI  HL 
# 23  14  33  16

genesets_hallmark <- read.gmt("~/Desktop/final_scripts/00_data/h.all.v2023.2.Hs.symbols.gmt")
genesets_hallmark <- genesets_hallmark[-c(38,43,47,48)]

# We excluded the following 4 hallmark signatures due to a
# lack of association with tumor processes or microenvironment in BC:
# pancreas beta cells, spermatogenesis, UV response up, UV response down.

# Create the GSVAParams object for the "gsva" method
params <- gsvaParam(as.matrix(ductal_pb_norm),  
                    kcdf="Gaussian", geneSets = genesets_hallmark)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)

table(rownames(module_sig) == ductal_meta$name)

ductal_meta <- cbind(ductal_meta, module_sig)

module_sig <- as.data.frame(module_sig)

module_sig$Subtypes = ductal_meta$ductal_subtype

data_plot <- as.data.frame(cbind(module_sig,
                                 matrix(rep(module_sig$Subtypes, length(table(module_sig$Subtypes)))
                                        ,ncol=length(table(module_sig$Subtypes)))))

data_plot[,c(47)] = NULL

colnames(data_plot)[47:50] <- c("EnR", "SF", "PI", "HL")
data_plot[,"EnR"] <- as.factor(ifelse(data_plot[,"EnR"]=="EnR", 1, 0))
data_plot[,"SF"] <- as.factor(ifelse(data_plot[,"SF"]=="SF", 1, 0))
data_plot[,"PI"] <- as.factor(ifelse(data_plot[,"PI"]=="PI", 1, 0))
data_plot[,"HL"] <- as.factor(ifelse(data_plot[,"HL"]=="HL", 1, 0))

colnames(data_plot)
data_plot[,1:46] = apply(data_plot[,1:46], 2, genefu::rescale, na.rm = T)

subtype = names(table(module_sig$Subtypes)) ## la colonna con i cluster, o i nomi dei cluster direttamente

effect <- matrix(0, nrow=length(colnames(data_plot[,1:46])), ncol=length( sort(unique(subtype)) ))
colnames(effect) <- sort(unique(subtype))
rownames(effect) <- colnames(data_plot[,1:46])

pvalue <- matrix(0, nrow=length(colnames(data_plot[,1:46])), ncol=length( sort(unique(subtype)) ))
colnames(pvalue) <- sort(unique(subtype))
rownames(pvalue) <- colnames(data_plot[,1:46])

ci <- matrix(0, nrow=length(colnames(data_plot[,1:46])), ncol=length( sort(unique(subtype)) ))
colnames(ci) <- sort(unique(subtype))
rownames(ci) <- colnames(data_plot[,1:46])

signatures = data_plot[,1:46]

colnames(data_plot) = gsub("-", "", colnames(data_plot))
colnames(data_plot) = gsub("\\+", "", colnames(data_plot))
colnames(data_plot) = gsub(" ", "_", colnames(data_plot))

colnames(signatures) = gsub("-", "", colnames(signatures))
colnames(signatures) = gsub("\\+", "", colnames(signatures))
colnames(signatures) = gsub(" ", "_", colnames(signatures))

for(j in 1:length(colnames(signatures)) ){
  print(paste("Signature n.", j))
  id <- sort(unique(subtype))
  eff <- NULL
  pval <- NULL
  low <- NULL
  up <- NULL
  
  for(i in 1:length(id)){
    print(id[i])
    data_plot[,colnames(signatures)[j]] <- as.numeric(as.character(data_plot[,colnames(signatures)[j]]))
    formula    <- as.formula(paste(id[i], " ~ ", colnames(signatures)[j]))
    res.logist <- glm(formula, data = data_plot, family=binomial)
    
    p = wilcox.test(as.formula(paste(colnames(signatures)[j], " ~ ", id[i])), data = data_plot)$p.value
    
    summary(res.logist)
    eff <- c(eff,round(exp(coef(res.logist))[2],2))
    pval <- c(pval,p)
    low <- c(low,round(exp(confint(res.logist,level=.95))[2,1],2))
    up <- c(up,round(exp(confint(res.logist,level=.95))[2,2],2))
  }
  names(pval) <- names(eff) <- id
  padjust <- p.adjust(pval,method = "fdr")
  
  effect[j,] <- eff
  pvalue[j,] <- pval
  ci[j,] <- paste(low,"-",up,sep="")
}
padjust <- matrix( p.adjust(pvalue,method="fdr"), ncol=ncol(pvalue), nrow=nrow(pvalue), dimnames=dimnames(pvalue))


effect_orig = effect
effect
effect[padjust>0.25] <- NA 
effect[effect<=.1] <- .1
effect[effect>=8 ] <- 10

effect <- t(effect)
e_ductal = effect

e_ductal <- e_ductal[c("EnR", "SF", "PI", "HL"),]
dim(e_ductal)

colnames(e_ductal) = gsub("_", " ", colnames(e_ductal))
colnames(e_ductal) = gsub("HALLMARK ", "", colnames(e_ductal))


library(gtools)
current_row_names <- rownames(e_ductal)
sorted_row_names <- current_row_names
e_ductal_sorted <- e_ductal[sorted_row_names, ]
rownames(e_ductal_sorted) <- sorted_row_names
print(rownames(e_ductal_sorted))
e_ductal = e_ductal_sorted

plot_legend_1 <- matrix(c("0.1",
                          "10"
) , nrow = 1, ncol = 2)
class(plot_legend_1) <- "numeric"
colnames(plot_legend_1) <- c(paste0("depletion"),
                             paste0("enrichment"))
rownames(plot_legend_1) <- "Direction"

plot_legend_2 <- matrix(c("0.4",
                          "0.2"
) , nrow = 1, ncol = 2)
class(plot_legend_2) <- "numeric"
colnames(plot_legend_2) <- c(paste0("large"),
                             paste0("small"))
rownames(plot_legend_2) <- "Effect size"

library(Cairo)
par(mar=c(3, 1, 6, 1))
add_space <- function(x) {x <- gsub("_", " ", x);x}
add_sign_pos <- function(x) {x <- gsub("_pos", "+", x);x}
add_sign_neg <- function(x) {x <- gsub("_neg", "-", x);x}

circle <- function(x,y,r,nsteps=100,...){  
  rs <- seq(0,2*pi,len=180)
  xc <- x + r * cos(rs+pi/2)
  yc <- y + r * sin(rs+pi/2)
  polygon(xc,yc,...)
}
plot(
  row(e_ductal)*1,
  col(e_ductal)*1,
  type="n",xlab="", ylab="",
  xlim=c(0.5,nrow(e_ductal)*1+0.5),
  ylim=c(0.5,ncol(e_ductal)*1+1),
  axes=FALSE,
  ann=FALSE, asp=1)

for(i in 1:nrow(e_ductal)){
  for(j in 1:ncol(e_ductal)){
    circle(i*1,j*1,abs(log2(e_ductal[i,j]))*0.13, # --> 0.05 is the size of the cirlces
           col=adjustcolor(ifelse(e_ductal[i,j]>1,"#009E73","#D55E00"), alpha.f = 0.8),lwd = .7, linecolor = "white", border = NA)
  }
}


axis(1,at=1:nrow(e_ductal)*1,labels=FALSE, pos = 0) # --> *1 is the distance between things on x axis
text(x=rep(.2,ncol(e_ductal)*1), y=1:ncol(e_ductal)*1, labels=add_space(add_sign_neg(add_sign_pos(colnames(e_ductal)))), cex = 0.9, cex.axis=1.5, las=2, srt = 0,pos = 2, xpd = T)
axis(2,at=1:ncol(e_ductal)*1,pos=0.5,labels=FALSE)
text(x=1:nrow(e_ductal)*1, y=rep(-1,nrow(e_ductal)*1), labels=add_space(add_sign_neg(add_sign_pos(rownames(e_ductal)))), cex.axis=1.5, srt = 0, las=2, xpd = T, srt = 90, adj = 1)


# Function to draw circles
circle <- function(x, y, r, nsteps=100, ...){  
  rs <- seq(0, 2 * pi, len = nsteps)
  xc <- x + r * cos(rs + pi/2)
  yc <- y + r * sin(rs + pi/2)
  polygon(xc, yc, ...)
}

# Legend 1
for(i in 1:nrow(plot_legend_1)){
  for(j in 1:ncol(plot_legend_1)){
    # Ensure circle is visible by adjusting coordinates and radius
    circle_x <- i * 9  # Adjusted to ensure it doesn't fall out of range
    circle_y <- j + 20
    circle(circle_x, circle_y, abs(log2(plot_legend_1[i, j])) * 0.1,
           col = adjustcolor(ifelse(plot_legend_1[i, j] > 1, "#009E73", "#D55E00"), alpha.f = 0.8),
           lwd = 0.6, border = NA, xpd = TRUE)
  }
}

text(x = rep(9, ncol(plot_legend_1) * 1), y = c(21, 22),  
     labels = add_space(add_sign_neg(add_sign_pos(colnames(plot_legend_1)))),
     cex.axis = 2, las = 2, srt = 0, pos = 4, xpd = TRUE)

text(x = rep(9.5, ncol(plot_legend_1) * 1), y = 22.5,
     labels = paste0("Direction"), cex = 1, las = 2, srt = 0, pos = 3, xpd = TRUE)

# Legend 2
for(i in 1:nrow(plot_legend_2)){
  for(j in 1:ncol(plot_legend_2)){
    circle_x <- i * 9  # Adjusted to ensure it doesn't fall out of range
    circle_y <- j + 15
    circle(circle_x, circle_y, plot_legend_2[i, j],
           col = adjustcolor("grey", alpha.f = 0.8), lwd = 0.6, border = NA, xpd = TRUE)
  }
}

text(x = rep(9, ncol(plot_legend_2) * 1), y = c(16, 17),  
     labels = add_space(add_sign_neg(add_sign_pos(colnames(plot_legend_2)))), cex.axis = 2, las = 2, srt = 0, pos = 4, xpd = TRUE)

text(x = rep(10, ncol(plot_legend_2) * 1), y = 17.5,
     labels = paste0("Effect size"), cex = 1, las = 2, srt = 0, pos = 3, xpd = TRUE)




####################################
## Fig 5c,d - ESR1 and ERBB2 expression across subtypes
####################################

library(ggplot2)
library(ggpubr)
library(dplyr)
library(tidyverse)

ductal_pb_norm <- readRDS("~/Desktop/final_scripts/00_data/ductal_pb_rpm_log.RDS")
ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")


####################################
## Fig 5d - ERBB2 expression across Spatial4HR+ subtypes
####################################

ERBB2_expr <- ductal_pb_norm["ERBB2", ]
ductal_meta$HER2_Score <- ERBB2_expr[match(ductal_meta$name, colnames(ductal_pb_norm))]

ductal_meta$ductal_subtype <- as.factor(ductal_meta$ductal_subtype)

ductal_subtype_comparisons <- combn(levels(ductal_meta$ductal_subtype), 2, simplify = FALSE)

subtype_colors <- c("SF"="#E46B71","EnR"="#8CBF40","PI" = "#38A7B0","HL"="#8C6BB1")

ductal_meta$ductal_subtype <- factor(ductal_meta$ductal_subtype,
                                     levels = c("EnR","SF","PI","HL"))

all_comps <- combn(levels(ductal_meta$ductal_subtype), 2, simplify = FALSE)

p_table <- compare_means(HER2_Score ~ ductal_subtype, data = ductal_meta,
                         method = "wilcox.test", p.adjust.method = "BH",
                         ref.group = NULL, paired = FALSE) %>%
  filter(group1 %in% levels(ductal_meta$ductal_subtype),
         group2 %in% levels(ductal_meta$ductal_subtype)) %>%
  # keep only significant (adjust if you want <= 0.1 etc.)
  filter(p.adj <= 0.05) %>%
  mutate(y.position = max(ductal_meta$HER2_Score, na.rm = TRUE) + 
           seq(0.2, 0.2* n(), by = 0.2),
         p.signif = case_when(
           p.adj < 0.001 ~ "***",
           p.adj < 0.01  ~ "**",
           p.adj < 0.05  ~ "*",
           TRUE ~ "ns"
         )) %>%
  dplyr::select(group1, group2, p.adj, p.signif, y.position) %>%
  rename(p.adj.signif = p.signif)

n_df <- ductal_meta %>% count(ductal_subtype)

p <- ggplot(ductal_meta, aes(x = ductal_subtype, y = HER2_Score, fill = ductal_subtype)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, size = 0.8, color = "black") +
  geom_jitter(aes(fill = ductal_subtype), shape = 21, color = "black",
              width = 0.15, size = 2.4, alpha = 0.85) +
  stat_pvalue_manual(p_table,
                     label = "p.adj.signif",
                     tip.length = 0.01, size = 5,
                     inherit.aes = FALSE) +      # <-- key fix  
  scale_fill_manual(values = subtype_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  labs(x = "Spatial Subtype", y = "ERBB2 Gene Expression") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x  = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 13),
    axis.title   = element_text(size = 16),
    panel.border = element_rect(color = "black", size = 0.6, fill = NA),
    legend.position = "none"
  )


p



####################################
## Fig 5c - ESR1 expression across Spatial4HR+ subtypes
####################################

ESR1_expr <- ductal_pb_norm["ESR1", ]
ductal_meta$ER_Score <- ESR1_expr[match(ductal_meta$name, colnames(ductal_pb_norm))]

ductal_meta$ductal_subtype <- as.factor(ductal_meta$ductal_subtype)

ductal_subtype_comparisons <- combn(levels(ductal_meta$ductal_subtype), 2, simplify = FALSE)

subtype_colors <- c("SF"="#E46B71","EnR"="#8CBF40","PI" = "#38A7B0","HL"="#8C6BB1")

ductal_meta$ductal_subtype <- factor(ductal_meta$ductal_subtype,
                                     levels = c("EnR","SF","PI","HL"))

all_comps <- combn(levels(ductal_meta$ductal_subtype), 2, simplify = FALSE)

p_table <- compare_means(ER_Score ~ ductal_subtype, data = ductal_meta,
                         method = "wilcox.test", p.adjust.method = "BH",
                         ref.group = NULL, paired = FALSE) %>%
  filter(group1 %in% levels(ductal_meta$ductal_subtype),
         group2 %in% levels(ductal_meta$ductal_subtype)) %>%
  # keep only significant (adjust if you want <= 0.1 etc.)
  filter(p.adj <= 0.05) %>%
  mutate(y.position = max(ductal_meta$ER_Score, na.rm = TRUE) + 
           seq(0.2, 0.2* n(), by = 0.2),
         p.signif = case_when(
           p.adj < 0.001 ~ "***",
           p.adj < 0.01  ~ "**",
           p.adj < 0.05  ~ "*",
           TRUE ~ "ns"
         )) %>%
  select(group1, group2, p.adj, p.signif, y.position) %>%
  rename(p.adj.signif = p.signif)

n_df <- ductal_meta %>% count(ductal_subtype)

p <- ggplot(ductal_meta, aes(x = ductal_subtype, y = ER_Score, fill = ductal_subtype)) +
  geom_boxplot(width = 0.58, outlier.shape = NA, size = 0.8, color = "black") +
  geom_jitter(aes(fill = ductal_subtype), shape = 21, color = "black",
              width = 0.15, size = 2.4, alpha = 0.85) +
  stat_pvalue_manual(p_table,
                     label = "p.adj.signif",
                     tip.length = 0.01, size = 5,
                     inherit.aes = FALSE) +      # <-- key fix  
  scale_fill_manual(values = subtype_colors) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  labs(x = "Spatial Subtype", y = "ESR1 Gene Expression") +
  theme_classic(base_size = 16) +
  theme(
    axis.text.x  = element_text(size = 14, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 13),
    axis.title   = element_text(size = 16),
    panel.border = element_rect(color = "black", size = 0.6, fill = NA),
    legend.position = "none"
  )

p



####################################
## Fig 5e - multivariate cox model
####################################

library(readxl)
library(dplyr)
library(survival)
library(survminer)
library(genefu)
library(org.Hs.eg.db)

ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")
ductal_pb_norm <- readRDS("~/Desktop/final_scripts/00_data/ductal_pb_rpm_log.RDS")

ductal_meta$ductal_subtype <- factor(ductal_meta$ductal_subtype,
                                     levels = c("EnR", "SF", "PI", "HL"))

head(ductal_meta)
table(ductal_meta$ductal_subtype)

clean_ki67 <- function(x) {
  # Extract the first number found in the string
  match <- regmatches(x, regexpr("\\d+(\\.\\d+)?", x))
  if (length(match) == 0 || is.na(match)) return(NA_real_)
  as.numeric(match)
}

ductal_meta$KI67_clean <- sapply(ductal_meta$KI67, clean_ki67)

hs <- org.Hs.eg.db
my.symbols <- rownames(ductal_pb_norm)
annotation = AnnotationDbi::select(hs,
                                   keys = my.symbols,
                                   columns = c("ENTREZID", "SYMBOL"),
                                   keytype = "SYMBOL")
colnames(annotation) = c("name", "EntrezGene.ID")
annotation = na.omit(annotation)
annotation = annotation[!duplicated(annotation$name), ]
rownames(annotation) = annotation$name
data("sig.ggi")

ggi_calc = ggi(data = t(as.matrix(ductal_pb_norm)), annot = annotation, do.mapping = T)
score_ggi_genefu = ggi_calc$score

table(names(score_ggi_genefu) == ductal_meta$name)
# all true

cox_data <- ductal_meta %>%
  dplyr::select(name, ductal_subtype, status, time, AGE, GRADE, 
                KI67_clean, N_STATUS, SIZE_MM, ADJUVANT_CHEMO) %>%
  mutate(
    subtype = factor(ductal_subtype, levels = c("EnR", "SF", "PI", "HL")),  # Set EnR as reference
    age = as.numeric(AGE),
    grade = as.factor(GRADE),
    ki67 = KI67_clean,
    nodal_status = as.factor(N_STATUS),
    size_mm = as.numeric(SIZE_MM),
    adjuvant_chemo = as.factor(ADJUVANT_CHEMO) 
  ) %>%
  filter(!is.na(subtype) & !is.na(status) & !is.na(time) &
           !is.na(age) & !is.na(grade) & !is.na(ki67) & !is.na(nodal_status) & !is.na(size_mm))

cox_data <- cbind(cox_data, score_ggi_genefu)

head(cox_data)

cox_model <- coxph(Surv(time, status) ~ subtype + age + grade + ki67 + 
                     nodal_status + size_mm + score_ggi_genefu + adjuvant_chemo, data = cox_data)

cox.zph(cox_model)
summary(cox_model)

# ============================================================
# Forest plot for multivariable Cox model
# ============================================================

basicForest = function(x, a, adj=NULL, xlim=NULL, xlog=FALSE, xlab="", col='black', cex.axis=.7,
                       titles=names(x), annotDir=NULL, lineHeight=1.5, cex.annot=0.7, colWidth=NULL, pch=16, beforePlot=NULL)
{ if (is.data.frame(x)) { nr = nrow(x); } else { nr = length(x[[1]]); }
  if (nr != nrow(a)) { stop("Not same number of rows in x and a"); }
  if (is.null(adj)) { adj = c("left", rep("right", length(x)-1)) }
  if (length(col)==1) { col = rep(col, nr); }
  
  plot.new(); plot.window(xlim=c(0,1), ylim=c(1, 0), xaxs='i', yaxs='i') 
  
  if (is.null(colWidth))  
  { cw = colMaxs(sapply(x, function(i) sapply(i, strwidth)));
  cw = pmax(cw, strwidth(names(x)));
  cw = cw + strwidth("M");
  }
  else { cw = colWidth[-length(colWidth)]/sum(colWidth); } 
  
  ch = abs(strheight("!"));
  if (is.numeric(lineHeight))
  { if (length(lineHeight)==1) { li = (1:(1+nr))*ch*lineHeight; }
    else { li = lineHeight*ch; }
  }
  else
  { if (lineHeight=="full")
  { li = seq(ch, 1-7*ch, len=nr+1); }
  }
  
  if (max(li)>1) { warning("Too many lines, does not fit in window"); }
  st = cumsum(c(0,cw));
  left = st[-length(st)]; right=st[-1]; center=(left+right)/2;
  cpos = left; cpos[adj=="right"]=right[adj=="right"]; cpos[adj=="center"]=center[adj=="center"];
  
  li2 = li-strheight("M")*.4;
  ad = c(left=0, right=1, center=.5)[adj]
  for (i in seq_along(x))
  { for (j in seq_along(x[[i]]))
  { text(labels=x[[i]][[j]], x=cpos[i], y=li2[j+1], adj=c(ad[i], 0)); }
    #rep(cpos[i], each=nr), y=li[-1], adj=ad[i]) }
    text(labels=titles[[i]], x=cpos[i], y=li2[1], font=2, adj=c(ad[i], 0))
  }
  
  if (is.null(xlim)) { xlim=range(a, na.rm=TRUE) };
  xx = grconvertX(c(st[length(st)]+1.5*strwidth("M"), 1), from='user', to='nfc')
  yy = grconvertY(c(li[length(li)],li[2]), from='user', to='nfc')
  old.par <- par( c('plt', 'usr', 'mgp', "cex.axis") )
  on.exit(par(old.par))
  
  par(new=TRUE, plt=c(xx,yy), mgp=c(1.5,.5,0), cex.axis=cex.axis);
  xl = range(a, na.rm=TRUE);
  plot.new(); plot.window(xlim=xlim, ylim=c(1,0), yaxs='i', xaxs='i', log=ifelse(xlog, "x", ""))
  if (!is.null(beforePlot)) { beforePlot(); }
  axis(1, line=1); title(xlab=xlab, line=2.5, cex.lab=cex.axis);
  lines(c(xlog, xlog), c(1-strheight("M")*1.5, 0+strheight("M")), col='lightgrey', xpd=TRUE);
  
  wp = which(!is.na(a[,2]));
  ty=rep(0, nrow(a));
  ty[which(a[,1]<xlim[1])]=ty[which(a[,1]<xlim[1])]+1;
  ty[which(a[,3]>xlim[2])]=ty[which(a[,3]>xlim[2])]+2;
  
  le=.05;
  w = which(a[,2]>=xlim[1] & a[,2]<=xlim[2]);
  a[which(a[,2]<xlim[1]),2] = (xlim[1]+a[which(a[,2]<xlim[1]),3])/2
  a[which(a[,2]>xlim[2]),2] = (xlim[2]+a[which(a[,2]>xlim[2]),1])/2
  
  li2 = li[-1]-li[2]; li2=li2/li2[length(li2)];
  
  for (i in wp)
  { suppressWarnings(arrows(a[i,2], li2[i], min(xlim[2], a[i,3]), li2[i], angle=c(30,90)[(a[i,3]<=xlim[2])+1],
                            length=le, col=col[i], xpd=NA, lwd=2*par('cex')))
    suppressWarnings(arrows(a[i,2], li2[i], max(xlim[1], a[i,1]), li2[i], angle=c(30,90)[(a[i,1]>=xlim[1])+1],
                            length=le, col=col[i], xpd=NA, lwd=2*par('cex')))
  }
  points(a[w,2], li2[w], pch=pch, col=col[w], xpd=NA, cex=2, family="Arial Unicode MS")
  
  if (!is.null(annotDir))
  { h = 1-strheight("M", cex=cex.axis)*5.5;
  v = strwidth(xlab, cex=cex.axis)*2;
  med = mean(log(xlim));
  arrows(exp(med-v), h, xlim[1], h, xpd=TRUE, length=le)
  arrows(exp(med+v), h, xlim[2], h, xpd=TRUE, length=le)
  text(annotDir, x= sqrt(c(exp(med-v), exp(med+v))*xlim), y=h-strheight("M"), xpd=NA, cex=cex.annot);
  }
  return(invisible(list(figInfo=list(xx=xx, yy=yy, xlim=xlim, ylim=c(1,0), li=li2), colWidth=cw,
                        linePos=li)))
}


# --- helper to print p-values nicely
formatNice <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("p<0.001")
  sprintf("p=%.3f", p)
}

# 1) Pull numbers from your fitted model
co <- summary(cox_model)$coefficients
ci <- summary(cox_model)$conf.int
# co rows = variables; columns include: "coef", "exp(coef)", "se(coef)", "Pr(>|z|)"
# ci columns include: "exp(coef)", "lower .95", "upper .95"

# 2) Build pretty labels for variables
pretty_names <- rownames(co)
pretty_names <- sub("^subtype",     "Subtype: ",      pretty_names)
pretty_names <- sub("^age",     "Age ",      pretty_names)
pretty_names <- sub("^grade",       "Grade: ",        pretty_names)
pretty_names <- sub("^nodal_status","Nodal: ",        pretty_names)
pretty_names <- sub("^size_mm$",    "Tumor size (mm)",pretty_names)
pretty_names <- sub("^ki67_10$",    "Ki67 (per 10%)", pretty_names)
pretty_names <- sub("^ki67$",       "Ki67 (%)",       pretty_names)
pretty_names <- sub("^score_ggi_genefu",       "GGI score",       pretty_names)
pretty_names <- sub("adjuvant_chemo1",       "Adjuvant Chemo",       pretty_names)

pretty_names

# 3) Text columns for the forest (left side)
r <- list(
  pretty_names,                             # Column 1: variable labels
  vapply(co[, "Pr(>|z|)"], formatNice, "")  # Column 2: p-values
)

# 4) HRs and 95% CIs (right side)
m <- cbind(
  ci[, "lower .95"],
  ci[, "exp(coef)"],
  ci[, "upper .95"]
)

# 5) Colors by significance (optional)
cols <- ifelse(co[, "Pr(>|z|)"] < 0.05, "#00AFBB", "slategray4")

# 6) Choose reasonable x-limits for HR axis
hr_range <- range(m, na.rm = TRUE)
xlim_hr <- c(max(0.4, hr_range[1]*0.9), min(5, hr_range[2]*1.1))

# 7) Plot using the provided basicForest()
## --- p-value formatting: plain numbers (sci for small) ---
fmt_p_plain <- function(p){
  if (is.na(p)) return(NA)
  if (p < 1e-3){
    s <- format(p, digits = 1, scientific = TRUE)        # e.g. "5e-04"
    s <- sub("e\\+?(-?\\d+)$", "\u00D710^{\\1}", s)      # "5×10^{-4}"
    return(parse(text = s))                               # pretty math text
  }
  sprintf("%.3f", p)
}

## --- left text columns ---
var_labels <- pretty_names
p_vals     <- vapply(co[, "Pr(>|z|)"], fmt_p_plain, character(1))
spacer_col <- rep("", length(var_labels))                 # thin blank column

## --- HR matrix (lower, HR, upper) ---
m <- cbind(ci[, "lower .95"], ci[, "exp(coef)"], ci[, "upper .95"])

## --- colors and limits (as you like) ---
cols    <- ifelse(co[, "Pr(>|z|)"] < 0.05, "#00AFBB", "slategray4")
xlim_hr <- range(m, na.rm = TRUE)

## --- optional: vertical 1.0 line & (no band here) ---
bf <- function(){
  abline(v = 1, col = "grey80", lwd = 1)
}


## --- columns: spacer first, then variable, then p ---
xcols <- list(
  rep("", length(var_labels)),        # thin spacer column
  var_labels,                         # variable names
  p_vals                              # p-values
)

## --- headers (note the blank for spacer) ---
headers <- c("", expression(bold("Variable")), expression(italic("p")))

par(mar = c(.5, .5, 2, 1))
basicForest(
  x = xcols,
  a = m,
  xlog = TRUE,
  xlab = "Hazard ratio",
  xlim = xlim_hr,
  titles = headers,
  annotDir = c("Better", "Worse"),
  col = cols,
  pch = 16,
  colWidth = c(.05, .70, .25, 0.5),
  cex.axis = 0.9,
  cex.annot = 0.9,
  beforePlot = bf
)




####################################
## Extended Data Fig. 9b - PAM50 alluvial plot
####################################

library(dplyr)
library(ggplot2)
library(ggalluvial)
library(scales)

# ensure factors and color palette
ductal_meta$ductal_subtype <- factor(
  ductal_meta$ductal_subtype,
  levels = c("EnR","SF","PI","HL")
)

ductal_meta$pam50 <- factor(
  ductal_meta$pam50,
  levels = c("Normal","LumA","LumB","Her2")
)

subtype_colors <- c("SF"="#E46B71","EnR"="#8CBF40","PI"="#38A7B0","HL"="#8C6BB1")

# summarize counts
flow_df <- ductal_meta %>%
  filter(!is.na(pam50), !is.na(ductal_subtype)) %>%
  dplyr::count(pam50, ductal_subtype, name = "n")

# clean and minimalist alluvial plot
p <- ggplot(flow_df,
            aes(y = n,
                axis1 = pam50,
                axis2 = ductal_subtype,
                fill  = ductal_subtype)) +
  geom_alluvium(width = 0, alpha = 0.9) +
  geom_stratum(width = 0.14, fill = "grey92", color = "grey50") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 3.2, color = "black", fontface = "plain") +
  scale_x_discrete(limits = c("PAM50", "Spatial subtype"),
                   expand = c(.05, .05)) +
  scale_fill_manual(values = subtype_colors, name = "Spatial subtype") +
  labs(x = NULL, y = "Samples") +
  theme_minimal(base_size = 11) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 10, face = "plain"),
    axis.title.y = element_text(size = 10),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 9),
    legend.key.size = unit(0.4, "cm"),
    plot.title = element_text(size = 11, face = "bold", hjust = 0),
    legend.position = "right"
  )

p

chisq_flow <- chisq.test(table(ductal_meta$pam50,
                               ductal_meta$ductal_subtype))
chisq_flow

# Pearson's Chi-squared test
# 
# data:  table(ductal_meta$pam50, ductal_meta$ductal_subtype)
# X-squared = 57.852, df = 9, p-value = 3.47e-09






####################################
## Extended Data Fig. 9c - Age
####################################

library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)

age_df <- ductal_meta %>% filter(!is.na(AGE), !is.na(ductal_subtype))

# Kruskal–Wallis p-value
kw_age <- kruskal.test(AGE ~ ductal_subtype, data = age_df)$p.value

# Publication boxplot (no title; larger fonts)
p_age <- ggplot(age_df, aes(x = ductal_subtype, y = AGE)) +
  geom_boxplot(aes(fill = ductal_subtype),
               width = 0.7, outlier.shape = NA,
               color = "grey25", linewidth = 0.3) +
  geom_jitter(aes(color = ductal_subtype),
              width = 0.12, size = 1.6, alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = subtype_colors, guide = "none") +
  scale_color_manual(values = subtype_colors, guide = "none") +
  labs(
    x = "Spatial Subtype",
    y = "Age at diagnosis",
    subtitle = paste0("Kruskal–Wallis p = ", formatC(kw_age, format = "e", digits = 2))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    axis.title = element_text(size = 13),
    axis.text  = element_text(size = 12, color = "black"),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 13, hjust = 0)
  )

p_age




####################################
## Extended Data Fig. 9d - Tumor size (mm) 
####################################

size_df <- ductal_meta %>% filter(!is.na(SIZE_MM), !is.na(ductal_subtype))

# Kruskal–Wallis p-value
kw_size <- kruskal.test(SIZE_MM ~ ductal_subtype, data = size_df)$p.value

# Publication boxplot (no main title; large fonts; minimal ink)
p_size <- ggplot(size_df, aes(x = ductal_subtype, y = SIZE_MM)) +
  geom_boxplot(aes(fill = ductal_subtype),
               width = 0.7, outlier.shape = NA,
               color = "grey25", linewidth = 0.3) +
  geom_jitter(aes(color = ductal_subtype),
              width = 0.12, size = 1.6, alpha = 0.7, show.legend = FALSE) +
  scale_fill_manual(values = subtype_colors, guide = "none") +
  scale_color_manual(values = subtype_colors, guide = "none") +
  labs(
    x = "Spatial Subtype",
    y = "Tumor size (mm)",
    subtitle = paste0("Kruskal–Wallis p = ", formatC(kw_size, format = "e", digits = 2))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    axis.title = element_text(size = 13),
    axis.text  = element_text(size = 12, color = "black"),
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 13, hjust = 0)
  )

p_size



####################################
# Extended Data Fig. 9e - Nodal status
####################################

ductal_meta$N_STATUS <- factor(ductal_meta$N_STATUS, levels = c("N0", "N+"))

# ---- compute Fisher’s exact p-value ----
contingency <- table(ductal_meta$ductal_subtype, ductal_meta$N_STATUS)
fisher_p <- fisher.test(contingency)$p.value

# ---- summarize for plotting ----
plot_df <- ductal_meta %>%
  filter(!is.na(ductal_subtype), !is.na(N_STATUS)) %>%
  count(ductal_subtype, N_STATUS, name = "n") %>%
  group_by(ductal_subtype) %>%
  mutate(total = sum(n),
         prop = n / total) %>%
  ungroup()

totals_df <- plot_df %>%
  distinct(ductal_subtype, total)

# ---- colors ----
nodal_cols <- c("N0" = "#C7C7C7", "N+" = "#009E73")

# ---- Nature-style stacked bar ----
p <- ggplot(plot_df, aes(x = ductal_subtype, y = prop, fill = N_STATUS)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.3) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = nodal_cols, name = "Nodal Status") +
  labs(
    x = "Spatial Subtype",
    y = "Percentage",
    subtitle = paste0("p = ",
                      formatC(fisher_p, format = "e", digits = 2))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.45, "cm"),
    legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 13, hjust = 0)
  )

p



####################################
## Extended Data Fig. 9f - Grade
####################################

ductal_meta$GRADE <- factor(ductal_meta$GRADE, levels = c("1","2","3"))

# Fisher’s exact p-value (4x3 table)
grade_tab <- table(ductal_meta$ductal_subtype, ductal_meta$GRADE)
fisher_p_grade <- fisher.test(grade_tab)$p.value

# summarize for plotting
grade_df <- ductal_meta %>%
  filter(!is.na(ductal_subtype), !is.na(GRADE)) %>%
  count(ductal_subtype, GRADE, name = "n") %>%
  group_by(ductal_subtype) %>%
  mutate(total = sum(n), prop = n / total) %>%
  ungroup()

## ----- palette -----
pal <- c("1" = "#1B9E77", 
         "2" = "#D95F02", 
         "3" = "#7570B3")

# plot
p_grade <- ggplot(grade_df, aes(x = ductal_subtype, y = prop, fill = GRADE)) +
  geom_col(width = 0.75, color = "white", linewidth = 0.3) +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.08))) +
  scale_fill_manual(values = pal, name = "Grade") +
  labs(
    x = "Spatial Subtype",
    y = "Percentage",
    subtitle = paste0("p = ", formatC(fisher_p_grade, format = "e", digits = 2))
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.3),
    axis.title = element_text(size = 13),
    axis.text = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    legend.key.size = unit(0.45, "cm"),
    legend.position = "right",
    plot.title = element_blank(),
    plot.subtitle = element_text(size = 13, hjust = 0)
  )

p_grade





############################################################
## EnR vs SF — Tumor and Stroma Pseudobulk Analysis
############################################################

############################################################
## Extended Data Fig. 9g 
## Tumor pseudobulk: no enrichment of cell-cycle programs
############################################################

library(dplyr)
library(stringr)
library(GSEABase)
library(GSVA)
library(limma)
library(pheatmap)
library(tibble)

## 1) Load data ---------------------------------------------
tumor_pb  <- readRDS("~/Desktop/final_scripts/00_data/50_percent_tumor_pb.RDS")
stroma_pb <- readRDS("~/Desktop/final_scripts/00_data/100_percent_stroma_pb.RDS")
ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")

## 2) Spatial subtype annotation ----------------------------
meta <- ductal_meta %>%
  mutate(
    sample = paste0("ST", orig.ident),
    spatial_subtype = recode(
      as.character(ductal_subtype),
      "1" = "SF",
      "2" = "EnR",
      "3" = "PI",
      "4" = "HL"
    ),
    spatial_subtype = factor(
      spatial_subtype,
      levels = c("EnR", "SF", "PI", "HL")
    )
  ) %>%
  dplyr::select(sample, spatial_subtype)

## Keep EnR vs SF only
meta_ES <- meta %>% filter(spatial_subtype %in% c("EnR", "SF"))

## 3) Load Hallmark gene sets 
library(qusage)
genesets <- read.gmt(file = "~/Desktop/final_scripts/00_data/h.all.v2023.2.Hs.symbols.gmt")

############################################################
## TUMOR PSEUDOBULKS
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)

## Prepare tumor matrix ------------------------------------
tumor_mat <- as.data.frame(tumor_pb)
colnames(tumor_mat) <- paste0("ST", colnames(tumor_mat))

tumor_samples <- intersect(colnames(tumor_mat), meta_ES$sample)
tumor_mat <- tumor_mat[, tumor_samples]

meta_tumor <- meta_ES %>%
  filter(sample %in% tumor_samples) %>%
  arrange(match(sample, tumor_samples))

meta_tumor$spatial_subtype <- droplevels(meta_tumor$spatial_subtype)

## GSVA (tumor) — YOUR STYLE --------------------------------
params_tumor <- gsvaParam(
  as.matrix(tumor_mat),
  kcdf = "Gaussian",
  geneSets = genesets
)

gsva_res_tumor <- gsva(params_tumor)

module_sig_tumor <- t(gsva_res_tumor)
module_sig_tumor_m <- as.matrix(module_sig_tumor)
colnames(module_sig_tumor_m) <- gsub("HALLMARK_", "", colnames(module_sig_tumor_m))

design_tumor <- model.matrix(~ spatial_subtype, data = meta_tumor)
fit_tumor <- lmFit(t(module_sig_tumor_m), design_tumor)
fit_tumor <- eBayes(fit_tumor)

stats_tumor <- topTable(
  fit_tumor,
  coef = "spatial_subtypeSF",
  number = Inf,
  adjust.method = "fdr"
) %>%
  rownames_to_column("Signature") %>%
  arrange(adj.P.Val)

cell_cycle_sigs <- c(
  "E2F_TARGETS",
  "G2M_CHECKPOINT",
  "MYC_TARGETS_V1",
  "MYC_TARGETS_V2",
  "MITOTIC_SPINDLE"
)

cell_cycle_sigs <- intersect(cell_cycle_sigs, colnames(module_sig_tumor_m))


df_cellcycle <- module_sig_tumor_m[
  meta_tumor$sample,
  cell_cycle_sigs,
  drop = FALSE
] %>%
  as.data.frame() %>%
  rownames_to_column("Sample") %>%
  pivot_longer(
    cols = all_of(cell_cycle_sigs),
    names_to = "Signature",
    values_to = "GSVA_score"
  ) %>%
  left_join(meta_tumor, by = c("Sample" = "sample"))


fdr_df <- stats_tumor %>%
  filter(Signature %in% cell_cycle_sigs) %>%
  dplyr::select(Signature, adj.P.Val) %>%
  mutate(
    FDR_label = paste0("FDR=", signif(adj.P.Val, 2))
  )

df_cellcycle <- df_cellcycle %>%
  left_join(fdr_df, by = "Signature")


pretty_names <- c(
  E2F_TARGETS    = "E2F targets",
  G2M_CHECKPOINT = "G2M checkpoint",
  MYC_TARGETS_V1 = "MYC targets V1",
  MYC_TARGETS_V2 = "MYC targets V2",
  MITOTIC_SPINDLE = "Mitotic spindle"
)

df_cellcycle$Signature_pretty <- pretty_names[df_cellcycle$Signature]

df_cellcycle$Signature_pretty <- factor(
  df_cellcycle$Signature_pretty,
  levels = pretty_names[cell_cycle_sigs])

subtype_cols <- c(
  EnR = "#8CBF40",
  SF  = "#E46B71"
)


label_df <- df_cellcycle %>%
  group_by(Signature, Signature_pretty) %>%
  summarise(
    y = max(GSVA_score, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  left_join(fdr_df, by = "Signature") %>%
  mutate(
    y = y + 0.12 * (max(df_cellcycle$GSVA_score, na.rm = TRUE) -
                      min(df_cellcycle$GSVA_score, na.rm = TRUE))
  )

p_tumor_cellcycle <- ggplot(
  df_cellcycle,
  aes(x = spatial_subtype, y = GSVA_score, fill = spatial_subtype)
) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.85) +
  geom_jitter(width = 0.15, size = 1.3, color = "black") +
  geom_text(
    data = label_df,
    aes(x = 1.5, y = y, label = FDR_label),
    inherit.aes = FALSE,
    size = 3.1
  ) +
  facet_wrap(~ Signature_pretty, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = subtype_cols) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", size = 10),
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(face = "bold")
  ) +
  labs(
    title = "Tumor pseudobulk: cell-cycle-associated programs",
    y = "GSVA score"
  )

print(p_tumor_cellcycle)


