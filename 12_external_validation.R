# 12_external_validation.R


##################################################
## External Validation Cohorts Data Preparation
##################################################

##################################################
## METABRIC
##################################################

library("pgirmess")
library(survival)
library(survminer)
library(patchwork)
library(dplyr)


# importing clinical data
samples = read.table("data_clinical_sample.txt", header = T, sep = "\t")
patients = read.table("data_clinical_patient.txt", header = T, sep = "\t")

# filtering patients for histology (ductal)
samples_lobular = samples[samples$CANCER_TYPE_DETAILED == "Breast Invasive Ductal Carcinoma", ]
# filtering patients for HER2 status and Horome Receptors status
ductal_her2neg = samples_lobular[samples_lobular$HER2_STATUS == "Negative", ]
ductal_her2neg_hrpos = ductal_her2neg[ductal_her2neg$ER_STATUS == "Positive" | ductal_her2neg$PR_STATUS == "Positive", ]
library(survival)
library(lubridate)
# transforming survival data in the correct format
patients$OS_STATUS = as.numeric(gsub(":.*", "", patients$OS_STATUS))
patients$RFS_STATUS = as.numeric(gsub(":.*", "", patients$RFS_STATUS))
# storing the patients ID
patients_filtered = patients[patients$PATIENT_ID %in% ductal_her2neg_hrpos$PATIENT_ID, ]

# opening gene expression data
metabric_expr = read.delim("data_mrna_agilent_microarray.txt", row.names = "Hugo_Symbol")
metabric_expr$Entrez_Gene_Id = NULL
colnames(metabric_expr) = sub("\\.", "-", colnames(metabric_expr))
# filtering gene expression data
expression_matrix = as.data.frame(metabric_expr[, colnames(metabric_expr) %in% ductal_her2neg_hrpos$SAMPLE_ID])

# final clean expression data:
expression_matrix = expression_matrix[,sort(colnames(expression_matrix))]
# final clean metadata:
df_analysis = cbind(patients_filtered, ductal_her2neg_hrpos[,2:13])

clinical_data_metabric <- df_analysis[,c(1,23,24)]

expression_metabric$NAME <- rownames(expression_metabric)

mart_export <- read.delim("/Users/bengisukarakose/Desktop/bc_cox_model/mart_export.txt")
colnames(mart_export) <- c("DESCRIPTION", "NAME")

expression_metabric <- left_join(expression_metabric, mart_export, by="NAME", multiple = "first" )

expression_metabric <- expression_metabric %>%
  select("NAME", "DESCRIPTION", everything())

expression_metabric$DESCRIPTION[is.na(expression_metabric$DESCRIPTION)] <- "info"
# saved as counts_metabric.txt

colnames(clinical_data_metabric) <- c("name", "status", "time")
clinical_data_metabric$time <- round(clinical_data_metabric$time * 30 , )
# saved as clin_rev_metabric.txt

## let's prepare final data
counts_metabric <- counts_metabric[,-2]
rownames(counts_metabric) <- counts_metabric$NAME
counts_metabric <- counts_metabric[,-1]

data_clinical_patient <- read.delim("~/data_clinical_patient.txt", comment.char="#")
data_clinical_sample <- read.delim("~/data_clinical_sample.txt", comment.char="#")

table(data_clinical_patient$PATIENT_ID == data_clinical_sample$PATIENT_ID)

metabric_meta <- cbind(data_clinical_patient, data_clinical_sample)
colnames(counts_metabric) = sub("\\.", "-", colnames(counts_metabric))
metabric_meta = metabric_meta[metabric_meta$PATIENT_ID %in% colnames(counts_metabric),]

# saved as metabric_meta.RDS and counts_metabric.RDS


##################################################
## SCAN-B
##################################################

library("pgirmess")
library(survival)
library(survminer)
library(patchwork)
library(dplyr)
library(openxlsx)
library(readxl)


metadata = as.data.frame(openxlsx::read.xlsx("Supplementary Data Table 1 - 2023-01-13.xlsx", 1))
# filtering for histology (ductal)
metadata_ductal = metadata[which(metadata$InvCa.type == "Ductal"), ]
# filtering for HER2 status
metadata_ductal_her2_neg = metadata_ductal[which(metadata_ductal$HER2 == "Negative"), ]
# filtering for Hormone Receptor status
metadata_ductal_her2_neg_hr_pos = metadata_ductal_her2_neg[which(metadata_ductal_her2_neg$ER == "Positive" | metadata_ductal_her2_neg$PR == "Positive"), ]
# taking just the follow-up cohort (where we have survival information)
metadata_ductal_her2_neg_hr_pos$Follow.up.cohort = as.character(metadata_ductal_her2_neg_hr_pos$Follow.up.cohort)
metadata_fu = metadata_ductal_her2_neg_hr_pos[which(metadata_ductal_her2_neg_hr_pos$Follow.up.cohort == TRUE), ]

# importing gene expression data
load("SCANB.9206.genematrix_noNeg.Rdata")
expression = as.data.frame(SCANB.9206.genematrix_noNeg)
gene_ids = read.table("Gene.ID.ann.txt", header = T)
names = gene_ids[gene_ids$Gene.ID %in% rownames(expression), "Gene.Name"]
entrez = gene_ids[gene_ids$Gene.ID %in% rownames(expression), "EntrezGene"]
expression$names = names
# removing duplicates (some gene names have duplicates - different isoforms - so a common thing to do is to choose the gene with the highest standard deviation [it should be the most informative])
duplicated = names(table(names)[table(names) > 1])
toremove = rep(0, times = length(duplicated))
i = 1
for (gene in duplicated) {
  ensg1 = rownames(expression[expression$names == gene,])[1]
  ensg2 = rownames(expression[expression$names == gene,])[2]
  sdv1 = sd(expression[rownames(expression) == ensg1, 1:9206])
  sdv2 = sd(expression[rownames(expression) == ensg2, 1:9206])
  v_name = c("sdv1","sdv2")
  maxi = v_name[which.max(c(sdv1,sdv2))]
  if (maxi == "sdv2") {
    toremove[i] = ensg1
  } else {
    toremove[i] = ensg2
  }
  i = i+1
}
expression_filter = expression[!rownames(expression) %in% toremove, ]
rownames(expression_filter) = expression_filter$names
expression_filter$names = NULL

# filtering gene expression with metadata information and log-transforming the data
expression_duc = expression_filter[,colnames(expression_filter) %in% metadata_fu$GEX.assay]

# final clean expression data
expression_duc = log2(expression_duc + 1)
# final clean metadata information
df_analysis = metadata_fu

expression_duc$NAME <- rownames(expression_duc)

mart_export <- read.delim("/Users/bengisukarakose/Desktop/bc_cox_model/mart_export.txt")
colnames(mart_export) <- c("DESCRIPTION", "NAME")

expression_duc <- left_join(expression_duc, mart_export, by="NAME", multiple = "first" )

expression_duc <- expression_duc %>%
  select("NAME", "DESCRIPTION", everything())

expression_duc$DESCRIPTION[is.na(expression_duc$DESCRIPTION)] <- "info"
## save as counts_scanb.txt

clinical_data_scanb <- df_analysis[,c(1,53,54)]
#DRFi days

colnames(clinical_data_scanb) <- c("name", "time", "status")

clinical_data_scanb <- clinical_data_scanb %>%
  select("name", "status", "time")
## save as clin_rev_scanb.txt


## let's prepare final data
counts_scanb <- counts_scanb[,-2]
rownames(counts_scanb) <- counts_scanb$NAME
counts_scanb <- counts_scanb[,-1]

clin_rev_scanb <- read.delim("~/clin_rev_scanb.txt")
clin_all_scanb <- read_excel("~/Supplementary Data Table 1 - 2023-01-13.xlsx")
clin_filtered <- clin_all_scanb %>% filter(GEX.assay %in% colnames(counts_scanb))

## saved as counts_scanb.RDS and clin_SCANB_filtered.RDS



##################################################
## Make gene signatures for Spatial4HR+ subtyping
##################################################

library(DESeq2)
library(dplyr)
library(pheatmap)

# Load data
ductal_pb <- readRDS("~/Desktop/bc_cox_model/data/pb_ductals_not_norm.RDS")
ductal_meta <- read.delim("~/Desktop/final_scripts/00_data/Spatial4HR_sample_metadata.txt")
ductal_meta$ductal_subtype <- as.factor(ductal_meta$ductal_subtype)

# Define DESeq2 comparison function
run_DESeq2_comparison <- function(count_data, metadata, group_column, target_level) {
  metadata$ductal_subtype_group <- ifelse(metadata[[group_column]] == target_level, target_level, "others")
  metadata$ductal_subtype_group <- factor(metadata$ductal_subtype_group, levels = c("others", target_level))
  dds <- DESeqDataSetFromMatrix(countData = count_data, colData = metadata, design = ~ ductal_subtype_group)
  dds <- DESeq(dds)
  results(dds, contrast = c("ductal_subtype_group", target_level, "others"))
}

# Run DESeq2 for each pattern
res_1_run <- run_DESeq2_comparison(ductal_pb, ductal_meta, "ductal_subtype", "SF")
res_2_run <- run_DESeq2_comparison(ductal_pb, ductal_meta, "ductal_subtype", "EnR")
res_3_run <- run_DESeq2_comparison(ductal_pb, ductal_meta, "ductal_subtype", "PI")
res_4_run <- run_DESeq2_comparison(ductal_pb, ductal_meta, "ductal_subtype", "HL")

# Filter significant DEGs
res_1 <- as.data.frame(res_1_run)
res_2 <- as.data.frame(res_2_run)
res_3 <- as.data.frame(res_3_run)
res_4 <- as.data.frame(res_4_run)

# Filter DEGs for each subtype: padj < 0.05 and log2FoldChange > threshold (e.g., >1)
threshold <- 0  # Adjust as needed
deg_1 <- res_1[!is.na(res_1$padj) & res_1$padj < 0.05 & res_1$log2FoldChange > threshold, ]
deg_2 <- res_2[!is.na(res_2$padj) & res_2$padj < 0.05 & res_2$log2FoldChange > threshold, ]
deg_3 <- res_3[!is.na(res_3$padj) & res_3$padj < 0.05 & res_3$log2FoldChange > threshold, ]
deg_4 <- res_4[!is.na(res_4$padj) & res_4$padj < 0.05 & res_4$log2FoldChange > threshold, ]


### filter for protein coding genes

library(biomaRt)
ensembl <- useEnsembl(
  biomart  = "genes",
  dataset  = "hsapiens_gene_ensembl"
)

deg_list <- list(
  deg_1 = deg_1,
  deg_2 = deg_2,
  deg_3 = deg_3,
  deg_4 = deg_4
)

all_de_genes <- unique(unlist(lapply(deg_list, rownames)))

gene_info <- getBM(
  attributes = c("hgnc_symbol", "gene_biotype"),
  filters    = "hgnc_symbol",
  values     = all_de_genes,
  mart       = ensembl
)

pc_symbols <- gene_info$hgnc_symbol[ gene_info$gene_biotype == "protein_coding" ]

deg_list_pc <- lapply(deg_list, function(df) {
  df[rownames(df) %in% pc_symbols, , drop = FALSE]
})

deg_1_pc <- deg_list_pc$deg_1
deg_2_pc <- deg_list_pc$deg_2
deg_3_pc <- deg_list_pc$deg_3
deg_4_pc <- deg_list_pc$deg_4

sapply(deg_list_pc, nrow)

deg_1 <- deg_1_pc
deg_2 <- deg_2_pc
deg_3 <- deg_3_pc
deg_4 <- deg_4_pc


# Select top marker genes for each subtype
top_n <- 30
top_genes_1 <- head(rownames(deg_1[order(-deg_1$log2FoldChange), ]), top_n)
top_genes_2 <- head(rownames(deg_2[order(-deg_2$log2FoldChange), ]), top_n)
top_genes_3 <- head(rownames(deg_3[order(-deg_3$log2FoldChange), ]), top_n)
top_genes_4 <- head(rownames(deg_4[order(-deg_4$log2FoldChange), ]), top_n)

# Compile gene signatures for each subtype
gene_signatures <- list(
  SF = top_genes_1,
  EnR = top_genes_2,
  PI = top_genes_3,
  HL = top_genes_4
)




##################################################
## Assign Spatial4HR+ subtypes to external cohorts - METABRIC
##################################################

library(dplyr)
library(survival)
library(survminer)

metabric_meta <- readRDS("~/Desktop/final_scripts/00_data/metabric_meta.RDS")
counts_metabric <- readRDS("~/Desktop/final_scripts/00_data/counts_metabric.RDS")
genesets_sig <- readRDS("~/Desktop/final_scripts/00_data/gene_signatures_protein_coding.rds")

# Create the GSVAParams object for the "gsva" method
params <- gsvaParam(as.matrix(counts_metabric),
                    kcdf="Gaussian", geneSets = genesets_sig)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)

table(rownames(module_sig) == metabric_meta$PATIENT_ID)

metabric_meta <- cbind(metabric_meta, module_sig)

metabric_meta$OS_STATUS = as.numeric(gsub(":.*", "", metabric_meta$OS_STATUS))
metabric_meta$RFS_STATUS = as.numeric(gsub(":.*", "", metabric_meta$RFS_STATUS))

# Now, assign every patient a subtype:
# For each patient, choose the subtype with the highest GSVA score
metabric_meta$predicted_subtype <- apply(module_sig, 1, function(x) names(x)[which.max(x)])
metabric_meta$predicted_subtype <- factor(metabric_meta$predicted_subtype, levels = c('EnR', 'SF', 'PI', 'HL'))

# View distribution of subtype assignments
print(table(metabric_meta$predicted_subtype))
# EnR  SF  PI  HL 
# 268 234 224 315 



###########################
# Fig. 7a - Survival Analysis Based on Predicted Subtypes in METABRIC
###########################

surv_obj_ext <- Surv(metabric_meta$RFS_MONTHS, metabric_meta$RFS_STATUS)
fit_ext <- survfit(Surv(RFS_MONTHS/12, RFS_STATUS) ~ predicted_subtype, data = metabric_meta)

# colors in the same order as the factor above (EnR, SF, PI, HL)
pal <- c("#8CBF40",  "#E46B71", "#38A7B0", "#8C6BB1")

ggsurvplot(
  fit_ext,
  data = metabric_meta,
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = TRUE,          # show row labels
  risk.table.y.text.col = TRUE,      # color row labels by strata
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (years)",
  ylab = "Survival Probability",
  title = "METABRIC Cohort Validation",
  legend.title = "Spatial Subtype",
  legend.labs = levels(metabric_meta$predicted_subtype),  # EnR, SF, PI, HL
  palette = pal,
  break.time.by = 4
)

# Perform pairwise log-rank tests using a formula defined directly:
pairwise_results <- pairwise_survdiff(Surv(RFS_MONTHS, RFS_STATUS) ~ predicted_subtype, 
                                      data = metabric_meta, 
                                      p.adjust.method = "BH")
print(pairwise_results)
# Pairwise comparisons using Log-Rank test  
# data:  metabric_meta and predicted_subtype 
#     EnR   SF    PI   
# SF 0.145 -     -    
# PI 0.023 0.329 -    
# HL 0.010 0.329 0.758
# P value adjustment method: BH 



###########################
# Fig. 7c - xCell inferred cell-type composition across Spatial4HR+ in METABRIC 
###########################

metabric_meta <- readRDS("~/Desktop/final_scripts/00_data/metabric_meta.RDS")
counts_metabric <- readRDS("~/Desktop/final_scripts/00_data/counts_metabric.RDS")

risk_groups <- as.data.frame(cbind(metabric_meta$PATIENT_ID, metabric_meta$predicted_subtype))
names(risk_groups) <- c("X1", "X2")
risk_groups$X2 <- factor(risk_groups$X2,
                         levels = c("EnR", "SF", "PI", "HL"))

sample_ids <- risk_groups$X1
counts_metabric <- counts_metabric[, colnames(counts_metabric) %in% sample_ids]

# computing xCell
library(xCell)
xcell_results = as.data.frame(t(xCellAnalysis(counts_metabric)))
cell_types = colnames(xcell_results)[c(2,4:14,18:23,25:26,31:34,38,41,44,48,46,49,50,54,57,61:64)]
table_xcell = as.data.frame(t(as.matrix(as.data.frame(xCellAnalysis(counts_metabric, cell.types.use = cell_types)))))
table_xcell$Subtypes = risk_groups$X2

data_plot <- as.data.frame(cbind(table_xcell,
                                 matrix(rep(table_xcell$Subtypes, length(table(table_xcell$Subtypes)))
                                        ,ncol=length(table(table_xcell$Subtypes)))))

data_plot[,c(38)] = NULL

colnames(data_plot)[38:41] <- c("EnR", "SF", "PI", "HL")
data_plot[,"EnR"] <- as.factor(ifelse(data_plot[,"EnR"]=="EnR", 1, 0))
data_plot[,"SF"] <- as.factor(ifelse(data_plot[,"SF"]=="SF", 1, 0))
data_plot[,"PI"] <- as.factor(ifelse(data_plot[,"PI"]=="PI", 1, 0))
data_plot[,"HL"] <- as.factor(ifelse(data_plot[,"HL"]=="HL", 1, 0))

colnames(data_plot)
data_plot[,1:37] = apply(data_plot[,1:37], 2, genefu::rescale, na.rm = T)

subtype = names(table(table_xcell$Subtypes))

effect <- matrix(0, nrow=length(colnames(data_plot[,1:37])), ncol=length( sort(unique(subtype)) ))
colnames(effect) <- sort(unique(subtype))
rownames(effect) <- colnames(data_plot[,1:37])

pvalue <- matrix(0, nrow=length(colnames(data_plot[,1:37])), ncol=length( sort(unique(subtype)) ))
colnames(pvalue) <- sort(unique(subtype))
rownames(pvalue) <- colnames(data_plot[,1:37])

ci <- matrix(0, nrow=length(colnames(data_plot[,1:37])), ncol=length( sort(unique(subtype)) ))
colnames(ci) <- sort(unique(subtype))
rownames(ci) <- colnames(data_plot[,1:37])

signatures = data_plot[,1:37]

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
effect[effect<=.25] <- .25 
effect[effect>=4] <- 4

effect <- t(effect)
#effect <- effect[, ncol(effect):1]

e_ductalST = effect

e_ductalST <- e_ductalST[c("EnR", "SF", "PI", "HL"),]

dim(e_ductalST)
rownames(e_ductalST) = gsub("_", " ", rownames(e_ductalST))
colnames(e_ductalST) = gsub("\\.", " ", colnames(e_ductalST))

library(gtools)
current_row_names <- rownames(e_ductalST)
sorted_row_names <- current_row_names
e_ductalST_sorted <- e_ductalST[sorted_row_names, ]
rownames(e_ductalST_sorted) <- sorted_row_names
print(rownames(e_ductalST_sorted))
e_ductalST = e_ductalST_sorted



plot_legend_1 <- matrix(c("0.25",
                          "4"
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
  row(e_ductalST)*1,
  col(e_ductalST)*1,
  type="n",xlab="", ylab="",
  xlim=c(0.5,nrow(e_ductalST)*1+0.5),
  ylim=c(0.5,ncol(e_ductalST)*1+1),
  axes=FALSE,
  ann=FALSE, asp=1)

for(i in 1:nrow(e_ductalST)){
  for(j in 1:ncol(e_ductalST)){
    circle(i*1,j*1,abs(log2(e_ductalST[i,j]))*0.2, # --> 0.05 is the size of the cirlces
           col=adjustcolor(ifelse(e_ductalST[i,j]>1,"#009E73","#D55E00"), alpha.f = 0.8),lwd = .7, linecolor = "white", border = NA)
  }
}


axis(1,at=1:nrow(e_ductalST)*1,labels=FALSE, pos = 0) # --> *1 is the distance between things on x axis
text(x=rep(.2,ncol(e_ductalST)*1), y=1:ncol(e_ductalST)*1, labels=add_space(add_sign_neg(add_sign_pos(colnames(e_ductalST)))), cex.axis=1.5, las=2, srt = 0,pos = 2, xpd = T)
axis(2,at=1:ncol(e_ductalST)*1,pos=0.5,labels=FALSE)
text(x=1:nrow(e_ductalST)*1, y=rep(-1,nrow(e_ductalST)*1), labels=add_space(add_sign_neg(add_sign_pos(rownames(e_ductalST)))), cex.axis=1.5, srt = 0, las=2, xpd = T, srt = 90, adj = 1)


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





###########################
# Fig. 7d - Hallmark pathway activity across Spatial4HR+ subtypes in METABRIC 
###########################

metabric_meta <- readRDS("~/Desktop/final_scripts/00_data/metabric_meta.RDS")
counts_metabric <- readRDS("~/Desktop/final_scripts/00_data/counts_metabric.RDS")

risk_groups <- as.data.frame(cbind(metabric_meta$PATIENT_ID, metabric_meta$predicted_subtype))
names(risk_groups) <- c("X1", "X2")
risk_groups$X2 <- factor(risk_groups$X2,
                         levels = c("EnR", "SF", "PI", "HL"))

sample_ids <- risk_groups$X1
counts_metabric <- counts_metabric[, colnames(counts_metabric) %in% sample_ids]

library(qusage)
library(GSVA)
library(dplyr)
library(survival)
library(survminer)

genesets_hallmark <- read.gmt("~/Desktop/final_scripts/00_data/h.all.v2023.2.Hs.symbols.gmt")
genesets_hallmark <- genesets_hallmark[-c(38,43,47,48)]

# We excluded the following 4 hallmark signatures due to a
# lack of association with tumor processes or microenvironment in BC:
# pancreas beta cells, spermatogenesis, UV response up, UV response down.

# Create the GSVAParams object for the "gsva" method
params <- gsvaParam(as.matrix(counts_metabric),
                    kcdf="Gaussian", geneSets = genesets_hallmark)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)
module_sig <- as.data.frame(module_sig)

table(rownames(module_sig) == risk_groups$X1)

module_sig$Subtypes = risk_groups$X2


#### hallmark ####
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

subtype = names(table(module_sig$Subtypes)) 

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

e_ductalST = effect

e_ductalST <- e_ductalST[c("EnR", "SF", "PI", "HL"),]

dim(e_ductalST)

colnames(e_ductalST) = gsub("_", " ", colnames(e_ductalST))
colnames(e_ductalST) = gsub("HALLMARK ", "", colnames(e_ductalST))


library(gtools)
current_row_names <- rownames(e_ductalST)
sorted_row_names <- current_row_names
e_ductalST_sorted <- e_ductalST[sorted_row_names, ]
rownames(e_ductalST_sorted) <- sorted_row_names
print(rownames(e_ductalST_sorted))
e_ductalST = e_ductalST_sorted



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
  row(e_ductalST)*1,
  col(e_ductalST)*1,
  type="n",xlab="", ylab="",
  xlim=c(0.5,nrow(e_ductalST)*1+0.5),
  ylim=c(0.5,ncol(e_ductalST)*1+1),
  axes=FALSE,
  ann=FALSE, asp=1)

for(i in 1:nrow(e_ductalST)){
  for(j in 1:ncol(e_ductalST)){
    circle(i*1,j*1,abs(log2(e_ductalST[i,j]))*0.13, # --> 0.05 is the size of the cirlces
           col=adjustcolor(ifelse(e_ductalST[i,j]>1,"#009E73","#D55E00"), alpha.f = 0.8),lwd = .7, linecolor = "white", border = NA)
  }
}


axis(1,at=1:nrow(e_ductalST)*1,labels=FALSE, pos = 0) # --> *1 is the distance between things on x axis
text(x=rep(.2,ncol(e_ductalST)*1), y=1:ncol(e_ductalST)*1, labels=add_space(add_sign_neg(add_sign_pos(colnames(e_ductalST)))), cex = 0.9, cex.axis=1.5, las=2, srt = 0,pos = 2, xpd = T)
axis(2,at=1:ncol(e_ductalST)*1,pos=0.5,labels=FALSE)
text(x=1:nrow(e_ductalST)*1, y=rep(-1,nrow(e_ductalST)*1), labels=add_space(add_sign_neg(add_sign_pos(rownames(e_ductalST)))), cex.axis=1.5, srt = 0, las=2, xpd = T, srt = 90, adj = 1)


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



##################################################
## Extended Data Fig. 10a  - Distribution of PAM50 subtypes across Spatial4HR+ 
## subtypes in the METABRIC cohort
##################################################

metabric_meta <- readRDS("~/Desktop/final_scripts/00_data/metabric_meta.RDS")

table(metabric_meta$CLAUDIN_SUBTYPE, metabric_meta$predicted_subtype)

metabric_meta$predicted_subtype <- factor(metabric_meta$predicted_subtype)
metabric_meta$pam50 <- factor(metabric_meta$CLAUDIN_SUBTYPE)

contingency <- table(metabric_meta$predicted_subtype, metabric_meta$pam50)
print(contingency)

chisq_res <- chisq.test(contingency)
print(chisq_res)

metabric_meta_filt <- metabric_meta[,c(1,42,43)]

library(ggplot2)
library(randomcoloR)

ggplot(metabric_meta_filt, aes(x = predicted_subtype, fill = pam50)) +
  geom_bar(position = "fill") +  # Proportional bars
  scale_y_continuous(labels = scales::percent_format()) +
  scale_fill_manual(values = distinctColorPalette(k = length(unique(metabric_meta$pam50)))) +
  labs(
    title = "Proportions of PAM50 classes by Ductal Subtype",
    x = "Ductal Subtype",
    y = "Percentage",
    fill = "PAM50"
  ) +
  theme_minimal(base_size = 14)



##################################################
## Fig. 7e  - ERBB2 copy-number status across Spatial4HR+ subtypes in METABRIC 
##################################################

colnames(metabric_meta) <- make.unique(colnames(metabric_meta))

metabric_meta$HER2_SNP6 <- factor(metabric_meta$HER2_SNP6)
table(metabric_meta$predicted_subtype, metabric_meta$HER2_SNP6)

metabric_meta$predicted_subtype <- factor(metabric_meta$predicted_subtype, levels = c('EnR', 'SF', 'PI', 'HL'))

contingency <- table(metabric_meta$predicted_subtype, metabric_meta$HER2_SNP6)
print(contingency)

chisq_res <- chisq.test(contingency)
print(chisq_res)
p_val <- paste0('p = ', round(chisq_res$p.value, 10))

ggplot(
  metabric_meta[!is.na(metabric_meta$HER2_SNP6), ],
  aes(x = predicted_subtype, fill = HER2_SNP6)
) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(
    subtitle = p_val,
    x = "Spatial4HR+ subtype",
    y = "Percentage",
    fill = "HER2 copy-number status"
  ) +
  theme_minimal(base_size = 14)





##################################################
## Assign Spatial4HR+ subtypes to external cohorts - SCAN-B
##################################################

library(dplyr)
library(survival)
library(survminer)

clin_filtered <- readRDS("~/Desktop/final_scripts/00_data/clin_SCANB_filtered.RDS")
counts_scanb <- readRDS("~/Desktop/final_scripts/00_data/counts_scanb.RDS")
genesets_sig <- readRDS("~/Desktop/final_scripts/00_data/gene_signatures_protein_coding.rds")

params <- gsvaParam(as.matrix(counts_scanb), kcdf="Gaussian", geneSets = genesets_sig)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)

table(rownames(module_sig) == clin_rev_scanb$name)
scanb_meta <- cbind(clin_rev_scanb, module_sig)

table(rownames(module_sig) == clin_filtered$GEX.assay)
scanb_big_meta <- cbind(clin_filtered, module_sig)

# Now, assign every patient a subtype:
# For each patient, choose the subtype with the highest GSVA score
scanb_big_meta$predicted_subtype <- apply(module_sig, 1, function(x) names(x)[which.max(x)])

print(table(scanb_big_meta$predicted_subtype))
# EnR   HL   PI   SF 
# 694 1063  752 1273 

# saveRDS(scanb_big_meta, 'scanb_big_meta.RDS')


###########################
# Fig. 7b - Survival Analysis Based on Predicted Subtypes in SCAN-B
###########################

library(dplyr)
library(stringr)
library(survival)
library(survminer)
library(gridExtra)

scanb_big_meta <- readRDS("~/Desktop/final_scripts/00_data/scanb_big_meta.RDS")
scanb_big_meta$predicted_subtype <- factor(scanb_big_meta$predicted_subtype, 
                                           levels = c('EnR', 'SF', 'PI', 'HL'))

surv_obj_ext <- Surv(scanb_big_meta$RFi_days, scanb_big_meta$RFi_event)
fit_ext <- survfit(surv_obj_ext ~ predicted_subtype, data = scanb_big_meta)

km_fit <- survfit(Surv(RFi_days/365, RFi_event) ~ predicted_subtype, data = scanb_big_meta)

pal <- c("#8CBF40", "#E46B71", "#38A7B0", "#8C6BB1")  # EnR, SF, PI, HL

ggsurvplot(
  km_fit,
  data = scanb_big_meta,
  risk.table = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text = TRUE,          # show row labels
  risk.table.y.text.col = TRUE,      # color row labels by strata
  pval = TRUE,
  conf.int = FALSE,
  xlab = "Time (years)",
  ylab = "Survival Probability",
  title = "SCAN-B Cohort Validation",
  legend.title = "Spatial Subtype",
  legend.labs = levels(scanb_big_meta$predicted_subtype),  # EnR, SF, PI, HL
  palette = pal,
  break.time.by = 2
)


pairwise_survdiff(Surv(RFi_days, RFi_event) ~ predicted_subtype, 
                                      data = scanb_big_meta, 
                                      p.adjust.method = "BH")


###########################
# Extended Data Fig. 10b - xCell inferred cell-type composition across Spatial4HR+ in SCAN-B 
###########################

scanb_big_meta <- readRDS("~/Desktop/final_scripts/00_data/scanb_big_meta.RDS")
counts_scanb <- readRDS("~/Desktop/final_scripts/00_data/counts_scanb.RDS")

risk_groups <- as.data.frame(cbind(scanb_big_meta$GEX.assay, scanb_big_meta$predicted_subtype))
names(risk_groups) <- c("X1", "X2")
# risk_groups$X1 <- gsub("-", ".", risk_groups$X1)
# 
# sample_ids <- risk_groups$X1
# counts_scanb <- counts_scanb[, colnames(counts_scanb) %in% sample_ids]


# computing xCell
library(xCell)
xcell_results = as.data.frame(t(xCellAnalysis(counts_scanb)))
cell_types = colnames(xcell_results)[c(2,4:14,18:23,25:26,31:34,38,41,44,48,46,49,50,54,57,61:64)]
table_xcell = as.data.frame(t(as.matrix(as.data.frame(xCellAnalysis(counts_scanb, cell.types.use = cell_types)))))
table_xcell$Subtypes = risk_groups$X2
table_xcell$Subtypes <- paste0(table_xcell$Subtypes)

data_plot <- as.data.frame(cbind(table_xcell,
                                 matrix(rep(table_xcell$Subtypes, length(table(table_xcell$Subtypes)))
                                        ,ncol=length(table(table_xcell$Subtypes)))))

data_plot[,c(38)] = NULL

colnames(data_plot)[38:41] <- names(table(table_xcell$Subtypes))
data_plot[,"EnR"] <- as.factor(ifelse(data_plot[,"EnR"]=="EnR", 1, 0))
data_plot[,"SF"] <- as.factor(ifelse(data_plot[,"SF"]=="SF", 1, 0))
data_plot[,"PI"] <- as.factor(ifelse(data_plot[,"PI"]=="PI", 1, 0))
data_plot[,"HL"] <- as.factor(ifelse(data_plot[,"HL"]=="HL", 1, 0))

colnames(data_plot)
data_plot[,1:37] = apply(data_plot[,1:37], 2, genefu::rescale, na.rm = T)

subtype = names(table(table_xcell$Subtypes)) ## la colonna con i cluster, o i nomi dei cluster direttamente

effect <- matrix(0, nrow=length(colnames(data_plot[,1:37])), ncol=length( sort(unique(subtype)) ))
colnames(effect) <- sort(unique(subtype))
rownames(effect) <- colnames(data_plot[,1:37])

pvalue <- matrix(0, nrow=length(colnames(data_plot[,1:37])), ncol=length( sort(unique(subtype)) ))
colnames(pvalue) <- sort(unique(subtype))
rownames(pvalue) <- colnames(data_plot[,1:37])

ci <- matrix(0, nrow=length(colnames(data_plot[,1:37])), ncol=length( sort(unique(subtype)) ))
colnames(ci) <- sort(unique(subtype))
rownames(ci) <- colnames(data_plot[,1:37])

signatures = data_plot[,1:37]

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
effect[effect<=.25] <- .25 
effect[effect>=4] <- 4

effect <- t(effect)
#effect <- effect[, ncol(effect):1]

e_ductalST = effect
e_ductalST <- e_ductalST[c("EnR", "SF", "PI", "HL"),]

dim(e_ductalST)

rownames(e_ductalST) = gsub("_", " ", rownames(e_ductalST))
colnames(e_ductalST) = gsub("\\.", " ", colnames(e_ductalST))

library(gtools)

plot_legend_1 <- matrix(c("0.25",
                          "4"
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
  row(e_ductalST)*1,
  col(e_ductalST)*1,
  type="n",xlab="", ylab="",
  xlim=c(0.5,nrow(e_ductalST)*1+0.5),
  ylim=c(0.5,ncol(e_ductalST)*1+1),
  axes=FALSE,
  ann=FALSE, asp=1)

for(i in 1:nrow(e_ductalST)){
  for(j in 1:ncol(e_ductalST)){
    circle(i*1,j*1,abs(log2(e_ductalST[i,j]))*0.2, # --> 0.05 is the size of the cirlces
           col=adjustcolor(ifelse(e_ductalST[i,j]>1,"#009E73","#D55E00"), alpha.f = 0.8),lwd = .7, linecolor = "white", border = NA)
  }
}


axis(1,at=1:nrow(e_ductalST)*1,labels=FALSE, pos = 0) # --> *1 is the distance between things on x axis
text(x=rep(.2,ncol(e_ductalST)*1), y=1:ncol(e_ductalST)*1, labels=add_space(add_sign_neg(add_sign_pos(colnames(e_ductalST)))), cex.axis=1.5, las=2, srt = 0,pos = 2, xpd = T)
axis(2,at=1:ncol(e_ductalST)*1,pos=0.5,labels=FALSE)
text(x=1:nrow(e_ductalST)*1, y=rep(-1,nrow(e_ductalST)*1), labels=add_space(add_sign_neg(add_sign_pos(rownames(e_ductalST)))), cex.axis=1.5, srt = 0, las=2, xpd = T, srt = 90, adj = 1)


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





###########################
# Extended Data Fig. 10c - Hallmark pathway activity across Spatial4HR+ subtypes in SCAN-B 
###########################

library(qusage)
library(GSVA)
library(dplyr)
library(survival)
library(survminer)

scanb_big_meta <- readRDS("~/Desktop/final_scripts/00_data/scanb_big_meta.RDS")
counts_scanb <- readRDS("~/Desktop/final_scripts/00_data/counts_scanb.RDS")

risk_groups <- as.data.frame(cbind(scanb_big_meta$GEX.assay, scanb_big_meta$predicted_subtype))
names(risk_groups) <- c("X1", "X2")

risk_groups$X2 <- factor(risk_groups$X2,
                         levels = c("EnR", "SF", "PI", "HL"))

genesets_hallmark <- read.gmt("~/Desktop/final_scripts/00_data/h.all.v2023.2.Hs.symbols.gmt")

genesets_hallmark <- genesets_hallmark[-c(38,43,47,48)]
# We excluded the following 4 hallmark signatures due to a
# lack of association with tumor processes or microenvironment in BC:
# pancreas beta cells, spermatogenesis, UV response up, UV response down.


# Create the GSVAParams object for the "gsva" method
params <- gsvaParam(as.matrix(counts_scanb),  
                    kcdf="Gaussian", geneSets = genesets_hallmark)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)
module_sig <- as.data.frame(module_sig)

table(rownames(module_sig) == risk_groups$X1)

module_sig$Subtypes = risk_groups$X2
# module_sig$Subtypes <- paste0('Subtype', module_sig$Subtypes)

#### hallmark ####
data_plot <- as.data.frame(cbind(module_sig,
                                 matrix(rep(module_sig$Subtypes, length(table(module_sig$Subtypes)))
                                        ,ncol=length(table(module_sig$Subtypes)))))

data_plot[,c(47)] = NULL

# cambia "group_1" etc con il nome dei clusters
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
e_ductalST = effect

e_ductalST <- e_ductalST[c("EnR", "SF", "PI", "HL"),]

dim(e_ductalST)

colnames(e_ductalST) = gsub("_", " ", colnames(e_ductalST))
colnames(e_ductalST) = gsub("HALLMARK ", "", colnames(e_ductalST))


library(gtools)
current_row_names <- rownames(e_ductalST)
sorted_row_names <- current_row_names
e_ductalST_sorted <- e_ductalST[sorted_row_names, ]
rownames(e_ductalST_sorted) <- sorted_row_names
print(rownames(e_ductalST_sorted))
e_ductalST = e_ductalST_sorted



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
  row(e_ductalST)*1,
  col(e_ductalST)*1,
  type="n",xlab="", ylab="",
  xlim=c(0.5,nrow(e_ductalST)*1+0.5),
  ylim=c(0.5,ncol(e_ductalST)*1+1),
  axes=FALSE,
  ann=FALSE, asp=1)

for(i in 1:nrow(e_ductalST)){
  for(j in 1:ncol(e_ductalST)){
    circle(i*1,j*1,abs(log2(e_ductalST[i,j]))*0.13, # --> 0.05 is the size of the cirlces
           col=adjustcolor(ifelse(e_ductalST[i,j]>1,"#009E73","#D55E00"), alpha.f = 0.8),lwd = .7, linecolor = "white", border = NA)
  }
}


axis(1,at=1:nrow(e_ductalST)*1,labels=FALSE, pos = 0) # --> *1 is the distance between things on x axis
text(x=rep(.2,ncol(e_ductalST)*1), y=1:ncol(e_ductalST)*1, labels=add_space(add_sign_neg(add_sign_pos(colnames(e_ductalST)))), cex = 0.9, cex.axis=1.5, las=2, srt = 0,pos = 2, xpd = T)
axis(2,at=1:ncol(e_ductalST)*1,pos=0.5,labels=FALSE)
text(x=1:nrow(e_ductalST)*1, y=rep(-1,nrow(e_ductalST)*1), labels=add_space(add_sign_neg(add_sign_pos(rownames(e_ductalST)))), cex.axis=1.5, srt = 0, las=2, xpd = T, srt = 90, adj = 1)


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



###########################
# Fig. 7h,i - Response to endocrine therapy combined with 
# CDK4/6 inhibition in the NeoRHEA cohort
###########################

library(DESeq2)
library(GSVA)
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(broom)

clinDf <- readRDS("~/Desktop/bc_cox_model/104_neoRhea/neoRhea/clinical/data/clinDf.RDS")
neoRhea_raw_counts <- readRDS("~/Desktop/bc_cox_model/104_neoRhea/neoRhea/bulk_rna/data/legacy/raw_cm.RDS")

# varianceStabilizingTransformation to make the dataset ready for GSVA
neoRhea_VST <- varianceStabilizingTransformation(neoRhea_raw_counts)

### Assign subtypes to neoRhea
genesets_sig <- readRDS("~/Desktop/final_scripts/00_data/gene_signatures_protein_coding.rds")
params <- gsvaParam(as.matrix(neoRhea_VST), kcdf="Gaussian", geneSets = genesets_sig)

gsva_res_filtered <- gsva(params)
module_sig <- t(gsva_res_filtered)

clinDf$`Patient ID`
rownames(module_sig)

ms <- module_sig %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(
    PatientID = sub("_[BS]$", "", sample),
    Timepoint = sub("^.*_([BS])$", "\\1", sample),
    Timepoint = recode(Timepoint, B = "Baseline", S = "Surgery")
  ) %>%
  dplyr::select(PatientID, Timepoint, everything(), -sample)


names(clinDf)[names(clinDf) == "Patient ID"] <- "PatientID"
clin_clean <- clinDf

merged_df <- ms %>%
  left_join(clin_clean, by = "PatientID")

merged_df <- merged_df %>%
  mutate(identifier = paste0(PatientID, "_", ifelse(Timepoint == "Baseline", "B", "S")))


# Now, assign every patient a subtype:
# For each patient, choose the subtype with the highest GSVA score
subtype_cols <- c("EnR","SF","PI","HL")

merged_df <- merged_df %>%
  rowwise() %>%
  mutate(
    Subtype_score = max(c_across(all_of(subtype_cols)), na.rm = TRUE),
    predicted_subtype  = subtype_cols[ which.max(c_across(all_of(subtype_cols))) ]
  ) %>%
  ungroup()

table(merged_df$predicted_subtype)
# EnR  HL  PI  SF 
# 18  34  17  35 

table(merged_df$PatientID, merged_df$predicted_subtype)

# We’ve got the predicted subtypes, and we want to test whether CDK4/6 inhibitor response 
# (based on either Ultrasound Response or Ki67 Response) differs between subtypes, 
# using only baseline samples.


######## filter for only ductals!!!!
merged_df_ductal <- merged_df[merged_df$`Histological Type` == 'Ductal',]
merged_df_ductal$predicted_subtype <- factor(merged_df_ductal$predicted_subtype, levels = c('EnR', 'SF', 'PI', 'HL'))


# --------------------------------------
## ---- USG Response 
# --------------------------------------

baseline_df <- merged_df_ductal %>%
  filter(Timepoint == "Baseline") %>%
  # keep only rows with a predicted_subtype
  filter(!is.na(predicted_subtype))

table(baseline_df$predicted_subtype)
# EnR  HL  PI  SF 
# 12  27   8   8 

run_assoc_test <- function(df, response_col) {
  df_use <- df %>% filter(!is.na(.data[[response_col]]))
  tab <- table(df_use$predicted_subtype, df_use[[response_col]])
  
  if (nrow(tab) < 1 || ncol(tab) < 2) {
    return(list(
      table = tab,
      test  = NA,
      method = "Insufficient categories to test",
      p.value = NA_real_
    ))
  }
  
  # Try chi-square; if any expected < 5, fall back to Fisher
  chi <- suppressWarnings(chisq.test(tab))
  if (any(chi$expected < 5)) {
    ft <- fisher.test(tab)
    return(list(table = tab, test = ft, method = "Fisher's exact test", p.value = ft$p.value))
  } else {
    return(list(table = tab, test = chi, method = "Chi-squared test", p.value = chi$p.value))
  }
}

## ---- Helper: response rate summary per subtype ----
calc_rates <- function(df, response_col) {
  df %>%
    filter(!is.na(.data[[response_col]])) %>%
    mutate(is_resp = .data[[response_col]] == "Responder") %>%
    group_by(predicted_subtype) %>%
    summarise(
      n               = n(),
      responders      = sum(is_resp, na.rm = TRUE),
      non_responders  = sum(!is_resp, na.rm = TRUE),
      response_rate   = responders / (responders + non_responders),
      .groups = "drop"
    ) %>%
    arrange(desc(response_rate))
}

## ---- ANALYSIS: Ultrasound Response ----
ultra_res <- run_assoc_test(baseline_df, "Ultrasound Response")
cat("\n=== Ultrasound Response x Subtype ===\n")
print(ultra_res$table)
cat("\nTest used:", ultra_res$method, " | p-value:", signif(ultra_res$p.value, 4), "\n")

ultra_rates <- calc_rates(baseline_df, "Ultrasound Response")
cat("\nUltrasound response rates by subtype:\n")
print(ultra_rates)

# Perform chi-square test
tbl <- table(baseline_df$predicted_subtype, baseline_df$`Ultrasound Response`)
pval <- chisq.test(tbl)$p.value
pval_text <- paste0("Global p = ", signif(pval, 3))

## plot
response_colors <- c(
  "Responder"     = "#0072B2",  # blue
  "Non responder" = "#D55E00"   # vermillion / orange
)

p_ultra <- ggplot(
  baseline_df %>% filter(!is.na(`Ultrasound Response`)),
  aes(x = predicted_subtype, fill = `Ultrasound Response`)
) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = response_colors) +
  labs(
    title = "Endocrine trx + CDK4/6inhibition\nUltrasound Response in NeoRHEA",
    x = "Spatial subtype",
    y = "Proportion",
    fill = "Ultrasound Response"
  ) +
  theme_minimal(base_size = 16) +
  theme(
    plot.title = element_text(size = 16),
    axis.title  = element_text(size = 16),
    axis.text   = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  ) +
  annotate(
    "text",
    x = Inf, y = 1.05, hjust = 1.1, vjust = 0,
    label = pval_text,
    size = 5
  )

print(p_ultra) ## Fig. 7h




# --------------------------------------
## ---- Ki67 Response 
# --------------------------------------

ki67_res <- run_assoc_test(baseline_df, "Ki67 Response")
cat("\n=== Ki67 Response x Subtype ===\n")
print(ki67_res$table)
cat("\nTest used:", ki67_res$method, " | p-value:", signif(ki67_res$p.value, 4), "\n")

ki67_rates <- calc_rates(baseline_df, "Ki67 Response")
cat("\nKi67 response rates by subtype:\n")
print(ki67_rates)

# Perform chi-square test
tbl <- table(baseline_df$predicted_subtype, baseline_df$`Ki67 Response`)
pval <- chisq.test(tbl)$p.value
pval_text <- paste0("Global p = ", signif(pval, 3))

## plot

ki67_colors <- c(
  "Responder"     = "#009E73",  # bluish green
  "Non responder" = "#CC79A7"   # reddish purple
)


p_ultra <- ggplot(
  baseline_df %>% filter(!is.na(`Ki67 Response`)),
  aes(x = predicted_subtype, fill = `Ki67 Response`)
) +
  scale_fill_manual(values = ki67_colors) +
  geom_bar(position = "fill") +
  labs(
    title = "Endocrine trx + CDK4/6inhibition\nKi67 Response in NeoRHEA",
    x = "Spatial subtype",
    y = "Proportion"
  ) +
  theme_minimal(base_size = 16) +  # 🔹 sets general font size
  theme(
    plot.title = element_text(size = 16),
    axis.title  = element_text(size = 16),
    axis.text   = element_text(size = 14),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  ) +
  annotate(
    "text",
    x = Inf, y = 1.05, hjust = 1.1, vjust = 0,
    label = pval_text,
    size = 5
  )

print(p_ultra) ## Fig. 7i

