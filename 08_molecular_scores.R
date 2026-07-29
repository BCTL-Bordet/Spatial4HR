# 08_molecular_scores.R


################################################################################
## Molecular subtyping and genomic risk scoring
################################################################################

library(AIMS)
library(genefu)
library(org.Hs.eg.db)
library(AnnotationDbi)

################################################################################
## 1. PAM50 intrinsic molecular subtyping using AIMS
################################################################################

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
  keytype = "SYMBOL"
)

## Remove duplicated gene-symbol mappings
annotation_not_dup <- annotation[
  !duplicated(annotation$SYMBOL),
]

## Run AIMS
matrice2 <- as.matrix(expression_matrix)

pam50_aims <- applyAIMS(
  matrice2,
  annotation_not_dup$ENTREZID
)

pam50class <- pam50_aims$cl

## Inspect subtype distribution
table(pam50class)

## Save PAM50 classifications




################################################################################
## 2. Genomic risk scoring using genefu
################################################################################

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

################################################################################
## 3. Genomic Grade Index
################################################################################

data("sig.ggi")

ggi_calc <- ggi(
  data = t(as.matrix(counts_ductal)),
  annot = annotation,
  do.mapping = TRUE
)

score_ggi_genefu <- ggi_calc$score

################################################################################
## 4. MammaPrint 70-gene score and risk classification
################################################################################

data("sig.gene70")

gene70_res <- gene70(
  data = t(as.matrix(counts_ductal)),
  annot = annotation,
  do.mapping = TRUE
)

score_gene70 <- gene70_res$score
risk_gene70 <- gene70_res$risk

################################################################################
## 5. Oncotype DX recurrence score and risk classification
################################################################################

data("sig.oncotypedx")

oncotype_res <- oncotypedx(
  data = t(as.matrix(counts_ductal)),
  annot = annotation,
  do.mapping = TRUE
)

score_oncotypedx <- oncotype_res$score
risk_oncotypedx <- oncotype_res$risk

################################################################################
## 6. Inspect outputs
################################################################################

score_ggi_genefu
score_gene70
risk_gene70
score_oncotypedx
risk_oncotypedx

################################################################################
## 7. Save score objects
################################################################################

