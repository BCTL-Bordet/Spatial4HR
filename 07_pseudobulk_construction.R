#07_pseudobulk_construction.R


################################################################################
## Whole-section, tumor and stromal pseudobulk construction
##
## Whole-section pseudobulk:
##   All retained Visium spots within each tumor
##
## Tumor pseudobulk:
##   Spots with >50% of their area annotated as invasive tumor
##
## Stromal pseudobulk:
##   Spots with 0% annotated tumor area
################################################################################

library(Seurat)
library(SeuratObject)
library(Matrix)
library(DESeq2)
library(dplyr)
library(tibble)

################################################################################
## 1. Settings
################################################################################

setwd("~/Desktop/final_scripts/00_data/deneme")

input_file <- "~/Desktop/final_scripts/00_data/object_merged_final.RDS"
output_dir <- "results/pseudobulk"

sample_column <- "orig.ident"
tumor_fraction_column <- "Tumor"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

################################################################################
## 2. Load merged spatial transcriptomics object
################################################################################

object_merged <- readRDS(input_file)

DefaultAssay(object_merged) <- "Spatial"

## Raw UMI count matrix: genes × spots
counts_mat <- GetAssayData(
  object = object_merged,
  assay = "Spatial",
  layer = "counts"
)

meta_data <- object_merged@meta.data

################################################################################
## 3. Confirm alignment and required metadata
################################################################################

required_columns <- c(
  sample_column,
  tumor_fraction_column
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta_data)
)

if (length(missing_columns) > 0) {
  stop(
    "Missing required metadata column(s): ",
    paste(missing_columns, collapse = ", ")
  )
}

if (!identical(
  colnames(counts_mat),
  rownames(meta_data)
)) {
  stop(
    "Spot order differs between the raw count matrix and metadata."
  )
}

if (anyNA(meta_data[[sample_column]])) {
  stop("Missing tumor/sample identifiers were detected.")
}

if (anyNA(meta_data[[tumor_fraction_column]])) {
  stop("Missing tumor-fraction values were detected.")
}

if (
  any(
    meta_data[[tumor_fraction_column]] < 0 |
    meta_data[[tumor_fraction_column]] > 1
  )
) {
  stop("Tumor fractions must lie between 0 and 1.")
}

meta_data$sample_id <- as.character(
  meta_data[[sample_column]]
)

meta_data$tumor_fraction <- meta_data[[tumor_fraction_column]]

samples <- unique(meta_data$sample_id)

################################################################################
## 4. Function to aggregate raw counts by tumor/sample
################################################################################

aggregate_pseudobulk <- function(
    count_matrix,
    metadata,
    selected_spots,
    sample_ids
) {
  
  selected_spots <- intersect(
    selected_spots,
    colnames(count_matrix)
  )
  
  selected_metadata <- metadata[
    selected_spots,
    ,
    drop = FALSE
  ]
  
  pseudobulk_matrix <- matrix(
    0,
    nrow = nrow(count_matrix),
    ncol = length(sample_ids),
    dimnames = list(
      rownames(count_matrix),
      sample_ids
    )
  )
  
  spot_summary <- data.frame(
    sample_id = sample_ids,
    n_spots = 0,
    has_pseudobulk = FALSE
  )
  
  for (sample_id in sample_ids) {
    
    sample_spots <- rownames(
      selected_metadata
    )[
      selected_metadata$sample_id == sample_id
    ]
    
    spot_summary$n_spots[
      spot_summary$sample_id == sample_id
    ] <- length(sample_spots)
    
    if (length(sample_spots) > 0) {
      
      pseudobulk_matrix[, sample_id] <- Matrix::rowSums(
        count_matrix[
          ,
          sample_spots,
          drop = FALSE
        ]
      )
      
      spot_summary$has_pseudobulk[
        spot_summary$sample_id == sample_id
      ] <- TRUE
    }
  }
  
  list(
    counts = pseudobulk_matrix,
    spot_summary = spot_summary
  )
}

################################################################################
## 5. Whole-section pseudobulk
################################################################################

whole_section_spots <- rownames(meta_data)

whole_section_result <- aggregate_pseudobulk(
  count_matrix = counts_mat,
  metadata = meta_data,
  selected_spots = whole_section_spots,
  sample_ids = samples
)

whole_section_counts <- whole_section_result$counts
whole_section_spot_summary <- whole_section_result$spot_summary

################################################################################
## 6. Tumor pseudobulk
##
## Tumor-rich spots contain >50% annotated invasive tumor area
################################################################################

tumor_spots <- rownames(
  meta_data
)[
  meta_data$tumor_fraction > 0.50
]

tumor_result <- aggregate_pseudobulk(
  count_matrix = counts_mat,
  metadata = meta_data,
  selected_spots = tumor_spots,
  sample_ids = samples
)

tumor_counts <- tumor_result$counts
tumor_spot_summary <- tumor_result$spot_summary

################################################################################
## 7. Stromal pseudobulk
##
## Stromal spots contain no annotated invasive tumor area
################################################################################

stromal_spots <- rownames(
  meta_data
)[
  meta_data$tumor_fraction == 0
]

stromal_result <- aggregate_pseudobulk(
  count_matrix = counts_mat,
  metadata = meta_data,
  selected_spots = stromal_spots,
  sample_ids = samples
)

stromal_counts <- stromal_result$counts
stromal_spot_summary <- stromal_result$spot_summary

################################################################################
## 8. Remove tumors without eligible spots from compartment matrices
################################################################################

tumor_samples_to_keep <- tumor_spot_summary$sample_id[
  tumor_spot_summary$has_pseudobulk
]

stromal_samples_to_keep <- stromal_spot_summary$sample_id[
  stromal_spot_summary$has_pseudobulk
]

tumor_counts_filtered <- tumor_counts[
  ,
  tumor_samples_to_keep,
  drop = FALSE
]

stromal_counts_filtered <- stromal_counts[
  ,
  stromal_samples_to_keep,
  drop = FALSE
]

################################################################################
## 9. Whole-section library-size normalization
##
## Counts per million followed by log2 transformation
################################################################################

whole_section_library_sizes <- colSums(
  whole_section_counts
)

if (any(whole_section_library_sizes == 0)) {
  stop(
    "At least one whole-section pseudobulk has zero total counts."
  )
}

whole_section_cpm <- sweep(
  whole_section_counts,
  2,
  whole_section_library_sizes,
  FUN = "/"
) * 1e6

whole_section_log2_cpm <- log2(
  whole_section_cpm + 1
)

################################################################################
## 10. Variance-stabilizing transformation of compartment pseudobulks
################################################################################

## Weighted or fractional counts were not used.
## Counts are integer sums of raw UMIs.

tumor_counts_filtered <- round(
  tumor_counts_filtered
)

stromal_counts_filtered <- round(
  stromal_counts_filtered
)

tumor_vst <- varianceStabilizingTransformation(
  tumor_counts_filtered,
  blind = TRUE
)

stromal_vst <- varianceStabilizingTransformation(
  stromal_counts_filtered,
  blind = TRUE
)

################################################################################
## 11. Summary table
################################################################################

pseudobulk_summary <- whole_section_spot_summary %>%
  select(
    sample_id,
    whole_section_spots = n_spots
  ) %>%
  left_join(
    tumor_spot_summary %>%
      select(
        sample_id,
        tumor_spots = n_spots,
        tumor_pseudobulk_available = has_pseudobulk
      ),
    by = "sample_id"
  ) %>%
  left_join(
    stromal_spot_summary %>%
      select(
        sample_id,
        stromal_spots = n_spots,
        stromal_pseudobulk_available = has_pseudobulk
      ),
    by = "sample_id"
  )

print(
  pseudobulk_summary[
    !pseudobulk_summary$stromal_pseudobulk_available,
  ]
)

