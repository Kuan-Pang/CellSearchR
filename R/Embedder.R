
#' Generate Seurat embeddings for both reference and query dataset
#'
#' A utility function to generate Seurat embeddings
#' for both reference and query dataset.
#'
#' Seurat embeddings uses PCA-based approach across both reference and query
#' dataset.
#'
#' @import Seurat
#' @param referenceDataset preprocessed reference dataset
#' @param queryDataset preprocessed query dataset
#'
#' @return Seurat embeddings
#'
seuratEmbedding <- function(referenceDataset, queryDataset){
  # return the computed PCA embedding with Seurat flavor
  return(FindTransferAnchors(
    reference = referenceDataset,
    query = queryDataset,
    dims = 1:30,
    reference.reduction = "pca",
    normalization.method = "SCT"
  ))
}
# [END]
