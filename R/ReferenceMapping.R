#' Perform reference mapping
#'
#' A function to perform reference mapping, which is to map the query dataset to
#' the reference dataset with specified embedding method.
#'
#' @param queryDataset A Seurat object of the query dataset
#' @param referenceDataset A Seurat object of the reference dataset;
#' @param embedMethod A string of the embedding method,
#'  which is one of the following:
#'    - "Seurat": the default embedding method using Seurat
#' 
#' @return predictions of reference mapping
#' @export
#'
#' @examples
#' predictions <- searchCell(covid.pbmc3k, cellxgene3k,
#'                            embedMethod = "Seurat")
#'
#' @references
#' Hao Y, Hao S, Andersen-Nissen E, III WMM, Zheng S, Butler A, Lee MJ, Wilk AJ,
#'  Darby C, Zagar M, Hoffman P, Stoeckius M, Papalexi E, Mimitou EP, Jain J,
#' Srivastava A, Stuart T, Fleming LB, Yeung B, Rogers AJ, McElrath JM,
#' Blish CA, Gottardo R, Smibert P, Satija R (2021). “Integrated analysis of
#' multimodal single-cell data.” Cell. doi:10.1016/j.cell.2021.04.048,
#' https://doi.org/10.1016/j.cell.2021.04.048.
#'
searchCell <- function(queryDataset,
                       referenceDataset,
                       embedMethod = "Seurat") {
  # check the user input
  if (class(queryDataset) != "Seurat") {
    stop("The queryDataset is not a Seurat object.")
  }
  if (class(referenceDataset) != "Seurat") {
    stop("The referenceDataset is not a Seurat object.")
  }
  
  # select the embedding method
  if (embedMethod == "Seurat") {
    predictions <-
      referenceMappingSeurat(queryDataset, referenceDataset)
  } else{
    stop("The embedMethod is not supported.")
  }
  return(predictions)
}

#' Perform reference mapping with Seurat
#' 
#' A function to perform reference mapping with Seurat.
#'
#' @param queryDataset 
#' @param referenceDataset 
#' @import Seurat
#'
#' @return predictions of reference mapping
#' 
#' @examples
#' predictions <- referenceMappingSeurat(covid.covid.pbmc3k, cellxgene3k)
#'
referenceMappingSeurat <- function(queryDataset,
                                   referenceDataset) {
  # dataset preprocessing
  print("Start reference mapping with Seurat")
  print("Preprocessing the reference dataset")
  referenceDataset <- NormalizeData(referenceDataset)
  referenceDataset <- FindVariableFeatures(referenceDataset)
  referenceDataset <- ScaleData(referenceDataset)
  referenceDataset <- RunPCA(referenceDataset)
  referenceDataset <- FindNeighbors(referenceDataset, dims = 1:30)
  referenceDataset <- FindClusters(referenceDataset)
  
  # dataset integration
  print("Integrating the reference dataset")
  referenceDataset <- SCTransform(referenceDataset)
  
  # query dataset
  print("Preprocessing the query dataset")
  queryDataset <- NormalizeData(queryDataset)
  print("Running reference mapping")
  queryAnchers <- seuratEmbedding(referenceDataset, queryDataset)
  print("Transfer annotation")
  predictions <- TransferData(
    anchorset = queryAnchers,
    refdata = referenceDataset$cell_type,
    dims = 1:30
  )
  class(predictions) <- "MapResult"
  return(predictions)
}

# [END]