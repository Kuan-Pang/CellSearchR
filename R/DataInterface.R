#' Build the reference dataset
#'
#' A function to build the reference dataset given the dataset name at cache directory
#'
#' The supported dataset names are:
#'  - "cellxgene3k": the default dataset with 10k cells from cellxgene blood collection
#'
#'
#' @param datasetName A string of the dataset name
#' @return Dataset Seurat object of the reference dataset
#'
#' @import Seurat
#' @import SeuratData
#' @export
#'
#'
#' @examples
#' # Example1: Install the default dataset (cellxgene3k)
#' InstallDataset("cellxgene3k")
#'
InstallDataset <- function(datasetName = "cellxgene3k") {
  if (datasetName == "cellxgene3k") {
    # this dataset is downloaded with the pkg
    dataset <- LoadH5Seurat("data/cellxgene3k.h5Seurat")
  } else if (datasetName == "pbmc3k") {
    # load the dataset using SeuratData
    InstallData("pbmc3k", force.reinstall = TRUE)
    dataset <- LoadData("pbmc3k")
  } else{
    stop("The dataset is not supported.")
  }
  return(dataset)
}
