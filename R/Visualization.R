#' Plot Annotation
#'
#' A function to plot the annotation of query dataset using UMAP
#' dimensionality reduction with Seurat.
#'
#' @param queryDataset query dataset
#' @param prediction reference mapping result
#'
#' @return a figure of annotation
#' @export
#'
#' @references McInnes, Leland,
#' John Healy, and James Melville. "Umap: Uniform manifold approximation
#' and projection for dimension reduction."
#' arXiv preprint arXiv:1802.03426 (2018).
#'
#' @examples
#' predictions <- searchCell(covid.pbmc3k, cellxgene3k)
#' plotAnnotation(covid.pbmc3k, predictions)
#'
plotAnnotation <- function(queryDataset, prediction){
  # check user input
  if (class(queryDataset) != "Seurat") {
    stop("The queryDataset is not a Seurat object.")
  }
  if (class(prediction) != "MapResult") {
    stop("The prediction is not a MapResult object.")
  }

  # prepare for UMAP
  print("Computing PCA for UMAP")
  queryDataset <- NormalizeData(queryDataset)
  queryDataset <- FindVariableFeatures(queryDataset)
  queryDataset <- ScaleData(queryDataset)
  queryDataset <- RunPCA(queryDataset)
  queryDataset <- FindNeighbors(queryDataset, dims = 1:30)
  # computing UMAP
  print("Computing UMAP")
  queryDataset <- RunUMAP(queryDataset, dims = 1:30,
                          return.model = TRUE)
  # assign predictions
  queryDataset$prediction <- prediction$predicted.id

  print("generating figure")
  plot <- DimPlot(queryDataset, reduction = "umap", group.by = "prediction",
          label = TRUE,
          label.size = 3,
          repel = TRUE) + NoLegend()
  return(plot)
}

# [END]
