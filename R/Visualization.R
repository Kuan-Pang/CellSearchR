

plotAnnotation <- function(queryDataset, prediction){
  queryDataset <- RunUMAP(queryDataset, dims = 1:30, reduction = "integrated.cca", return.model = TRUE)
  queryDataset$prediction <- prediction
  DimPlot(queryDataset, reduction = "umap", group.by = "prediction", 
          label = TRUE, 
          label.size = 3,
          repel = TRUE) 
          + NoLegend() 
           + ggtitle("Query Annotations")
}