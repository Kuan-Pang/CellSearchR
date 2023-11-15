library(CellSearchR)


test_that("searchCell for identifcal dataset", {
  
  query <- cellxgene3k
  reference <- cellxgene3k
  
  predictions <- searchCell(query, reference, "Seurat")
  
  expect_type(predictions, "list")
  expect_s3_class(predictions, "MapResult")
  expect_length(predictions, 82)
})

test_that("searchCell error upon invalid  dataset input", {
  
  query <- cellxgene3k
  reference <- cellxgene3k
  
  expect_error(searchCell(query, reference, "A"))
  expect_error(searchCell(1, reference, "Seurat"))
  expect_error(searchCell(query, 1, "Seurat"))
})


test_that("searchCell for different dataset", {
  
  query <- covid.pbmc3k
  reference <- cellxgene3k
  
  predictions <- searchCell(query, reference, "Seurat")
  
  expect_type(predictions, "list")
  expect_s3_class(predictions, "MapResult")
  expect_length(predictions, 82)
})

# [END]