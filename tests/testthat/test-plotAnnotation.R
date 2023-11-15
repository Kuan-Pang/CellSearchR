library(CellSearchR)


test_that("plotAnnotation input checks", {
  
  query <- cellxgene3k
  predictions <- searchCell(query, query, "Seurat")
  
  expect_error(plotAnnotation(1, predictions))
  expect_error(plotAnnotation(query, 1))
})

test_that("plotAnnotation output checks", {
  
  query <- cellxgene3k
  predictions <- searchCell(query, query, "Seurat")
  
  expect_type(plotAnnotation(query, predictions), "list")
})

# [END]