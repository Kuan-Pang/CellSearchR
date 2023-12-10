#' Launch Shiny App for CellSearchR
#'
#' A function that launches the Shiny app for CellSearchR
#' This app allows users to upload a query dataset and choose a reference
#' datset for reference mapping.
#'
#' @return No return value but open up a Shiny page.
#'
#' @examples
#' \dontrun{
#'
#' CellSearchR::runCellSearchR()
#' }
#'
#' @references
#' Grolemund, G. (2015). Learn Shiny - Video Tutorials.
#' \href{https://shiny.rstudio.com/tutorial/}{Link}
#'
#' @export
#' @importFrom shiny runApp
runCellSearchR <- function() {
  appDir <- system.file("shiny-scripts",
                        package = "cellSearchR")
  actionShiny <- shiny::runApp(appDir, display.mode = "normal")
  return(actionShiny)
}
# [END]
