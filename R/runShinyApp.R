#' @title Run Shiny app
#' @description Runs the event prediction Shiny app.
#'
#' @author Kaifeng Lu, \email{kaifenglu@@gmail.com}
#'
#' @export
runShinyApp <- function() {
  shiny::shinyAppDir(system.file("shinyApp", package = "eventPred"))
}
