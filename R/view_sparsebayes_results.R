view_sparsebayes_results <- function(model_estimates, anat) {

  # check class of input

  uifl <- system.file("shiny", "ui.R", package = "hcicaR")
  serverfl <- system.file("shiny", "server.R", package = "hcicaR")
  appDir <- system.file("shiny", "app.R", package = "hcicaR")

  if (appDir == "") {
    stop("Could not find shiny app directory. Try re-installing `hcicaR`.", call. = FALSE)
  }


  shinyOptions(image_data = model_estimates,
               anat = anat)

  #source(system.file("app.R", package = "my_pkg", local = TRUE, chdir = TRUE))$value
  source(uifl)
  source(serverfl)
  source(appDir)$value

  #shiny::runApp(appDir, display.mode = "normal")

}
