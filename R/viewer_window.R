# See https://stackoverflow.com/questions/49470474/saving-r-shiny-app-as-a-function-with-arguments-passed-to-the-shiny-app
# for why launch the app this way
# allows inclusion of arguments and keeps global env clean
launch_app <- function(image_data, anat = NULL){
  shinyOptions(image_data = image_data,
               anat = anat)
  source(system.file("app.R", package = "my_pkg", local = TRUE, chdir = TRUE))$value
}





if (0 == 1){

  image_data <- getShinyOption("image_data")


  ui <- fluidPage(

    title = "Brain Network Editor",

    fluidRow(
      splitLayout(cellWidths = c("33%", "33%", "33%"),
                  plotOutput("sag_plot", click = "plot_click_sag"),
                  plotOutput("cor_plot", click = "plot_click_cor"),
                  plotOutput("axi_plot", click = "plot_click_axi"))
    ),

    fluidRow(

      column(4,
             sliderInput("pos_cutoff_slider", "Positive Cutoff:",
                         min = 0, max = 1.0,
                         value = 0.0, step = 0.01),
             sliderInput("neg_cutoff_slider", "Negative Cutoff:",
                         min = -1.0, max = 0.0,
                         value = 0.0, step = 0.01),
      ),
      column(4,
             selectInput("type","Select Component", choices = c(1:2)),
      ),
      column(4,
             verbatimTextOutput("info"),
      )
    )


  ) # end of ui defn







  #
  # Server
  #


  server <- function(input, output) {

    selected_slices <- reactiveValues(cor_slice = 1,
                                      axi_slice = 1,
                                      sag_slice = 1)

    brain_dims <- reactiveValues(sag_dim = 1,
                                 cor_dim = 1,
                                 axi_dim = 1)

    # Clicks within the sagittal viewer window
    observeEvent(input$plot_click_sag$y, {
      selected_slices$axi_slice <- round(input$plot_click_sag$y)
    })
    observeEvent(input$plot_click_sag$x, {
      selected_slices$cor_slice <- round(input$plot_click_sag$x)
    })

    # Clicks within the coronal viewer window
    observeEvent(input$plot_click_cor$y, {
      selected_slices$axi_slice <- round(input$plot_click_cor$y)
    })
    observeEvent(input$plot_click_cor$x, {
      selected_slices$sag_slice <- round(input$plot_click_cor$x)
    })

    # Clicks within the axial viewer window
    observeEvent(input$plot_click_axi$y, {
      selected_slices$cor_slice <- round(input$plot_click_axi$y)
    })
    observeEvent(input$plot_click_axi$x, {
      selected_slices$sag_slice <- round(input$plot_click_axi$x)
    })

    # Cutoffs for the overlay - control how much of overlay shows up
    pos_cutoff <- reactive({
      input$pos_cutoff_slider
    })
    neg_cutoff <- reactive({
      input$neg_cutoff_slider
    })

    # Load the underlay
    underlay <- reactive({
      file  <- input$load_underlay
      ext   <- tools::file_ext(file$datapath)
      req(file)
      validate(need(ext == "nii", "Please select a Nifti file"))
      anat = readNifti(file$datapath)
      anat = anat / max(anat, na.rm=TRUE)

      brain_dims$sag_dim = dim(anat)[1]
      brain_dims$cor_dim = dim(anat)[2]
      brain_dims$axi_dim = dim(anat)[3]

      selected_slices$sag_slice <- floor(dim(anat)[1] / 2)
      selected_slices$cor_slice <- floor(dim(anat)[2] / 2)
      selected_slices$axi_slice <- floor(dim(anat)[3] / 2)

      melt(anat)
    })




    # # Load the overlay
    # overlay <- reactive({
    #
    #   file  <- input$load_overlay
    #   ext   <- tools::file_ext(file$datapath)
    #   req(file)
    #   validate(need(ext == "nii", "Please select a Nifti file"))
    #   over = readNifti(file$datapath)
    #
    #   if (input$scale_unit == TRUE){
    #     over = over / max(abs(over))
    #   }
    #
    #   melt(over)
    # })


    # TODO - masking rule
    # Create mask from overlay
    overlay_mask <- reactive({

      mask_inds = which(((overlay()$value > pos_cutoff() & overlay()$value >= 0.0) |
                           (overlay()$value < neg_cutoff() & overlay()$value <= 0.0)))

      mask <- overlay()[mask_inds, ]

      # return indices of non-zero values
      mask  %>% rename("overlay_value" = "value")
    })



    sag_underlay <-  reactive({
      underlay() %>% dplyr::filter(Var1 == selected_slices$sag_slice)
    })
    sag_overlay <- reactive({
      overlay() %>% dplyr::filter(Var1 == selected_slices$sag_slice,
                                  (value > pos_cutoff() & value >= 0.0) |
                                  (value < neg_cutoff() & value <= 0.0))
    })

    cor_underlay <-  reactive({
      underlay() %>% dplyr::filter(Var2 == selected_slices$cor_slice )
    })
    cor_overlay <- reactive({
      overlay() %>% dplyr::filter(Var2 == selected_slices$cor_slice, value > pos_cutoff())
    })

    axi_underlay <-  reactive({
      underlay() %>% dplyr::filter(Var3 == selected_slices$axi_slice )
    })
    axi_overlay <- reactive({
      overlay() %>% dplyr::filter(Var3 == selected_slices$axi_slice, value > pos_cutoff())
    })


    # Axial Plot Update
    output$sag_plot <- renderPlot({
      req(sag_underlay)
      plt <- ggplot(sag_underlay(), aes(x = Var2, y = Var3, fill = value )) +
        geom_tile() +
        theme(legend.position="none") +
        scale_fill_gradient(low = "black", high = 'white', na.value = "black")

      req(sag_overlay)
      plt <- plt +
        new_scale("fill") +
        geom_tile(data = sag_overlay(), aes(x = Var2, y = Var3, fill = value), alpha=0.5) +
        scale_fill_gradient(low = "white", high = "red")

      # Highlighting the selected point - only needed if the current ROI
      # is not being drawn, otherwise will be included in ROI
      if (input$draw_ROI == FALSE){
        new_data <- data.frame(x = selected_slices$cor_slice, y = selected_slices$axi_slice, value = 1)
        plt <- plt +
          new_scale("fill") +
          geom_tile(data = new_data, aes(x = x, y = y, fill = value), color = "Yellow") +
          scale_fill_gradient(low = "green", high = "green")
      }

      if (input$draw_ROI == TRUE){
        new_data <- current_ROI()
        plt <- plt +
          new_scale("fill") +
          geom_tile(data = new_data, aes(x = cor, y = axi, fill = value), color = "Yellow") +
          scale_fill_gradient(low = "green", high = "yellow")
      }

      plt
    })

    output$cor_plot <- renderPlot({
      req(cor_underlay)
      plt <- ggplot(cor_underlay(), aes(x = Var1, y =Var3, fill = value )) +
        geom_tile() +
        theme(legend.position="none") +
        new_scale("fill") +
        geom_tile(data = cor_overlay(), aes(x = Var1, y = Var3, fill = value), alpha=0.5) +
        scale_fill_gradient(low = "white", high = "red")

      # Highlighting the selected point - only needed if the current ROI
      # is not being drawn, otherwise will be included in ROI
      if (input$draw_ROI == FALSE){
        new_data <- data.frame(x = selected_slices$sag_slice, y = selected_slices$axi_slice, value = 1)
        plt <- plt +
          new_scale("fill") +
          geom_tile(data = new_data, aes(x = x, y = y, fill = value), color = "Yellow") +
          scale_fill_gradient(low = "green", high = "green")
      }

      if (input$draw_ROI == TRUE){
        new_data <- current_ROI()
        plt <- plt +
          new_scale("fill") +
          geom_tile(data = new_data, aes(x = sag, y = axi, fill = value), color = "Yellow") +
          scale_fill_gradient(low = "green", high = "yellow")
      }

      plt
    })

    output$axi_plot <- renderPlot({
      req(axi_underlay)
      plt <- ggplot(axi_underlay(), aes(x = Var1, y = Var2, fill = value )) + geom_tile() +
        theme(legend.position="none") +
        new_scale("fill") +
        geom_tile(data = axi_overlay(), aes(x = Var1, y = Var2, fill = value), alpha=0.5) +
        scale_fill_gradient(low = "white", high = "red")

      # Highlighting the selected point - only needed if the current ROI
      # is not being drawn, otherwise will be included in ROI
      if (input$draw_ROI == FALSE){
        new_data <- data.frame(x = selected_slices$sag_slice, y = selected_slices$cor_slice, value = 1)
        plt <- plt +
          new_scale("fill") +
          geom_tile(data = new_data, aes(x = x, y = y, fill = value), color = "Yellow") +
          scale_fill_gradient(low = "green", high = "green")
      }

      if (input$draw_ROI == TRUE){
        new_data <- current_ROI()
        plt <- plt +
          new_scale("fill") +
          geom_tile(data = new_data, aes(x = sag, y = cor, fill = value), color = "Yellow") +
          scale_fill_gradient(low = "green", high = "yellow")
      }

      plt
    })

    output$info <- renderText({
      c(
        paste0("Voxel index: (", selected_slices$sag_slice, ", ",
               selected_slices$cor_slice, ", ",
               selected_slices$axi_slice, ")\n"),
        paste0("Overlay value at voxel index: ", overlay() %>% filter(Var1 == selected_slices$sag_slice,
                                                                      Var2 == selected_slices$cor_slice,
                                                                      Var3 == selected_slices$axi_slice) %>%
                 pull(value))
      )
    })

    # Save the output
    observeEvent(input$save_brain_map, {
      if(nrow(ROI_list$data) > 0){

        # Combine all ROIs
        all_vox <- ROI_list$data %>% unnest(coords)

        new_img <- array(0, dim= c(brain_dims$sag_dim, brain_dims$cor_dim, brain_dims$axi_dim ))
        for (i in 1:nrow(all_vox)){
          new_img[all_vox$sag[i], all_vox$cor[i], all_vox$axi[i]] =
            new_img[all_vox$sag[i], all_vox$cor[i], all_vox$axi[i]] + all_vox$value[i]
        }

        new_nii = asNifti(new_img)
        writeNifti(new_nii, input$save_name)

      }
    })

  }

  shinyApp(ui, server)
















}
