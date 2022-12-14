server <- function(input, output) {

  # input
  image_data <- getShinyOption("image_data")
  anat       <- getShinyOption("anat")

  color_range_1 <- colorRampPalette(c("blue", "white", "red"))

  # Format the covariate data
  all_colnames <- colnames(image_data)
  P <- (length(all_colnames) - 5) / 2
  cov_colnames <- all_colnames[6:(6 + P - 1)]

  sel <- rep(0, length(cov_colnames))
  tab <- tibble(cov_colnames, sel)

  output$covariateList <- DT::renderDataTable(tab,
                                        options = list(scrollX = TRUE),
                                        rownames = FALSE)

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
    brain_dims$sag_dim = max(anat[,1])
    brain_dims$cor_dim = max(anat[,2])
    brain_dims$axi_dim = max(anat[,3])

    selected_slices$sag_slice <- floor(brain_dims$sag_dim / 2)
    selected_slices$cor_slice <- floor(brain_dims$cor_dim / 2)
    selected_slices$axi_slice <- floor(brain_dims$axi_dim / 2)

    anat
  })




  selected_S0_colormap <- reactive({
    scale_fill_gradient(low = "black", high = 'white', na.value = "black")
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

  overlay <- reactive({

    image_data %>%
      filter(Q == as.numeric(input$select_component_dropdown)) %>%
      mutate(value = S0)

  })

  # additional_overlays <- reactive({
  #
  #   image_data %>%
  #     filter(Q == 1) %>%
  #     mutate(value = S0)
  #
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
    overlay() %>% dplyr::filter(sag == selected_slices$sag_slice,
                                (value > pos_cutoff() & value >= 0.0) |
                                  (value < neg_cutoff() & value <= 0.0))
  })

  cor_underlay <-  reactive({
    underlay() %>% dplyr::filter(Var2 == selected_slices$cor_slice )
  })
  cor_overlay <- reactive({
    overlay() %>% dplyr::filter(cor == selected_slices$cor_slice,
                                (value > pos_cutoff() & value >= 0.0) |
                                  (value < neg_cutoff() & value <= 0.0))
  })

  axi_underlay <-  reactive({
    underlay() %>% dplyr::filter(Var3 == selected_slices$axi_slice )
  })
  axi_overlay <- reactive({
    overlay() %>% dplyr::filter(axi == selected_slices$axi_slice,
                                (value > pos_cutoff() & value >= 0.0) |
                                  (value < neg_cutoff() & value <= 0.0))
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
      geom_tile(data = sag_overlay(), aes(x = cor, y = axi, fill = value), alpha=0.5) +
      scale_fill_gradientn(colours=color_range_1(100),
                           breaks=c(-1.0, 0.0, 1.0))

    # Highlighting the selected point - only needed if the current ROI
    # is not being drawn, otherwise will be included in ROI
    # if (input$draw_ROI == FALSE){
    #   new_data <- data.frame(x = selected_slices$cor_slice, y = selected_slices$axi_slice, value = 1)
    #   plt <- plt +
    #     new_scale("fill") +
    #     geom_tile(data = new_data, aes(x = x, y = y, fill = value), color = "Yellow") +
    #     scale_fill_gradient(low = "green", high = "green")
    # }

    plt
  })

  output$cor_plot <- renderPlot({
    req(cor_underlay)
    plt <- ggplot(cor_underlay(), aes(x = Var1, y =Var3, fill = value )) +
      geom_tile() +
      theme(legend.position="none") +
      scale_fill_gradient(low = "black", high = 'white', na.value = "black")

    plt <- plt +
      new_scale("fill") +
      geom_tile(data = cor_overlay(), aes(x = sag, y = axi, fill = value), alpha=0.5) +
      scale_fill_gradient(low = "white", high = "red")

    # Highlighting the selected point - only needed if the current ROI
    # is not being drawn, otherwise will be included in ROI
    # if (input$draw_ROI == FALSE){
    #   new_data <- data.frame(x = selected_slices$sag_slice, y = selected_slices$axi_slice, value = 1)
    #   plt <- plt +
    #     new_scale("fill") +
    #     geom_tile(data = new_data, aes(x = x, y = y, fill = value), color = "Yellow") +
    #     scale_fill_gradient(low = "green", high = "green")
    # }

    plt
  })

  output$axi_plot <- renderPlot({
    req(axi_underlay)
    plt <- ggplot(axi_underlay(), aes(x = Var1, y = Var2, fill = value )) +
      geom_tile() +
      theme(legend.position="none") +
      scale_fill_gradient(low = "black", high = 'white', na.value = "black")

    plt <- plt +
      new_scale("fill") +
      geom_tile(data = axi_overlay(), aes(x = sag, y = cor, fill = value), alpha=0.5) +
      scale_fill_gradient(low = "white", high = "red")

    # Highlighting the selected point - only needed if the current ROI
    # is not being drawn, otherwise will be included in ROI
    # if (input$draw_ROI == FALSE){
    #   new_data <- data.frame(x = selected_slices$sag_slice, y = selected_slices$cor_slice, value = 1)
    #   plt <- plt +
    #     new_scale("fill") +
    #     geom_tile(data = new_data, aes(x = x, y = y, fill = value), color = "Yellow") +
    #     scale_fill_gradient(low = "green", high = "green")
    # }

    plt
  })



}
