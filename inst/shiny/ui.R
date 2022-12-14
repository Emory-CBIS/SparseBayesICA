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
           selectInput("select_component_dropdown", "Select Component", choices = c(1:Q)),
           DT::dataTableOutput("covariateList")
    ),
    column(4,
           verbatimTextOutput("info"),
    )
  )


) # end of ui defn
