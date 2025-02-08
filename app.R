library(shiny)
library(tidyverse)

# Ensure your backend functions are loaded
source("R/backend_2.R")

ui <- fluidPage(
  # Overall title and description
  titlePanel("Immune Stimulation Analysis"),

  # Overall description displayed above the tabs:
  tags$div(
    style = "padding: 10px; background-color: #f0f0f0;",
    tags$p("This Shiny app is designed to analyze immune stimulation data by generating interactive plots. Whole blood samples were stimulated using a variety of agents and degree of stimulation was quantified via several reagents/readouts using phospho-Flow cytometry."),
    tags$p("Use the tabs below to explore stimulation vs. variance, compare stimulation across reagents, and identify correlations between reagents."),
    tags$p("The data was derived from this study by Bjornson-Hooper et al.: ",
      tags$a(href = "https://doi.org/10.3389/fimmu.2022.867015", "https://doi.org/10.3389/fimmu.2022.867015", target = "_blank")
    )
  ),

  # Tabset for different types of plots:
  tabsetPanel(
    # Tab 1: Volcano Plot
    tabPanel("Stimulation vs. Variance",
             sidebarLayout(
               sidebarPanel(
                 # Tab-specific description
                 tags$p("This tab displays a scatterplot that shows the relationship between the median difference from basal (non-stimulated) and median variance for the selected stimulus. This is helpful for identifying stimuli that affect multiple cell types."),
                 selectInput("vol_stimulus", "Select Stimulus:",
                             choices = unique(data_obj$Condition)),
                 actionButton("go_volcano", "Generate Plot"),
                 downloadButton("download_volcano", "Download Volcano Plot (PNG)")
               ),
               mainPanel(
                 plotOutput("volcanoPlot")
               )
             )
    ),
    # Tab 2: Boxplot
    tabPanel("Stimulation Across Reagents",
             sidebarLayout(
               sidebarPanel(
                 tags$p("This tab generates boxplots that compare the stimulation across different reagents for a given cell type. This is useful for determining which reagents show the greatest degree of stimulation and variation across the cohort. If 'Split by Gender' is selected, statistical significance is annotated on the plot."),
                 selectInput("box_cell", "Select Cell Type:",
                             choices = unique(data_obj$population)),
                 selectInput("box_stimulus", "Select Stimulus:",
                             choices = unique(data_obj$Condition)),
                 checkboxInput("box_gender", "Split by Gender", value = TRUE),
                 actionButton("go_box", "Generate Plot"),
                 downloadButton("download_box", "Download Boxplot (PNG)")
               ),
               mainPanel(
                 plotOutput("boxPlot")
               )
             )
    ),
    # Tab 3: Correlation Plot
    tabPanel("Identifying Correlated Reagents",
             sidebarLayout(
               sidebarPanel(
                 tags$p("This tab creates correlation plots between two reagents for a specific cell type and stimulus. Use the selectors below to choose the reagents and see the relationship. Use this to determine if multiple reagents are highly correlated and redundant in this context."),
                 selectInput("corr_cell", "Select Cell Type:",
                             choices = unique(data_obj$population)),
                 selectInput("corr_stimulus", "Select Stimulus:",
                             choices = unique(data_obj$Condition)),
                 selectInput("corr_read1", "Select Reagent 1:",
                             choices = unique(data_obj$reagent)),
                 selectInput("corr_read2", "Select Reagent 2:",
                             choices = unique(data_obj$reagent)),
                 actionButton("go_corr", "Generate Plot"),
                 downloadButton("download_corr", "Download Correlation Plot (PNG)")
               ),
               mainPanel(
                 plotOutput("corrPlot")
               )
             )
    )
  )
)

server <- function(input, output, session) {

  # Volcano Plot: render when the user clicks "Generate Plot"
  observeEvent(input$go_volcano, {
    output$volcanoPlot <- renderPlot({
      p <- plot_median_variance(
        stimulus = input$vol_stimulus,
        data_obj = data_obj,
        save_plot = FALSE
      )
      print(p)
    })
  })

  output$download_volcano <- downloadHandler(
    filename = function() {
      paste0("volcano_", input$vol_stimulus, ".png")
    },
    content = function(file) {
      p <- plot_median_variance(
        stimulus = input$vol_stimulus,
        data_obj = data_obj,
        save_plot = FALSE
      )
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300, device = "png")
    }
  )

  # Boxplot: render when the user clicks "Generate Plot"
  observeEvent(input$go_box, {
    output$boxPlot <- renderPlot({
      p <- generate_boxplot(
        cell_type = input$box_cell,
        stimulus = input$box_stimulus,
        gender = input$box_gender,
        save_plot = FALSE
      )
      print(p)
    })
  })

  output$download_box <- downloadHandler(
    filename = function() {
      paste0("boxplot_", input$box_stimulus, ".png")
    },
    content = function(file) {
      p <- generate_boxplot(
        cell_type = input$box_cell,
        stimulus = input$box_stimulus,
        gender = input$box_gender,
        save_plot = FALSE
      )
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300, device = "png")
    }
  )

  # Correlation Plot: render when the user clicks "Generate Plot"
  observeEvent(input$go_corr, {
    output$corrPlot <- renderPlot({
      p <- correlation_plot(
        cell_type = input$corr_cell,
        stimulus = input$corr_stimulus,
        read1 = input$corr_read1,
        read2 = input$corr_read2,
        save_plot = FALSE
      )
      print(p)
    })
  })

  output$download_corr <- downloadHandler(
    filename = function() {
      paste0("correlation_", input$corr_cell, "_", input$corr_stimulus, ".png")
    },
    content = function(file) {
      p <- correlation_plot(
        cell_type = input$corr_cell,
        stimulus = input$corr_stimulus,
        read1 = input$corr_read1,
        read2 = input$corr_read2,
        save_plot = FALSE
      )
      ggsave(file, plot = p, width = 6, height = 4, dpi = 300, device = "png")
    }
  )
}

shinyApp(ui = ui, server = server)
