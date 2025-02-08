library(shiny)
library(tidyverse)

# Ensure your backend functions are loaded
source("R/backend_2.R")

ui <- fluidPage(
  titlePanel("Immune Stimulation Analysis"),
  tabsetPanel(
    # Tab 1: Volcano Plot
    tabPanel("Stimulation vs. variance",
             sidebarLayout(
               sidebarPanel(
                 selectInput("vol_stimulus", "Select Stimulus:",
                             choices = unique(data_obj$Condition)),
                 actionButton("go_volcano", "Generate Plot"),
                 # Download button for volcano plot
                 downloadButton("download_volcano", "Save")
               ),
               mainPanel(
                 plotOutput("volcanoPlot")
               )
             )
    ),
    # Tab 2: Boxplot
    tabPanel("Stimulation across reagents",
             sidebarLayout(
               sidebarPanel(
                 selectInput("box_cell", "Select Cell Type:",
                             choices = unique(data_obj$population)),
                 selectInput("box_stimulus", "Select Stimulus:",
                             choices = unique(data_obj$Condition)),
                 checkboxInput("box_gender", "Split by Gender", value = TRUE),
                 actionButton("go_box", "Generate Plot"),
                 # Download button for boxplot
                 downloadButton("download_box", "Save")
               ),
               mainPanel(
                 plotOutput("boxPlot")
               )
             )
    ),
    # Tab 3: Correlation Plot
    tabPanel("Identifying correlated reagents",
             sidebarLayout(
               sidebarPanel(
                 selectInput("corr_cell", "Select Cell Type:",
                             choices = unique(data_obj$population)),
                 selectInput("corr_stimulus", "Select Stimulus:",
                             choices = unique(data_obj$Condition)),
                 selectInput("corr_read1", "Select Reagent 1:",
                             choices = unique(data_obj$reagent)),
                 selectInput("corr_read2", "Select Reagent 2:",
                             choices = unique(data_obj$reagent)),
                 actionButton("go_corr", "Generate Plot"),
                 # Download button for correlation plot
                 downloadButton("download_corr", "Save")
               ),
               mainPanel(
                 plotOutput("corrPlot")
               )
             )
    )
  )
)

server <- function(input, output, session) {

  # Volcano Plot generation and download
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
      paste0("volcano_", input$vol_stimulus, ".pdf")
    },
    content = function(file) {
      # Generate the plot (without saving to disk in the function)
      p <- plot_median_variance(
        stimulus = input$vol_stimulus,
        data_obj = data_obj,
        save_plot = FALSE
      )
      # Save the plot to the temporary file in pdf format;
      # adjust width, height, and dpi as desired.
      ggsave(file, plot = p, width = 6, height = 4, device = "pdf")
    }
  )

  # Boxplot generation and download
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
      paste0("boxplot_", input$box_stimulus, ".pdf")
    },
    content = function(file) {
      p <- generate_boxplot(
        cell_type = input$box_cell,
        stimulus = input$box_stimulus,
        gender = input$box_gender,
        save_plot = FALSE
      )
      ggsave(file, plot = p, width = 10, height = 8, device = "pdf")
    }
  )

  # Correlation Plot generation and download
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
      paste0("correlation_", input$corr_cell, "_", input$corr_stimulus, ".pdf")
    },
    content = function(file) {
      p <- correlation_plot(
        cell_type = input$corr_cell,
        stimulus = input$corr_stimulus,
        read1 = input$corr_read1,
        read2 = input$corr_read2,
        save_plot = FALSE
      )
      ggsave(file, plot = p, width = 6, height = 4, device = "pdf")
    }
  )
}

shinyApp(ui = ui, server = server)
