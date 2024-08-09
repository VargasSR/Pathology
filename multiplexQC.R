library(shiny)
library(dplyr)
library(ggplot2)
library(data.table)
library(openxlsx)

# Define UI
ui <- fluidPage(
  titlePanel("Multiplex Data Analysis Tool"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload File"),
      actionButton("processData", "Read File"),
      selectInput("selectBiomarker", "Select Biomarker:", choices = NULL),
      numericInput("pickThreshold","Threshold", value = 200, min = 0, max = 36000),
      actionButton("updateResults", "Update Results"),
      textOutput("subheadingtext")
    
    ),
    mainPanel(
      plotOutput("positive_intensity_plot"),
      textOutput("thresholdText"),
      textOutput("median_pos_intensity"),
      textOutput("mean_pos_intensity")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 200 * 1024^2)
  
  output$subheadingtext <- renderText("Export cleaned and intensity thresholded data for the selected biomarker")
  
  # Read data from uploaded file
  data <- reactive({
    req(input$file)
    read.table(input$file$datapath, sep = "\t", header = TRUE)
  })
  
  processedData <- eventReactive(input$processData, {
    req(data())
    
    ignore_headers <- c("RED", "GREEN", "BLUE", "DAPI")
    channels_to_ignore <- sapply(ignore_headers, paste, collapse = "|")
    
    # Select relevant means and parent columns
    selected_columns <- data() %>%
      select(matches("mean|parent|object.ID", ignore.case = TRUE)) %>%
      select(-matches(channels_to_ignore))
    
    cell_means <- selected_columns %>%
      select(matches("cell|parent|object.ID", ignore.case = TRUE))
    
    colnames(cell_means) <- sub("\\.\\..", "", colnames(cell_means))
    
    # Extract list of unique column names (biomarkers), excluding parent
    unique_headers <- colnames(cell_means)
    biomarkers <- unique_headers[-grep("parent|object.ID", unique_headers, ignore.case = TRUE)]
    
    list(data = data(), biomarkers = biomarkers, cell_means = cell_means)
  })
  
  # Update the selectInput choices when data is processed
  observeEvent(processedData(), {
    updateSelectInput(
      session,
      "selectBiomarker",
      choices = processedData()$biomarkers
    )
  })
  
  selectedBiomarkerData <- eventReactive(input$updateResults, {
    req(processedData()$cell_means, input$pickThreshold, input$selectBiomarker)
    
    cell_means <- processedData()$cell_means
    
    # Create data frame for selected biomarker + Parent column
    selected_biomarker_means <- cell_means[, c("Object.ID", "Parent", input$selectBiomarker), drop = FALSE]
    
    # Calculate 95% quantile
    quantile95 <- quantile(selected_biomarker_means[[input$selectBiomarker]], probs = 0.95)
    
    # Make quantile bindex
    quantile_bindex <- selected_biomarker_means[[input$selectBiomarker]] < quantile95
    
    # Filter by quantile bindex
    filtered_means <- subset(selected_biomarker_means, quantile_bindex)
    
    # Make threshold bindex
    threshold_bindex <- selected_biomarker_means[[input$selectBiomarker]] > input$pickThreshold
    
    # Filter data by threshold
    thresholded_means <- subset(selected_biomarker_means, threshold_bindex)
    
    # Calculate average of data above threshold
    mean_positive_intensity <- round(mean(thresholded_means[[input$selectBiomarker]]), 2)
    
    # Calculate median of data above threshold
    median_positive_intensity <- round(median(thresholded_means[[input$selectBiomarker]]), 2)
    
    list(
      filtered_means = filtered_means,
      mean_positive_intensity = mean_positive_intensity,
      median_positive_intensity = median_positive_intensity
    )
  })
  
  output$positive_intensity_plot <- renderPlot({
    req(input$updateResults)
    
    plot <- ggplot(selectedBiomarkerData()$filtered_means, aes(x = .data[[input$selectBiomarker]], fill = Parent, colour = Parent)) +
      geom_histogram(aes(x = .data[[input$selectBiomarker]]), alpha = 0.5, position = "identity", bins = 50) +
      geom_vline(aes(xintercept = input$pickThreshold), color = "black") +
      geom_vline(aes(xintercept = selectedBiomarkerData()$mean_positive_intensity), color = "Red") +
      geom_vline(aes(xintercept = selectedBiomarkerData()$median_positive_intensity), color = "Blue") +
      labs(x = "Mean Intensity per cell", title = input$selectBiomarker)
    
    plot
  })
  
  output$thresholdText <- renderText({
    paste("Threshold: ", input$pickThreshold)
  })
  
  output$median_pos_intensity <- renderText({
    paste("Median: ", selectedBiomarkerData()$median_positive_intensity)
  })
  
  output$mean_pos_intensity <- renderText({
    paste("Mean: ", selectedBiomarkerData()$mean_positive_intensity)
  })
  
  
}

# Run the app
shinyApp(ui = ui, server = server)
