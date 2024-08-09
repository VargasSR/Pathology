library(shiny)
library(dplyr)
library(ggplot2)
library(data.table)
library(openxlsx)

linebreaks <- function(n){HTML(strrep(br(), n))}

# Define UI
ui <- fluidPage(
  titlePanel("Multiplex Data QC and Filtering Tool"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Select File"),
      actionButton("processData", "Read Data"),
      selectInput("selectBiomarker", "Select Biomarker:", choices = NULL),
      numericInput("pickThreshold","Select Positive Intensity Threshold", value = 0, min = 0, max = 36000),
      sliderInput("pickQuantile", "Data View Cut-off", min = .90, max = 1, value = 1, step = 0.01),
      textOutput("subheadingtext1"),
      tags$style("#subheadingtext1{color: black;
font-size: 9px;
}"
      ),
      p( linebreaks(1)),
      actionButton("updateResults", "Update and View Plot"),
    ),
    mainPanel(
      plotOutput("positive_intensity_plot"),
      textOutput("thresholdText"),
      tags$style("#thresholdText{color: black;
font-size: 20px;
}"
      ),
      textOutput("medianText"),
      tags$style("#medianText{color: blue;
font-size: 20px;
}"
      ),
      textOutput("meanText"),
      tags$style("#meanText{color: red;

font-size: 20px;
}"
      ),
      tableOutput("stats_table")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  options(shiny.maxRequestSize = 200 * 1024^2)
  
  output$subheadingtext1 <- renderText("Range:90%-100%")
  

  
  # Read data from uploaded file
  data <- reactive({
    req(input$file)
    read.table(input$file$datapath, sep = "\t", header = TRUE)
  })
  
  processedData <- eventReactive(input$processData, {
    req(data())
    
    ignore_headers <- c("Study", "Name", "Image", "LayerData", "Object")
    channels_to_ignore <- sapply(ignore_headers, paste, collapse = "|")
    
    # Select relevant means 
    cell_means <- data() %>%
      select(matches("entire.cell", ignore.case = TRUE)) %>%
      select(-matches(channels_to_ignore))
    
    colnames(cell_means) <- 
      sub("Biomarker.intensity..Entire.cell....", "", colnames(cell_means))
    
    # Extract list of unique column names (biomarkers)
    unique_headers <- colnames(cell_means)
    biomarkers <- unique_headers
    
    list(data = data(), 
         biomarkers = biomarkers, 
         cell_means = cell_means)
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

    
    # Calculate quantile based on the selected value
    quantile_value <- quantile(cell_means[[input$selectBiomarker]], probs = input$pickQuantile)
    
    # Make quantile bindex
    quantile_bindex <- cell_means[[input$selectBiomarker]] < quantile_value
    
    # Filter by quantile bindex
    filtered_means <- subset(cell_means, quantile_bindex)
    
    # Make threshold bindex
    threshold_bindex <- cell_means[[input$selectBiomarker]] > input$pickThreshold
    
    # Cutoff data by threshold
    thresholded_means <- subset(cell_means, threshold_bindex)
    
    # Calculate number of cells above threshold
    num_cells_above_threshold <- nrow(thresholded_means)
    
    # Calculate total number of cells  
    total_cells <- nrow(cell_means)
    
    # Calculate percent of cells above threshold
    percent_cells_above_threshold <- round(num_cells_above_threshold / total_cells * 100, 2)
    
    # Calculate average of data above threshold
    mean_positive_intensity <- round(mean(thresholded_means[[input$selectBiomarker]]), 2)
    
    # Calculate median of data above threshold
    median_positive_intensity <- round(median(thresholded_means[[input$selectBiomarker]]), 2)
    
    
    
    list(
      filtered_means = filtered_means,
      thresholded_means = thresholded_means,
      mean_positive_intensity = mean_positive_intensity,
      median_positive_intensity = median_positive_intensity,
      percent_cells_above_threshold = percent_cells_above_threshold
    )
  })
  
  # Generate Pos Intensity plot
  output$positive_intensity_plot <- renderPlot({
    req(input$updateResults)
    
    plot <- ggplot(selectedBiomarkerData()$filtered_means, aes(x = .data[[input$selectBiomarker]])) +
      geom_density(color="darkblue", fill="lightblue") +
      geom_vline(aes(xintercept = input$pickThreshold), color = "black") +
      geom_vline(aes(xintercept = selectedBiomarkerData()$mean_positive_intensity), color = "Red") +
      geom_vline(aes(xintercept = selectedBiomarkerData()$median_positive_intensity), color = "Blue") +
      labs(x = "Mean Intensity per cell", title = input$selectBiomarker)
    
    plot
  })
  
  
  # Text outputs of vertical lines
  output$thresholdText <- renderText({
    paste("Threshold: ", input$pickThreshold)
  })
  
  output$medianText <- renderText({
    paste("Median above threshold: ", selectedBiomarkerData()$median_positive_intensity)
  })
  
  output$meanText <- renderText({
    paste("Mean above threshold: ", selectedBiomarkerData()$mean_positive_intensity)
  })
  
  # Display the group info table
  output$stats_table <- renderTable({
    selectedBiomarkerData()$group_stats
  }
  )
  }


# Run the app
shinyApp(ui = ui, server = server)
