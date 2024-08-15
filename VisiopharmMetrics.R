library(shiny)
library(dplyr)

# Define UI
ui <- fluidPage(
  titlePanel("Visiopharm Validation Calculator V3"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload .tsv File"),
      actionButton("processData", "Process Validation Data")
    ),
    mainPanel(
      "Aggregate Metrics", tableOutput("Metrics")),
  ),
  p("Sensitivity(Recall,true positive rate) is the probability of a positive test result actually being positive." ),
  p( "Specificity (true negative rate) is the probability of a negative test result actually being negative."),
  p("Accuracy (experimental bias) is how close a given set of measurements are to their true value."),
  p("Precision (positive predictive value, spread of data) is how close the measurements are to each other."),
  p("F1 Score is the harmonic mean of Sensitivity and Precision, represents both biases in one metric."),
  p(" Values range from 0-1, >.70 is good.")
)



# Define server logic
server <- function(input, output) {
  
  # Read data from uploaded file
  loaded_data <- reactive({
    req(input$file)
    read.delim(input$file$datapath, sep = "\t", header = TRUE)
  })
  
  # Performance metric Calculation
  results <- eventReactive(input$processData, {
    
    data <- loaded_data()
    
    #select relevant columns
    
    orginal_headers <- names(data)
    
    ignore_headers <-c("Study.level.1|Study.level.2|Image")
    
    channels_to_ignore <- sapply(list(ignore_headers), paste, collapse="|")
    
    #select relevant columns
    selected_columns <-
      data %>% select(matches("positives|negatives|AreaAll", ignore.case = TRUE))
    selected_columns <-
      selected_columns %>% select(-matches(channels_to_ignore))
    
    #remove zero columns
    selected_columns <-
     selected_columns[, colSums(selected_columns != 0) > 0]
    
    #extract relevant labels
    relevent_columns  <-
      selected_columns %>% select(matches("positives|negatives", ignore.case = TRUE))
    
    relevent_headers <- names(relevent_columns)
    
    #extract first word in headers
    first_words <- sapply(strsplit(relevent_headers, "\\."), function(x) x[1])
    
    #unique first words
    unique_headers <-unique(first_words)
    
    #Aggregate calculation
    raw_headers <- colnames(selected_columns)
    
    # Calculate column sums
    sums <- colSums(selected_columns)
    
    # Create new dataframe with headers and sums
    aggregate_data <- data.frame(sums)
    
    #Extract Total area
    TotalArea <- aggregate_data["AreaAll", 1]
    
    #initialize results
    results <- data.frame(class = unique_headers,
                          accuracy = NA,
                          precision = NA,
                          sensitivity = NA, 
                          specificity = NA,
                          f1 = NA)
    
    # Loop through each class 
    for(header in unique_headers) {
      
  
      # Extract true positives, false positives, false negatives
      tp <- aggregate_data[paste0(header, ".True.Positives"),]
      if(is.null(tp)) {
        tp <- 0
      }
      
      fp <- aggregate_data[paste0(header, ".False.Positives"),] 
      if(is.null(fp)) {
        fp <- 0
      }
      
      fn <- aggregate_data[paste0(header, ".False.Negatives"),]
      if(is.null(fn)) {
        fn <- 0
      }
      
      tn <- TotalArea - tp - fp - fn
      
      
      
      # Calculate metrics
      accuracy <- (tp + tn)/(tp + fp + fn + tn)
      precision <- tp / (tp + fp)
      sensitivity <- tp / (tp + fn)
      specificity <- tn / (tn + fp)
      f1 <- (2 * (precision * sensitivity) / (precision + sensitivity))
      
      # Store in results
      results$accuracy[results$class == header] <-accuracy
      results$precision[results$class == header] <- precision
      results$sensitivity[results$class == header] <- sensitivity
      results$specificity[results$class == header] <- specificity
      results$f1[results$class == header] <- f1
      results[is.na(results)] <- 0.0
    
    }
    
    return(results)
    
  })
  
  # Display Results
  output$Metrics <- renderTable({
    results()
  })
  
}

# Run the app
shinyApp(ui = ui, server = server)

