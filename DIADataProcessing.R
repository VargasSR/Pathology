library(shiny)
library(ggplot2)
library(nortest)
library(dplyr)

# Define UI
ui <- fluidPage(
  titlePanel("Image Analysis Data Processing and Statistics"),
  sidebarLayout(
    sidebarPanel(
      fileInput("file", "Upload .tsv File"),
      actionButton("processData", "Read File"),
      selectInput("group", "Select Group:", choices = NULL),
      textInput("groupName", "Edit Group Name:"),
      actionButton("saveInfo", "Update Group Name"),
      selectInput("targetVariable", "Select Target Variable:", choices = NULL),
      actionButton("analyzeData", "Analyze Data"),
      downloadButton("exportData", "Export Dataset")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Group Information", tableOutput("groupTable")),
        tabPanel("Uploaded Data", tableOutput("resultsTable")),
        tabPanel("Scatterplot", plotOutput("scatterplot")),
        tabPanel("Normality Test", tableOutput("normalityTable")),
        tabPanel("Group Statistics", tableOutput("groupStatsTable")),
        tabPanel("Columnized Data for Prism", tableOutput("columnsTable"))
      )
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Read data from uploaded file
  data <- reactive({
    req(input$file)
    read.delim(input$file$datapath, sep = "\t")
  })
  
  # Process data and split column "Name"
  processedData <- eventReactive(input$processData, {
    data <- data()
    
    splitNames <- strsplit(data$Name, " ")
    
    # Extract the split values into separate columns
    data$StudyID <- sapply(splitNames, function(x) x[1])
    data$AnimalID <- sapply(splitNames, function(x) x[2])
    data$Tissue <- sapply(splitNames, function(x) x[3])
    data$Group <- sapply(splitNames, function(x) x[4])
    
    # Specify the columns to exclude
    exclude_columns <- c("Study.level.1", "Study.level.2", "Study.level.3", "Image", "LayerData", "Tissue", "StudyID")
    
    # Select all columns except the excluded ones
    data <- data[, !colnames(data) %in% exclude_columns]
    data
  })
  
  # Display processed data in a table
  output$resultsTable <- renderTable({
    processedData()
  })
  
  # Update group selection choices based on processed data
  observe({
    updateSelectInput(inputId = "group", choices = unique(processedData()$Group))
    saveGroupInfo() 
  })
  
  # Update target variable choices based on processed data
  observe({
    updateSelectInput(inputId = "targetVariable", choices = colnames(processedData()))
  })
  
  # Store edited group information
  groupInfo <- reactiveValues()
  
  # Save edited group information
  saveGroupInfo <- eventReactive(input$saveInfo, {
    group <- input$group
    groupName <- input$groupName
    groupInfo[[group]] <- list(groupName = groupName)
    processedData()  
  })
  
  # Create a new column "Group_Name" in the processed data table
  analyzedData <- eventReactive(input$analyzeData, {
    data <- processedData()
    
    # Get the group name for each row based on the "Group" column value
    data$Group_Name <- sapply(data$Group, function(x) {
      if (is.null(groupInfo[[x]])) {
        ""
      } else {
        groupInfo[[x]]$groupName
      }
    })
    
    data
  })
  
  # Display analyzed data in a table
  output$resultsTable <- renderTable({
    analyzedData()
  })
  
  # Display edited group information in a table
  output$groupTable <- renderTable({
    data.frame(Group = names(groupInfo), 
               GroupName = sapply(names(groupInfo), function(x) groupInfo[[x]]$groupName))
  })
  
  # Generate scatterplot
  output$scatterplot <- renderPlot({
    data <- analyzedData()
    targetVariable <- input$targetVariable
    
    if (!is.numeric(data[[targetVariable]])) {
      return("Error: Selected target variable is not numeric.")
    }
    
    ggplot(data, aes(x = Group_Name, y = !!sym(targetVariable))) +
      geom_point() +
      labs(x = "Group Name", y = targetVariable)
  })
  
  # Perform normality tests
  output$normalityTable <- renderTable({
    data <- analyzedData()
    targetVariable <- input$targetVariable
    
    if (!is.numeric(data[[targetVariable]])) {
      return("Error: Selected target variable is not numeric.")
    }
    
    # Shapiro-Wilk test
    shapiro_result <- shapiro.test(data[[targetVariable]])
    shapiro_p_value <- format(shapiro_result$p.value, digits = 5)
    shapiro_statistic <- format(shapiro_result$statistic, digits = 5)
    
    # Anderson-Darling test
    ad_result <- ad.test(data[[targetVariable]])
    ad_p_value <- format(ad_result$p.value, digits = 5)
    ad_statistic <- format(ad_result$statistic, digits = 5)
    
    # Create a data frame with the test results
    result_df <- data.frame(Test = c("Shapiro-Wilk", "Anderson-Darling"),
                            Variable = c(targetVariable, targetVariable),
                            p_value = c(shapiro_p_value, ad_p_value),
                            statistic = c(shapiro_statistic, ad_statistic),
                            Normal = c(ifelse(as.numeric(shapiro_p_value) > 0.05, "Yes", "No"),
                                       ifelse(as.numeric(ad_p_value) > 0.05, "Yes", "No")))
    
    result_df
  })
  
  # Calculate group statistics
  output$groupStatsTable <- renderTable({
    data <- analyzedData()
    targetVariable <- input$targetVariable
    
    if (!is.numeric(data[[targetVariable]])) {
      return("Error: Selected target variable is not numeric.")
    }
    
    group_stats <- data %>%
      group_by(Group_Name) %>%
      summarise(
        Mean = format(mean(!!sym(targetVariable)), digits = 5),
        Median = format(median(!!sym(targetVariable)), digits = 5),
        SD = format(sd(!!sym(targetVariable)), digits = 5),
        Semi_IQR = format(IQR(!!sym(targetVariable))/2, digits = 5)
      )
    
    group_stats
  })
  # Organize into columns
  output$columnsTable <- renderTable({
    
    data <- analyzedData()
    targetVariable <- input$targetVariable
    
    if(!is.numeric(data[[targetVariable]])) {
      return("Error: Selected target variable is not numeric.") 
    }
    
    group_names <- unique(data$Group_Name)
    
    columns <- sapply(group_names, function(x) {
      data[data$Group_Name == x, targetVariable]
    })
    
    colnames(columns) <- group_names
    
    as.data.frame(columns)
    
  })
  
  
  # Export dataset
  output$exportData <- downloadHandler(
    filename = function() {
      "processed_data.tsv"
    },
    content = function(file) {
      data <- analyzedData()
      targetVariable <- input$targetVariable
      
      if (!is.numeric(data[[targetVariable]])) {
        return("Error: Selected target variable is not numeric.")
      }
      
      transformedColumnName <- paste0("transformed_", targetVariable)
      
      if (transformedColumnName %in% colnames(data)) {
        data[[transformedColumnName]] <- NULL
      }
      
      data[[transformedColumnName]] <- log(data[[targetVariable]])
      
      saveData <- data
      
      write.table(saveData, file = file, sep = "\t", row.names = FALSE, quote = FALSE)
    }
  )
  
}

# Run the app
shinyApp(ui = ui, server = server)
