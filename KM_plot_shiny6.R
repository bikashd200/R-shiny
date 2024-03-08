require(FirebrowseR)
library(shiny)
library(survminer)
library(survival)
library(tidyverse)
library(tidyr)

# Define Shiny app options
options(shiny.maxRequestSize = 1000*1024^2)  # Set maximum request size to 30 MB

#https://github.com/bikashd200/R-shiny/tree/main/Data

##penguins_csv <- "https://raw.githubusercontent.com/jcheng5/simplepenguins.R/main/penguins.csv"
# Define UI
ui <- fluidPage(
  titlePanel("Kaplan-Meier Plot"),
  
  sidebarLayout(
    sidebarPanel(
      # Gene list input
      textInput("gene", "Gene Name:", value = ""),
      # Select from dropdown list
      selectInput(inputId = "cohort",
                  label = "Disease Cohort:",
                  choices = c("ACC","BLCA","BRCA","CESC","CHOL","COAD","COADREAD","DLBC","ESCA","FPPP","GBM","GBMLGG","HNSC","KICH","KIPAN","KIRC","KIRP","LAML","LGG","LIHC","LUAD","LUSC","MESO","OV","PAAD","PCPG","PRAD","READ","SARC","SKCM","STAD","STES","TGCT","THCA","THYM","UCEC","UCS","UVM"), # Provide your codes here
                  selected = "ACC"), # Default selection
      
      
      # Dropdown menu to choose between median_value and quantile_value75
      selectInput("threshold_type", "Choose Threshold Type:",
                  choices = c("Median Value" = "median", "Quantile Value 75" = "quantile"),
                  selected = "median")
    ),
    
    mainPanel(
      # KM plot output
      plotOutput("km_plot", width = "80%", height = "600px")
    )
  )
)

# Define server logic
server <- function(input, output) {
  
  # Load gene expression data
  gene_exp <- reactive({
    req(input$cohort, input$gene)
    all.Received = F
    page.Counter = 1
    page.size = 150
    gene_exp = list()
    #cancer.Type
    while(all.Received == F){
      gene_exp[[page.Counter]] = Samples.mRNASeq(format = "csv",
                                                 cohort = input$cohort,
                                                 protocol = "RSEM",
                                                 gene = input$gene,
                                                 sort_by = "cohort",
                                                 page_size = page.size,
                                                 page = page.Counter)
      if(page.Counter > 1)
        colnames(gene_exp[[page.Counter]]) = colnames(gene_exp[[page.Counter-1]])
      
      if(nrow(gene_exp[[page.Counter]]) < page.size){
        all.Received = T
      } else{
        page.Counter = page.Counter + 1
      }
    }
    
    gene_exp = do.call(rbind, gene_exp)
    dim(gene_exp)
    
    gene_exp = gene_exp[, c("tcga_participant_barcode", "expression_log2", "z.score", "cohort", "sample_type")]
    
    return(gene_exp)
  })
  
  # Load clinical data
  clinical_data <- reactive({
    req(input$cohort)
    
    all.Received = F
    page.Counter = 1
    page.size = 150
    clinical_data = list()
    
    while(all.Received == F){
      clinical_data[[page.Counter]] = Samples.Clinical(format = "csv",
                                                       cohort = input$cohort,
                                                       page_size = page.size,
                                                       page = page.Counter)
      if(page.Counter > 1)
        colnames(clinical_data[[page.Counter]]) = colnames(clinical_data[[page.Counter-1]])
      
      if(nrow(clinical_data[[page.Counter]]) < page.size){
        all.Received = T
      } else{
        page.Counter = page.Counter + 1
      }
    }
    
    clinical_data = do.call(rbind, clinical_data)
    dim(clinical_data)
    #clinical_data <- Samples.Clinical(format = "csv", cohort = cohort)
    
    # change certain values the way they are encoded
    clinical_data$deceased <- ifelse(clinical_data$vital_status == "dead", TRUE, FALSE)
    
    clinical_data$overall_survival <- ifelse(is.na(clinical_data$days_to_death) & is.na(clinical_data$days_to_last_followup), NA,
                                             ifelse(clinical_data$vital_status == 'alive', clinical_data$days_to_last_followup, clinical_data$days_to_death))
    
    return(clinical_data)
  })
  
  # KM plot
  output$km_plot <- renderPlot({
      gene_exp <- gene_exp()
      clinical_data <- clinical_data()
      datagene1 <- merge(gene_exp, clinical_data, by = 'tcga_participant_barcode')
      
      median_value <- median(datagene1$expression_log2)
      quantile_value <- quantile(datagene1$expression_log2, .75)
      quantile_value75 <- quantile_value[[1]]
      
      # Choose threshold value based on the selected type
      if (input$threshold_type == "median") {
        threshold_value <- median_value
      } else if (input$threshold_type == "quantile") {
        threshold_value <- quantile_value75
      }
      # Denote which cases have higher or lower expression than median count
      datagene1$strata <- ifelse(datagene1$expression_log2 >= threshold_value, "HIGH", "LOW")
      table(datagene1$strata)
      # Add clinical information
      #datagene$RowNames <- gsub('\\.01.*', '', datagene$RowNames)
      #datagene <- merge(datagene, clinical, by.x = 'RowNames', by.y = 'RowNames')
      
      # Convert values
      datagene1$overall_survival <- as.numeric(datagene1$overall_survival)
      
      # Create title
      title <- paste(input$gene, "_", input$cohort, "_survival", sep = "")
      
      # Fitting survival curve
      fit <- survfit(Surv(overall_survival, deceased) ~ strata, data = datagene1)
      
      # Create KM plot
      surv_plot <- ggsurvplot(fit,
                              data = datagene1,
                              pval = TRUE,
                              risk.table = TRUE,
                              title = title,
                              ggtheme = theme_classic2(base_size=16, base_family = "Arial"),
                              font.family = "Arial")
      
      
      # Output KM plot
      print(surv_plot)
      
  
  })
}

# Run the application
shinyApp(ui = ui, server = server)

