# app.R #
#---------------------------------------------------------------------------
#  This application is governed by the GNU General Public License version 3.
#  You can  use, modify and/ or redistribute this code under the terms
#  of the GNU General Public License: http://www.gnu.org/licenses/
#
#  Goknur Giner, The Walter and Eliza Hall Institute of Medical Research
#  March 16th, 2017
# Local path: Google Drive\Projects\Github\fry-shiny-dashboard
#---------------------------------------------------------------------------

library(shiny)
library(shinydashboard)
source("global.R")
options(shiny.maxRequestSize = 100*1024^2) # increase the size of the file can be uploaded

# define header and its appearance
header <- dashboardHeader(
            title = "FRY PATHWAY ANALYSIS TOOL",
            titleWidth = 300
          )

# define sidebar and its appearance
sidebar <- dashboardSidebar(
            sidebarMenu(id = "tabs",
              menuItem("HOME", icon = icon("home"), tabName = "home"),
              menuItem("WATCH TUTORIAL ", icon = icon("tv"), tabName = "video"),
              hr(),
              menuItem("UPLOAD DATA", icon = icon("upload"), tabName = "upload"),
              menuItem("CHOOSE PATHWAYS", icon = icon("database"), tabName = "pathway"),
              menuItem("RESULTS", icon = icon("table"), tabName = "results"),
              menuItem("CHARTS", icon = icon("bar-chart"), tabName = "charts"),
              menuItem("GENE EXPLORER", icon = icon("list"), tabName = "genes"),
              hr(),
              menuItem("Twitter", icon = icon("twitter"), tabName = "twitter"),
              menuItem("Github", icon = icon("github"), tabName = "github"),
              menuItem("Contact", icon = icon("envelope"), tabName = "email")
            ),
            width = 300
          )

# define body and its appearance
body <- dashboardBody(
          tags$head(
            # link the tags to design.css file and make the changes in appearances in design.css
            tags$link(rel = "stylesheet", type = "text/css", href = "design.css") 
          ),
          tabItems(
            tabItem(tabName = "home",
              h1("Home")
            ),
            tabItem(tabName = "video",
              h2("Youtube link is coming soon")
            ),
            tabItem(
              tabName = "upload",
              box(
                title = "UPLOAD DATA", width = 12, solidHeader = TRUE, status = "primary",
                fileInput("counts", label = h5("EXPRESSION MATRIX")),
                fileInput("design", label = h5("DESIGN MATRIX")),
                conditionalPanel("output.fileUploaded", uiOutput("contrastcomp")),
                checkboxInput("example", label = "USE EXAMPLE DATASET AND CONTINUE")
              )
              
            ),
            tabItem(
              tabName = "pathway",
              box(
                title = "SELECT DATABASE AND PATHWAY", width = 12, solidHeader = TRUE, status = "primary",
              selectInput("database", label = h5("SELECT DATABASE"),
                choices = list("GO", "KEGG", "MSigDB_HALLMARK", "REACTOME"), selected = "KEGG"),
              p(h5("SELECT PATHWAYS")),
              conditionalPanel(condition = "input.database == 'GO'",
                selectizeInput("goSelected", label = NULL,
                  choices = NULL, multiple = TRUE)),
              
              conditionalPanel(condition = "input.database == 'KEGG'",
                selectizeInput("keggSelected", label = h5("Type to search for KEGG pathways"),
                  choices = keggAll, multiple = TRUE)),
              
              conditionalPanel(condition = "input.database == 'MSigDB_HALLMARK'",
                selectizeInput("msigdbSelected", label = h5("Type to search for MSigDB_HALLMARK pathways"),
                  choices = msigdbAll, multiple = TRUE)),
              
              conditionalPanel(condition = "input.database == 'REACTOME'",
                selectizeInput("reactomeSelected", label = h5("Type to search for REACTOME pathways"),
                  choices = reactomeAll, multiple = TRUE)),
              
              checkboxInput("allGeneSets", label = "Select all gene sets"),
              
              # p("or"),
              
              # fileInput("geneList", label = h5("Upload the gene list of interest")),
              
              actionButton("run", label = h4("APPLY GENE SET TEST")),
              h3("Continue next tab to explore the results")
              )
            ),
            tabItem(
              tabName = "results",
              conditionalPanel("input.run",
                tabBox(width = "500px", height = "5000px",
                  tabPanel(h4("RESULTS"),
                    fluidRow(
                      box(
                        title = textOutput("contrastselected"),
                          #, textOutput("contrastselected")), 
                        width = 12, solidHeader = TRUE, status = "primary",
                        dataTableOutput('fryTable')
                      )
                    ),
                    fluidRow(
                      box(
                        title = "SAVE TABLE", width = 2, solidHeader = TRUE, status = "primary",
                          fluidRow(
                            column(12, wellPanel(
                              radioButtons('saving_type', h5("Download"), 
                                choices = c("All", "Selected", "Filtered")),
                              radioButtons("filetype", h5("File type"),
                                choices = c(".csv", ".txt", ".xlsx")),
                              downloadButton('downloadData', 'Save')
                            )
                              # p(downloadButton('pval_dl', 'Download'))
                            )
                            )
                      )
                    )
                  ),
                  tabPanel(h4("PATHWAY NETWORK"))
                )
              )
            ),
            tabItem(
              tabName = "charts"
            ),
            tabItem(
              tabName = "genes"
            )
          )
        )

server <- function(input, output, session) {
  updateSelectizeInput(session, "goSelected", choices = goAll, server = TRUE)
  
  output$fileUploaded <- reactive({
    if (is.null(input$design))
      0
    else
      1
  })
  
  outputOptions(output, 'fileUploaded', suspendWhenHidden = FALSE)
  
  output$contrastcomp <- renderUI({
    data <- read.table(input$design$datapath)
    selectInput("contrast", 
      "Choose the contrast of interest", choices = as.list(colnames(data)[-1]))
  })
  
  output$contrastselected <- renderText({
    paste("PATHWAY SIGNIFICANCE TABLE FOR CONTRAST ", input$contrast)
  })
  # 
  # outputOptions(output, 'conselected', suspendWhenHidden = FALSE)
  
  applyGS <- eventReactive(input$run, {
  if (input$database == 'GO') {
    if (!input$allGeneSets & !is.null(input$goSelected)) input$goSelected
    else if (input$allGeneSets) goAll
  }
  
  else if (input$database == 'MSigDB_HALLMARK') {
    if (!input$allGeneSets & !is.null(input$msigdbSelected)) input$msigdbSelected
    else if (input$allGeneSets) msigdbAll
  }
  
  else if (input$database == 'REACTOME') {
    if (!input$allGeneSets & !is.null(input$reactomeSelected)) input$reactomeSelected
    else if (input$allGeneSets) reactomeAll
  }
  
  else {
    if (!input$allGeneSets & !is.null(input$keggSelected)) input$keggSelected
    else if (input$allGeneSets) keggAll
  }
})
  
  databaseInput <- reactive({
    
    # Take a dependency on input$run so the reactive object fryTable now
    # depends on action button "Apply gene set test"
    input$run
    
    # We use isolate to avoid the dependency on database object. Because we
    # do not want gene set test to be applied prior to the selection of 
    # new gene sets
    
    database <- isolate(input$database)
    
    if (database == 'GO') gene.sets <- as.list(GO)[applyGS()]
    
    else if (database == 'KEGG') gene.sets <- KEGG[applyGS()]
    
    else if (database == 'REACTOME') gene.sets <- REACTOME[applyGS()]
    
    else gene.sets <- Hs.H[applyGS()]
    
    if (input$example) {
      cnt <- read.table("../../Google Drive/Files/datasets/dgelist.txt")
      des <- read.table("../../Google Drive/Files/datasets/design.txt")
    }
    else{
      cnt <- read.table(input$counts$datapath)
      des <- read.table(input$design$datapath) 
    }
    
    # is the following step necessary or null sets have already been discarded
    # PathwayName <- names(goAll[goAll %in% names(gene.sets)])
    # gene.sets <- gene.sets[!sapply(gene.sets, is.null)]
    
    idx <- ids2indices(gene.sets, rownames(cnt))
    idx <- idx[!sapply(idx, is.null)]
    
    if (input$example) fry <- fry(cnt, design = des, index = idx, sort = "directional")
    else fry <- fry(cnt, design = des, index = idx, sort = "directional", contrast = input$contrast)
    PathwayID <- rownames(fry)
    rownames(fry) <- NULL
    fry <- format(fry, digits = 2)
    
    if (database == 'GO') {
      m <- match(PathwayID, goAll)
      PathwayName <- names(goAll[m])
      fry.table <- data.frame(PathwayID = PathwayID,
        PathwayName = PathwayName, fry)
    }
    
    else if (database == 'KEGG') { m <- match(PathwayID, keggAll)
    PathwayName <- names(keggAll[m])
    fry.table <- data.frame(PathwayID = PathwayID,
      PathwayName = PathwayName, fry)
    }
    
    else if (database == 'REACTOME') 
      fry.table <- data.frame(PathwayName = PathwayID, fry)
    
    else fry.table <- data.frame(PathwayName = PathwayID, fry)
  })
  
  output$fryTable <- DT::renderDataTable({
    databaseInput()
  }, options = list(orderClasses = TRUE, pageLength = 5, autoWidth = TRUE),
    filter = 'top'
  )
  
  output$downloadData <- downloadHandler(
    
    # This function returns a string which tells the client
    # browser what name to use when saving the file.
      filename = function() {
      paste0(Sys.Date(), "-", input$database, "-", input$saving_type, input$filetype)
    },
    
    # This function should write data to a file given to it by
    # the argument 'file'.
    content = function(file) {
      if (input$filetype %in% c(".csv", ".txt")) {
        sep <- switch(input$filetype, ".csv" = ",", ".txt" = "\t")
        if (input$saving_type == "Filtered")
          s = input$fryTable_rows_all
        else if (input$saving_type == "Selected") 
          s = input$fryTable_rows_selected
        # Write to a file specified by the 'file' argument
        if (input$saving_type == "All") 
          write.table(format(databaseInput(), scientific = TRUE, digits = 3),
            file, sep = sep, row.names = FALSE)
        else  write.table(format(databaseInput()[s, , drop = FALSE], scientific = TRUE, digits = 3),
          file, sep = sep, row.names = FALSE)
      }
      else {
        if (input$saving_type == "Filtered")
          s = input$fryTable_rows_all
        else if (input$saving_type == "Selected") 
          s = input$fryTable_rows_selected
        if (input$saving_type == "All") 
          write.xlsx(format(databaseInput(), scientific = TRUE, digits = 3),
            file, row.names = FALSE)
        else  write.xlsx(format(databaseInput()[s, , drop = FALSE], scientific = TRUE, digits = 3),
          file, row.names = FALSE)
      }
    }
  )
}

ui <- dashboardPage(header, sidebar, body, skin = "blue")
shinyApp(ui, server)
