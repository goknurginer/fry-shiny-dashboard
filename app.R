# app.R #
#---------------------------------------------------------------------------
#  This application is governed by the GNU General Public License version 3.
#  You can  use, modify and/ or redistribute this code under the terms
#  of the GNU General Public License: http://www.gnu.org/licenses/
#
#  Goknur Giner, The Walter and Eliza Hall Institute of Medical Research
#  March 16th, 2017
#---------------------------------------------------------------------------

library(shiny)
library(shinydashboard)
source("global.R")
header <- dashboardHeader(
            title = "Fry Pathway Analysis Tool",
            titleWidth = 270
          )

sidebar <- dashboardSidebar(
            sidebarMenu(id = "tabs",
              menuItem("Home", icon = icon("home"), tabName = "home"),
              menuItem("Watch Introduction Video", icon = icon("tv"), tabName = "video"),
              hr(),
              menuItem("Upload Data", icon = icon("upload"), tabName = "upload"),
              menuItem("Choose Pathways", icon = icon("database"), tabName = "pathway"),
              menuItem("Results", icon = icon("table"), tabName = "results"),
              menuItem("Charts", icon = icon("bar-chart"), tabName = "charts"),
              menuItem("Gene Explorer", icon = icon("list"), tabName = "genes"),
              hr(),
              menuItem("Twitter", icon = icon("twitter"), tabName = "twitter"),
              menuItem("Github", icon = icon("github"), tabName = "github"),
              menuItem("Email", icon = icon("envelope"), tabName = "email")
            ),
            width = 270
          )
body <- dashboardBody(
          tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "bootstrap.css"),
            tags$style(HTML("
              .box.box-solid.box-primary>.box-header {
              color:#fff;
              background:#666666
              }
              
              .box.box-solid.box-warning>.box-header {
              color:#fff;
              background:#00bdc0
              }
              
              .box.box-solid.box-primary{
              border-bottom-color:#00bdc0;
              border-left-color:#00bdc0;
              border-right-color:#00bdc0;
              border-top-color:#00bdc0;
              }

              .skin-blue .main-header .logo {
               background-color: #00bbbf;
              }
              .skin-blue .main-header .logo:hover {
              background-color: #00bbbf;
              }
              
              .skin-blue .main-header .navbar {
              background-color: #00bbbf;
              } 

              .skin-blue .main-sidebar {
              background-color: #a1e8ed;
              }

              .skin-blue .main-sidebar .sidebar .sidebar-menu .active a{
                              background-color: #ffa477;
              }

              .skin-blue .main-sidebar .sidebar .sidebar-menu a{
                              background-color: #a1e8ed;
                              color: #ffffff;
              }

              .skin-blue .main-sidebar .sidebar .sidebar-menu a:hover{
                              background-color: #ff9503;
              }

              .skin-blue .main-header .navbar .sidebar-toggle:hover{
                              background-color: #ff9503;
              }

              "))
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
              fileInput("counts", label = "Upload expression values"),
              fileInput("design", label = "Upload design of the experiment"),
              checkboxInput("example", label = "Use example dataset and continue")
            ),
            tabItem(
              tabName = "pathway",
              selectInput("database", label = h3("Choose a database"),
                choices = list("GO", "KEGG", "MSig_HALLMARK", "REACTOME")),
              br(),
              p(h3("Choose pathways")),
              conditionalPanel(condition = "input.database == 'GO'",
                selectizeInput("goSelected", label = h5("Type in the box to search for GO pathways"),
                  choices = NULL, multiple = TRUE)),
              
              conditionalPanel(condition = "input.database == 'KEGG'",
                selectizeInput("keggSelected", label = h5("Type to search for KEGG pathways"),
                  choices = keggAll, multiple = TRUE)),
              
              conditionalPanel(condition = "input.database == 'MSig_HALLMARK'",
                selectizeInput("MsigSets", label = h5("Type to search for MSig_HALLMARK pathways"),
                  choices = msigAll_hall, multiple = TRUE)),
              
              conditionalPanel(condition = "input.database == 'REACTOME'",
                selectizeInput("reactomeSelected", label = h5("Type to search for REACTOME pathways"),
                  choices = reactomeAll, multiple = TRUE)),
              
              checkboxInput("allGeneSets", label = "Select all gene sets"),
              
              # p("or"),
              
              # fileInput("geneList", label = h5("Upload the gene list of interest")),
              
              actionButton("run", label = h4(strong("Apply gene set test"))),
              h3("Continue next tab to explore the results")
            ),
            tabItem(
              tabName = "results",
              tabBox(width = "500px", height = "5000px",
                tabPanel("Results",
                  fluidRow(
                    box(
                      title = "Pathway Significance Table", width = 12, solidHeader = TRUE, status = "primary",
                      dataTableOutput('fryTable')
                    )
                  ),
                  fluidRow(
                    box(
                      title = "Save table", width = 2, solidHeader = TRUE, status = "warning", 
                      conditionalPanel(condition = "input.run",
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
                  )
                ),
                tabPanel("Pathway Network")
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



options(shiny.maxRequestSize=100*1024^2)
msigAll_hall <- names(Hs.H)
server <- function(input, output, session) {
  updateSelectizeInput(session, "goSelected", choices = goAll, server = TRUE)
  
  applyGS <- eventReactive(input$run, {
    
    if(input$database == 'GO') {
      if (!input$allGeneSets & !is.null(input$goSelected)) input$goSelected
      else if (input$allGeneSets) goAll
    }
    
    else if (input$database == 'MSig_HALLMARK') {
      if (!input$allGeneSets & !is.null(input$MsigSets)) input$MsigSets
      else if (input$allGeneSets) msigAll_hall
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
    # do not want gene set test to be applied when we pick a different database
    # before we select the gene sets
    database <- isolate(input$database)
    
    if(database == 'GO') gene.sets <- as.list(GO)[applyGS()]
    
    else if(database == 'KEGG') gene.sets <- KEGG[applyGS()]
    
    else if(database == 'REACTOME') gene.sets <- REACTOME[applyGS()]
    
    else gene.sets <- Hs.H[applyGS()]
    
    if (input$example) {
      cnt <- read.table("dgelist.txt")
      des <- read.table("design.txt")
    }
    else{
      cnt <- read.table(input$counts$name)
      des <- read.table(input$design$name) 
    }
    
    # is the following step necessary or null sets have already been discarded
    #PathwayName <- names(goAll[goAll %in% names(gene.sets)])
    #gene.sets <- gene.sets[!sapply(gene.sets, is.null)]
    
    idx <- ids2indices(gene.sets, rownames(cnt))
    idx <- idx[!sapply(idx, is.null)]
    
    fry <- fry(cnt, design =des, index = idx, sort="directional")
    PathwayID <- rownames(fry)
    
    if(database == 'GO') {
      m <- match(PathwayID, goAll)
      PathwayName <- names(goAll[m])
      fry.table <- data.frame(PathwayID = PathwayID,
        PathwayName = PathwayName, fry)
    }
    
    else if(database == 'KEGG') { m <- match(PathwayID, keggAll)
    PathwayName <- names(keggAll[m])
    fry.table <- data.frame(PathwayID = PathwayID,
      PathwayName = PathwayName, fry)
    }
    
    else if(database == 'REACTOME') 
      fry.table <- data.frame(PathwayID = PathwayID, fry)
    
    else fry.table <- data.frame(PathwayID = PathwayID, fry)
  })
  
  output$fryTable <- DT::renderDataTable({
    format(databaseInput(), scientific = TRUE, digits = 3)
  }, options = list(orderClasses = TRUE), filter = 'top'
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
      if(input$filetype %in% c(".csv", ".txt")) {
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
