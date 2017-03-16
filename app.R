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
#source("global.R")
header <- dashboardHeader(
            title = "Fry Pathway Analysis Tool",
            titleWidth = 350
          )

sidebar <- dashboardSidebar(
            sidebarMenu(id = "tabs",
              menuItem("Read Instructions", icon = icon("book"), tabName = "instructions"),
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
            width = 350
          )

body <- dashboardBody(
          tags$head(
            tags$link(rel = "stylesheet", type = "text/css", href = "custom.css")
          ),
          tabItems(
            tabItem(tabName = "instructions",
              h2("Instructions are coming soon")
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
              tabBox(
                tabPanel("Tables",
                  fluidRow(
                    box(
                      title = "Title 1", width = 12, solidHeader = TRUE, status = "primary",
                      "Box content"
                    )
                  ),
                  fluidRow(
                    box(
                      title = "Save the table", width = 6, solidHeader = TRUE, status = "primary",
                      conditionalPanel(condition = "input.run",
                        fluidRow(
                          column(6, wellPanel(
                            radioButtons('saving_type', h5("Download"), 
                              choices = c("All rows", "Selected rows", "Filtered rows"))
                            
                            # p(downloadButton('pval_dl', 'Download'))
                          )
                          ),
                          column(6, wellPanel(
                            radioButtons("filetype", h5("File type"),
                              choices = c(".csv", ".txt", ".xlsx"))
                                    )
                          )
                        ),
                        fluidRow(
                          column(12, wellPanel(
                            textInput('filename', label = h5("Name the file")),
                            downloadButton('downloadData', 'Download table')
                          )
                          )
                        )
                        
                      )
                    )
                  )
                ),
                tabPanel("Network Maps")
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

ui <- dashboardPage(header, sidebar, body, skin = "blue")
server <- function(input, output) { }

shinyApp(ui, server)
