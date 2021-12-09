library(shiny)
library(shinythemes)
library(pryr)

plnUI <- shinyUI({
  navbarPage("PLaNet",
             theme = shinytheme("cosmo"),
             tabPanel("Using PLaNet",
                      includeHTML("introduction.html")
             ),
             tabPanel("Example Dataset",
             # Row 1: Title of the application (example tab) --------------------------------
             fluidRow(
               # Application title
               titlePanel("PLaNet - Phenotypic Linkage Network Web App - Example Data"),
               hr()
             ),
             # Row 2: Example data input, plot and statistical data -------------------------
             fluidRow(width=12,
                              # Tabset containing the Example Data and User Data panel
                              column(width=3, # Example data panel
                                     sidebarPanel(width=12,
                                                  h4("Gene Visualisation"),
                                                  radioButtons(
                                                    "geneSelect",
                                                    "Set of Genes to Visualise:",
                                                    c("Gene Search"="textbox",
                                                      "Select Gene Set (GWAS)"="dropmenu",
                                                      "Select Gene Set (HGMD)"="dropmenuHGMD",
                                                      "Select Multiple Gene Sets"="doubleGeneset",
                                                      "Upload Gene List"="upload"),
                                                    selected="dropmenu"
                                                  ),
                                                  uiOutput("genesetUI"),
                                                  verbatimTextOutput("test"),
                                                  uiOutput("errors")
                                     ),
                                     actionButton(inputId="seedShuffle",label="Shuffle Network",icon("random")),
                                     checkboxGroupInput(
                                       inputId="advancedOptions",
                                       label="Advanced Options",
                                       choices=c("Select"="select")
                                     ),
                                     uiOutput("advancedOutput"),
                                     uiOutput("exportui"),
                                     plotOutput("gnetClusterBox"),
#                                     actionButton("store_position", "Store positions!"),
                                     downloadLink('downloadNetwork', 'Download network',
                                                  icon=shiny::icon("download")),
                                     actionButton("store_position","Store Network Positions"),
                                     tableOutput("sessionInfo"),
                                     class="col-md-3"), # Define the column as bootstrap grid system width 3 (medium pixelation =>768)
                      column(width=6, # Column containing the visNetwork PLN gene network
                            mainPanel(
                                       h4("Phenotypic Linkage Network"),
                                       textOutput("userGenes"),
                                       visNetworkOutput("visNet",
                                                        height="700px"), # Output the visNetwork PLN
                                       uiOutput("nodeSliderOut"),
                                       width="100%"),
                                     class="col-md-6"),
                      column(width=3, # Column containing the data checkbox
                             selectInput("plnSelect",
                                         h3("Select Phenotypic Linkage Network:"),
                                         selected=NULL,
                                         choices=c("General"="General",
                                                   "Metabolic"="T2D",
                                                   "Nervous"="Nervous"),
                                         multiple=FALSE,
                                         width="auto"
                                         ),
#                             actionButton("plotButton","Plot Graph"),
                             uiOutput("dataCheckbox"),
                             verbatimTextOutput("seedNumber"),
                             tableOutput("dataInfo"),
                             plotOutput("clusterBox"),
                             h3("iGraph Gene Distances"),
                             h5("Gene 1"),
                              textInput(
                                "GeneDist1", 
                                label = NULL, 
                                placeholder = "e.g. FANCL, BRD7, GABRA2"),
                             h5("Gene 2"),
                              textInput(
                                "GeneDist2", 
                                label = NULL, 
                                placeholder = "e.g. FANCL, BRD7, GABRA2"),
#                             plotOutput("biclusteringBox"),
                             class="col-md-3")# Define the column as bootstrap grid system width 6 (medium pixelation =>768)
             )
             )
  )
})