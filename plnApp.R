library(shiny)
library(dplyr)
library(igraph)
library(network)
library(visNetwork)
library(shinythemes)
library(ggplot2)

source("PLaNetUI.R",local=TRUE)
source("PLaNetServer2.R",local=TRUE)

# Run the application 
shinyApp(
  ui = plnUI, 
  server = plnServer
)