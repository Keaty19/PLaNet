library(shiny)
library(dplyr)
library(igraph)
library(network)
library(visNetwork)
library(shinythemes)
library(ggplot2)
library(data.table)
library(scales)
library(pryr)
library(readr)
library(tidyverse)
library(data.table)

plnServer <- function(input,output,session) {
  # Loading In Data ------------------------------------------------------------
  
#  load("sig_data.RData")
#  nodes=fread("sigNodes.csv")
#  edges=fread("sigEdges.csv")
#  load("/Users/samkeat/Desktop/Neo4j/Barebones/Full_Data2.Rdata")

  nodes=fread("Nodes.csv")
  edges=fread("Edges.csv")
  geneList=fread("geneList.csv")
  
#  edges=edges[-c(as.numeric(row.names(edges[edges$from == edges$to,]))),]
#  edges=edges[!edges$dataset=="d10",]
#  plnClusterGen=readRDS("General_GWAScat_Dataset_Clust.rds")
  plnClusterGen=fread("General_GWAScat_Dataset_Clust2.csv")
  plnClusterT2d=fread("T2d_GWAScat_Dataset_Clust2.csv")
  plnClusterNer=fread("Nervous_GWAScat_Dataset_Clust2.csv")
  
  biclustering=fread("Biclustering2.csv")
  clustering=fread("clustering_protein.csv")
  geneList=fread("geneList2_p.csv")
  hgmd=fread("HGMD_genes.csv")
  
  # Defining Variables ---------------------------------------------------------
  
  choiceCol=list(
    tags$span("BrainSpan RNAseq",style = "color: red;font-weight: bold;"),
    tags$span("Gene Coexpression",style = "color: green;font-weight: bold;"),
    tags$span("Cotranslation Data",style = "color: olive;font-weight: bold;"),
    tags$span("Gene Ontology Annotation Database",style = "color: blue;font-weight: bold;"),
    tags$span("GTEx",style = "color: orange;font-weight: bold;"),
    tags$span("InterPro Protein Families",style = "color: grey;font-weight: bold;"),
    tags$span("Kegg Reactome",style = "color: pink;font-weight: bold;"),
    tags$span("Text Literature",style = "color: purple;font-weight: bold;"),
    tags$span("CDR Protein-Protein Interactions",style = "color: brown;font-weight: bold;")
  )
    fullSel=c("1","2","3","4","5",
           "6","7","8","9")
   
  # Functions --------------------------------------------------------------------
  
  # genesetUI - Select Gene Set ------------------------------------------------  
  output$genesetUI <- renderUI(
    {
      if (is.null(input$geneSelect))
        return()
      
      switch(input$geneSelect,
             # Render TextArea -------------------------------------------------
             "textbox" = div(
               h5("Enter a List of Genes to Visualise (Separated by Commas):"),
               textAreaInput(
                 "geneList", 
                 label = NULL, 
                 placeholder = "e.g. FANCL, BRD7, GABRA2",
                 resize = "vertical"),
               # Gene Neighbour Slider
               sliderInput(
                 "nodeNeighbour",
                 h3("Gene Neighbours:"),
                 min = 0,
                 max = 5,
                 value = 1,
                 width = "auto"
               )
#               actionButton("userSearch", "Update")
             ),
             
             # Render SelectInput ----------------------------------------------
             "dropmenu" = div(
               selectInput( 
                 "geneExample",
                 label = "Select Geneset:",
                 selected = NULL,
                 choices = c("None"="None",
                             "Parkinson's Disease (PD)"="Parkinson's disease",
                             "Alzheimer's Disease (AD)"="Alzheimer's disease",
                             "Multiple sclerosis (MS)"="Multiple sclerosis",
                             "Schizophrenia"="Schizophrenia",
                             "Prostate Cancer"="Prostate cancer",
                             "Breast Cancer"="Breast cancer",
                             "Asthma"="Asthma",
                             "Coronary Artery Disease (CAD)"="Coronary artery disease",
                             "Type 2 Diabetes (T2D)"="Type 2 diabetes"
                             ),
                 multiple = FALSE,
                 width = "auto"
               ),
               renderTable(
                 data.frame(
                   "NumberOfGenes"=length(userGenes())
                   )
               ),
               # Gene Neighbour Slider
               sliderInput(
                 "nodeNeighbour",
                 h3("Gene Neighbours:"),
                 min = 0,
                 max = 5,
                 value = 1,
                 width = "auto"
               )
             ),
            "dropmenuHGMD" = div(
              selectInput( 
                "geneExampleHGMD",
                label = "Select HGMD Geneset:",
                selected = NULL,
                choices = c("None",
                            "Amyotrophic lateral sclerosis",
                            "Cognitive decline",
                            "Dementia",
                            "Frontotemporal dementia",
                            "Multiple system atrophy",
                            "Neuronal ceroid lipofuscinosis",
                            "Parkinson's disease",
                            "Spinal muscular atrophy",
                            "Spinocerebellar ataxia"
                            ),
                multiple = FALSE,
                width = "auto"
              ),
              # Gene Neighbour Slider
              sliderInput(
                "nodeNeighbour",
                h3("Gene Neighbours:"),
                min = 0,
                max = 5,
                value = 1,
                width = "auto"
              )
            ),
            "doubleGeneset"=div(
              selectInput(
                "doubleGeneset1",
                label = 
#                  "Select Geneset 1:",
                list(
                  tags$span("Select Geneset 1:",style = "color: darkgreen;font-weight: bold;")
                ),
                selected = NULL,
                choices = c("None"="None",
                            "Parkinson's Disease (PD)"="Parkinson's disease",
                            "Alzheimer's Disease (AD)"="Alzheimer's disease",
                            "Multiple sclerosis (MS)"="Multiple sclerosis",
                            "Schizophrenia"="Schizophrenia",
                            "Prostate Cancer"="Prostate cancer",
                            "Breast Cancer"="Breast cancer",
                            "Asthma"="Asthma",
                            "Coronary Artery Disease (CAD)"="Coronary artery disease",
                            "Type 2 Diabetes (T2D)"="Type 2 diabetes"
                            ),
                multiple = FALSE,
                width = "auto"
              ),
              selectInput(
                "doubleGeneset2",
                label = 
#                  "Select Geneset 2:",
                  list(
                    tags$span("Select Geneset 2:",style = "color: orangered;font-weight: bold;")
                  ),
                  selected = NULL,
                choices = c("None"="None",
                            "Parkinson's Disease (PD)"="Parkinson's disease",
                            "Alzheimer's Disease (AD)"="Alzheimer's disease",
                            "Multiple sclerosis (MS)"="Multiple sclerosis",
                            "Schizophrenia"="Schizophrenia",
                            "Prostate Cancer"="Prostate cancer",
                            "Breast Cancer"="Breast cancer",
                            "Asthma"="Asthma",
                            "Coronary Artery Disease (CAD)"="Coronary artery disease",
                            "Type 2 Diabetes (T2D)"="Type 2 diabetes"
                            ),
                multiple = FALSE,
                width = "auto"
              ),
              # Gene Neighbour Slider
              sliderInput(
                "nodeNeighbour",
                h3("Gene Neighbours:"),
                min = 0,
                max = 5,
                value = 1,
                width = "auto"
              )
            ),
            "upload" = div(
              fileInput("geneUpload", 
                        label="Upload List of Genes:",
                        multiple = FALSE, # Only accept one singular upload
                        accept = c("text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")),
              # Gene Neighbour Slider
              sliderInput(
                "nodeNeighbour",
                h3("Gene Neighbours:"),
                min = 0,
                max = 5,
                value = 1,
                width = "auto"
              )
            )
      )
    }
  )
  
  # advancedOutput: Advanced Options Output ------------------------------------
  
  output$advancedOutput <- renderUI(
    {
      if (is.null(input$advancedOptions))
        return()
      switch(input$advancedOptions,
             "select"=div(
               numericInput(
                 inputId = "seedNumInput", 
                 label = "Enter Seed Number", 
                 width = "100%",
                 value=1)
             )
      )
    })
    
  # dataCheckbox: Select Dataset Sources Checkboxes ----------------------------
  
  output$dataCheckbox <- renderUI(
    {
      if (is.null(input$geneSelect))
        return()
      
      switch(input$plnSelect,
             # Dataset Checkboxes for General PLN
             "General"=div(
               checkboxGroupInput(
                 "datasetCheckbox", 
                 label = "Select Dataset:", 
                  choiceNames = choiceCol,
                  choiceValues=fullSel
               )
             ),
             # Dataset Checkboxes for Metabolic (T2D) PLN
             "T2D"=div(
               checkboxGroupInput(
                 "datasetCheckbox", 
                 label = "Select Dataset:", 
                 choiceNames = choiceCol,
                 choiceValues=fullSel
                 )
               ),
             # Dataset Checkboxes for Nervous System PLN
             "Nervous"=div(
               checkboxGroupInput(
                 "datasetCheckbox", 
                 label = "Select Dataset:", 
                 choiceNames = choiceCol,
                 choiceValues=fullSel
                 )
               )
      )
      }
  )
  
  # selectDataset: Reactive Output of the datasetCheckbox Function --------------
  selectDataset <- reactive(
    {
      input$datasetCheckbox
    }
  )
  
  # Search String Split  ----------------------------------------------------
  gene <- reactive({
    # Select all Dataset Sources by Default
    updateCheckboxGroupInput(session,
                             "datasetCheckbox",
                             selected=fullSel
                             )
    # Stringsplit the Gene Inputs
      unlist(
        strsplit(
          as.character(input$geneList), 
          split = "[;, [:space:]]+", 
          fixed = FALSE
        )
      )
    })
  
  # nodeSlider: Gene Slider ----------------------------------------------------
  
  output$nodeSliderOut <- renderUI(
    sliderInput("nodeSlider",
              h3("Number of Genes to Display"), # Slider to select the number of nodes to display
              min=0,
#              max=500,
              max=dim(nodes[pln==as.character(input$plnSelect)])[1],
              value=50,
              width="100%")
  )
  
  nodeSlider <- reactive({ # Code for the slider that selects number of genes
    input$nodeSlider
  })
  
  # dataInfo: PLN Statistics Table ---------------------------------------------
  
  output$dataInfo <- renderTable({
    if (is.null(input$pln))
      return()
    data.frame(
      "Genes"=dim(nodes[pln==input$plnSelect])[1],
      "Relationships"=dim(edges[pln==input$plnSelect])[1],
      "GenesDisplayed"=as.character(nodeSlider())
    )
  })
  
  range01 <- function(x){(x-min(x))/(max(x)-min(x))}
  
  # clusterBox: Clustering Test PLN Comparative Boxplots -----------------------
          
  output$clusterBox <- renderPlot({
    if (is.null(input$geneExample))
      return()
    
    if (input$geneExample == "Parkinson's disease"){
      geneList2=clustering[trait=="Parkinson's disease"]
    } else if (input$geneExample == "Alzheimer's disease"){
      geneList2=clustering[trait=="Alzheimer's disease"]
    } else if (input$geneExample == "Type 2 diabetes"){
      geneList2=clustering[trait=="Type 2 diabetes"]
    } else if (input$geneExample == "Multiple sclerosis"){
      geneList2=clustering[trait=="Multiple sclerosis"]
    } else if (input$geneExample == "Schizophrenia"){
      geneList2=clustering[trait=="Schizophrenia"]
    } else if (input$geneExample == "Breast cancer"){
      geneList2=clustering[trait=="Breast cancer"]
    } else if (input$geneExample == "Prostate cancer"){
      geneList2=clustering[trait=="Prostate cancer"]
    } else if (input$geneExample == "Asthma"){
      geneList2=clustering[trait=="Asthma"]
    } else if (input$geneExample == "Coronary artery disease"){
      geneList2=clustering[trait=="Coronary artery disease"]
    }
    
#    if (input$geneExample!="None"){
#      geneList2=clustering[clustering$trait==input$geneExample,]
#    } else {
#      geneList2=NULL
#    }
    if (input$geneExample!="None"){
    if (geneList2[1,4]==0){
      titleClustGen="pval < 0.001"
    } else {
      titleClustGen=paste("pval =",round(geneList2[1,4],4),sep=" ")
    }
    if (geneList2[1,7]==0){
      titleClustT2d="pval < 0.001"
    } else {
      titleClustT2d=paste("pval =",round(geneList2[1,7],4),sep=" ")
    }
    if (geneList2[1,10]==0){
      titleClustNer="pval < 0.001"
    } else {
      titleClustNer=paste("pval =",round(geneList2[1,10],4),sep=" ")
    }
    
    list=data.frame("General"=range01(c(as.data.frame(geneList2)[1,3],
                                        as.data.frame(geneList2)[,2])),
                     "Metabolic"=range01(c(as.data.frame(geneList2)[1,6],
                                           as.data.frame(geneList2)[,5])),
                     "Nervous"=range01(c(as.data.frame(geneList2)[1,9],
                                         as.data.frame(geneList2)[,8]))
    )
    
    boxplot(list$General,
            list$Metabolic,
            list$Nervous,
            main=paste(geneList2[1,1],"Clustering",sep=" "),
            ylab="Sum of Weighted Links Within Gene Network",
            col="darkblue",
            xaxt="n")
    axis(side = 1,
         at = 1:3,
         labels=c(
           paste("General",titleClustGen,sep="\n"),
           paste("Metabolic",titleClustT2d,sep="\n"),
           paste("Nervous",titleClustNer,sep="\n")
         ),
         lwd.ticks = FALSE)
    stripchart(list[1,],
               vertical=T,
               method="overplot",
               add=T,
               col="red",
               pch=20,
               cex=3)
    }
    })

  # gnetClusterBox: Clustering Test Gene Network (PLN Data) Comparative Boxplots ------
  output$gnetClusterBox <- renderPlot({
    if (is.null(input$geneExample))
      return()
    
    if (input$plnSelect=="General"){
      plnCluster=plnClusterGen
    } else if (input$plnSelect=="T2D"){
      plnCluster=plnClusterT2d
    } else if (input$plnSelect=="Nervous"){
      plnCluster=plnClusterNer
    }
    
    if (input$geneExample == "Parkinson's disease"){
      geneList3=plnCluster[trait=="Parkinson's disease"]
    } else if (input$geneExample == "Alzheimer's disease"){
      geneList3=plnCluster[trait=="Alzheimer's disease"]
    } else if (input$geneExample == "Type 2 diabetes"){
      geneList3=plnCluster[trait=="Type 2 diabetes"]
    } else if (input$geneExample == "Multiple sclerosis"){
      geneList3=plnCluster[trait=="Multiple sclerosis"]
    } else if (input$geneExample == "Schizophrenia"){
      geneList3=plnCluster[trait=="Schizophrenia"]
    } else if (input$geneExample == "Breast cancer"){
      geneList3=plnCluster[trait=="Breast cancer"]
    } else if (input$geneExample == "Prostate cancer"){
      geneList3=plnCluster[trait=="Prostate cancer"]
    } else if (input$geneExample == "Asthma"){
      geneList3=plnCluster[trait=="Asthma"]
    } else if (input$geneExample == "Coronary artery disease"){
      geneList3=plnCluster[trait=="Coronary artery disease"]
    }
    
    if (input$geneExample!="None"){
      if (geneList3[1,4]==0){
        titleBrain="pval < 0.001"
      } else {
        titleBrain=paste("pval =",round(geneList3[1,4],4),sep=" ")
      }
      if (geneList3[1,7]==0){
        titleCoex="pval < 0.001"
      } else {
        titleCoex=paste("pval =",round(geneList3[1,7],4),sep=" ")
      }
      if (geneList3[1,10]==0){
        titleCotrans="pval < 0.001"
      } else {
        titleCotrans=paste("pval =",round(geneList3[1,10],4),sep=" ")
      }
      if (geneList3[1,13]==0){
        titleGO="pval < 0.001"
      } else {
        titleGO=paste("pval =",round(geneList3[1,13],4),sep=" ")
      }
      if (geneList3[1,16]==0){
        titleGtex="pval < 0.001"
      } else {
        titleGtex=paste("pval =",round(geneList3[1,16],4),sep=" ")
      }
      if (geneList3[1,19]==0){
        titleInter="pval < 0.001"
      } else {
        titleInter=paste("pval =",round(geneList3[1,19],4),sep=" ")
      }
      if (geneList3[1,22]==0){
        titleKegg="pval < 0.001"
      } else {
        titleKegg=paste("pval =",round(geneList3[1,22],4),sep=" ")
      }
      if (geneList3[1,25]==0){
        titleLit="pval < 0.001"
      } else {
        titleLit=paste("pval =",round(geneList3[1,25],4),sep=" ")
      }
      if (geneList3[1,28]==0){
        titlePPI="pval < 0.001"
      } else {
        titlePPI=paste("pval =",round(geneList3[1,28],4),sep=" ")
      }
      
      list2=data.frame("Brainspan"=range01(c(as.data.frame(geneList3)[1,3],
                                             as.data.frame(geneList3)[,2])),
                       "Coexpression"=range01(c(as.data.frame(geneList3)[1,6],
                                                as.data.frame(geneList3)[,5])),
                       "Cotrans"=range01(c(as.data.frame(geneList3)[1,9],
                                           as.data.frame(geneList3)[,8])),
                       "GO"=range01(c(as.data.frame(geneList3)[1,12],
                                      as.data.frame(geneList3)[,11])),
                       "GTEx"=range01(c(as.data.frame(geneList3)[1,15],
                                        as.data.frame(geneList3)[,14])),
                       "Interpro"=range01(c(as.data.frame(geneList3)[1,18],
                                            as.data.frame(geneList3)[,17])),
                       "Kegg"=range01(c(as.data.frame(geneList3)[1,21],
                                        as.data.frame(geneList3)[,20])),
                       "Literature"=range01(c(as.data.frame(geneList3)[1,24],
                                              as.data.frame(geneList3)[,23])),
                       "PPI"=range01(c(as.data.frame(geneList3)[1,27],
                                       as.data.frame(geneList3)[,26]))
      )
      
      boxplot(list2$Brainspan,
              list2$Coexpression,
              list2$Cotrans,
              list2$GO,
              list2$GTEx,
              list2$Interpro,
              list2$Kegg,
              list2$Literature,
              list2$PPI,
              main=paste(geneList3[1,1],"Clustering",sep=" "),
              ylab="Sum of Weighted Links Within Gene Network",
              col=c(
                "red",
                "green3",
                "yellow4",
                "blue",
                "orange",
                "grey",
                "pink",
                "purple",
                "brown"
              ),
              xaxt="n")
      axis(side = 1,
           at = 0:10,
           mgp = c(0, 0.15, 0),
           labels=c(
             "",
             paste("Brainspan",titleBrain,sep="\n"),
             paste("Coexpression",titleCoex,sep="\n"),
             paste("Cotrans",titleCotrans,sep="\n"),
             paste("GO",titleGO,sep="\n"),
             paste("GTEx",titleGtex,sep="\n"),
             paste("Interpro",titleInter,sep="\n"),
             paste("Kegg",titleKegg,sep="\n"),
             paste("Literature",titleLit,sep="\n"),
             paste("PPI",titlePPI,sep="\n"),
             ""
           ),
           cex.axis=0.9,
           las=3,
           lwd.ticks = FALSE)
      stripchart(list2[1,],
                 vertical=T,
                 method="overplot",
                 add=T,
                 col="red",
                 pch=20,
                 cex=3)
    }
  })
  
  
  # biclusteringBox: Biclustering Test on Gene Sets Comparative Boxplots -------
  output$biclusteringBox <- renderPlot({
    if (input$geneExample == "Parkinson's Disease (PD)"){
      geneList4=biclustering[trait=="Parkinson's disease"]
    } else if (input$geneExample == "Alzheimer's Disease (AD)"){
      geneList4=biclustering[trait=="Alzheimer's disease"]
    } else if (input$geneExample == "Type 2 Diabetes"){
      geneList4=biclustering[trait=="Type 2 diabetes"]
    } else if (input$geneExample == "Multiple Sclerosis (MS)"){
      geneList4=biclustering[trait=="Multiple sclerosis"]
    } else if (input$geneExample == "Depression"){
      geneList4=biclustering[trait=="Depression"]
    } else if (input$geneExample == "Lewy Body Dementia (LBD)"){
      geneList4=biclustering[trait=="Dementia with Lewy bodies"]
    } else if (input$geneExample == "Schizophrenia"){
      geneList4=biclustering[trait=="Schizophrenia"]
    } else if (input$geneExample == "Progressive Supranuclear Palsy (PSP)"){
      geneList4=biclustering[trait=="Progressive supranuclear palsy"]
    }
      
    if (geneList4[1,4]==0){
      titleBiclustGen="pval < 0.001"
    } else {
      titleBiclustGen=paste("pval =",round(geneList4[1,4],4),sep=" ")
    }
    if (geneList4[1,7]==0){
      titleBiclustT2d="pval < 0.001"
    } else {
      titleBiclustT2d=paste("pval =",round(geneList4[1,7],4),sep=" ")
    }
    if (geneList4[1,10]==0){
      titleBiclustNer="pval < 0.001"
    } else {
      titleBiclustNer=paste("pval =",round(geneList4[1,10],4),sep=" ")
    }
      
    list3=data.frame("General"=range01(c(as.data.frame(geneList4)[1,3],
                                         as.data.frame(geneList4)[,2])),
                     "Metabolic"=range01(c(as.data.frame(geneList4)[1,6],
                                           as.data.frame(geneList4)[,5])),
                     "Nervous"=range01(c(as.data.frame(geneList4)[1,9],
                                         as.data.frame(geneList4)[,8]))
                     )
      
    boxplot(list3$General,
            list3$Metabolic,
            list3$Nervous,
            main=paste(geneList4[1,1],"Biclustering",sep=" "),
            ylab="Sum of Weighted Links Within Gene Network",
            col="darkblue",
              xaxt="n")
      axis(side = 1,
           at = 1:3,
           labels=c(
             paste("General",titleBiclustGen,sep="\n"),
             paste("Metabolic",titleBiclustT2d,sep="\n"),
             paste("Nervous",titleBiclustNer,sep="\n")
           ),
           lwd.ticks = FALSE)
      stripchart(list3[1,],
                 vertical=T,
                 method="overplot",
                 add=T,
                 col="red",
                 pch=20,
                 cex=3)
  })
  # exportui: Capture and Download PLN Plot ------------------------------------
  output$exportui <-  renderUI(
    
    #Validate#
#    {
#      if (input$gene_box == "textbox") {
#        validate(
#          need(nrow(appData()) > 0 && gene() >= 1, "")
#          , need(c( input$humancells, input$mousecells), '' )
#        )
#      }
#      else {
#        validate(
#          need(input$selectgene > 0 , ''), 
#          need(c( input$humancells, input$mousecells), '' )
#        )
#      }
      
      div(    
        # UI Render Download Options ########################################
        
        
        hr(),
        fluidRow(
            h4("Export Plot:")
        ),
        div(
          id = "downloader",
          fluidRow(
            # UI Render File Type #############################################
              "File Type:",
              selectInput(
                'devices',
                label = NULL,
                choices = c( "png", "pdf", "jpeg" ),
                selected = "png",
                width = "300px"
            ),
            
            # UI Render Download Button ------------------------------------------------------
              HTML("<br>"),
              downloadButton('downloadplot', 'Download Plot')
          )
        )
      )
#    }
  )
  
  # Render Download Plot ----------------------------------------------------
  output$downloadplot <- downloadHandler(
    filename = function(){
      paste(
        'ClusterBot.',
        Sys.Date(),
        '.',
        input$devices,
        sep = '')
    },
    content = function(file){
      req(visNet())
#      ggsave(
#        file,
#        plot = visNet(),
#        device = input$devices,
#        width = (as.numeric(input$dimension[1]) * 25.4)/96,
#        height = (as.numeric(input$dimension[2]) * 25.4)/96,
#        scale = 1.5,
#        units = "mm",
#        dpi = 96
#      )
    }
  )
  
  # userGenes: Separate Good Data ------------------------------------------------------
  userGenes <- reactive(
    {
      if(input$geneSelect == "textbox") {
        unique(toupper(
          (as.list(gene())[
            (as.list(gene()) %in% nodes$node | toupper(as.list(gene())) %in% toupper(nodes$node)) 
        ])
        )
        )
      }
      else if(input$geneSelect == "dropmenu") {
        if (input$geneExample == "Parkinson's disease") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
            as.vector(geneList[trait == "Parkinson's disease"]$gene)
        }
        else if (input$geneExample == "Alzheimer's disease") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(geneList[trait == "Alzheimer's disease"]$gene)
        }
        else if (input$geneExample == "Type 2 diabetes") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(geneList[trait == "Type 2 diabetes"]$gene)
        }
        else if (input$geneExample == "Multiple sclerosis") {
                    updateCheckboxGroupInput(session,
                                             "datasetCheckbox",
                                             selected=fullSel
                    )
          
          as.vector(geneList[geneList$trait == "Multiple sclerosis",]$gene)
        }
        else if (input$geneExample == "Schizophrenia") {
                    updateCheckboxGroupInput(session,
                                             "datasetCheckbox",
                                             selected=fullSel
                    )
          
          as.vector(geneList[trait == "Schizophrenia"]$gene)
        }
        else if (input$geneExample == "Breast cancer") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(geneList[trait == "Breast cancer"]$gene)
        }
        else if (input$geneExample == "Prostate cancer") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(geneList[trait == "Prostate cancer"]$gene)
        }
        else if (input$geneExample == "Asthma") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(geneList[trait == "Asthma"]$gene)
        }
        else if (input$geneExample == "Coronary artery disease") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(geneList[trait == "Coronary artery disease"]$gene)
        }
      }
      else if(input$geneSelect == "dropmenuHGMD") {
        if (input$geneExampleHGMD == "Amyotrophic lateral sclerosis") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(hgmd[trait == "Amyotrophic lateral sclerosis"]$gene)
        }
        else if (input$geneExampleHGMD == "Cognitive decline") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(hgmd[trait == "Cognitive decline"]$gene)
        }
        else if (input$geneExampleHGMD == "Dementia") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(hgmd[trait == "Dementia"]$gene)
        }
        else if (input$geneExample == "Frontotemporal dementia") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(hgmd[trait == "Frontotemporal dementia"]$gene)
        }
        else if (input$geneExampleHGMD == "Multiple system atrophy") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(hgmd[trait == "Multiple system atrophy"]$gene)
        }
        else if (input$geneExampleHGMD == "Neuronal ceroid lipofuscinosis") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(hgmd[trait == "Neuronal ceroid lipofuscinosis"]$gene)
        }
        else if (input$geneExampleHGMD == "Parkinson's disease") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(hgmd[trait == "Parkinson's disease"]$gene)
        }
        else if (input$geneExampleHGMD == "Spinal muscular atrophy") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(hgmd[trait == "Spinal muscular atrophy"]$gene)
        }
        else if (input$geneExampleHGMD == "Spinocerebellar ataxia") {
          updateCheckboxGroupInput(session,
                                   "datasetCheckbox",
                                   selected=fullSel
          )
          
          as.vector(hgmd[trait == "Spinocerebellar ataxia"]$gene)
        }
      }
      else if (input$geneSelect == "doubleGeneset") {
       updateCheckboxGroupInput(session,
                                 "datasetCheckbox",
                                 selected=fullSel
        )
        as.vector(geneList[trait == input$doubleGeneset1 | trait == input$doubleGeneset2]$gene)
        }
      else if (input$geneSelect == "upload") {
        req(input$geneUpload)
        
        updateCheckboxGroupInput(session,
                                 "datasetCheckbox",
                                 selected=fullSel
        )
        
        tryCatch(
          {
            df <- data.table::fread(
              input$geneUpload$datapath,
              data.table=TRUE
                                    )
          },
          error = function(e) {
            # return a safeError if a parsing error occurs
            stop(safeError(e))
          }
        )
        
        if(dim(df)[1]>0) {
          as.data.frame(df)[,1]
        }
        else if(dim(df)[1]==0) {
          as.vector(colnames(df))
        }
      }
    }
  )
  
  # Seed Selection -------------------------------------------------------------
  
  seedNum <- eventReactive(input$seedShuffle,{
    return(sample(1:100,1,replace=TRUE))
  },ignoreNULL=FALSE)
  
  seed <- reactive({
      if (is.null(input$advancedOptions)){
        return(seedNum())
        } else if (!is.null(input$advancedOptions)){
          return(input$seedNumInput)
          }
    })
  
  output$seedNumber <- renderText({
    #    if (is.na(seedNumInput()) & is.na(seedInput())){
    return(paste("Seed Number: ",seed(),sep=""))
  })
   
  # sessionInfo: Show Session Information --------------------------------------
  
  output$sessionInfo <- renderTable(
    data.frame("GenesMemory"=paste(sprintf(as.numeric(object_size(nodes))/1e+06,fmt='%#.3f'),"MB",sep=" "),
               "RelationshipsMemory"=paste(sprintf(as.numeric(object_size(edges))/1e+06,fmt='%#.3f'),"MB",sep=" "),
               "TotalMemory"=paste(sprintf(as.numeric(mem_used())/1e+06,fmt='%#.3f'),"MB",sep=" "))
#               "GraphMemory"=as.character(object_size(visNet())))
    )
  
  # VisNetwork Example Plot Section ----------------------------------------------
  
  output$visNet <- renderVisNetwork({ # Output for creating the visNetwork PLN
    # Changing Input Datasets and Number of Nodes Displayed ====================
#    req(input$datasetCheckbox)
    if (is.null(input$datasetCheckbox))
      return()
    
    if (input$plnSelect=="General"){
      nodes=nodes[pln=="General"]
      edges=edges[pln=="General"]
    }
    else if (input$plnSelect=="T2D"){
      nodes=nodes[pln=="T2D"]
      edges=edges[pln=="T2D"]
    }
    else if (input$plnSelect=="Nervous"){
      nodes=nodes[pln=="Nervous"]
      edges=edges[pln=="Nervous"]
    }
    
    #nodes=as.data.frame(nodes) # Ensure input nodes are in data.frame format
    #edges=as.data.frame(edges) # Ensure input edges are in data.frame format
    
    dataIn=NULL
    dataIn=as.vector(selectDataset())
    
    genes=as.vector(userGenes())
    
    inputPLN=input$plnSelect
    
    if (seed()=="") {
      seed=1
      } else {
        seed=seed()
      }
    
    genes2=geneList[geneList$trait == input$doubleGeneset2,]$gene
    
    # Customise Nodes and Edges based on Selections ============================
    # Completely subset neighbour pick
    
    if (length(genes)==0) {
      vis.nodes=nodes[pln==inputPLN]
      vis.edges=edges[pln==inputPLN]
      vis.edges=vis.edges[dataset.id %in% dataIn]
    }
    else if (length(genes)==1 & input$nodeNeighbour==0) {
      vis.nodes=nodes[pln==inputPLN]
      vis.edges=edges[pln==inputPLN]
      genesID=vis.nodes[toupper(node) %in% toupper(genes)]$id
      genesID2=vis.nodes[toupper(vis.nodes$node) %in% toupper(nodes[nodes$node %in% genes2,]$node),]$id
      vis.nodes=vis.nodes[id %in% genesID]
      vis.nodes[id %in% genesID]$node.type=2
      vis.nodes[vis.nodes$id %in% genesID2,]$node.type=3
      vis.nodes[!(vis.nodes$id %in% genesID | vis.nodes$id %in% genesID2),]$node.type=1
      vis.nodes[!(id %in% genesID)]$node.type=1
      vis.edges=vis.edges[!(from %in% genesID)]
      vis.edges=vis.edges[!(to %in% genesID)]
      vis.edges=vis.edges[dataset.id %in% dataIn]
    }
    else {
      if (input$nodeNeighbour==0) {
        vis.nodes=nodes[pln==inputPLN]
        vis.edges=edges[pln==inputPLN]
        genesID=vis.nodes[toupper(node) %in% toupper(genes)]$id
#        genesID2=vis.nodes[toupper(vis.nodes$node) %in% toupper(nodes[nodes$node %in% genes2,]$node),]$id
        vis.nodes[id %in% genesID]$node.type=2
#        vis.nodes[vis.nodes$id %in% genesID2,]$node.type=3
#        vis.nodes[!(vis.nodes$id %in% genesID | vis.nodes$id %in% genesID2),]$node.type=1
        vis.nodes[!(id %in% genesID)]$node.type=1
        vis.nodes=vis.nodes[id %in% genesID]
        vis.edges=vis.edges[from %in% genesID]
        vis.edges=vis.edges[to %in% genesID]
        vis.edges=vis.edges[dataset.id %in% dataIn]
      }
      else {
        vis.nodes=nodes[pln==inputPLN]
        vis.edges=edges[pln==inputPLN]
        genesID=vis.nodes[toupper(node) %in% toupper(genes)]$id
#        genesID2=vis.nodes[toupper(vis.nodes$node) %in% toupper(nodes[nodes$node %in% genes2,]$node),]$id
        vis.nodes[id %in% genesID]$node.type=2
#        vis.nodes[vis.nodes$id %in% genesID2,]$node.type=3
#        vis.nodes[!(vis.nodes$id %in% genesID | vis.nodes$id %in% genesID2),]$node.type=1
        vis.nodes[!(id %in% genesID)]$node.type=1
        vis.edges=vis.edges[dataset.id %in% dataIn]
        vis.edges=vis.edges[from %in% genesID | to %in% genesID]
        vis.nodes=vis.nodes[id %in% vis.edges$from | id %in% vis.edges$to]
        if (as.numeric(input$nodeNeighbour) > 1) {
          rep=1
          repeat {
            vis.edges=edges[edges[pln==inputPLN]$from %in% vis.edges$from | edges[pln==inputPLN]$to %in% vis.edges$to]
            vis.nodes=nodes[nodes[pln==inputPLN]$id %in% vis.edges$from | nodes[pln==input$plnSelect]$id %in% vis.edges$to]
            vis.nodes[id %in% genesID]$node.type=2
#            vis.nodes[vis.nodes$id %in% genesID2,]$node.type=3
#            vis.nodes[!(vis.nodes$id %in% genesID | vis.nodes$id %in% genesID2),]$node.type=1
            vis.nodes[!(id %in% genesID)]$node.type=1
            vis.edges=vis.edges[dataset.id %in% dataIn]
            rep=rep+1
            if (rep==as.numeric(input$nodeNeighbour)) {
              break
            }
          }
        }
      }
    }
    
    if (length(dataIn)<1){
      vis.nodes=NULL
      vis.edges=NULL
    }
    else {
      
#      vis.nodes=vis.nodes[vis.nodes$id %in% vis.edges$from | 
#                            vis.nodes$id %in% vis.edges$to,]
    #set.seed(1)
#    rows=sample(nrow(vis.nodes))
#    vis.nodes=vis.nodes[rows,]
    
    vis.nodes$color.background=c(rgb(30,144,255,180,maxColorValue=255),
#                                 rgb(11,102,35,255,maxColorValue=255),
                                 rgb(255,102,0,255,maxColorValue=255)
    )[vis.nodes$node.type]
    
    vis.nodes$borderWidth=c(1,3,3)[vis.nodes$node.type]
    
#    vis.edges$width=1
#    vis.edges[vis.edges$weight >= 0.5,]$width=2
    
    vis.edges$width=rescale(vis.edges$weight,to=c(3,7),from=c(0,1))
    
    vis.edges$color=c("red",
                      "green",
                      "olive",
                      "blue",
                      "orange",
                      "grey",
                      "pink",
                      "purple",
                      "brown"
#                      "black"
                      )[vis.edges$dataset.id]
    
    #vis.nodes$size=sum(vis.nodes$node %in% vis.edges)
     
    #freqSize=as.data.frame(table(data.frame(stack(vis.edges[,1:2]))))
  
    
    #vis.nodes$size=
       
    # Shape of node depends on node type 
        vis.nodes$title  <- vis.nodes$node # Text on click
        vis.nodes$label  <- vis.nodes$node # Node label # Node size
        vis.nodes$font.size <- 25
#        vis.nodes$font.vadjust <- 10
        
        
        vis.nodes=vis.nodes[order(vis.nodes$node.type,decreasing=TRUE),]
        vis.nodes=head(vis.nodes,n=nodeSlider())

        vis.edges=vis.edges[vis.edges$from %in% vis.nodes$id | vis.edges$to %in% vis.nodes$id,]
        vis.nodes=vis.nodes[vis.nodes$id %in% vis.edges$from | vis.nodes$id %in% vis.edges$to,]
        
        datasets=c("Brainspan",
                   "Coexpression_Array",
                   "Cotrans",
                   "GO",
                   "GTEx",
                   "Interpro",
                   "Kegg_Reactome",
                   "Literature",
                   "PPI")
        
    # Colour of node depends on node type
        vis.nodes$color.border <- c("black",
                                    #"yellow",
                                    "yellow")[vis.nodes$node.type] # Nodes will have black border
        vis.nodes$color.highlight.background <- c("orange")[vis.nodes$node.type] # Nodes will appear orange when highlighted
        vis.nodes$color.highlight.border <- c("darkred")[vis.nodes$node.type] # Node border will appear dark red when highlighted
#        vis.nodes$group=nodes$vis.nodes.label # Nodes are grouped by node label (GENE, PROTEIN etc.)
    
    if (nrow(vis.edges)>1) {
      vis.edges$title <- paste(paste("Weight: ",vis.edges$weight,sep=""),
                             paste("Dataset: ",datasets[vis.edges$dataset.id],sep=""),
                             sep="<br />") # what the line displays on click
    }
    
    # Plotting pln ##############################################################
        
    visNetwork(nodes=vis.nodes,edges=vis.edges,height="200%",width="100%") %>%
      visOptions(autoResize=TRUE,highlightNearest=TRUE) %>% # Select by node label (e.g. GENE)
      visLayout(randomSeed=seed) %>%
      visInteraction(dragNodes=FALSE) %>%
      visPhysics(stabilization=TRUE,solver = "barnesHut") %>%
#                 forceAtlas2Based = list(gravitationalConstant = -100)) %>%
      visEdges(
        physics=FALSE,
        smooth=
          list(
          enabled=TRUE,
          type="discrete",
          roundness="1"),
        length=350,
        shadow=FALSE,
        arrows=NULL
               ) %>%
      visNodes(shadow=TRUE,
               shape="dot")
    }
  })
  
  
  
  observeEvent(input$store_position, {
    
    
    visNetworkProxy("visNet") %>% 
      visOptions(autoResize=TRUE,highlightNearest=TRUE) %>% # Select by node label (e.g. GENE)
      visLayout(randomSeed=seed) %>%
      visInteraction(dragNodes=FALSE) %>%
      visPhysics(stabilization=TRUE,solver = "barnesHut") %>%
      visEdges(
        physics=FALSE,
        smooth=
          list(
            enabled=TRUE,
            type="discrete",
            roundness="1"),
        length=350,
        shadow=FALSE,
        arrows=NULL
      ) %>%
      visNodes(shadow=TRUE,
               shape="dot") %>%
      visGetPositions()
  })
  
  nodes_positions <- reactive({
    positions <- input$visNet_positions
    if (!is.null(positions)){
      nodes_positions <- do.call("rbind", lapply(positions, function(x){ data.frame(x = x$x, y = x$y)}))
      nodes_positions$id <- names(positions)
      nodes_positions
    } else {
      NULL
    }
  })
  
  output$downloadNetwork <- downloadHandler(
    filename = function() {
      paste('network-', Sys.Date(), '.html', sep='')
    },
    content = function(con) {
      nodes_positions <- nodes_positions()
      if(!is.null(nodes_positions)){
        nodes_save <- merge(nodes, nodes_positions, by = "id", all = T)
      } else  {
        nodes_save <- nodes
      }
      
      visNetwork(nodes = nodes_save, edges = edges, height = "800px") %>%
        visOptions(highlightNearest = TRUE) %>% visExport() %>%
        visPhysics(enabled = FALSE) %>% visEdges(smooth = FALSE) %>% visSave(con)
    }
  )

# Render Download Plot ----------------------------------------------------
output$downloadplot <- downloadHandler(
  filename = function(){
    paste(
      '`PLaNet.'
      , Sys.Date()
      , '.'
      , input$devices
      , sep = '')
  }
  , content = function(file){
    req(myplot())
    ggsave(
      file
      , plot = myplot()
      , device = input$devices
      , width = (as.numeric(input$dimension[1]) * 25.4)/96
      , height = (as.numeric(input$dimension[2]) * 25.4)/96
      , scale = 1.5
      , units = "mm"
      , dpi = 96
    )
  }
)
}