#
#
packages <- c('shiny', 'BiocManager', 'ggforce', 'yaml', 'scatterD3', 'EDASeq', 'DT', 'EnhancedVolcano', 'tools',
              'ComplexHeatmap', 'gridExtra', 'circlize', 'openxlsx', 'reshape2', 'kableExtra', 'rlang', 'shinydashboard',
              'DESeq2', 'tidyverse', 'fgsea', 'shinyFiles', 'cowplot')
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

library(shiny)
library(BiocManager)
options(repos = BiocManager::repositories())
library(ggforce)
#library(grid)
#library(ggpubr)
library(yaml)
library(scatterD3)
library(EDASeq)
library(DT)
library(EnhancedVolcano)
library(tools)
library(ComplexHeatmap)
library(gridExtra)
#library(gtable)
library(circlize)
library(openxlsx)
library(reshape2)
library(kableExtra)
library(rlang)
library(shinydashboard)
#library(clusterProfiler)
library(DESeq2)
library(tidyverse)
library(fgsea)
library(shinyFiles)
library(cowplot)

####


#library(plyr)
#source("/yerkes-cifs/runs/home/gktharp/gitrepos/Coding_Club/BulkRNASeq/gktharp/GKTDESeq2Extentions.R")
source("functions.R")


#library(shiny)

options(shiny.maxRequestSize = 99*1024^2)
ui <- fluidPage(
  
  # App title ----
  #titlePanel("Uploading Files"),
  
  navbarPage("Bulk RNA-seq analysis",
             ## Setup ----
             navbarMenu("Upload",
                        # Load in Page ----
                        tabPanel("Analysis Data",
                                 fileInput("analysis_obj", "Upload Analysis.RData object",
                                           multiple=FALSE, accept=c(".RData", ".rdata", ".Rdata"))
                        ),
                        tabPanel("Gene list",
                                 fileInput("genelist_file", "Upload genelist",
                                           multiple = FALSE,
                                           accept = c("text/csv",
                                                      "text/comma-separated-values,text/plain",
                                                      ".csv")),
                                 #checkboxInput("header", "Header", TRUE),
                                 radioButtons("genelist_sep", "Separator",
                                              choices = c(Comma = ",",
                                                          Semicolon = ";",
                                                          Tab = "\t",
                                                          Space = " ",
                                                          Newline = "\n"),
                                              selected = " "),
                                 actionButton("load_genelist_file", "Load genelist file")
                        )
             ),
             #     sidebarLayout(
             #       sidebarPanel(
             #         tags$h5("Metadata"),
             #         tags$hr(),
             #         fileInput("file2", "Choose config yaml",
             #              multiple = TRUE,
             #              accept = c(".yaml",
             #                         "text/comma-separated-values,text/plain",
             #                         ".yml")),
             #         tags$hr(),
             #         fileInput("file1", "Choose metadata CSV",
             #              multiple = TRUE,
             #              accept = c("text/csv",
             #                         "text/comma-separated-values,text/plain",
             #                         ".csv")),
             #         checkboxInput("header", "Header", TRUE),
             #         radioButtons("sep", "Separator",
             #                 choices = c(Comma = ",",
             #                             Semicolon = ";",
             #                             Tab = "\t",
             #                             Space = " "),
             #                 selected = " "),
             #         tags$hr(),
             #         radioButtons("disp", "Display",
             #                 choices = c(Head = "head",
             #                             All = "all"),
             #                 selected = "head")
             #         
             #       ),
             #       mainPanel(
             #         dataTableOutput("exp.design")
             #       )
             #     )  
             #   ),
             #   tabPanel("Config",
             #     sidebarLayout(
             #       sidebarPanel(
             #         tags$h5("Config"),
             #         tags$hr(),
             #         fileInput("file2", "Choose yaml",
             #              multiple = TRUE,
             #              accept = c(".yaml",
             #                         "text/comma-separated-values,text/plain",
             #                         ".yml")),
             #         tags$hr()
             # 
             # 
             #       ),
             #       mainPanel(
             #         uiOutput("config")
             #       )
             #     )
             #   )
             # ),
             ## Samples ----
             tabPanel("QC",
                      navlistPanel(
                        tabPanel("Read counts",
                                 plotOutput("mapping")
                        ),
                        tabPanel('RLE normalization',
                                 plotOutput('RLE_raw'),
                                 plotOutput('RLE_norm')
                        ),
                        tabPanel("PCA",
                                 plotOutput("PCA.plot")
                        ),
                        widths=c(2,10)
                      )
             ),
             tabPanel("DGE Analysis",
                      fluidRow(
                        sidebarLayout(
                          sidebarPanel(
                            uiOutput("SelectGroup1"),
                            # Horizontal line ----
                            
                            uiOutput("SelectGroup2"),
                            tags$hr(),
                            actionButton("pairwisecomp", "Run Analysis"),
                            width=3
                          ),
                          mainPanel(
                            verbatimTextOutput("summary_out"),
                            width=9
                          )
                        )
                      ),
                      fluidRow(
                        column(12,
                               scatterD3Output("volcano_plot")
                        )
                        #                         column(4,
                        #                                # fluidRow(
                        #                                #   tableOutput("click_selected")
                        #                                # ),
                        #                                fluidRow(
                        #                                  # DTOutput("lassod.points"),
                        #                                  #tags$head(tags$style("#lassod.points{color:red; font-size:12px; font-style:italic; 
                        # #overflow-y:scroll; max-height: 50px; background: ghostwhite;}"))
                        #                                ),
                        #                                fluidRow(
                        #                                  # uiOutput("add.lassod.button")#actionButton("add.lassod.points", "Add selected genes to genelist")
                        #                                )
                        
                      ),
                      fluidRow(
                        checkboxInput("view_lassod_genes", "View lassod genes in table"),
                        # actionButton("reset_comptable", "Reset table")
                      ),
                      fluidRow(
                        DTOutput('comptable.rendered')
                      ),
                      fluidRow(
                        verbatimTextOutput("selected.genes"),
                        uiOutput("add.selected.button")
                      )
             ),
             tabPanel("Genes of interest",
                      tabsetPanel(
                        tabPanel("Boxplot",
                                 sidebarLayout(
                                   sidebarPanel(
                                     textInput("genename", "Gene name", value="Gdf15"),
                                     actionButton("renderGenePlot", "Go")
                                   ),
                                   mainPanel(
                                     plotOutput("gene.comparison")
                                   )
                                 )
                        ),
                        tabPanel("Heatmap",
                                 sidebarLayout(
                                   sidebarPanel(
                                     fileInput("heatmap.genes", "Upload a genelist", multiple=FALSE,
                                               accept=c(".txt","text/comma-separated-values,text/plain")),
                                     tags$p("or"),
                                     checkboxInput("use.genelist", "Use my genelist"),
                                     tags$hr(),
                                     actionButton("heatmap.render", "Render"),
                                     downloadButton("download.heatmap", "Save heatmap")
                                   ),
                                   mainPanel(
                                     plotOutput("heatmap.plot")
                                   )
                                 )
                        ),
                        tabPanel("My genelist",
                                 mainPanel(
                                   fluidRow(
                                     column(5,
                                            tags$h3("Your genelist"),
                                            tags$h6("Upload genes from files in the upload tab")
                                     ),
                                     column(2),
                                     column(5,
                                            tags$h4("Add to genelist")
                                     )
                                   ),
                                   fluidRow(
                                     column(5,
                                            dataTableOutput("genelisttable"),
                                            uiOutput('remove.selected.button'),
                                            downloadButton("download.genelist", "Save genelist")
                                     ),
                                     column(2),
                                     column(5,
                                            textInput('manual_genelist_addition', "Search genes to add"),
                                            actionButton('search_genelist', 'Search'),
                                            tags$hr(),
                                            dataTableOutput("searched_genes_table"),
                                            actionButton("add_selected_searched_genes", "Add selected genes")
                                     )
                                   )
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                   
                                 )
                        )
                      )
             ),
             tabPanel("GSEA",
                      tabsetPanel(
                        tabPanel("Overview",
                                 sidebarLayout(
                                   sidebarPanel(
                                     uiOutput('gsea.collections.UI'),
                                     uiOutput('gsea.group1.UI'),
                                     uiOutput('gsea.group2.UI'),
                                     tags$hr(),
                                     
                                     actionButton("gsea.render", "Render")
                                   ),
                                   mainPanel(
                                     textOutput("GSEA_textsum"),
                                     plotOutput("gsea.toppaths"),
                                     plotOutput("GSEA_overview_dotplot")
                                   )
                                 )
                        ),
                        tabPanel(
                          "Collection", 
                          
                          fluidRow(
                            sidebarLayout(
                              sidebarPanel(
                                uiOutput("select.gsea.geneset"),
                                textOutput('gsea.groups')
                              ),
                              mainPanel(
                                plotOutput("gsea.geneset.plot")
                              )
                            )
                          ),
                          fluidRow(
                            column(8, align='center', offset=2,
                                   h3("Leading edge genes"),
                                   h6("Expression data for the leading genes for this comparison")
                            )
                          ),
                          fluidRow(
                            uiOutput("gsea.genes.button")
                          ),
                          fluidRow(
                            dataTableOutput("gsea.LE_table")
                          )
                        ),
                        tabPanel(
                          'x-compare',
                          sidebarLayout(
                            sidebarPanel(
                              uiOutput("gsea.xcomp.inputs.UI"),
                              uiOutput("gsea.xcomp.collection.UI"),
                              actionButton("gsea.xcomp.render", "Render")
                            ),
                            mainPanel(
                              plotOutput("gsea.xcomp.plot")
                            )
                          )
                        ),
                        tabPanel("Heatmaps",
                                 sidebarPanel(
                                   uiOutput('gsea.heatmap.collection.UI'),
                                   #checkboxInput('gsea.heatmap.prune', 'Only show most variable genes'),
                                   numericInput('gsea.heatmap.genes.N', 'Only show N most variable genes',30),
                                   actionButton('gsea.heatmap.render', 'Render')
                                 ),
                                 mainPanel(
                                   plotOutput('gsea.heatmap.plot')
                                 )
                                 
                                 
                                 #selected.gsea.geneset.heatmap
                        )
                      )
             )
  )
)


