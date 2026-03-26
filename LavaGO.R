#########################
#        LavaGO         #
#      Sept 2024        #
#########################

####################################################
#enter the path to CSV files containing FC and Pval#
####################################################
path2files<-("C:/path/to/your/files/")  #enter the path of a folder containing your expression data


##############################################################
#uncomment and run only once to install Bioconductor packages#
#delete or comment afterward                                 #
##############################################################
# if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# install.packages("remotes")
# remotes::install_github("PNNL-Comp-Mass-Spec/MSnSet.utils")
# BiocManager::install("clusterProfiler")
# BiocManager::install("MSnID")
# BiocManager::install("org.Hs.eg.db")

#in Rstudio (https://posit.co/download/rstudio-desktop/), all the following packages should install automatically. 
#you might need to install Rtools as well for some of them  https://cran.r-project.org/bin/windows/Rtools/


#BiocManager::install("ReactomePA")

R.utils::setOption("clusterProfiler.download.method","auto")
library(plotly)
library(ggiraph)
library(shiny)
library(shinyWidgets)
library(data.table)
library(shinyalert)
library(shinycssloaders)
require(htmlwidgets)
library(dplyr)
library(DT)
library(colourpicker)
library(ggridges)
library(org.Hs.eg.db)
library(MSnID)
library(dplyr)
library(clusterProfiler)
library(ggplot2)
library(calibrate)
library(tibble)
library(ggupset)  #??
library(enrichplot)
library(DOSE)


library(enrichR)
library(stringr)
# library('biomaRt')
# library(ReactomePA)


#set list of Enricher libraries to build the interface
websiteLive <- getOption("enrichR.live")
allEnricherDBS<-NULL
if (websiteLive) {
  # listEnrichrSites()
  setEnrichrSite("Enrichr") # Human genes 
  dbs <- listEnrichrDbs()
  allEnricherDBS<-str_sort(dbs$libraryName)
}


ui <- fluidPage(
  
  tabsetPanel(
    
    tabPanel("Volcano plot",
  
      sidebarLayout(
        sidebarPanel(
          br(),
          tags$head(tags$style("#DataPlot{width:100 !important;}")),
          
          selectInput("CCfile",
                      label = "File",
                      choices = list.files(path2files, pattern = "\\.csv$"),
                      multiple = F
          ),
          
          numericInput(inputId = "pvalue_threshold",
                       label = "Set pvalue threshold",
                       min = 0,
                       max= 1,
                       value = 0.05,
                       step = 0.01,
                       ),
          
          
          # set logfc threshold
          uiOutput("logfc_slider"),
          
          
          #show/hide up or down expressed genes
          checkboxInput("show_up_genes",
                        "Show up-regulated genes",
                        value = TRUE),
          checkboxInput("show_down_genes",
                        "Show down-regulated genes",
                        value = TRUE),
          
          #show/hide annotations
          checkboxInput("show_labels",
                        "Show labels",
                        value = TRUE),
          
          
          # show/hide logfc lines
          checkboxInput("show_logfc_threshold",
                        "Show fold change threshold line",
                        value = TRUE),
          
          # show/hide logfc and pval line
          checkboxInput("show_pvalue_threshold",
                        "Show significance threshold line",
                        value = TRUE),
          
          checkboxInput("color_by_de",
                        "Color significantly different features",
                        TRUE),
          
          
          sliderInput("labelsize","Label size", 0, 10, 3, step = 0.5),
          sliderInput("ptsize","Dot size", 0, 10, 1, step = 0.5),
          sliderInput("transp","Dot transparency",0,100,50),
          
          # gene selector menu
          uiOutput("gene_selector"),
          
          div(style = "display: flex;",
              colourInput("color1", "enriched", "#45FA2D"),
              colourInput("color2", "not enriched", "orangered"),
          ),
          actionButton("resetcol","reset", style = "height: 30px"),
          
          
          
          width = 2
        ),
        mainPanel(
          
          plotlyOutput(outputId = "VolcanoPlot", height = 900, width = "100%")%>% withSpinner(color="lightblue"),
          
          checkboxInput("filtervolcanotable",label = "show significant entries only", value = TRUE),
          DTOutput(outputId = "VolcanoTable"),
          downloadButton(outputId = 'downloadTableVolcano',"Download table as CSV"),

          width = 10
          
        )
      )
    ), # end of volcano plot

    tabPanel("GSEA",
      sidebarLayout(
        sidebarPanel(
          br(),
          
          selectInput("plottype", 
                      label = "Plot type",
                      choices = c(
                                  "ridge plot",
                                  "dot plot",
                                  "Enrichment Map",
                                  "Category Netplot",
                                  # "heat plot",
                                  "upset plot"
                      ),
                      multiple = F,
                      selected = "ridge plot"
          ),
          
          
          selectInput("ontology", 
                      label = "Select ontology",
                      choices = c("All GO"            ="ALL",
                                  "Biological process"="BP",
                                  "Molecular function"="MF",
                                  "Cellular component"="CC",
                                  "KEGG"              ="KEGG",
                                  # "KEGG module"       ="mKEGG",
                                  "WikiPathways"      ="WP"
                      ),
                      multiple = F,
                      selected = "BP"
          ),
          
          selectInput("ranking", 
                      label = "Gene ranking method",
                      choices = c("-log10(p-value)*sign(logFC)"="log",
                                  "logFC"="logFC",
                                  "T test"="T"
                                  
                      ),
                      multiple = F
          ),
          
          
          numericInput(inputId = "pvalue_threshold_GSEA",
                       label = "Set pvalue threshold",
                       min = 0,
                       max= 1,
                       value = 0.05,
                       step = 0.01,
          ),
          
          
          
          
          selectInput("pAdjustMethodGSEA", 
                      label = "Pvalue Adjust. Method",
                      choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                      selected = "BH",
                      multiple = F
          ),
          
          
          numericInput("minGsize", label = "minimum gene set size", value = 10),
          numericInput("maxGsize", label = "maximum gene set size", value = 500),
          numericInput("perm", label = "number of permutations", value = 10000),   #50000

          numericInput(inputId = "showcat", label = "Show categories", value = 20),
          
          textAreaInput("gocat", "show categories in plot",rows = 5, value = NULL ),

          width = 2
        ),
        
        mainPanel(
          
          htmlOutput(outputId = "commentGSEA", inline = TRUE),
          htmlOutput(outputId = "nomatchGSEA", inline = TRUE),

          conditionalPanel(
            condition = "input.plottype != 'ridge plot'",
            plotOutput(outputId = "gseaPlot", height = 1000, width = "100%")%>% withSpinner(color="lightblue"),          #height = 900, 
          ),
          conditionalPanel(
            condition = "input.plottype == 'ridge plot'",
            girafeOutput(outputId = "gseaPlotridge", height = 1000, width = "100%") %>% withSpinner(color="lightblue"),  #height = 900, 
          ),
         
         plotOutput(outputId = "gseaPlot2", height = 300, width = "100%"),
         
         downloadButton(outputId = 'downloadTablegsea',"Download table as CSV"),
         downloadButton(outputId = 'downloadranklistgsea',"Download pre-ranked gene list as CSV"),
         
         DTOutput(outputId = "gseaTable", width = 2000),
           
         width = 10
        )
      )
    ),# end of GSEA
    
    
    
    ##################################################################ORA#######################################################
    tabPanel("ORA",
             sidebarLayout(
               sidebarPanel(
                 br(),
                 
                 selectInput("plottype2", 
                             label = "Plot type",
                             choices = c(
                               "dot plot",
                               "Enrichment Map",
                               "Category Netplot",
                               "upset plot"
                             ),
                             multiple = F,
                             selected = "dot plot"
                 ),
                 
                 
                 selectInput(inputId = "ontology2", 
                             label = "Select ontology",
                             choices = c("All GO"            ="ALL",
                                         "Biological process"="BP",
                                         "Molecular function"="MF",
                                         "Cellular component"="CC",
                                         "KEGG"              ="KEGG",
                                         "Pathway Commons"   ="PC"
                             ),
                             multiple = F,
                             selected = "BP"
                 ),
                 
                 
                 numericInput(inputId = "pvalue_threshold_ORA",
                              label = "Set pvalue threshold",
                              min = 0,
                              max= 1,
                              value = 0.05,
                              step = 0.01,
                 ),
                 
                 
                 selectInput("pAdjustMethodORA", 
                             label = "Pvalue Adjust. Method",
                             choices = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"),
                             selected = "BH",
                             multiple = F
                 ),
                 
                 
                 numericInput("minGSizeORA", label = "minimum gene set size", value = 10),
                 numericInput("maxGSizeORA", label = "maximum gene set size", value = 500),
                 
                 
                 numericInput(inputId = "showcat2", label = "Show categories", value = 20),
                 
                 textAreaInput("gocat2", "show categories in plot",rows = 5, value = NULL ),
                 
                 
                 actionButton("refresh", "Refresh"),
                 
                 width = 2
               ),
               
               mainPanel(
                 htmlOutput(outputId = "commentORA"),
                 htmlOutput(outputId = "nomatchORA", inline = TRUE),
                 
                 plotOutput(outputId = "ORAPlot", height = 1000, width = "100%")%>% withSpinner(color="lightblue"),            #height = 900, 
                 

                 downloadButton(outputId = 'downloadTableora',"Download table as CSV"),
                 DTOutput(outputId = "ORAtable", width = 2000),
                 
                 
                 width = 10
               )

             )
    ),#################################################################### end of ORA##########################################################
    
  )
)

server <- function(input, output, session) {
  
  adata     <- reactiveVal(NULL)
  gseadata  <- reactiveVal(NULL)
  ranges    <- reactiveValues(x = NULL, y = NULL)
  tabledatavolcano <- reactiveVal(NULL)
  tabledatagsea    <- reactiveVal(NULL)
  tabledataora     <- reactiveVal(NULL)
  geneList         <- reactiveVal(NULL)

  observeEvent(input$resetcol, {
    updateColourInput(session, "color1", value = "#45FA2D")
    updateColourInput(session, "color2", value = "orangered")
  })
  
  #to dowload the table content
  output$downloadTableVolcano <- downloadHandler(
    filename = function() {
      paste0( sub('\\..*$', '', basename(input$CCfile)),"-Volcano-Selection-table-",Sys.time(), ".csv")
    },
    
    content = function(fname){
      if(! is.null(tabledatavolcano()) )
        write.csv(tabledatavolcano(), fname, row.names = F)
    }
  )
  #to dowload the table content
  output$downloadTablegsea <- downloadHandler(
    filename = function() {
      paste0( sub('\\..*$', '', basename(input$CCfile)),"-GSEA-table-",Sys.time(), ".csv")
    },
    
    content = function(fname){
      if(! is.null(tabledatagsea()) )
        write.csv(tabledatagsea(), fname, row.names = F)
    }
  )
  
  
  #to dowload the table content
  output$downloadTableora <- downloadHandler(
    filename = function() {
      paste0( sub('\\..*$', '', basename(input$CCfile)),"-ORA-table-",Sys.time(), ".csv")
    },
    
    content = function(fname){
      if(! is.null(tabledataora()) )
        write.csv(tabledataora(), fname, row.names = F)
    }
  )
  
  
  #to download the pre-ranked gene list
  output$downloadranklistgsea <- downloadHandler(
    filename = function() {
      paste0( sub('\\..*$', '', basename(input$CCfile)),"-GSEA-ranked-gene-list-",Sys.time(), ".txt")
    },
    
    content = function(fname){
      if(! is.null(geneList()) )
        write.table(x = geneList(), file = fname, row.names = T, sep = "\t", quote = F, col.names = F)
    }
  )
  
  # select genes to highlight
  output$gene_selector <- renderUI({
    data<-adata()
    selectInput("selected_genes",
                "Select gene(s) to label",
                sort(data$GeneSymbol),
                multiple = TRUE,
                selectize= TRUE)
    
  })
  
  
  output$logfc_slider <- renderUI({
    
    data<-adata()

    if(is.null(data)){  return(NULL)  }
    
    
    
    numericInput("logfc_threshold",
                 "Set log fold change threshold",
                 min = 0,
                 max = max(data$logFC), #round(max(data[[input$logfc_col]])),
                 value = 1,
                 step = 1
    )
    
  })
  
  
  observeEvent(input$CCfile, {
    #update Menu
    file<-paste0(path2files,input$CCfile)
    
    #reset values
    adata(NULL)
    gseadata(NULL)
    tabledatavolcano(NULL)
    tabledatagsea(NULL)
    tabledataora(NULL)
    geneList(NULL)
    

    data <- read.csv(file, header=TRUE,stringsAsFactors = F)
    
    
    
    
    #colnames(data)[which(names(data) == "log2FoldChange")] <- "logFC"
    
    
    ###############for RNAseq files, col names are different####################
    colnames(data)[which(names(data) == "stat")]           <- "t"
    colnames(data)[which(names(data) == "pvalue")]         <- "P.Value"
    colnames(data)[which(names(data) == "padj")]           <- "adj.P.Val"
    colnames(data)[which(names(data) == "log2FoldChange")] <- "logFC"
    colnames(data)[which(names(data) == "X")]              <- "GeneSymbol"
    ############################################################################
    
    
    
    #check if we have gene symbols. find them from ENSEMBL IDs otherwise
    if(!"GeneSymbol" %in% colnames(data)){
      if("Gene" %in% colnames(data)){
        if(grepl('ENS', data$Gene[1])  ){
          showNotification("Gene names are missing. Converting ENSEMBL IDs to gene symbols.")

          genes<-data$Gene
          genes<-sub("\\..*", "", genes) #removes gene id version
          data$Gene2<-genes              #add col with new ENSEMBL IDs

          conv_tbl <- AnnotationDbi::select(org.Hs.eg.db, columns = c('SYMBOL'), keytype = 'ENSEMBL', keys = data$Gene2) #converts IDs to names

          data <- data %>%
            mutate(ENSEMBL = Gene2) %>%
            left_join(conv_tbl, by = "ENSEMBL", multiple = "first")

          colnames(data)[which(names(data) == "SYMBOL")] <- "GeneSymbol" #rename col for compatibility
          data$GeneSymbol[is.na(data$GeneSymbol)] <- ''  #removes NAs
        }
      }
    }
      


    # convert pval to -log10(pval)
    data <- mutate(data, log_pval = -log10(data$P.Value))
    
    adata(data)
  })
  
  
  output$VolcanoPlot <- renderPlotly({
      data<-adata()
      
      
      highlight_genes<-input$selected_genes
      show_labels      <- input$show_labels
      logfc_thresh     <- as.numeric(input$logfc_threshold)
      pvalue_thresh    <- as.numeric(input$pvalue_threshold)
      color_by_de      <- input$color_by_de
      
      show_pvalue_thresh  <- input$show_pvalue_threshold
      show_logfc_thresh   <- input$show_logfc_threshold
      
      show_up_genes   <- input$show_up_genes
      show_down_genes <- input$show_down_genes
      
      x_label             <- "log FC"
      y_label             <- "-log10 Pval"
      legend_title        <-"Differentially Expressed"
      xlim = ranges$x
      ylim = ranges$y
      
      color1<-input$color1
      color2<-input$color2
      ptsize<-input$ptsize
      transp<-(100-input$transp)/100
      labelsize<-input$labelsize

      
      #update show_up_genes
      de_vec <- NULL
      if(show_up_genes == TRUE & show_down_genes==TRUE)
        de_vec <- abs(data$logFC) >= logfc_thresh & data$P.Value <= pvalue_thresh 
      else if(show_up_genes == TRUE)
        de_vec <- abs(data$logFC) >= logfc_thresh & data$P.Value <= pvalue_thresh & data$logFC>0
      else if(show_down_genes == TRUE)
        de_vec <- abs(data$logFC) >= logfc_thresh & data$P.Value <= pvalue_thresh & data$logFC<0
      else #FALSE FALSE
        de_vec<-rep(FALSE,nrow(data))
      
      
      
      

      
      # Reorder de_vec factors so that TRUE is first
      de_vec <- factor(de_vec, levels=c("TRUE", "FALSE"))
      
      # build base of plot
      volcano <- ggplot(data, aes(x = .data$logFC, y = log_pval, text=paste0("Gene:",.data$GeneSymbol) ))
      
      


      
      # if show_devec is true color by DE genes
      if (color_by_de) {

        volcano <- volcano + geom_point(alpha = transp, size = ptsize, color = color1, data = ~subset(., de_vec == TRUE)) + 
                           geom_point(alpha = transp, size = ptsize, color = color2, data = ~subset(., de_vec == FALSE))


      } else {
        volcano <- volcano +
          geom_point(alpha = .6)
      }
      
      if (any(!is.null(c(xlim, ylim)))) {
        volcano <- volcano + coord_cartesian(xlim = xlim, ylim = ylim, expand = FALSE)
      }
      
      # if show_pvalue_thresh = true add vline layer
      if (show_pvalue_thresh) {
        volcano <- volcano + geom_hline(yintercept = -log10(pvalue_thresh), linetype = "dashed", col = "grey", size = 1)
      }
      
      # if show_logfc_thresh = true add hline layer
      if (show_logfc_thresh) {
        volcano <- volcano + geom_vline(xintercept = c(logfc_thresh, -logfc_thresh), linetype = "dashed", col = "grey", size = 1)
      }
      
      
      
      #text
      if(show_labels==TRUE){
        if(length(highlight_genes)>0){
          volcano <- volcano +
            geom_text(
              label=ifelse(data$GeneSymbol %in% highlight_genes,data$GeneSymbol,NA), # to display genes to render
              na.rm = T,
              size = labelsize,
              nudge_y = 0.1
            )
        }else{
          volcano <- volcano +
            geom_text(
              label=ifelse(de_vec==TRUE,data$GeneSymbol,NA),   #to display all relevant genes
              na.rm = T,
              size = labelsize,
              nudge_y = 0.1
            )
        }
      }
      
      # add finishing touches to plot
      volcanoPlot <- volcano + labs(x = x_label, y = y_label, color = legend_title) + theme_classic(base_size = 12)
      
      # display plot
      tryCatch(
        expr = {
          ggplotly(volcanoPlot)
        },
        error = function(e){ 
          NULL
        }
      )
  })
  
  
  
  output$VolcanoTable <-  renderDT({
    data<-adata()
    dataRAW<-data

    filter<-input$filtervolcanotable
    show_up_genes   <- input$show_up_genes
    show_down_genes <- input$show_down_genes
    
    preurl<-"<a  target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    posturl<-"</a>"
    data$GeneSymbol<-paste0("/",data$GeneSymbol)
    data$GeneSymbol<-stringr::str_replace_all(data$GeneSymbol,"/([^/]+)",paste0(preurl,"\\1\">\\1",posturl))
    
    
    if(filter==TRUE){
      
      logfc_thresh  <- as.numeric(input$logfc_threshold)
      pvalue_thresh <- as.numeric(input$pvalue_threshold)
      
      
      #update show_up_genes
      de_vec <- NULL
      if(show_up_genes == TRUE & show_down_genes==TRUE)
        de_vec <- abs(data$logFC) >= logfc_thresh & data$P.Value <= pvalue_thresh
      else if(show_up_genes == TRUE)
        de_vec <- abs(data$logFC) >= logfc_thresh & data$P.Value <= pvalue_thresh & data$logFC>0
      else if(show_down_genes == TRUE)
        de_vec <- abs(data$logFC) >= logfc_thresh & data$P.Value <= pvalue_thresh & data$logFC<0
      else #FALSE FALSE
        de_vec<-rep(FALSE,nrow(data))
      

      data<-data[de_vec,]
      
    }
    tabledatavolcano(dataRAW)
    data
    
  }, escape = FALSE, filter = "top", selection = 'none')
  
  
  
  
  
  
  computeGSEA<- eventReactive(list(input$ontology,input$ranking,input$perm,input$minGsize,input$maxGsize, input$CCfile, input$pvalue_threshold_GSEA, input$pAdjustMethodGSEA),{
    
    output$nomatchGSEA <- renderUI({ HTML("") })
    
    ontology<-input$ontology
    ranking<-input$ranking
    perm     <- as.numeric(input$perm)
    minGsize <- as.numeric(input$minGsize)
    maxGsize <- as.numeric(input$maxGsize)
    pvalue_threshold_GSEA<-as.numeric(input$pvalue_threshold_GSEA)
    pAdjustMethodGSEA    <-input$pAdjustMethodGSEA
    
    res<-adata()
    
    
    #see tuto: https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/gsea.html
    
    ## Fetch GENE SYMBOLS to ENTREZID conversion table
    conv_tbl <- AnnotationDbi::select(org.Hs.eg.db, columns = c('ENTREZID'), keytype = 'SYMBOL', keys = res$GeneSymbol)
    
    
    # Add ENTREZID column to dataset
    res <- res %>% 
      mutate(SYMBOL = GeneSymbol) %>% 
      left_join(conv_tbl, by = "SYMBOL", multiple = "first")
    
    
    
    ################################################
    

     nbrow1<-nrow(res)
     #removes rows with no ENTREZ gene ID
     #res2<-na.omit(res$ENTREZID)
     #nbrow2<-nrow(res2)
     nbrow2<-sum(! is.na( res$ENTREZID ) )

     output$commentGSEA <- renderUI({ HTML(paste0((nbrow1-nbrow2)," (",round(100-(nbrow2/nbrow1*100),2),"%) rows are missing a gene symbol or don\'t match an ENTREZ gene ID and have been removed from the analysis.")) })
    
    ################################################
    
    
    #why use Pval (not adjusted): https://www.biostars.org/p/9608060/

    #score and rank genes
    if(ranking=="log"){
      geneList <- res %>%
        filter(!is.na(ENTREZID), !is.na(logFC)) %>% 
        mutate(ranking_metric = -log10(P.Value)*sign(logFC)) %>% 
        group_by(ENTREZID) %>% 
        summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
        arrange(-ranking_metric) %>% # sort descending (important!)
        tibble::deframe() # convert to named vector
    }else if(ranking=="T"){
      geneList <- res %>%
        filter(!is.na(ENTREZID), !is.na(t)) %>% 
        mutate(ranking_metric = t) %>% 
        group_by(ENTREZID) %>% 
        summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
        arrange(-ranking_metric) %>% # sort descending (important!)
        tibble::deframe() # convert to named vector
    }else if(ranking=="logFC"){
      geneList <- res %>%
        filter(!is.na(ENTREZID), !is.na(logFC)) %>% 
        mutate(ranking_metric = logFC) %>% 
        group_by(ENTREZID) %>% 
        summarise(ranking_metric = mean(ranking_metric, na.rm = TRUE)) %>% 
        arrange(-ranking_metric) %>% # sort descending (important!)
        tibble::deframe() # convert to named vector
      }
    
     
     ########export raw inputs and results !!!!!!!!!!!!!!
     write.table(file = paste0(path2files,"/LavaGO-GSEA-genelist.txt"),     
                 x = geneList, row.names = T, quote = F, col.names = F, sep = " ",  append = F)
     write.table(file = paste0(path2files,"/LavaGO-GSEA-data.txt"), 
                 x = res, row.names = F, append = F)
     ##########################
     
     
     
    
    set.seed(123) #to get consistent results
    
    #compute GSEA
    if(ontology %in% c("ALL","BP","MF","CC")){
      cgsea_res <- clusterProfiler::gseGO(geneList = geneList, 
                         ont = ontology,    #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three
                         OrgDb="org.Hs.eg.db",
                         minGSSize = minGsize, #10,  #minimal size of each geneSet for analyzing
                         maxGSSize = maxGsize, #500, #maximal size of genes annotated for testing
                         eps = 0,#1e-10,              #This parameter sets the boundary for calculating the p value  0=more precise Pval
                         nPermSimple = perm, #100000, #10000
                         seed = TRUE,
                         pvalueCutoff = pvalue_threshold_GSEA,
                         pAdjustMethod = pAdjustMethodGSEA,
                         keyType = 'ENTREZID'
                         
      )
    }
    else if(ontology=="KEGG"){
      cgsea_res <- clusterProfiler::gseKEGG(geneList = geneList, 
                           organism = "hsa",
                           minGSSize = minGsize,
                           maxGSSize = maxGsize,
                           eps = 0,
                           nPerm = perm,
                           seed = TRUE,
                           pvalueCutoff = pvalue_threshold_GSEA,
                           pAdjustMethod = pAdjustMethodGSEA
      )
    }
    else if(ontology=="mKEGG"){  #not working
      cgsea_res <- clusterProfiler::gseMKEGG(geneList = geneList, 
                            organism = "hsa",
                            minGSSize = minGsize, 
                            maxGSSize = maxGsize,
                            eps = 0,
                            nPerm = perm,
                            seed = TRUE,
                            pvalueCutoff = pvalue_threshold_GSEA,
                            pAdjustMethod = pAdjustMethodGSEA
      )
    }
    else if(ontology=="WP"){
      cgsea_res <- clusterProfiler::gseWP(geneList = geneList, 
                         organism = "Homo sapiens"
      )
    }
    
    
    
    
    gseadata(cgsea_res) # to save the results and show in table
    geneList(geneList)
    
    
    
    
    
    ########export raw inputs and results !!!!!!!!!!!!!!
    write.table(file = paste0(path2files,"/LavaGO-GSEA-results.txt"), 
                x = cgsea_res@result, row.names = F, append = F)
    #################
    
    
    return(cgsea_res)
    
    
  })
  
  
  output$gseaPlot <- renderPlot({
    
    
    gocat <- input$gocat
    gocat <- unlist(strsplit(gocat,"\n",fixed = T))
    
    showcat  <- as.numeric(input$showcat)
    plottype <- input$plottype
    
    
    cgsea_res<-computeGSEA()
    geneList<-geneList()         
    
    
    
    
    tryCatch(
      expr = {


############################
    
    
    if(plottype=="dot plot"){
      
      if( length(gocat)==0 ){
        dotplot(cgsea_res, showCategory=showcat, split=".sign") + facet_grid(.~.sign)
      }else{
        dotplot(cgsea_res, showCategory=gocat, split=".sign") + facet_grid(.~.sign)
      }
    }  #not working
    # else if(plottype=="bar plot"){
    #   barplot(cgsea_res, showCategory=showcat) #+ ggtitle("barplot for GSEA")
    #   
    # }
    
    else if(plottype=="Enrichment Map"){
      pwt <- pairwise_termsim(cgsea_res)
      
      
      if( length(gocat)==0 ){
        emapplot(pwt, showCategory = showcat)
      }else{
        emapplot(pwt, showCategory = gocat)
      }
    
    }else if(plottype=="Category Netplot"){
      ## convert gene ID to Symbol
      cgsea_res <- setReadable(cgsea_res, 'org.Hs.eg.db', 'ENTREZID')
      
      pwt <- pairwise_termsim(cgsea_res)
      
      
      if( length(gocat)==0 ){
        cnetplot(pwt, categorySize="pvalue", foldChange=geneList, showCategory = showcat, colorEdge = TRUE)
      }else{
        cnetplot(pwt, categorySize="pvalue", foldChange=geneList, showCategory = gocat, colorEdge = TRUE)
      }

    }else if(plottype=="upset plot"){
      upsetplot(cgsea_res)
      
    }
    

        
###############################
        
      },
      error = function(e){ 
        NULL
      }
    )
    
  })
  
  
  output$gseaPlotridge <- renderGirafe({

    showcat <- as.numeric(input$showcat)
    gocat   <- input$gocat
    gocat   <- unlist(strsplit(gocat,"\n",fixed = T))
    
    # plottype<-input$plottype
    cgsea_res <- computeGSEA()
    
    
    
    tryCatch(
      expr = {
    
############################

        if( length(gocat)==0 ){
          toto<-enrichplot::ridgeplot(cgsea_res, showCategory=showcat, label_format = 30)+theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
          x <- girafe(ggobj = toto, height_svg = showcat*0.3)
        }else{
          toto<-enrichplot::ridgeplot(cgsea_res, showCategory=gocat,   label_format = 30)+theme(axis.text.x = element_text(size = 8), axis.text.y = element_text(size = 8))
          x <- girafe(ggobj = toto, height_svg = length(gocat)*0.3  )
        }
        
        x <- girafe_options(x, opts_zoom(min = .7, max = 10), opts_hover(css = "fill:wheat;stroke:orange;r:5pt;"), opts_selection_key() )
        
##############################
    
      },
      error = function(e){ 
        output$nomatchGSEA <- renderUI({ HTML("<br><b>No result</b>") })
        NULL
      }
    )
    
        
  })
  
  
  
  
  
observe({
  output$gseaTable <-  renderDT({
    
    enriched<-gseadata()
    ## convert gene ID to Symbol
    enriched <- setReadable(enriched, 'org.Hs.eg.db', 'ENTREZID') 
    enriched <- enriched@result
    
    enrichedRAW<-enriched
    
    #make link for genes
    preurl<-"<a  target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene="
    posturl<-"</a>"
    enriched$core_enrichment<-paste0("/",enriched$core_enrichment)
    enriched$core_enrichment<-stringr::str_replace_all(enriched$core_enrichment,"/([^/]+)",paste0(preurl,"\\1\">\\1",posturl))
    
    #make link for GO terms
    preurl<-"<a  target=\"_blank\" href=\"https://amigo.geneontology.org/amigo/term/"
    enriched$ID<-stringr::str_replace_all(enriched$ID,"(GO:[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))

    
    #make link for Wikipathway
    preurl<-"<a  target=\"_blank\" href=\"https://www.wikipathways.org/pathways/"
    enriched$ID<-stringr::str_replace_all(enriched$ID,"(WP[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))
    
    
    #make link for KEGG
    preurl<-"<a  target=\"_blank\" href=\"https://www.genome.jp/dbget-bin/www_bget?pathway+"
    enriched$ID<-stringr::str_replace_all(enriched$ID,"(hsa[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))
    
    
    
    
    
    
    
    
    enriched<-enriched 
    
    tabledatagsea(enrichedRAW) 
    
    if(nrow(enriched)>0)
      enriched
    
  }, escape = FALSE, filter = "top", selection = list(mode="single", target = "row"))
})
  
  
  
  observeEvent(input$gseaTable_rows_selected, {

    output$gseaPlot2 <- renderPlot({

      req(input$gseaTable_rows_selected)
      enriched<-gseadata()
      enriched<-enriched@result
      
      goid<-enriched[input$gseaTable_rows_selected,1]
      goname<-enriched[input$gseaTable_rows_selected,2]
      
      cgsea_res<-computeGSEA()
      gseaplot(cgsea_res, by = "all", title = goname, geneSetID = goid)
      
    })
  
  })
  
  
  
################################################################################
#############################ORA################################################
################################################################################
  
  # enriched       <- reactiveVal(NULL)
  ORAdata        <- reactiveVal(NULL)
  
  
  computeORA<- eventReactive(list(input$ontology2, input$CCfile, input$gseaTable, input$refresh, input$pAdjustMethodORA, input$pvalue_threshold_ORA, input$minGSizeORA, input$maxGSizeORA),{  #, tabledatagsea
    ontology<-input$ontology2
    
    
    pAdjustMethodORA    <- input$pAdjustMethodORA
    pvalue_threshold_ORA<- as.numeric(input$pvalue_threshold_ORA)
    minGSizeORA         <- as.numeric(input$minGSizeORA)
    maxGSizeORA         <- as.numeric(input$maxGSizeORA)
    
    output$nomatchORA <- renderUI({ HTML("") })
    
    data<-adata()
    filter<-input$filtervolcanotable
    if(filter==TRUE){
      logfc_thresh  <- as.numeric(input$logfc_threshold)
      pvalue_thresh <- as.numeric(input$pvalue_threshold)
      show_up_genes   <- input$show_up_genes
      show_down_genes <- input$show_down_genes
      
      
      de_vec <- NULL
      if(show_up_genes == TRUE & show_down_genes==TRUE)
        de_vec <- abs(data$logFC) >= logfc_thresh & data$P.Value <= pvalue_thresh
      else if(show_up_genes == TRUE)
        de_vec <- abs(data$logFC) >= logfc_thresh & data$P.Value <= pvalue_thresh & data$logFC>0
      else if(show_down_genes == TRUE)
        de_vec <- abs(data$logFC) >= logfc_thresh & data$P.Value <= pvalue_thresh & data$logFC<0
      else #FALSE FALSE
        de_vec<-rep(FALSE,nrow(data))
      
      data<-data[de_vec,]
    }
    selection<-data$GeneSymbol
    
    print ( paste0( length(selection), " genes to analyze" ) )
    
    ## Fetch GENE SYMBOLS to ENTREZID conversion table
    conv_tbl <- AnnotationDbi::select(org.Hs.eg.db, columns = c('ENTREZID'), keytype = 'SYMBOL', keys = selection)
    
    
    
    
    
    ######update
    res <- data %>% 
      mutate(SYMBOL = GeneSymbol) %>% 
      left_join(conv_tbl, by = "SYMBOL", multiple = "first")
    
    nbrow1<-nrow(res)
    nbrow2<-sum(! is.na( res$ENTREZID ) )

    selectionENTREZ<-unique(na.omit(res$ENTREZID))
    #######
    
    
    
    
    
    ########export raw inputs and results !!!!!!!!!!!!!!
    write.table(file = paste0(path2files,"/LavaGO-ORA-genelist.txt"), 
                x = selectionENTREZ, row.names = F, quote = F, col.names = F, append = F)
    write.table(file = paste0(path2files,"/LavaGO-ORA-data.txt"), 
                x = res, row.names = F, append = F)
    ##########################
    
    
    
    
    

    
    
    output$commentORA <- renderUI({ HTML(paste0("ORA is performed on the gene list selection displayed after the volcano plot(",length(selection)," genes). Set up the fold change and the P Value threshold accordingly and make sure to tick \"show significant entries only\".
                                    <br>", (nbrow1-nbrow2)," (",round(100-(nbrow2/nbrow1*100),2),"%) rows are missing a gene symbol or don\'t match an ENTREZ gene ID and have been removed from the analysis.
                                    <br><b>Gene selection parameters:</b> logfc>=",logfc_thresh," P.Value<=",pvalue_thresh," up genes=",show_up_genes," down genes=",show_down_genes )) })
    

    
    #compute ORA
    if(ontology %in% c("ALL","BP","MF","CC")){

      
      yy <- clusterProfiler::enrichGO(gene = selectionENTREZ, 
                                      OrgDb = 'org.Hs.eg.db', 
                                      ont=ontology,                       #one of "BP", "MF", and "CC" subontologies, or "ALL" for all three
                                      pvalueCutoff = pvalue_threshold_ORA,
                                      pAdjustMethod = pAdjustMethodORA,
                                      minGSSize = minGSizeORA, 
                                      maxGSSize = maxGSizeORA,
                                      keyType = 'ENTREZID'
                                      )
      
    }
    else if(ontology=="KEGG"){
      yy <- enrichKEGG(gene = selectionENTREZ,
                       pvalueCutoff = pvalue_threshold_ORA,
                       pAdjustMethod = pAdjustMethodORA,
                       minGSSize = minGSizeORA,
                       maxGSSize = maxGSizeORA
            )
    }
    
    else if(ontology=="WP"){
      yy <- enrichWP(gene = selectionENTREZ,
                     organism = "Homo sapiens"
      )
    }
    
    
    
    
    else if(ontology=="PC"){
      yy <- enrichPC( unique(na.omit(selection)) )
    }
    
    
    ########export raw inputs and results !!!!!!!!!!!!!!
    write.table(file = paste0(path2files,"/LavaGO-ORA-results.txt"), 
                x = yy@result, append = F)
    ###############
    
    
    # ORAdata(yy) # to save the results and show in table
    ORAdata(yy)
    return(yy)
    
  })
  
  
  

  
  
  output$ORAPlot <-  renderPlot({

    ORAlines<-as.numeric(input$showcat2)

    gocat   <-input$gocat2
    gocat   <-unlist(strsplit(gocat,"\n",fixed = T))
    showcat <-as.numeric(input$showcat2)
    plottype<-input$plottype2

    ORAres<-computeORA()
    
    
    
    print(ORAres)
    
    
    
    if(plottype=="dot plot"){

      if( length(gocat)==0 ){
        toto<-dotplot(ORAres, showCategory=showcat)
      }else{
        toto<-dotplot(ORAres, showCategory=gocat)
      }
    }
    else if(plottype=="Enrichment Map"){
      pwt <- pairwise_termsim(ORAres)
      
      
      if( length(gocat)==0 ){
        toto<-emapplot(pwt, showCategory = showcat)
      }else{
        toto<-emapplot(pwt, showCategory = gocat)
      }
      
    }else if(plottype=="Category Netplot"){
      ## convert gene ID to Symbol
      ORAres <- setReadable(ORAres, 'org.Hs.eg.db', 'ENTREZID')
      
      pwt <- pairwise_termsim(ORAres)
      
      
      if( length(gocat)==0 ){
        toto<-cnetplot(pwt, categorySize="pvalue", showCategory = showcat, colorEdge = TRUE)
      }else{
        toto<-cnetplot(pwt, categorySize="pvalue", showCategory = gocat, colorEdge = TRUE)  
      }
      
    }else if(plottype=="upset plot"){
      toto<-upsetplot(ORAres)
      
    }
    
    
    tryCatch(
      expr = {
        plot(toto)
      },
      error = function(e){ 
        output$nomatchORA <- renderUI({ HTML("<br><b>No result</b>") })
        NULL
      }
    )
    
  })
  
  
  observe({
    output$ORAtable <-  renderDT({
      
      enriched<-ORAdata() #gseadata()
      
      if(input$ontology2 != "PC"){
        ## convert gene ID to Symbol
        enriched <- setReadable(enriched, 'org.Hs.eg.db', 'ENTREZID') 
      }
      enriched <- enriched@result
      enrichedRAW<-enriched
      
      #make link for genes
      preurl<-"<a  target=\"_blank\" class=\"btn btn-primary\" href=\"https://www.genecards.org/cgi-bin/carddisp.pl?gene="
      posturl<-"</a>"
      enriched$geneID<-paste0("/",enriched$geneID) #core_enrichment
      enriched$geneID<-stringr::str_replace_all(enriched$geneID,"/([^/]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      #make link for GO terms
      preurl<-"<a  target=\"_blank\" href=\"https://amigo.geneontology.org/amigo/term/"
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(GO:[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      #make link for Wikipathway
      preurl<-"<a  target=\"_blank\" href=\"https://www.wikipathways.org/pathways/"
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(WP[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      #make link for KEGG
      preurl<-"<a  target=\"_blank\" href=\"https://www.genome.jp/dbget-bin/www_bget?pathway+"
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(hsa[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      
      
      
      #make link for REACTOME
      preurl<-"<a  target=\"_blank\" href=\"https://reactome.org/content/detail/"
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(R-HSA-[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      #make link for PathBank
      preurl<-"<a  target=\"_blank\" href=\"https://pubchem.ncbi.nlm.nih.gov/pathway/"
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(pathbank:SMP[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))

      
      #make link for HumanCyc
      preurl<-"<a  target=\"_blank\" href=\"https://apps.pathwaycommons.org/pathways?uri="
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(humancyc:Pathway[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      
      #make link for HumanCyc
      preurl<-"<a  target=\"_blank\" href=\"https://apps.pathwaycommons.org/pathways?uri="
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(pid:pid_[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      #make link for HumanCyc
      preurl<-"<a  target=\"_blank\" href=\"https://apps.pathwaycommons.org/pathways?uri="
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(inoh:id[^ ]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      #make link for HumanCyc
      preurl<-"<a  target=\"_blank\" href=\"https://apps.pathwaycommons.org/pathways?uri="
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(netpath:Pathway_[^ ]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      #make link for HumanCyc
      preurl<-"<a  target=\"_blank\" href=\"https://apps.pathwaycommons.org/pathways?uri="
      enriched$ID<-stringr::str_replace_all(enriched$ID,"(biofactoid:[^ ]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      #make link for Panther
      preurl<-"<a  target=\"_blank\" href=\"https://pantherdb.org/pathway/pathwayDiagram.jsp?catAccession="
      enriched$ID<-stringr::str_replace_all(enriched$ID,"panther.pathway:(P[0-9]+)",paste0(preurl,"\\1\">\\1",posturl))
      
      
      #save for download
      tabledataora(enrichedRAW)
      
      
      

      
      enriched$zScore        <-signif(enriched$zScore,3)
      enriched$FoldEnrichment<-signif(enriched$FoldEnrichment,3)
      enriched$pvalue        <-signif(enriched$pvalue,3)
      enriched$qvalue        <-signif(enriched$qvalue,3)
      enriched$p.adjust      <-signif(enriched$p.adjust,3)
      enriched$RichFactor    <-signif(enriched$RichFactor,3)
      enriched
      
      
      
    }, escape = FALSE, filter = "top", selection = list(mode="single", target = "row"))
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)

