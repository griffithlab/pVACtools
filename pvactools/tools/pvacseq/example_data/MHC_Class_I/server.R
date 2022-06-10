# load libraries
library(shiny)
library(ggplot2)
library(DT)
library(reshape2)
library(jsonlite)
library(tibble)
library(tidyr)
library(plyr)
library(dplyr)

source("anchor_and_helper_functions.R")
source("styling.R")

#specify max shiny app upload size (currently 300MB)
options(shiny.maxRequestSize=300*1024^2)

server <- shinyServer(function(input, output, session) {
  
  ## helper function defined for generating shinyInputs in mainTable (Evaluation dropdown menus)
  shinyInput = function(data, FUN, len, id, ...) { 
    inputs = character(len) 
    for (i in seq_len(len)) { 
      inputs[i] = as.character(FUN(paste0(id, i), label = NULL, ...,selected=data[i,'Evaluation']))
    } 
    inputs
  } 
  
  ## helper function defined for generating shinyInputs in mainTable (Investigate button)
  shinyInputSelect = function(FUN, len, id, ...) {
    inputs = character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  ## helper function defined for getting values of shinyInputs in mainTable (Evaluation dropdown menus)
  shinyValue = function(id, len, data) { 
    unlist(lapply(seq_len(len), function(i) { 
      value = input[[paste0(id, i)]] 
      if (is.null(value)) {
        data[i,'Evaluation']
      }
      else { 
        value
      }
    })) 
  } 

  ##############################DATA UPLOAD TAB################################### 
  
  #reactive values defined for row selection, main table, metrics data, additional data, and dna cutoff 
  df <- reactiveValues(
    selectedRow = 1,
    mainTable = NULL,
    dna_cutoff = 0.5,
    metricsData = NULL,
    additionalData = NULL,
    gene_list = NULL,
    allele_expr_high = 3,
    allele_expr_low = 1,
    comments = data.frame("N/A")
  )
 
  #Option 1: User uploaded main aggregate report file
  observeEvent(input$mainDataInput$datapath,{
    mainData <- read.table(input$mainDataInput$datapath, sep = '\t',  header = FALSE, stringsAsFactors = FALSE, check.names=FALSE)
    colnames(mainData) <- mainData[1,]
    mainData <- mainData[-1,]
    row.names(mainData) <- NULL
    mainData$`Eval` <- shinyInput(mainData, selectInput,nrow(mainData),"selecter_", choices=c("Pending", "Accept", "Reject", "Review"), width="60px")
    mainData$Select <- shinyInputSelect(actionButton, nrow(mainData), "button_" , label = "Investigate", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
    mainData$`IC50 MT` <- as.numeric(mainData$`IC50 MT`)
    mainData$`%ile MT` <- as.numeric(mainData$`%ile MT`)
    mainData$`RNA Depth` <- as.integer(mainData$`RNA Depth`)
    df$mainTable <- mainData
    dna_vaf <- as.numeric(as.character(unlist(df$mainTable['DNA VAF'])))
    df$dna_cutoff <- max(dna_vaf[dna_vaf < 0.6])
    df$mainTable$`Tier Count` <- apply(df$mainTable, 1, function(x) tier_numbers(x, input$anchor_contribution, df$dna_cutoff, df$allele_expr_high, df$allele_expr_low, unlist(x["Pos"]), anchor_mode="default"))
    df$mainTable$`Gene of Interest` <- apply(df$mainTable,1, function(x) {any(x['Gene'] == df$gene_list)})
    if("Comments" %in% colnames(df$mainTable))
    {
      df$comments <- data.frame(data = df$mainTable$`Comments`,nrow=nrow(df$mainTable),ncol=1)
    }
    else{
      df$comments <- data.frame(matrix("No comments",nrow=nrow(df$mainTable)),ncol=1)
    }
    rownames(df$comments) <- df$mainTable$ID
  })
  
  #Option 1: User uploaded metrics file
  observeEvent(input$metricsDataInput,{
    df$metricsData <- fromJSON(input$metricsDataInput$datapath)
  })
  
  #Option 1: User uploaded additional data file 
  observeEvent(input$additionalDataInput,{
    addData <- read.table(input$additionalDataInput$datapath, sep = '\t',  header = FALSE, stringsAsFactors = FALSE, check.names=FALSE)
    colnames(addData) <- addData[1,]
    addData <- addData[-1,]
    row.names(addData) <- NULL
    df$additionalData <- addData
  })
  
  #Option 1: User uploaded additional gene list 
  observeEvent(input$gene_list,{
    gene_list <- read.table(input$gene_list$datapath, sep = '\t',  header = FALSE, stringsAsFactors = FALSE, check.names=FALSE)
    df$gene_list <- gene_list
    df$mainTable$`Gene of Interest` <- apply(df$mainTable,1, function(x) {any(x['Gene'] == df$gene_list)})
  })
  

  #Option 2: Load from default (relative) file path for aggregate report file 
   observeEvent(input$loadDefaultmain,{
     data <- getURL("https://raw.githubusercontent.com/griffithlab/pVACtools/835a7e4ae8b660a362c0c0b54140e26830d72bf2/tools/pvacview/data/H_NJ-HCC1395-HCC1395.all_epitopes.aggregated.tsv")
     mainData <- read.table(text=data, sep = '\t', header = FALSE, stringsAsFactors = FALSE, check.names=FALSE)
     colnames(mainData) <- mainData[1,]
     mainData <- mainData[-1,]
     row.names(mainData) <- NULL
     mainData$`Eval` <- shinyInput(mainData, selectInput,nrow(mainData),"selecter_", choices=c("Pending", "Accept", "Reject", "Review"), width="60px")
     mainData$Select <- shinyInputSelect(actionButton, nrow(mainData), "button_" , label = "Investigate", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
     mainData$`IC50 MT` <- as.numeric(mainData$`IC50 MT`)
     mainData$`%ile MT` <- as.numeric(mainData$`%ile MT`)
     mainData$`RNA Depth` <- as.integer(mainData$`RNA Depth`)
     df$mainTable <- mainData
     dna_vaf <- as.numeric(as.character(unlist(df$mainTable['DNA VAF'])))
     df$dna_cutoff <- max(dna_vaf[dna_vaf < 0.6])
     df$mainTable$`Tier Count` <- apply(df$mainTable, 1, function(x) tier_numbers(x, input$anchor_contribution, df$dna_cutoff, df$allele_expr_high, df$allele_expr_low, unlist(x["Pos"]), anchor_mode="default"))
     df$mainTable$`Gene of Interest` <- apply(df$mainTable,1, function(x) {any(x['Gene'] == df$gene_list)})
     if("Comments" %in% colnames(df$mainTable))
     {
       df$comments <- data.frame(data = df$mainTable$`Comments`,nrow=nrow(df$mainTable),ncol=1)
     }
     else{
       df$comments <- data.frame(matrix("No comments",nrow=nrow(df$mainTable)),ncol=1)
     }
     rownames(df$comments) <- df$mainTable$ID
     metricsdata <- getURL("https://raw.githubusercontent.com/griffithlab/pVACtools/6423b8b65f2e3f5cc2979f33b86c4650a6aa4570/tools/pvacview/data/H_NJ-HCC1395-HCC1395.all_epitopes.aggregated.metrics.json")
     df$metricsData <- fromJSON(txt = metricsdata)
     updateTabItems(session, "tabs", "explore")
   })

  observeEvent(input$visualize,{
    updateTabItems(session, "tabs", "explore")
  })
  
  #reactions for once "regenerate table" command is submitted
  observeEvent(input$submit,{
      session$sendCustomMessage('unbind-DT', 'mainTable')
      df$dna_cutoff <- input$dna_cutoff
      df$mainTable$`Evaluation` <- shinyValue("selecter_",nrow(df$mainTable), df$mainTable)
      df$mainTable$`Mutated Positions` <- apply(df$mainTable, 1, function(x) calculate_mutation_info(df$metricsData[[x[["ID"]]]]))
      df$mainTable$`Best HLA allele` <- apply(df$mainTable, 1, function(x) df$metricsData[[x[["ID"]]]]$best_hla_allele)
      if (input$use_anchor){
        df$mainTable$`Tier` <- apply(df$mainTable, 1, function(x) tier(x, input$anchor_contribution, input$dna_cutoff, input$allele_expr_high, input$allele_expr_low, unlist(x["Mutated Positions"]), x["Best HLA allele"]))
        df$mainTable$`Tier Count` <- apply(df$mainTable, 1, function(x) tier_numbers(x, input$anchor_contribution, input$dna_cutoff, input$allele_expr_high, input$allele_expr_low, unlist(x["Pos"]), hla_allele = x["Best HLA allele"]))
      }else{
        df$mainTable$`Tier` <- apply(df$mainTable, 1, function(x) tier(x, input$anchor_contribution, input$dna_cutoff, input$allele_expr_high, input$allele_expr_low, unlist(x["Mutated Positions"]), x["Best HLA allele"], anchor_mode="default"))
        df$mainTable$`Tier Count` <- apply(df$mainTable, 1, function(x) tier_numbers(x, input$anchor_contribution, input$dna_cutoff, input$allele_expr_high, input$allele_expr_low, unlist(x["Pos"]), hla_allele = x["Best HLA allele"], anchor_mode="default"))
      }
      df$mainTable$`Gene of Interest` <- apply(df$mainTable,1, function(x) {any(x['Gene'] == df$gene_list)})
      tier_sorter <- c("Pass", "Relaxed", "LowExpr", "Anchor", "Subclonal", "Poor", "NoExpr")
      df$mainTable$`Rank_ic50` <- NA
      df$mainTable$`Rank_expr` <- NA
      df$mainTable$`Rank_ic50` <- rank(as.numeric(df$mainTable$`IC50 MT`), ties.method = "first")
      df$mainTable$`Rank_expr` <- rank(desc(as.numeric(df$mainTable$`Allele Expr`)), ties.method = "first")
      df$mainTable$`Rank` <- df$mainTable$`Rank_ic50` + df$mainTable$`Rank_expr`
      df$mainTable <- df$mainTable %>%
        arrange(factor(Tier, levels = tier_sorter), Rank)
      df$mainTable$`Rank` <- NULL
      df$mainTable$`Rank_ic50` <- NULL
      df$mainTable$`Rank_expr` <- NULL
      df$mainTable$Select <- shinyInputSelect(actionButton, nrow(df$mainTable), "button_" , label = "Investigate", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
      df$mainTable$`Eval` <- shinyInput(df$mainTable, selectInput,nrow(df$mainTable),"selecter_", choices=c("Pending", "Accept", "Reject", "Review"), width="60px")
  })
  
  #determine hla allele count in order to generate column tooltip locations correctly 
  hla_count <- reactive({
    which(colnames(df$mainTable)=="Gene" )-1
  })
  
  #class type of user-provided additional file 
  type <- reactive({
    switch(input$hla_class,
           class_i = 1,
           class_ii = 2)
  })
  
  output$type_text <- renderText({
    input$add_file_label
  })
  
  output$max_dna <- renderText({
    if (is.null(df$mainTable)){
      return ("N/A")
    }
    dna_vaf <- as.numeric(as.character(unlist(df$mainTable['DNA VAF'])))
    dna_cutoff <- max(dna_vaf[dna_vaf < 0.6])
    dna_cutoff 
  })
  
  output$comment_text <- renderText({
    if (is.null(df$mainTable)){
      return ("N/A")
    }
    df$comments[selectedID(),1]
  })
  
  ##############################PEPTIDE EXPLORATION TAB###########################################                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
  
  ##main table display with color/background/font/border configurations 
  output$mainTable = DT::renderDataTable(
    if (is.null(df$mainTable)){
      return ()
    }
    else{
      datatable(df$mainTable[, !(colnames(df$mainTable) == "ID") & !(colnames(df$mainTable) == "Evaluation") & !(colnames(df$mainTable) == "Mutated Positions") & !(colnames(df$mainTable) == "Best HLA allele") & !(colnames(df$mainTable) == "Comments")]
    , escape = FALSE, callback = JS(callBack(hla_count())), class = 'stripe',
          options=list(lengthChange = FALSE, dom = 'Bfrtip', 
                   columnDefs = list(list(className = 'dt-center', targets =c(0:hla_count()-1)), list(visible=FALSE, targets=c(25-(7-hla_count()),26-(7-hla_count()))),
                                     list(orderable=TRUE, targets=0)),
                   buttons = list(I('colvis')), 
                   initComplete = htmlwidgets::JS(
                     "function(settings, json) {",
                     paste("$(this.api().table().header()).css({'font-size': '", '10pt', "'});"), 
                     "}"),
                   rowCallback = JS(rowCallback(hla_count(), df$selectedRow-1)),
                   preDrawCallback = JS('function() { 
                                        Shiny.unbindAll(this.api().table().node()); }'), 
                   drawCallback = JS('function() { 
                                     Shiny.bindAll(this.api().table().node()); } ')
                   ),
      selection = 'none', 
      extensions = c("Buttons"))}
    %>% formatStyle('IC50 MT', backgroundColor = styleInterval(c(50,100,200,300,400,500,600, 700, 800, 900, 1000), 
                                                                c("#00FF00", "#00EE00","#00D500","#00BC00","#00A300", "#008B00", "#FFFF00", "#FFEB00", "#FFD800","#FFC500","#FFB100", "#FF9999"))
                     ,fontWeight = styleInterval(c(1000), c('normal', 'bold')), border= styleInterval(c(1000), c('normal','2px solid red')))
    %>% formatStyle('%ile MT', backgroundColor = styleInterval(c(0.1,0.2,0.3,0.4,0.5,0.75,1,1.5,2),
                                                                c("#00EE00","#00D500","#00BC00","#00A300", "#008B00", "#FFFF00", "#FFEB00", "#FFD800","#FFC500", "#FF9999")))
    %>% formatStyle( 'Tier', color = styleEqual(c('Pass', 'Relaxed', 'Poor','Anchor','Subclonal','LowExpr', 'NoExpr'), c('green','lightgreen', 'orange', '#b0b002', '#D4AC0D', 'salmon', 'red')) )
    %>% formatStyle(c('RNA VAF'),background = styleColorBar(range(0,1), 'lightblue'), backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat', backgroundPosition = 'right')
    %>% formatStyle(c('DNA VAF'),background = styleColorBar(range(0,1), 'lightblue'), backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat', backgroundPosition = 'right') 
    %>% formatStyle(c('RNA Expr'),background = styleColorBar(range(0,50), 'lightblue'), backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat', backgroundPosition = 'right')
    %>% formatStyle(c('RNA Depth'),background = styleColorBar(range(0,200), 'lightblue'), backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat', backgroundPosition = 'right')
    %>% formatStyle(c('Allele Expr'),background = styleColorBar(range(0,(max(as.numeric(as.character(unlist(df$mainTable['RNA VAF'])))*50))), 'lightblue'), backgroundSize = '98% 88%', backgroundRepeat = 'no-repeat', backgroundPosition = 'right')
    %>% formatStyle(c('Allele Expr'),'Tier Count', fontWeight = styleEqual(c('2'), c('bold')), border= styleEqual(c('2'), c('2px solid red')))
    %>% formatStyle(c('IC50 MT','Allele Expr'),'Tier Count', fontWeight = styleEqual(c('3'), c('bold')), border= styleEqual(c('3'), c('2px solid red')))
    %>% formatStyle(c('IC50 MT'),'Tier Count', fontWeight = styleEqual(c('4'), c('bold')), border= styleEqual(c('4'), c('2px solid red')))
    %>% formatStyle(c('IC50 WT','Pos'),'Tier Count', fontWeight = styleEqual(c('5'), c('bold')), border= styleEqual(c('5'), c('2px solid red')))
    %>% formatStyle(c('DNA VAF'),'Tier Count', fontWeight = styleEqual(c('6'), c('bold')), border= styleEqual(c('6'), c('2px solid red')))
    %>% formatStyle(c('Allele Expr'),'Tier Count', fontWeight = styleEqual(c('7'), c('bold')), border= styleEqual(c('7'), c('2px solid red')))
    %>% formatStyle(c('Gene Expression'),'Tier Count', fontWeight = styleEqual(c('8'), c('bold')), border= styleEqual(c('8'), c('2px solid red')))
    %>% formatStyle(c('RNA VAF','RNA Depth'),'Tier Count', fontWeight = styleEqual(c('8'), c('bold')), border= styleEqual(c('8'), c('2px solid green')))
    %>% formatStyle(c('RNA Expr'),'Tier Count', fontWeight = styleEqual(c('9'), c('bold')), border= styleEqual(c('9'), c('2px solid red')))
    %>% formatStyle(c('RNA VAF'),'Tier Count', fontWeight = styleEqual(c('10'), c('bold')), border= styleEqual(c('10'), c('2px solid red')))
    %>% formatStyle(c('RNA VAF','RNA Expr'),'Tier Count', fontWeight = styleEqual(c('11'), c('bold')), border= styleEqual(c('11'), c('2px solid red')))
    %>% formatStyle(c('IC50 WT','Pos'),'Tier Count', fontWeight = styleEqual(c('13'), c('bold')), border= styleEqual(c('13'), c('2px solid red')))
    %>% formatStyle(c('DNA VAF'),'Tier Count', fontWeight = styleEqual(c('14'), c('bold')), border= styleEqual(c('14'), c('2px solid red')))
    %>% formatStyle(c('DNA VAF','IC50 WT','Pos'),'Tier Count', fontWeight = styleEqual(c('15'), c('bold')), border= styleEqual(c('15'), c('2px solid red')))
    %>% formatStyle(c('IC50 WT','Pos', 'DNA VAF', 'Allele Expr'),'Tier Count', fontWeight = styleEqual(c('23'), c('bold')), border= styleEqual(c('23'), c('2px solid red')))
    %>% formatStyle(c('DNA VAF', 'Allele Expr'),'Tier Count', fontWeight = styleEqual(c('22'), c('bold')), border= styleEqual(c('22'), c('2px solid red')))
    %>% formatStyle(c('IC50 WT','Pos', 'Allele Expr'),'Tier Count', fontWeight = styleEqual(c('21'), c('bold')), border= styleEqual(c('21'), c('2px solid red')))
    %>% formatStyle(c('Allele Expr'),'Tier Count', fontWeight = styleEqual(c('20'), c('bold')), border= styleEqual(c('20'), c('2px solid red')))
    %>% formatStyle(c('Gene'),'Gene of Interest', fontWeight = styleEqual(c(TRUE), c('bold')), border= styleEqual(c(TRUE), c('2px solid green')))
    , server=FALSE)
  
  #help menu for main table
  observeEvent(input$help, {
    showModal(modalDialog(
      title = "Aggregate Report of Best Candidates by Mutation",
      h5("* Hover over individual column names to see further description of specific columns. (HLA allele columns excluded)"),
      h4(" HLA specific columns:", style="font-weight: bold"),
      h5(" Number of good binding peptides for each specific HLA-allele.", br(),
         " The same peptide could be counted in multiple columns if it was predicted to be a good binder for multiple HLA alleles."),
      h4(" Color scale for IC50 MT column:", style="font-weight: bold"),
      h5(" lightgreen to darkgreen (0nM to 500nM); ", br(),"yellow to orange (500nM to 1000nM);", br()," red (> 1000nM) "),
      h4(" Color scale for %ile MT column:", style="font-weight: bold"),
      h5(" lightgreen to darkgreen (0-0.5%);", br()," yellow to orange (0.5% to 2 %);", br()," red (> 2%) "),
      h4(" Bar backgrounds:", style="font-weight: bold"), 
      h5(" RNA VAF and DNA VAF: Bar graphs range from 0 to 1", br(),
         " RNA Depth: Bar graph ranging from 0 to maximum value of RNA depth values across variants", br(),
         " RNA Expr: Bar graph ranging from 0 to 50 (this is meant to highlight variants with lower expression values for closer inspection)", br(),
         " Allele Expr: Bar graph ranging from 0 to (50 * maximum value of RNA VAF values across variants) "),
      h4(" Tier Types:", style="font-weight: bold"),
      h5(" Variants are ordered by their Tiers in the following way: Pass, Relaxed, LowExpr, Anchor, Subclonal, Poor, NoExpr.
           Within the same tier, variants are ordered by the sum of their ranking in binding affinity and allele expression (i.e. lower binding 
           affinity and higher allele expression is prioritized.)"),
      h5(" NoExpr: Mutant allele is not expressed ", br(),
         " LowExpr: Mutant allele has low expression (Allele Expr < 1)", br(),
         " Subclonal: Likely not in the founding clone of the tumor (DNA VAF > max(DNA VAF)/2)", br(),
         " Anchor: Mutation is at an anchor residue in the shown peptide, and the WT allele has good binding (WT IC50 <1000)", br(),
         " Poor: Fails two or more of the above criteria", br(),
         " Relaxed: Passes the above criteria (1 < Allele Expr < 3), has decent MT binding (IC50 < 1000)", br(),
         " Pass: Passes the above criteria, has strong MT binding (IC50 < 500) and strong expression (Allele Expr > 3)"
      ),
    ))
  })
  
  ##update table upon selecting to investigate each individual row 
  observeEvent(input$select_button, {
    if (is.null(df$mainTable)){
      return ()
    }
    df$selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    session$sendCustomMessage('unbind-DT', 'mainTable')
    df$mainTable$`Evaluation` <- shinyValue("selecter_",nrow(df$mainTable), df$mainTable)
    df$mainTable$`Eval` <- shinyInput(df$mainTable, selectInput,nrow(df$mainTable),"selecter_", choices=c("Pending", "Accept", "Reject", "Review"), width="60px")
    dataTableProxy("mainTable") %>% 
      selectPage((df$selectedRow-1) %/% 10 + 1)
  })
  
  ##selected row text box
  output$selected <- renderText({
    if (is.null(df$mainTable)){
      return ()
    }
    df$selectedRow
  })
  
  ##selected id update 
  selectedID <- reactive({
    if (is.null(df$selectedRow)) {
      df$mainTable$ID[1]
    }
    else{
      df$mainTable$ID[df$selectedRow]
    }
  })
  
  ## Update comments section based on selected row
  observeEvent(input$comment, {
    if (is.null(df$mainTable)){
      return ()
    }
    df$comments[selectedID(), 1] <- input$comments
  })
  
  ##display of genomic information 
  output$metricsTextGenomicCoord = renderText({
    if (is.null(df$metricsData)){
      return ()
    }
    selectedID()
  })
  
  ##display of openCRAVAT link for variant
  output$url <- renderUI({
    if (is.null(df$mainTable)){
      return ()
    }
    id <- strsplit(selectedID(), "-")
    chromosome <- id[[1]][1]
    start <- id[[1]][2]
    stop <- id[[1]][3]
    ref <- id[[1]][4]
    alt <- id[[1]][5]
    url <- a("OpenCRAVAT variant report", href=paste("https://run.opencravat.org/webapps/variantreport/index.html?chrom=",chromosome,"&pos=",stop,"&ref_base=",ref,"&alt_base=",alt, sep=""), target="_blank")
    HTML(paste(url))
  })
  
  ##display of RNA VAF
  output$metricsTextRNA = renderText({
    if (is.null(df$metricsData)){
      return ()
    }
    df$metricsData[[selectedID()]]$`RNA VAF`
  })
  
  ##display of DNA VAF
  output$metricsTextDNA = renderText({
    if (is.null(df$metricsData)){
      return ()
    }
    df$metricsData[[selectedID()]]$`DNA VAF`
  })
  
  ##display of MT IC50 from additional data file
  output$addData_IC50 = renderText({
    if (is.null(df$additionalData)){
      return ()
    }
    df$additionalData[df$additionalData$ID == selectedID(),]$`IC50 MT`
  })
  
  ##display of MT percentile from additional data file 
  output$addData_percentile = renderText({
    df$additionalData[df$additionalData$ID == selectedID(),]$`%ile MT`
  })
  
  ##transcripts table displaying transcript id and transcript expression values 
  output$transcriptsTable = renderDT(
    {
      withProgress(message = 'Loading Transcripts Table', value = 0, {
      GB_transcripts <- data.frame()
      if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0){
        GB_transcripts <- data.frame("Transcripts" = df$metricsData[[selectedID()]]$good_binders_transcript,"Expression" = df$metricsData[[selectedID()]]$transcript_expr)
      }
      else{
        GB_transcripts <- data.frame("Transcripts" = character(), "Expression" = character())
      }
      incProgress(0.5)
      names(GB_transcripts) <- c("Transcripts producing good binding peptides", "Transcript Expression")
      incProgress(0.5)
      GB_transcripts 
      })
    },
    selection = list(mode='single', selected = '1')
  )
  
  ##update selected transcript id
  selectedTranscript <- reactive({
    selection <- input$transcriptsTable_rows_selected
    if (is.null(selection)){
      selection <- 1
    }
    df$metricsData[[selectedID()]]$good_binders_transcripts[selection]
  })
  
  ##display transcript expression 
  output$metricsTextTranscript = renderText({
    if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0){
      df$metricsData[[selectedID()]]$good_binders[[selectedTranscript()]]$`transcript_expr`
    }
    else {
      "N/A"
    }
  })
  
  ##display gene expression 
  output$metricsTextGene = renderText({
    if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0){
      df$metricsData[[selectedID()]]$`gene_expr`
    }
    else {
      "N/A"
    }
  })
  
  ##display peptide table with coloring 
  output$peptideTable = renderDT({
      withProgress(message = 'Loading Peptide Table', value = 0, {
        if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0 & !is.null(df$metricsData)){
          peptide_data <- df$metricsData[[selectedID()]]$good_binders[[selectedTranscript()]]
          peptide_names <- names(peptide_data)
          for(i in 1:length(peptide_names)){
            peptide_data[[peptide_names[[i]]]]$individual_ic50_calls <- NULL
            peptide_data[[peptide_names[[i]]]]$individual_percentile_calls <- NULL
          }
          incProgress(0.5)
          peptide_data <- as.data.frame(peptide_data)
          incProgress(0.5)
          dtable <- datatable(do.call("rbind",lapply(peptide_names, table_formatting, peptide_data)), options =list(
            pageLength = 10,
            columnDefs = list(list(defaultContent="X",
              targets = c(2:hla_count()+1)),
              list(orderable=TRUE, targets=0)),
            rowCallback = JS('function(row, data, index, rowId) {',
                             'console.log(rowId)','if(((rowId+1) % 4) == 3 || ((rowId+1) % 4) == 0) {',
                             'row.style.backgroundColor = "#E0E0E0";','}','}')
          ), selection = list(mode='single', selected = '1')) %>% formatStyle('Type', fontWeight = styleEqual('MT','bold'), color = styleEqual('MT', '#E74C3C'))
          dtable$x$data[[1]] <- as.numeric(dtable$x$data[[1]])
          dtable
        }
        else {
          incProgress(1)
          datatable(data.frame("Peptide Datatable"=character()), selection = list(mode='single', selected = '1'))
        }
      })
    }
  )

  ##update selected peptide data 
  selectedPeptideData <- reactive({
    selection <- input$peptideTable_rows_selected
    if (is.null(selection)){
      selection <- 1
    }
    peptide_names <- names(df$metricsData[[selectedID()]]$good_binders[[selectedTranscript()]])
    index <- floor((as.numeric(selection)+1)/2)
    df$metricsData[[selectedID()]]$good_binders[[selectedTranscript()]][[peptide_names[index]]]
  })
  
  ##Add legend for anchor heatmap 
  output$peptideFigureLegend<- renderPlot({
    colors <- colorRampPalette(c("lightblue", "blue"))(99)[seq(1,99,7)]
    color_pos = data.frame(d=as.character(seq(1,99,7)), x1=seq(0.1,1.5,0.1), x2=seq(0.2,1.6,0.1), y1=rep(1,15), y2=rep(1.1,15), colors=colors)
    color_label = data.frame(x=c(0.1,0.8,1.6), y=rep(0.95,3), score=c(0,0.5,1))
    p1 <- ggplot() + 
      scale_y_continuous(limits=c(0.90, 1.2), name="y") + scale_x_continuous(limits=c(0,1.7), name="x") +
      geom_rect(data=color_pos, aes(xmin=x1, xmax=x2, ymin=y1, ymax=y2, fill=colors), color="black", alpha=1) +
      scale_fill_identity()
    p1 <- p1 + geom_text(data=color_label, aes(x=x, y=y, label=score), size=4, fontface =2) +
      annotate(geom="text", x = 0.5, y= 1.18, label="Normalized Anchor Score", size=4, fontface =2) +
      coord_fixed() +
      theme_void() + theme(legend.position = "none", panel.border = element_blank(), plot.margin=margin(0, 0, 0, 0, "cm"))
    print(p1)
  })
  
  ##Anchor Heatmap overlayed on selected peptide sequences 
  output$anchorPlot <- renderPlot({
    if (is.null(df$metricsData)){
      return ()
    }
    withProgress(message = 'Loading Anchor Heatmap', value = 0, {
      if (type() == 2){
        p1 <- ggplot() + annotate(geom="text", x = 10, y = 20, label = "No data available for Class II HLA alleles", size=6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p1)
      }
      else if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0) {
        peptide_data <- df$metricsData[[selectedID()]]$good_binders[[selectedTranscript()]]
        peptide_names <- names(peptide_data)
        for(i in 1:length(peptide_names)){
          peptide_data[[peptide_names[[i]]]]$individual_ic50_calls <- NULL
          peptide_data[[peptide_names[[i]]]]$individual_percentile_calls <- NULL
        }
        peptide_data <- as.data.frame(peptide_data)
        p1 <- ggplot() + scale_x_continuous(limits=c(0,80)) + scale_y_continuous(limits=c(-31, 1)) 
        all_peptides <- list()
        incProgress(0.1)
        for(i in 1:length(peptide_names)){
          mutation_pos <- as.numeric(df$metricsData[[selectedID()]]$good_binders[[selectedTranscript()]][[peptide_names[i]]]$`mutation_position`)
          wt_peptide <- as.character(df$metricsData[[selectedID()]]$good_binders[[selectedTranscript()]][[peptide_names[i]]]$`wt_peptide`)
          df_mt_peptide <- data.frame("aa"=unlist(strsplit(peptide_names[i],"", fixed = TRUE)), "x_pos" = c(1:nchar(peptide_names[i])))
          df_mt_peptide$mutation <- 'not_mutated'
          df_mt_peptide$type <- 'mt'
          df_mt_peptide$y_pos <- (i*2-1)*-1
          df_mt_peptide$length <- nchar(peptide_names[i])
          df_mt_peptide[mutation_pos, 'mutation'] <- 'mutated'
          df_wt_peptide <- data.frame("aa"=unlist(strsplit(wt_peptide,"", fixed = TRUE)), "x_pos" = c(1:nchar(wt_peptide)))
          df_wt_peptide$mutation <- 'not_mutated'
          df_wt_peptide$type <- 'wt'
          df_wt_peptide$y_pos <- (i*2)*-1
          df_wt_peptide$length <- nchar(wt_peptide)
          all_peptides[[i]] <- rbind(df_mt_peptide, df_wt_peptide)
        }
        incProgress(0.4)
        all_peptides <- do.call(rbind, all_peptides)
        peptide_table <- do.call("rbind",lapply(peptide_names, table_formatting, peptide_data))
        peptide_table_filtered <- Filter(function(x) length(unique(x))!=1, peptide_table)
        peptide_table_names <- names(peptide_table_filtered)
        hla_list <- peptide_table_names[grepl("^HLA-*", peptide_table_names)]
        
        hla_data <- data.frame(hla = hla_list)
        hla_sep <- max(nchar(peptide_table$`Peptide Sequence`))
        hla_data$y_pos <- 1
        hla_data$x_pos <- hla_sep/2
        pad <- 3
        all_peptides_multiple_hla <- list()
        incProgress(0.1)
        for(i in 1:length(hla_list)){
          hla_data$x_pos[i] <- hla_data$x_pos[i]+(hla_sep+pad)*(i-1)
          omit_rows <- which(is.na(peptide_table_filtered[names(peptide_table_filtered) == hla_list[[i]]]))*-1
          all_peptides_multiple_hla[[i]] <- all_peptides[!(all_peptides$y_pos %in% omit_rows),]
          all_peptides_multiple_hla[[i]]$color_value <- apply(all_peptides_multiple_hla[[i]], 1, function(x) peptide_coloring(hla_list[[i]],x))
          all_peptides_multiple_hla[[i]]$x_pos <-  all_peptides_multiple_hla[[i]]$x_pos+(hla_sep+pad)*(i-1)
        }
        incProgress(0.2)
        all_peptides_multiple_hla <- do.call(rbind, all_peptides_multiple_hla)
        
        h_line_pos <- data.frame(y_pos = seq(min(all_peptides_multiple_hla['y_pos'])-0.5,max(all_peptides_multiple_hla['y_pos'])-1.5, 2), x_pos=c(min(all_peptides_multiple_hla['x_pos'])-1))
        h_line_pos <- rbind(h_line_pos,data.frame(x_pos = max(all_peptides_multiple_hla['x_pos'])+1, y_pos = seq(min(all_peptides_multiple_hla['y_pos'])-0.5,max(all_peptides_multiple_hla['y_pos'])-1.5, 2)))
        incProgress(0.2)
        p1 <- p1 + 
          geom_rect(data=all_peptides_multiple_hla, aes(xmin=x_pos-0.5, xmax=1+x_pos-0.5, ymin=.5+y_pos, ymax=-.5+y_pos), fill=all_peptides_multiple_hla$color_value) +
          geom_text(data=all_peptides_multiple_hla, aes(x=x_pos, y=y_pos, label = aa, color=mutation), size=5) +
          geom_text(data=hla_data, aes(x=x_pos, y=y_pos, label=hla), size=5, fontface="bold") +
          geom_line(data = h_line_pos, (aes(x=x_pos, y=y_pos, group=y_pos)), linetype = "dashed")
        
        p1 <- p1 + scale_color_manual("mutation", values=c("not_mutated" = "#000000", "mutated" = "#e74c3c"))
        p1 <- p1 + theme_void() + theme(legend.position = "none", panel.border = element_blank())
        print(p1)
      }
      else {
        p1 <- ggplot() + annotate(geom="text", x = 10, y = 20, label = "No data available", size=6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p1)
      }
    })
  }, height = 500, width = 1000)
  
  ##updating IC50 binding score for selected peptide pair 
  bindingScoreDataIC50 <- reactive({
    if (is.null(df$metricsData)){
      return ()
    }
    if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0){
      algorithm_names <- data.frame(algorithms=selectedPeptideData()$individual_ic50_calls$algorithms)
      wt_data <- as.data.frame(selectedPeptideData()$individual_ic50_calls$WT, check.names=FALSE)
      colnames(wt_data) <- paste(colnames(wt_data),"_WT_Score", sep="")
      mt_data <- as.data.frame(selectedPeptideData()$individual_ic50_calls$MT, check.names=FALSE)
      colnames(mt_data) <- paste(colnames(mt_data),"_MT_Score", sep="")
      full_data <- cbind(algorithm_names,mt_data,wt_data) %>%
        gather('col', 'val', colnames(mt_data)[1]:tail(colnames(wt_data), n=1)) %>%
        separate(col, c('HLA_allele','Mutant','Score'), sep='\\_') %>%
        spread('Score', val)
      full_data
    }
    else{
      return()
    }
  })
  
  ##plotting IC5 binding score violin plot 
  output$bindingData_IC50 <- renderPlot({
    withProgress(message = 'Loading Binding Score Plot (IC50)', value = 0, {
      if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0){
        line.data <- data.frame(yintercept = c(500,1000), Cutoffs = c("500nM", "1000nM"), color=c("#28B463","#EC7063"))
        hla_allele_count <- length(unique(bindingScoreDataIC50()$HLA_allele))
        incProgress(0.5)
        p <- ggplot(data=bindingScoreDataIC50() , aes(x=Mutant, y=Score, color=Mutant), trim=FALSE) + geom_violin() + facet_grid(cols = vars(HLA_allele)) + scale_y_continuous(trans="log10") + #coord_trans(y = "log10") +
          stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)) + 
          geom_jitter(data=bindingScoreDataIC50(), aes(shape=algorithms), size=5, stroke = 1, position=position_jitter(0.3)) + 
          scale_shape_manual(values = 0:8) +
          geom_hline(aes(yintercept = yintercept, linetype = Cutoffs), line.data, color=rep(line.data$color, hla_allele_count)) + 
          scale_color_manual(values = rep(c("MT"="#D2B4DE", "WT"="#F7DC6F"),hla_allele_count)) + 
          theme(strip.text = element_text(size=15), axis.text = element_text(size=10), axis.title = element_text(size=15), axis.ticks = element_line(size=3), legend.text = element_text(size=15), legend.title = element_text(size=15))
        incProgress(0.5)
        print(p)}
      else{
        p <- ggplot() + annotate(geom="text", x = 10, y = 20, label = "No data available", size=6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p)
      }
    })
  })
  
  ##updating percentile binding score for selected peptide pair 
  bindingScoreDataPercentile <- reactive({
    if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0){
      algorithm_names <- data.frame(algorithms=selectedPeptideData()$individual_percentile_calls$algorithms)
      wt_data <- as.data.frame(selectedPeptideData()$individual_percentile_calls$WT, check.names=FALSE)
      colnames(wt_data) <- paste(colnames(wt_data),"_WT_Score", sep="")
      mt_data <- as.data.frame(selectedPeptideData()$individual_percentile_calls$MT, check.names=FALSE)
      colnames(mt_data) <- paste(colnames(mt_data),"_MT_Score", sep="")
      full_data <- cbind(algorithm_names,mt_data,wt_data) %>%
        gather('col', 'val', colnames(mt_data)[1]:tail(colnames(wt_data), n=1)) %>%
        separate(col, c('HLA_allele','Mutant','Score'), sep='\\_') %>%
        spread('Score', val)
      full_data
    }
    else{
      return()
    }
  })
  
  ##plotting percentile binding score violin plot 
  output$bindingData_percentile <- renderPlot({
    withProgress(message = 'Loading Binding Score Plot (Percentile)', value = 0, {
      if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0){
        line.data <- data.frame(yintercept = c(0.5,2), Cutoffs = c("0.5%", "2%"), color=c("#28B463","#EC7063"))
        hla_allele_count <- length(unique(bindingScoreDataPercentile()$HLA_allele))
        incProgress(0.5)
        p <- ggplot(data=bindingScoreDataPercentile() , aes(x=Mutant, y=Score, color=Mutant), trim=FALSE) + geom_violin() + facet_grid(cols = vars(HLA_allele)) + scale_y_continuous(trans="log10") + #coord_trans(y = "log10") +
          stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)) + 
          geom_jitter(data=bindingScoreDataPercentile(), aes(shape=algorithms), size=5, stroke = 1, position=position_jitter(0.3)) + 
          scale_shape_manual(values = 0:8) +
          geom_hline(aes(yintercept = yintercept, linetype = Cutoffs), line.data, color=rep(line.data$color, hla_allele_count)) + 
          scale_color_manual(values = rep(c("MT"="#D2B4DE", "WT"="#F7DC6F"),hla_allele_count)) + 
          theme(strip.text = element_text(size=15), axis.text = element_text(size=10), axis.title = element_text(size=15), axis.ticks = element_line(size=3), legend.text = element_text(size=15), legend.title = element_text(size=15))
        incProgress(0.5)
        print(p)}
      else{
        p <- ggplot() + annotate(geom="text", x = 10, y = 20, label = "No data available", size=6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p)
      }
    })
  })
  
  ##plotting binding data table with IC50 and percentile values
  output$bindingDatatable <- renderDT({
    withProgress(message = 'Loading binding datatable', value = 0, {
      if (length(df$metricsData[[selectedID()]]$good_binders_transcript) != 0){
      binding_data <- bindingScoreDataIC50()
      names(binding_data)[names(binding_data) == 'Score'] <- 'IC50 Score'
      binding_data['% Score'] <- bindingScoreDataPercentile()['Score']
      binding_data['Score'] <- paste(round(as.numeric(binding_data$`IC50 Score`),2)," (%: ",round(as.numeric(binding_data$`% Score`),2),")", sep="")
      binding_data['IC50 Score'] <- NULL
      binding_data['% Score'] <- NULL
      binding_reformat <- dcast(binding_data, HLA_allele + Mutant ~ algorithms, value.var = "Score")
      incProgress(1)
      dtable <- datatable(binding_reformat, options =list(
          pageLength = 10,
          lengthMenu = c(10),
          rowCallback = JS('function(row, data, index, rowId) {',
                           'console.log(rowId)','if(((rowId+1) % 4) == 3 || ((rowId+1) % 4) == 0) {',
                           'row.style.backgroundColor = "#E0E0E0";','}','}')
        )) %>% formatStyle('Mutant', fontWeight = styleEqual('MT','bold'), color = styleEqual('MT', '#E74C3C'))
        dtable
      }
      else {
        incProgress(1)
        datatable(data.frame("Binding Predictions Datatable"=character()))
      }
    })
  })

  
  ##############################EXPORT TAB##############################################

  #evalutation overview table
  output$checked <- renderTable({
    if (is.null(df$mainTable)){
      return ()
    }
    Evaluation <- data.frame(selected=shinyValue("selecter_",nrow(df$mainTable), df$mainTable))
    data <- as.data.frame(table(Evaluation))
    data$Count <- data$Freq
    data$Freq <- NULL
    data
  })
  
  #export table display with options to download 
  output$ExportTable = renderDataTable({
    if (is.null(df$mainTable)){
      return ()
    }
    data <- df$mainTable[, !(colnames(df$mainTable) == "Evaluation") & !(colnames(df$mainTable) == "Eval") & !(colnames(df$mainTable) == "Select") & !(colnames(df$mainTable) == "Tier Count") & !(colnames(df$mainTable) == "Comments") & !(colnames(df$mainTable) == "Gene of Interest")  & !(colnames(df$mainTable) == "Mutated Positions") & !(colnames(df$mainTable) == "Best HLA allele")]
    col_names <- colnames(data)
    data <- data.frame(data, Evaluation=shinyValue("selecter_",nrow(df$mainTable), df$mainTable))
    colnames(data) <- c(col_names,"Evaluation")
    comments <- data.frame("ID" = row.names(df$comments), Comments = df$comments[,1])
    data <- join(data, comments)
    data
    }, escape = FALSE, server = FALSE, rownames = FALSE,
    options=list(dom = 'Bfrtip', 
                 buttons = list(list(
                   extend = 'csvHtml5', 
                   filename = input$exportFileName,
                   fieldSeparator = '\t',
                   text = 'Download as TSV',
                   extension = '.tsv'), list(
                     extend = 'excel', 
                     filename = input$exportFileName,
                     text = 'Download as excel'
                   )
                 ),
                 initComplete = htmlwidgets::JS(
                   "function(settings, json) {",
                   paste0("$(this.api().table().header()).css({'font-size': '", '8pt', "'});"), 
                   "}")
    ),
    selection = 'none',
    extensions = c("Buttons"))
  })
