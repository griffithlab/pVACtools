library(shiny)
library(ggplot2)
library(DT)
library(reshape2)
library(jsonlite)
library(tibble)
library(tidyr)
library(plyr)
library(dplyr)
library("stringr")
library(plotly)
library(shinyWidgets)
library(colourpicker)

source("anchor_and_helper_functions.R", local = TRUE)
source("styling.R")

#specify max shiny app upload size (currently 300MB)
options(shiny.maxRequestSize = 300 * 1024^2)
options(shiny.host = '0.0.0.0')
options(shiny.port = 3333)

server <- shinyServer(function(input, output, session) {
  ## pVACtools version

  output$version <- renderText({"pVACtools version 4.2.0"})

  ##############################DATA UPLOAD TAB###################################
  ## helper function defined for generating shinyInputs in mainTable (Evaluation dropdown menus)
  shinyInput <- function(data, FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), label = NULL, ..., selected = data[i, "Evaluation"]))
    }
    inputs
  }
  ## helper function defined for generating shinyInputs in mainTable (Investigate button)
  shinyInputSelect <- function(FUN, len, id, ...) {
    inputs <- character(len)
    for (i in seq_len(len)) {
      inputs[i] <- as.character(FUN(paste0(id, i), ...))
    }
    inputs
  }
  ## helper function defined for getting values of shinyInputs in mainTable (Evaluation dropdown menus)
  shinyValue <- function(id, len, data) {
    unlist(lapply(seq_len(len), function(i) {
      value <- input[[paste0(id, i)]]
      if (is.null(value)) {
        data[i, "Evaluation"]
      } else {
        value
      }
    }))
  }
  #reactive values defined for row selection, main table, metrics data, additional data, and dna cutoff
  df <- reactiveValues(
    selectedRow = 1,
    mainTable = NULL,
    dna_cutoff = NULL,
    metricsData = NULL,
    additionalData = NULL,
    gene_list = NULL,
    binding_threshold = NULL,
    use_allele_specific_binding_thresholds = NULL,
    aggregate_inclusion_binding_threshold = NULL,
    percentile_threshold = NULL,
    allele_specific_binding_thresholds = NULL,
    allele_expr = NULL,
    anchor_mode = NULL,
    anchor_contribution = NULL,
    comments = data.frame("N/A"),
    pageLength = 10
  )
  #Option 1: User uploaded main aggregate report file
  observeEvent(input$mainDataInput$datapath, {
    #session$sendCustomMessage("unbind-DT", "mainTable")
    mainData <- read.table(input$mainDataInput$datapath, sep = "\t",  header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(mainData) <- mainData[1, ]
    mainData <- mainData[-1, ]
    row.names(mainData) <- NULL
    mainData$`Eval` <- shinyInput(mainData, selectInput, nrow(mainData), "selecter_", choices = c("Pending", "Accept", "Reject", "Review"), width = "90px")
    mainData$Select <- shinyInputSelect(actionButton, nrow(mainData), "button_", label = "Investigate", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
    mainData$`IC50 MT` <- as.numeric(mainData$`IC50 MT`)
    mainData$`%ile MT` <- as.numeric(mainData$`%ile MT`)
    mainData$`RNA Depth` <- as.integer(mainData$`RNA Depth`)
    mainData$`TSL`[is.na(mainData$`TSL`)] <- "NA"
    df$mainTable <- mainData
    df$metricsData <- NULL
  })
  #Option 1: User uploaded metrics file
  observeEvent(input$metricsDataInput, {
    df$metricsData <- fromJSON(input$metricsDataInput$datapath)
    df$binding_threshold <- df$metricsData$`binding_threshold`
    df$use_allele_specific_binding_thresholds <- df$metricsData$`use_allele_specific_binding_thresholds`
    df$allele_specific_binding_thresholds <- df$metricsData$`allele_specific_binding_thresholds`
    df$aggregate_inclusion_binding_threshold <- df$metricsData$`aggregate_inclusion_binding_threshold`
    df$percentile_threshold <- df$metricsData$`percentile_threshold`
    df$dna_cutoff <- df$metricsData$vaf_clonal
    df$allele_expr <- df$metricsData$allele_expr_threshold
    df$anchor_mode <- ifelse(df$metricsData$`allele_specific_anchors`, "allele-specific", "default")
    df$allele_specific_anchors <- df$metricsData$`allele_specific_anchors`
    df$anchor_contribution <- df$metricsData$`anchor_contribution_threshold`
    hla <- df$metricsData$alleles
    converted_hla_names <- unlist(lapply(hla, function(x) {
      if (grepl("HLA-", x)) {
        strsplit(x, "HLA-")[[1]][2]
      } else {
        x
      }
    }))
    if (!("Ref Match" %in% colnames(df$mainTable))) {
      df$mainTable$`Ref Match` <- "Not Run"
    }
    columns_needed <- c("ID", converted_hla_names, "Gene", "AA Change", "Num Passing Transcripts", "Best Peptide", "Best Transcript", "TSL",	"Allele",
                        "Pos", "Prob Pos", "Num Passing Peptides", "IC50 MT",	"IC50 WT", "%ile MT",	"%ile WT", "RNA Expr", "RNA VAF",
                        "Allele Expr", "RNA Depth", "DNA VAF",	"Tier",	"Ref Match", "Evaluation", "Eval", "Select")
    if ("Comments" %in% colnames(df$mainTable)) {
      columns_needed <- c(columns_needed, "Comments")
      df$comments <- data.frame(data = df$mainTable$`Comments`, nrow = nrow(df$mainTable), ncol = 1)
    }else {
      df$comments <- data.frame(matrix("No comments", nrow = nrow(df$mainTable)), ncol = 1)
    }
    df$mainTable <- df$mainTable[, columns_needed]
    df$mainTable$`Tier Count` <- apply(df$mainTable, 1, function(x) tier_numbers(x, df$anchor_contribution, df$dna_cutoff, df$allele_expr, x["Pos"], x["Allele"], x["TSL"], df$metricsData[1:15], df$anchor_mode, df$allele_specific_binding_thresholds, df$use_allele_specific_binding_thresholds, df$binding_threshold))
    df$mainTable$`Gene of Interest` <- apply(df$mainTable, 1, function(x) {any(x["Gene"] == df$gene_list)})
    rownames(df$comments) <- df$mainTable$ID
    df$mainTable$`Scaled BA` <- apply(df$mainTable, 1, function(x) scale_binding_affinity(df$allele_specific_binding_thresholds, df$use_allele_specific_binding_thresholds, df$binding_threshold, x["Allele"], x["IC50 MT"]))
    df$mainTable$`Scaled percentile` <- apply(df$mainTable, 1, function(x) {ifelse(is.null(df$percentile_threshold), as.numeric(x["%ile MT"]), as.numeric(x["%ile MT"]) / (df$percentile_threshold))})
    df$mainTable$`Bad TSL` <- apply(df$mainTable, 1, function(x) {x["TSL"] == "NA" | (x["TSL"] != "NA" & x["TSL"] != "Not Supported" & x["TSL"] > df$metricsData$maximum_transcript_support_level)})
    df$mainTable$`Col RNA Expr` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["RNA Expr"]), 0, x["RNA Expr"])})
    df$mainTable$`Col RNA VAF` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["RNA VAF"]), 0, x["RNA VAF"])})
    df$mainTable$`Col Allele Expr` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["Allele Expr"]), 0, x["Allele Expr"])})
    df$mainTable$`Col RNA Depth` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["RNA Depth"]), 0, x["RNA Depth"])})
    df$mainTable$`Col DNA VAF` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["DNA VAF"]), 0, x["DNA VAF"])})
    if (is.null(df$percentile_threshold)) {
      df$mainTable$`Percentile Fail` <- apply(df$mainTable, 1, function(x) {FALSE})
    }else {
      df$mainTable$`Percentile Fail` <- apply(df$mainTable, 1, function(x) {ifelse(as.numeric(x["%ile MT"]) > as.numeric(df$percentile_threshold), TRUE, FALSE)})
    }
    df$mainTable$`Has Prob Pos` <- apply(df$mainTable, 1, function(x) {ifelse(x["Prob Pos"] != "None", TRUE, FALSE)})
  })
  #Option 1: User uploaded additional data file
  observeEvent(input$additionalDataInput, {
    addData <- read.table(input$additionalDataInput$datapath, sep = "\t",  header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(addData) <- addData[1, ]
    addData <- addData[-1, ]
    row.names(addData) <- NULL
    df$additionalData <- addData
  })
  #Option 1: User uploaded additional gene list
  observeEvent(input$gene_list, {
    gene_list <- read.table(input$gene_list$datapath, sep = "\t",  header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    df$gene_list <- gene_list
    df$mainTable$`Gene of Interest` <- apply(df$mainTable, 1, function(x) {any(x["Gene"] == df$gene_list)})
  })
  #Option 2: Load from HCC1395 demo data from github
   observeEvent(input$loadDefaultmain, {
     ## Class I demo aggregate report
     #session$sendCustomMessage("unbind-DT", "mainTable")
     withProgress(message = "Loading Demo Data", value = 0, {
       load(url("https://github.com/griffithlab/pVACtools/raw/52ced64ad04bf627ef900fa6ade38f5d366a783f/pvactools/tools/pvacview/HCC1395_demo_data.rda"))
       incProgress(0.3)
       #data <- getURL("https://raw.githubusercontent.com/griffithlab/pVACtools/0359d15c/pvactools/tools/pvacview/data/H_NJ-HCC1395-HCC1395.Class_I.all_epitopes.aggregated.tsv")
       #mainData <- read.table(text = data, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
       colnames(mainData) <- mainData[1, ]
       mainData <- mainData[-1, ]
       row.names(mainData) <- NULL
       mainData$`Eval` <- shinyInput(mainData, selectInput, nrow(mainData), "selecter_", choices = c("Pending", "Accept", "Reject", "Review"), width = "90px")
       mainData$Select <- shinyInputSelect(actionButton, nrow(mainData), "button_", label = "Investigate", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
       mainData$`IC50 MT` <- as.numeric(mainData$`IC50 MT`)
       mainData$`%ile MT` <- as.numeric(mainData$`%ile MT`)
       mainData$`RNA Depth` <- as.integer(mainData$`RNA Depth`)
       mainData$`TSL`[is.na(mainData$`TSL`)] <- "NA"
       df$mainTable <- mainData
       incProgress(0.1)
       ## Class I demo metrics file
       #metricsdata <- getURL("https://raw.githubusercontent.com/griffithlab/pVACtools/0359d15c/pvactools/tools/pvacview/data/H_NJ-HCC1395-HCC1395.Class_I.all_epitopes.aggregated.metrics.json")
       #df$metricsData <- fromJSON(txt = metricsdata)
       df$metricsData <- metricsData
       df$binding_threshold <- df$metricsData$`binding_threshold`
       df$allele_specific_binding_thresholds <- df$metricsData$`allele_specific_binding_thresholds`
       df$use_allele_specific_binding_thresholds <- df$metricsData$`use_allele_specific_binding_thresholds`
       df$aggregate_inclusion_binding_threshold <- df$metricsData$`aggregate_inclusion_binding_threshold`
       df$percentile_threshold <- df$metricsData$`percentile_threshold`
       df$dna_cutoff <- df$metricsData$vaf_clonal
       df$allele_expr <- df$metricsData$allele_expr_threshold
       df$anchor_mode <- ifelse(df$metricsData$`allele_specific_anchors`, "allele-specific", "default")
       df$allele_specific_anchors <- df$metricsData$`allele_specific_anchors`
       df$anchor_contribution <- df$metricsData$`anchor_contribution_threshold`
       hla <- df$metricsData$alleles
       incProgress(0.1)
       converted_hla_names <- unlist(lapply(hla, function(x) {
         if (grepl("HLA-", x)) {
           strsplit(x, "HLA-")[[1]][2]
         } else {
           x
         }
       }))
       if (!("Ref Match" %in% colnames(df$mainTable))) {
         df$mainTable$`Ref Match` <- "Not Run"
       }
       columns_needed <- c("ID", converted_hla_names, "Gene", "AA Change", "Num Passing Transcripts", "Best Peptide", "Best Transcript", "TSL",	"Allele",
                           "Pos", "Prob Pos", "Num Passing Peptides", "IC50 MT",	"IC50 WT", "%ile MT",	"%ile WT", "RNA Expr", "RNA VAF",
                           "Allele Expr", "RNA Depth", "DNA VAF",	"Tier",	"Ref Match", "Evaluation", "Eval", "Select")
       if ("Comments" %in% colnames(df$mainTable)) {
         columns_needed <- c(columns_needed, "Comments")
         df$comments <- data.frame(data = df$mainTable$`Comments`, nrow = nrow(df$mainTable), ncol = 1)
       }else {
         df$comments <- data.frame(matrix("No comments", nrow = nrow(df$mainTable)), ncol = 1)
       }
       df$mainTable <- df$mainTable[, columns_needed]
       df$mainTable$`Tier Count` <- apply(df$mainTable, 1, function(x) tier_numbers(x, df$anchor_contribution, df$dna_cutoff, df$allele_expr, x["Pos"], x["Allele"], x["TSL"], df$metricsData[1:15], df$anchor_mode, df$allele_specific_binding_thresholds, df$use_allele_specific_binding_thresholds, df$binding_threshold))
       df$mainTable$`Gene of Interest` <- apply(df$mainTable, 1, function(x) {any(x["Gene"] == df$gene_list)})
       if ("Comments" %in% colnames(df$mainTable)) {
         df$comments <- data.frame(data = df$mainTable$`Comments`, nrow = nrow(df$mainTable), ncol = 1)
       }else {
         df$comments <- data.frame(matrix("No comments", nrow = nrow(df$mainTable)), ncol = 1)
       }
       rownames(df$comments) <- df$mainTable$ID
       incProgress(0.2)
       ## Class II additional demo aggregate report
       add_data <- getURL("https://raw.githubusercontent.com/griffithlab/pVACtools/0359d15c/pvactools/tools/pvacview/data/H_NJ-HCC1395-HCC1395.Class_II.all_epitopes.aggregated.tsv")
       addData <- read.table(text = add_data, sep = "\t",  header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
       colnames(addData) <- addData[1, ]
       addData <- addData[-1, ]
       row.names(addData) <- NULL
       df$additionalData <- addData
       incProgress(0.1)
       ## Hotspot gene list autoload
       gene_data <- getURL("https://raw.githubusercontent.com/griffithlab/pVACtools/0359d15c/pvactools/tools/pvacview/data/cancer_census_hotspot_gene_list.tsv")
       gene_list <- read.table(text = gene_data, sep = "\t",  header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
       df$gene_list <- gene_list
       df$mainTable$`Gene of Interest` <- apply(df$mainTable, 1, function(x) {any(x["Gene"] == df$gene_list)})
       df$mainTable$`Scaled BA` <- apply(df$mainTable, 1, function(x) scale_binding_affinity(df$allele_specific_binding_thresholds, df$use_allele_specific_binding_thresholds, df$binding_threshold, x["Allele"], x["IC50 MT"]))
       df$mainTable$`Scaled percentile` <- apply(df$mainTable, 1, function(x) {ifelse(is.null(df$percentile_threshold), as.numeric(x["%ile MT"]), as.numeric(x["%ile MT"]) / (df$percentile_threshold))})
       df$mainTable$`Bad TSL` <- apply(df$mainTable, 1, function(x) {x["TSL"] == "NA" | (x["TSL"] != "NA" & x["TSL"] != "Not Supported" & x["TSL"] > df$metricsData$maximum_transcript_support_level)})
       df$mainTable$`Col RNA Expr` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["RNA Expr"]), 0, x["RNA Expr"])})
       df$mainTable$`Col RNA VAF` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["RNA VAF"]), 0, x["RNA VAF"])})
       df$mainTable$`Col Allele Expr` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["Allele Expr"]), 0, x["Allele Expr"])})
       df$mainTable$`Col RNA Depth` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["RNA Depth"]), 0, x["RNA Depth"])})
       df$mainTable$`Col DNA VAF` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["DNA VAF"]), 0, x["DNA VAF"])})
       incProgress(0.1)
       if (is.null(df$percentile_threshold)) {
         df$mainTable$`Percentile Fail` <- apply(df$mainTable, 1, function(x) {FALSE})
       }else {
         df$mainTable$`Percentile Fail` <- apply(df$mainTable, 1, function(x) {ifelse(as.numeric(x["%ile MT"]) > as.numeric(df$percentile_threshold), TRUE, FALSE)})
       }
       df$mainTable$`Has Prob Pos` <- apply(df$mainTable, 1, function(x) {ifelse(x["Prob Pos"] != "None", TRUE, FALSE)})
      updateTabItems(session, "tabs", "explore")
      incProgress(0.1)
     })
   })
   ##Clear file inputs if demo data load button is clicked
   output$aggregate_report_ui <- renderUI({
     input$loadDefaultmain
     fileInput(inputId = "mainDataInput", label = "1. Neoantigen Candidate Aggregate Report (tsv required)",
               accept =  c("text/tsv", "text/tab-separated-values,text/plain", ".tsv"))
   })
   output$metrics_ui <- renderUI({
     input$loadDefaultmain
     fileInput(inputId = "metricsDataInput", label = "2. Neoantigen Candidate Metrics file (json required)",
               accept = c("application/json", ".json"))
   })
   output$add_file_ui <- renderUI({
     input$loadDefaultmain
     fileInput(inputId = "additionalDataInput", label = "3. Additional Neoantigen Candidate Aggregate Report (tsv required)",
               accept =  c("text/tsv", "text/tab-separated-values,text/plain", ".tsv"))
   })

  ##Visualize button
  observeEvent(input$visualize, {
    updateTabItems(session, "tabs", "explore")
  })
  ##Parameter UIs
  output$allele_specific_anchors_ui <- renderUI({
    current_is_allele_specific_anchors_set <- df$allele_specific_anchors
    checkboxInput(
      "use_anchor",
      "If you want to use allele-specific anchor calculations, please check this box. Otherwise anchors will be calculated as 1,2 and n-1,n for n-mer peptides.",
      value = current_is_allele_specific_anchors_set,
      width = NULL
    )
  })
  output$anchor_contribution_ui <- renderUI({
    current_anchor_contribution_threshold <- df$anchor_contribution
    sliderInput("anchor_contribution", "Contribution cutoff for determining anchor locations", 0.5, 0.9, current_anchor_contribution_threshold, step = 0.1, width = 400)
  })
  output$binding_threshold_ui <- renderUI({
    current_binding <- df$binding_threshold
    max_cutoff <- df$aggregate_inclusion_binding_threshold
    numericInput("binding_threshold", "Binding Threshold", current_binding, min = 0, max = max_cutoff, step = 10, width = 500)
  })
  output$allele_specific_binding_ui <- renderUI({
    current_is_allele_specific_binding_set <- df$use_allele_specific_binding_thresholds
    checkboxInput(
      "allele_specific_binding",
      "If you want to use allele-specific binding thresholds for tiering purposes please check this box.",
      value = current_is_allele_specific_binding_set,
      width = NULL
    )
  })
  output$percentile_threshold_ui <- renderUI({
    current_percentile <- df$percentile_threshold
    numericInput("percentile_threshold", "Percentile Threshold", current_percentile, min = 0, max = 100, step = 0.01, width = 500)
  })
  output$dna_cutoff_ui <- renderUI({
    current_dna_cutoff <- df$dna_cutoff
    numericInput("dna_cutoff", "Clonal DNA VAF (Anything lower than 1/2 of chosen VAF level will be considered subclonal)", current_dna_cutoff, min = 0, max = 1, step = 0.01, width = 500)
  })
  output$allele_expr_ui <- renderUI({
    current_allele_expr <- df$allele_expr
    numericInput("allele_expr", "Allele Expression cutoff to be considered a Pass variant. Note that this criteria is also used in determining Anchor and Subclonal variants.", current_allele_expr, min = 0, max = 100, step = 0.1, width = 500)
  })
  #reactions for once "regenerate table" command is submitted
  observeEvent(input$submit, {
    session$sendCustomMessage("unbind-DT", "mainTable")
    df$binding_threshold <- input$binding_threshold
    df$use_allele_specific_binding_thresholds <- input$allele_specific_binding
    df$percentile_threshold <- input$percentile_threshold
    df$dna_cutoff <- input$dna_cutoff
    df$allele_expr <- input$allele_expr
    df$allele_specific_anchors <- input$use_anchor
    df$anchor_contribution <- input$anchor_contribution
    df$mainTable$`Evaluation` <- shinyValue("selecter_", nrow(df$mainTable), df$mainTable)
    if (input$use_anchor) {
      df$anchor_mode <- "allele-specific"
      df$anchor_contribution <- input$anchor_contribution
    }else {
      df$anchor_mode <- "default"
    }
    df$mainTable$`Tier` <- apply(df$mainTable, 1, function(x) tier(x, df$anchor_contribution, input$dna_cutoff, input$allele_expr, x["Pos"], x["Allele"], x["TSL"], df$metricsData[1:15], df$anchor_mode, df$use_allele_specific_binding_thresholds, df$binding_threshold))
    df$mainTable$`Tier Count` <- apply(df$mainTable, 1, function(x) tier_numbers(x, df$anchor_contribution, input$dna_cutoff, input$allele_expr, x["Pos"], x["Allele"], x["TSL"], df$metricsData[1:15], df$anchor_mode, df$allele_specific_binding_thresholds, df$use_allele_specific_binding_thresholds, df$binding_threshold))
    df$mainTable$`Scaled BA` <- apply(df$mainTable, 1, function(x) scale_binding_affinity(df$allele_specific_binding_thresholds, df$use_allele_specific_binding_thresholds, df$binding_threshold, x["Allele"], x["IC50 MT"]))
    df$mainTable$`Scaled percentile` <- apply(df$mainTable, 1, function(x) {ifelse((is.null(df$percentile_threshold) || is.na(df$percentile_threshold)), as.numeric(x["%ile MT"]), as.numeric(x["%ile MT"]) / (df$percentile_threshold))})
    if (is.null(df$percentile_threshold) || is.na(df$percentile_threshold)) {
      df$mainTable$`Percentile Fail` <- apply(df$mainTable, 1, function(x) {FALSE})
    }else {
      df$mainTable$`Percentile Fail` <- apply(df$mainTable, 1, function(x) {ifelse(as.numeric(x["%ile MT"]) > as.numeric(df$percentile_threshold), TRUE, FALSE)})
    }
    tier_sorter <- c("Pass", "LowExpr", "Anchor", "Subclonal", "Poor", "NoExpr")
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
    df$mainTable$Select <- shinyInputSelect(actionButton, nrow(df$mainTable), "button_", label = "Investigate", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
    df$mainTable$`Eval` <- shinyInput(df$mainTable, selectInput, nrow(df$mainTable), "selecter_", choices = c("Pending", "Accept", "Reject", "Review"), width = "90px")
  })
  #reset tier-ing with original parameters
  observeEvent(input$reset_params, {
    session$sendCustomMessage("unbind-DT", "mainTable")
    df$binding_threshold <- df$metricsData$`binding_threshold`
    df$allele_specific_binding_thresholds <- df$metricsData$`allele_specific_binding_thresholds`
    df$use_allele_specific_binding_thresholds <- df$metricsData$`use_allele_specific_binding_thresholds`
    df$percentile_threshold <- df$metricsData$`percentile_threshold`
    df$dna_cutoff <- df$metricsData$`vaf_clonal`
    df$allele_expr <- df$metricsData$`allele_expr`
    df$anchor_mode <- ifelse(df$metricsData$`allele_specific_anchors`, "allele-specific", "default")
    df$allele_specific_anchors <- df$metricsData$`allele_specific_anchors`
    df$anchor_contribution <- df$metricsData$`anchor_contribution_threshold`
    df$mainTable$`Evaluation` <- shinyValue("selecter_", nrow(df$mainTable), df$mainTable)
    df$mainTable$`Tier` <- apply(df$mainTable, 1, function(x) tier(x, df$anchor_contribution, df$dna_cutoff, df$allele_expr, x["Pos"], x["Allele"], x["TSL"], df$metricsData[1:15], df$anchor_mode, df$use_allele_specific_binding_thresholds, df$binding_threshold))
    df$mainTable$`Tier Count` <- apply(df$mainTable, 1, function(x) tier_numbers(x, df$anchor_contribution, df$dna_cutoff, df$allele_expr, x["Pos"], x["Allele"], x["TSL"], df$metricsData[1:15], df$anchor_mode, df$allele_specific_binding_thresholds, df$use_allele_specific_binding_thresholds, df$binding_threshold))
    df$mainTable$`Scaled BA` <- apply(df$mainTable, 1, function(x) scale_binding_affinity(df$allele_specific_binding_thresholds, df$use_allele_specific_binding_thresholds, df$binding_threshold, x["Allele"], x["IC50 MT"]))
    df$mainTable$`Scaled percentile` <- apply(df$mainTable, 1, function(x) {ifelse(is.null(df$percentile_threshold), as.numeric(x["%ile MT"]), as.numeric(x["%ile MT"]) / (df$percentile_threshold))})
    if (is.null(df$percentile_threshold)) {
      df$mainTable$`Percentile Fail` <- apply(df$mainTable, 1, function(x) {FALSE})
    }else {
      df$mainTable$`Percentile Fail` <- apply(df$mainTable, 1, function(x) {ifelse(as.numeric(x["%ile MT"]) > as.numeric(df$percentile_threshold), TRUE, FALSE)})
    }
    tier_sorter <- c("Pass", "LowExpr", "Anchor", "Subclonal", "Poor", "NoExpr")
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
    df$mainTable$Select <- shinyInputSelect(actionButton, nrow(df$mainTable), "button_", label = "Investigate", onclick = 'Shiny.onInputChange(\"select_button\",  this.id)')
    df$mainTable$`Eval` <- shinyInput(df$mainTable, selectInput, nrow(df$mainTable), "selecter_", choices = c("Pending", "Accept", "Reject", "Review"), width = "90px")
  })
  #determine hla allele count in order to generate column tooltip locations correctly
  hla_count <- reactive({
    which(colnames(df$mainTable) == "Gene") - 1
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
  output$paramTable <- renderTable(
    data <- data.frame(
      "Parameter" = c("Tumor Purity", "VAF Clonal", "VAF Subclonal", "Allele Expression for Passing Variants",
                      "Binding Threshold", "Binding Threshold for Inclusion into Metrics File", "Maximum TSL",
                      "Percentile Threshold", "Allele Specific Binding Thresholds",
                      "MT Top Score Metric", "WT Top Score Metric",
                      "Allele Specific Anchors Used", "Anchor Contribution Threshold"),
      "Value" = c(if (is.null(df$metricsData$tumor_purity)) {"NULL"}else {df$metricsData$tumor_purity},
                  df$metricsData$`vaf_clonal`, df$metricsData$`vaf_subclonal`, df$metricsData$`allele_expr_threshold`,
                  df$metricsData$binding_threshold, df$metricsData$`aggregate_inclusion_binding_threshold`,
                  df$metricsData$maximum_transcript_support_level,
                  if (is.null(df$metricsData$percentile_threshold)) {"NULL"}else { df$metricsData$percentile_threshold},
                  df$metricsData$use_allele_specific_binding_thresholds,
                  df$metricsData$mt_top_score_metric, df$metricsData$wt_top_score_metric,
                  df$metricsData$allele_specific_anchors, df$metricsData$anchor_contribution_threshold)
    ), digits = 3
  )
  output$bindingParamTable <- renderTable(
    if (df$metricsData$use_allele_specific_binding_thresholds) {
      data <- data.frame(
        "HLA Alleles" = df$metricsData$alleles,
        "Binding Cutoffs" = unlist(lapply(df$metricsData$alleles, function(x) {
          if (x %in% names(df$metricsData$allele_specific_binding_thresholds)) {
            df$metricsData$allele_specific_binding_thresholds[[x]]
          } else {
            df$metricsData$binding_threshold
          }
        }
        )))
    } else {
      data <- data.frame(
        "HLA Alleles" = df$metricsData$alleles,
        "Binding Cutoffs" = unlist(lapply(df$metricsData$alleles, function(x) df$metricsData$binding_threshold))
      )
    }
  )
  output$currentParamTable <- renderTable(
    data <- data.frame(
      "Parameter" = c("VAF Clonal", "VAF Subclonal", "Allele Expression for Passing Variants",
                      "Binding Threshold", "Binding Threshold for Inclusion into Metrics File", "Maximum TSL",
                      "Percentile Threshold", "Allele Specific Binding Thresholds",
                      "MT Top Score Metric", "WT Top Score Metric",
                      "Allele Specific Anchors Used", "Anchor Contribution Threshold"),
      "Value" = c(
        df$dna_cutoff,
        df$dna_cutoff / 2,
        df$allele_expr,
        df$binding_threshold,
        df$metricsData$`aggregate_inclusion_binding_threshold`,
        df$metricsData$maximum_transcript_support_level,
        if (is.null(df$percentile_threshold) || is.na(df$percentile_threshold)) {"NULL"}else { df$percentile_threshold},
        df$use_allele_specific_binding_thresholds,
        df$metricsData$mt_top_score_metric,
        df$metricsData$wt_top_score_metric,
        df$allele_specific_anchors, df$anchor_contribution)
    ), digits = 3
  )
  output$currentBindingParamTable <- renderTable(
    if (df$use_allele_specific_binding_thresholds) {
      data <- data.frame(
        "HLA Alleles" = df$metricsData$alleles,
        "Binding Cutoffs" = unlist(lapply(df$metricsData$alleles, function(x) {
          if (x %in% names(df$metricsData$allele_specific_binding_thresholds)) {
            df$metricsData$allele_specific_binding_thresholds[[x]]
          } else {
            df$binding_threshold
          }
        }
        )))
    } else {
      data <- data.frame(
        "HLA Alleles" = df$metricsData$alleles,
        "Binding Cutoffs" = unlist(lapply(df$metricsData$alleles, function(x) df$binding_threshold))
      )
    }
  )
  output$comment_text <- renderUI({
    if (is.null(df$mainTable)) {
      return(HTML("N/A"))
    }
    HTML(paste(df$comments[selectedID(), 1]))
  })
  observeEvent(input$page_length, {
    if (is.null(df$mainTable)) {
      return()
    }
    df$pageLength <- as.numeric(input$page_length)
    session$sendCustomMessage("unbind-DT", "mainTable")
    df$mainTable$`Evaluation` <- shinyValue("selecter_", nrow(df$mainTable), df$mainTable)
    df$mainTable$`Eval` <- shinyInput(df$mainTable, selectInput, nrow(df$mainTable), "selecter_", choices = c("Pending", "Accept", "Reject", "Review"), width = "90px")
  })
  output$filesUploaded <- reactive({
    val <- !(is.null(df$mainTable) | is.null(df$metricsData))
    print(val)
  })
  outputOptions(output, "filesUploaded", suspendWhenHidden = FALSE)
  ##############################PEPTIDE EXPLORATION TAB################################
  ##main table display with color/background/font/border configurations
  output$mainTable <- DT::renderDataTable(
    if (is.null(df$mainTable) | is.null(df$metricsData)) {
      return(datatable(data.frame("Aggregate Report" = character())))
    }else {
      datatable(df$mainTable[, !(colnames(df$mainTable) == "ID") & !(colnames(df$mainTable) == "Evaluation") & !(colnames(df$mainTable) == "Comments")],
                escape = FALSE, callback = JS(callback(hla_count(), df$metricsData$mt_top_score_metric)), class = "stripe",
                options = list(lengthChange = FALSE, dom = "Bfrtip", pageLength = df$pageLength,
                               columnDefs = list(list(defaultContent = "NA", targets = c(hla_count() + 10, (hla_count() + 12):(hla_count() + 17))),
                                                 list(className = "dt-center", targets = c(0:hla_count() - 1)), list(visible = FALSE, targets = c(1:(hla_count()-1), (hla_count()+2), (hla_count()+4), -1:-12)),
                                                 list(orderable = TRUE, targets = 0)), buttons = list(I("colvis")),
                               initComplete = htmlwidgets::JS(
                                 "function(settings, json) {",
                                 paste("$(this.api().table().header()).css({'font-size': '", "10pt", "'});"),
                                 "}"),
                               rowCallback = JS(rowcallback(hla_count(), df$selectedRow - 1)),
                               preDrawCallback = JS("function() {
                                        Shiny.unbindAll(this.api().table().node()); }"),
                               drawCallback = JS("function() { 
                                     Shiny.bindAll(this.api().table().node()); } ")),
                selection = "none",
                extensions = c("Buttons"))
    }
    %>% formatStyle("IC50 MT", "Scaled BA", backgroundColor = styleInterval(c(0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6, 1.8, 2),
                                                                            c("#68F784", "#60E47A", "#58D16F", "#4FBD65", "#47AA5A", "#3F9750", "#F3F171", "#F3E770", "#F3DD6F", "#F0CD5B", "#F1C664", "#FF9999"))
                    , fontWeight = styleInterval(c(1000), c("normal", "bold")), border = styleInterval(c(1000), c("normal", "2px solid red")))
    %>% formatStyle("%ile MT", "Scaled percentile", backgroundColor = styleInterval(c(0.2, 0.4, 0.6, 0.8, 1, 1.25, 1.5, 1.75, 2),
                                                                                    c("#68F784", "#60E47A", "#58D16F", "#4FBD65", "#47AA5A", "#F3F171", "#F3E770", "#F3DD6F", "#F1C664", "#FF9999")))
    %>% formatStyle("Tier", color = styleEqual(c("Pass", "Poor", "Anchor", "Subclonal", "LowExpr", "NoExpr"), c("green", "orange", "#b0b002", "#D4AC0D", "salmon", "red")))
    %>% formatStyle(c("RNA VAF"), "Col RNA VAF", background = styleColorBar(range(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "right")
    %>% formatStyle(c("DNA VAF"), "Col DNA VAF", background = styleColorBar(range(0, 1), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "right")
    %>% formatStyle(c("RNA Expr"), "Col RNA Expr", background = styleColorBar(range(0, 50), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "right")
    %>% formatStyle(c("RNA Depth"), "Col RNA Depth", background = styleColorBar(range(0, 200), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "right")
    %>% formatStyle(c("Allele Expr"), "Col Allele Expr", background = styleColorBar(range(0, (max(as.numeric(as.character(unlist(df$mainTable["Col RNA VAF"]))) * 50))), "lightblue"), backgroundSize = "98% 88%", backgroundRepeat = "no-repeat", backgroundPosition = "right")
    %>% formatStyle(c("Allele Expr"), "Tier Count", fontWeight = styleEqual(c("2"), c("bold")), border = styleEqual(c("2"), c("2px solid red")))
    %>% formatStyle(c("IC50 MT", "Allele Expr"), "Tier Count", fontWeight = styleEqual(c("3"), c("bold")), border = styleEqual(c("3"), c("2px solid red")))
    %>% formatStyle(c("IC50 MT"), "Tier Count", fontWeight = styleEqual(c("4"), c("bold")), border = styleEqual(c("4"), c("2px solid red")))
    %>% formatStyle(c("Allele Expr", "%ile MT"), "Tier Count", fontWeight = styleEqual(c("102"), c("bold")), border = styleEqual(c("102"), c("2px solid red")))
    %>% formatStyle(c("IC50 MT", "Allele Expr", "%ile MT"), "Tier Count", fontWeight = styleEqual(c("103"), c("bold")), border = styleEqual(c("103"), c("2px solid red")))
    %>% formatStyle(c("%ile MT"), "Tier Count", fontWeight = styleEqual(c("104"), c("bold")), border = styleEqual(c("104"), c("2px solid red")))
    %>% formatStyle(c("IC50 MT", "%ile MT"), "Tier Count", fontWeight = styleEqual(c("105"), c("bold")), border = styleEqual(c("105"), c("2px solid red")))
    %>% formatStyle(c("IC50 WT", "Pos"), "Tier Count", fontWeight = styleEqual(c("5"), c("bold")), border = styleEqual(c("5"), c("2px solid red")))
    %>% formatStyle(c("DNA VAF"), "Tier Count", fontWeight = styleEqual(c("6"), c("bold")), border = styleEqual(c("6"), c("2px solid red")))
    %>% formatStyle(c("Allele Expr"), "Tier Count", fontWeight = styleEqual(c("7"), c("bold")), border = styleEqual(c("7"), c("2px solid red")))
    %>% formatStyle(c("Gene Expression"), "Tier Count", fontWeight = styleEqual(c("8"), c("bold")), border = styleEqual(c("8"), c("2px solid red")))
    %>% formatStyle(c("RNA VAF", "RNA Depth"), "Tier Count", fontWeight = styleEqual(c("8"), c("bold")), border = styleEqual(c("8"), c("2px solid green")))
    %>% formatStyle(c("RNA Expr", "Tier Count"), fontWeight = styleEqual(c("9"), c("bold")), border = styleEqual(c("9"), c("2px solid red")))
    %>% formatStyle(c("RNA VAF"), "Tier Count", fontWeight = styleEqual(c("10"), c("bold")), border = styleEqual(c("10"), c("2px solid red")))
    %>% formatStyle(c("RNA VAF", "RNA Expr"), "Tier Count", fontWeight = styleEqual(c("11"), c("bold")), border = styleEqual(c("11"), c("2px solid red")))
    %>% formatStyle(c("IC50 WT", "Pos"), "Tier Count", fontWeight = styleEqual(c("13"), c("bold")), border = styleEqual(c("13"), c("2px solid red")))
    %>% formatStyle(c("DNA VAF"), "Tier Count", fontWeight = styleEqual(c("14"), c("bold")), border = styleEqual(c("14"), c("2px solid red")))
    %>% formatStyle(c("DNA VAF", "IC50 WT", "Pos"), "Tier Count", fontWeight = styleEqual(c("15"), c("bold")), border = styleEqual(c("15"), c("2px solid red")))
    %>% formatStyle(c("IC50 WT", "Pos", "DNA VAF", "Allele Expr"), "Tier Count", fontWeight = styleEqual(c("23"), c("bold")), border = styleEqual(c("23"), c("2px solid red")))
    %>% formatStyle(c("DNA VAF", "Allele Expr"), "Tier Count", fontWeight = styleEqual(c("22"), c("bold")), border = styleEqual(c("22"), c("2px solid red")))
    %>% formatStyle(c("IC50 WT", "Pos", "Allele Expr"), "Tier Count", fontWeight = styleEqual(c("21"), c("bold")), border = styleEqual(c("21"), c("2px solid red")))
    %>% formatStyle(c("Allele Expr"), "Tier Count", fontWeight = styleEqual(c("20"), c("bold")), border = styleEqual(c("20"), c("2px solid red")))
    %>% formatStyle(c("Gene"), "Gene of Interest", fontWeight = styleEqual(c(TRUE), c("bold")), border = styleEqual(c(TRUE), c("2px solid green")))
    %>% formatStyle(c("TSL"), "Bad TSL", fontWeight = styleEqual(c(TRUE), c("bold")), border = styleEqual(c(TRUE), c("2px solid red")))
    %>% formatStyle(c("%ile MT"), "Percentile Fail", border = styleEqual(c(TRUE), c("2px solid red")))
    %>% formatStyle(c("Prob Pos"), "Has Prob Pos", fontWeight = styleEqual(c(TRUE), c("bold")), border = styleEqual(c(TRUE), c("2px solid red")))
    %>% formatStyle(c("Ref Match"), "Ref Match", fontWeight = styleEqual(c("True"), c("bold")), border = styleEqual(c("True"), c("2px solid red")))
    %>% formatStyle("Best Peptide", fontFamily="monospace")
    , server = FALSE)
  #help menu for main table
  observeEvent(input$help, {
    showModal(modalDialog(
      title = "Aggregate Report of Best Candidates by Mutation",
      h5("* Hover over individual column names to see further description of specific columns. (HLA allele columns excluded)"),
      h4(" HLA specific columns:", style = "font-weight: bold"),
      h5(" Number of good binding peptides for each specific HLA-allele.", br(),
         " The same peptide could be counted in multiple columns if it was predicted to be a good binder for multiple HLA alleles."),
      h4(" Color scale for IC50 MT column:", style = "font-weight: bold"),
      h5(" lightgreen to darkgreen (0nM to 500nM); ", br(), "yellow to orange (500nM to 1000nM);", br(), " red (> 1000nM) "),
      h4(" Color scale for %ile MT column:", style = "font-weight: bold"),
      h5(" lightgreen to darkgreen (0-0.5%);", br(), " yellow to orange (0.5% to 2 %);", br(), " red (> 2%) "),
      h4(" Bar backgrounds:", style = "font-weight: bold"),
      h5(" RNA VAF and DNA VAF: Bar graphs range from 0 to 1", br(),
         " RNA Depth: Bar graph ranging from 0 to maximum value of RNA depth values across variants", br(),
         " RNA Expr: Bar graph ranging from 0 to 50 (this is meant to highlight variants with lower expression values for closer inspection)", br(),
         " Allele Expr: Bar graph ranging from 0 to (50 * maximum value of RNA VAF values across variants) "),
      h4(" Tier Types:", style = "font-weight: bold"),
      h5(" Variants are ordered by their Tiers in the following way: Pass, LowExpr, Anchor, Subclonal, Poor, NoExpr.
           Within the same tier, variants are ordered by the sum of their ranking in binding affinity and allele expression (i.e. lower binding
           affinity and higher allele expression is prioritized.)"),
      h5(" NoExpr: Mutant allele is not expressed ", br(),
         " LowExpr: Mutant allele has low expression (Allele Expr < allele expression threshold)", br(),
         " Subclonal: Likely not in the founding clone of the tumor (DNA VAF > max(DNA VAF)/2)", br(),
         " Anchor: Mutation is at an anchor residue in the shown peptide, and the WT allele has good binding (WT IC50 < binding threshold)", br(),
         " Poor: Fails two or more of the above criteria", br(),
         " Pass: Passes the above criteria, has strong MT binding (IC50 < 500) and strong expression (Allele Expr > allele expression threshold)"
      ),
    ))
  })
  ##update table upon selecting to investigate each individual row
  observeEvent(input$select_button, {
    if (is.null(df$mainTable)) {
      return()
    }
    df$selectedRow <- as.numeric(strsplit(input$select_button, "_")[[1]][2])
    session$sendCustomMessage("unbind-DT", "mainTable")
    df$mainTable$`Evaluation` <- shinyValue("selecter_", nrow(df$mainTable), df$mainTable)
    df$mainTable$`Eval` <- shinyInput(df$mainTable, selectInput, nrow(df$mainTable), "selecter_", choices = c("Pending", "Accept", "Reject", "Review"), width = "90px")
    dataTableProxy("mainTable") %>%
      selectPage((df$selectedRow - 1) %/% df$pageLength + 1)
  })
  ##selected row text box
  output$selected <- renderText({
    if (is.null(df$mainTable)) {
      return()
    }
    df$selectedRow
  })
  ##selected id update
  selectedID <- reactive({
    if (is.null(df$selectedRow)) {
      df$mainTable$ID[1]
    }else {
      df$mainTable$ID[df$selectedRow]
    }
  })
  output$selectedPeptide <- reactive({
    if (is.null(df$selectedRow)) {
      df$mainTable$`Best Peptide`[1]
    }else {
      df$mainTable$`Best Peptide`[df$selectedRow]
    }
  })
  output$selectedAAChange <- reactive({
    if (is.null(df$selectedRow)) {
      df$mainTable$`AA Change`[1]
    }else {
      df$mainTable$`AA Change`[df$selectedRow]
    }
  })
  output$selectedPos <- reactive({
    if (is.null(df$selectedRow)) {
      df$mainTable$`Pos`[1]
    }else {
      df$mainTable$`Pos`[df$selectedRow]
    }
  })
  output$selectedGene <- reactive({
    if (is.null(df$selectedRow)) {
      df$mainTable$`Gene`[1]
    }else {
      df$mainTable$`Gene`[df$selectedRow]
    }
  })
  ## Update comments section based on selected row
  observeEvent(input$comment, {
    if (is.null(df$mainTable)) {
      return()
    }
    df$comments[selectedID(), 1] <- input$comments
  })
  ##display of genomic information
  output$metricsTextGenomicCoord <- renderText({
    if (is.null(df$metricsData)) {
      return()
    }
    selectedID()
  })
  ##display of openCRAVAT link for variant
  output$url <- renderUI({
    if (is.null(df$mainTable)) {
      return()
    }
    id <- strsplit(selectedID(), "-")
    chromosome <- id[[1]][1]
    start <- id[[1]][2]
    stop <- id[[1]][3]
    ref <- id[[1]][4]
    alt <- id[[1]][5]
    url <- a("OpenCRAVAT variant report", href = paste("https://run.opencravat.org/webapps/variantreport/index.html?chrom=", chromosome, "&pos=", stop, "&ref_base=", ref, "&alt_base=", alt, sep = ""), target = "_blank")
    HTML(paste(url))
  })
  ##display of RNA VAF
  output$metricsTextRNA <- renderText({
    if (is.null(df$metricsData)) {
      return()
    }
    df$metricsData[[selectedID()]]$`RNA VAF`
  })
  ##display of DNA VAF
  output$metricsTextDNA <- renderText({
    if (is.null(df$metricsData)) {
      return()
    }
    df$metricsData[[selectedID()]]$`DNA VAF`
  })
  ##display of MT IC50 from additional data file
  output$addData_IC50 <- renderText({
    if (is.null(df$additionalData)) {
      return()
    }
    df$additionalData[df$additionalData$ID == selectedID(), ]$`IC50 MT`
  })
  ##display of MT percentile from additional data file
  output$addData_percentile <- renderText({
    df$additionalData[df$additionalData$ID == selectedID(), ]$`%ile MT`
  })
  ##display of Best Peptide from additional data file
  output$addData_peptide <- renderText({
    df$additionalData[df$additionalData$ID == selectedID(), ]$`Best Peptide`
  })
  ##display of Corresponding HLA allele from additional data file
  output$addData_allele <- renderText({
    df$additionalData[df$additionalData$ID == selectedID(), ]$`Allele`
  })
  ##display of Best Transcript from additional data file
  output$addData_transcript <- renderText({
    df$additionalData[df$additionalData$ID == selectedID(), ]$`Best Transcript`
  })
  ##transcripts sets table displaying sets of transcripts with the same consequence
  output$transcriptSetsTable <- renderDT({
    withProgress(message = "Loading Transcript Sets Table", value = 0, {
      GB_transcripts <- data.frame()
      best_transcript <- df$mainTable[df$mainTable$ID == selectedID(), ]$`Best Transcript`
      if (length(df$metricsData[[selectedID()]]$sets) != 0) {
        GB_transcripts <- data.frame(
          "Transcript Sets" = df$metricsData[[selectedID()]]$sets,
          "# Transcripts" = df$metricsData[[selectedID()]]$transcript_counts,
          "# Peptides" = df$metricsData[[selectedID()]]$peptide_counts,
          "Total Expr" = df$metricsData[[selectedID()]]$set_expr
        )
        names(GB_transcripts) <- c("Transcripts Sets", "#Transcripts", "# Peptides", "Total Expr")
        best_transcript_set <- NULL
        incProgress(0.5)
        for (i in 1:length(df$metricsData[[selectedID()]]$sets)){
          transcript_set <- df$metricsData[[selectedID()]]$good_binders[[df$metricsData[[selectedID()]]$sets[i]]]$`transcripts`
          transcript_set <- lapply(transcript_set, function(x) strsplit(x, "-")[[1]][1])
          if (best_transcript %in% transcript_set) {
            best_transcript_set <- df$metricsData[[selectedID()]]$sets[i]
            best_transcript_set_id <- i
          }
        }
        incProgress(0.5)
        datatable(GB_transcripts, selection = list(mode = "single", selected = best_transcript_set_id), style="bootstrap") %>%
          formatStyle("Transcripts Sets", backgroundColor = styleEqual(c(best_transcript_set), c("#98FF98")))
      }else {
        GB_transcripts <- data.frame("Transcript Sets" = character(), "# Transcripts" = character(), "# Peptides" = character(), "Total Expr" = character())
        names(GB_transcripts) <- c("Transcripts Sets", "#Transcripts", "# Peptides", "Total Expr")
        incProgress(0.5)
        datatable(GB_transcripts)
        incProgress(0.5)
      }
    })
  })
  ##update selected transcript set id
  selectedTranscriptSet <- reactive({
    selection <- input$transcriptSetsTable_rows_selected
    if (is.null(selection)) {
      selection <- 1
    }
    df$metricsData[[selectedID()]]$sets[selection]
  })
  
  ##transcripts table displaying transcript id and transcript expression values
  output$transcriptsTable <- renderDT({
    withProgress(message = "Loading Transcripts Table", value = 0, {
      GB_transcripts <- data.frame()
      best_transcript <- df$mainTable[df$mainTable$ID == selectedID(), ]$`Best Transcript`
      if (length(df$metricsData[[selectedID()]]$sets) != 0) {
        GB_transcripts <- data.frame("Transcripts" = df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$`transcripts`,
                                     "Expression" = df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$`transcript_expr`,
                                     "TSL" = df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$`tsl`,
                                     "Biotype" = df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$`biotype`,
                                     "Length" = df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$`transcript_length`)
        GB_transcripts$`Best Transcript` <- apply(GB_transcripts, 1, function(x) grepl(best_transcript, x["Transcripts"], fixed = TRUE))
        incProgress(0.5)
        names(GB_transcripts) <- c("Transcripts in Selected Set", "Expression", "Transcript Support Level", "Biotype", "Transcript Length (#AA)", "Best Transcript")
        incProgress(0.5)
        datatable(GB_transcripts, options = list(columnDefs = list(list(defaultContent = "N/A", targets = c(3)), list(visible = FALSE, targets = c(-1))))) %>%
          formatStyle(c("Transcripts in Selected Set"), "Best Transcript", backgroundColor = styleEqual(c(TRUE), c("#98FF98")))
      }else {
        GB_transcripts <- data.frame("Transcript" = character(), "Expression" = character(), "TSL" = character(), "Biotype" = character(), "Transcript Length (#AA)"= character(), "Length" = character())
        incProgress(0.5)
        names(GB_transcripts) <- c("Transcripts in Selected Set", "Expression", "Transcript Support Level", "Biotype", "Transcript Length (#AA)", "Best Transcript")
        incProgress(0.5)
        datatable(GB_transcripts)
      }
    })
  })
  
  ##display transcript expression
  output$metricsTextTranscript <- renderText({
    if (length(df$metricsData[[selectedID()]]$sets) != 0) {
      df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$`transcript_expr`
    }else {
      "N/A"
    }
  })
  ##display gene expression
  output$metricsTextGene <- renderText({
    if (length(df$metricsData[[selectedID()]]$sets) != 0) {
      df$metricsData[[selectedID()]]$`gene_expr`
    }else {
      "N/A"
    }
  })
  ##display peptide table with coloring
  output$peptideTable<- renderDT({
    withProgress(message = "Loading Peptide Table", value = 0, {
      if (length(df$metricsData[[selectedID()]]$sets) != 0 & !is.null(df$metricsData)) {
        peptide_data <- df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$`peptides`
        best_peptide <- df$mainTable[df$mainTable$ID == selectedID(), ]$`Best Peptide`
        peptide_names <- names(peptide_data)
        for (i in 1:length(peptide_names)) {
          peptide_data[[peptide_names[[i]]]]$individual_ic50_calls <- NULL
          peptide_data[[peptide_names[[i]]]]$individual_percentile_calls <- NULL
          peptide_data[[peptide_names[[i]]]]$individual_el_calls <- NULL
          peptide_data[[peptide_names[[i]]]]$individual_el_percentile_calls <- NULL
        }
        incProgress(0.5)
        peptide_data <- as.data.frame(peptide_data)
        incProgress(0.5)
        dtable <- datatable(do.call("rbind", lapply(peptide_names, table_formatting, peptide_data)), options = list(
          pageLength = 10,
          columnDefs = list(list(defaultContent = "X",
                                 targets = c(2:hla_count() + 1)),
                            list(orderable = TRUE, targets = 0),
                            list(visible = FALSE, targets = c(-1, -2))),
          rowCallback = JS("function(row, data, index, rowId) {",
                           "if(((rowId+1) % 4) == 3 || ((rowId+1) % 4) == 0) {",
                           'row.style.backgroundColor = "#E0E0E0";', "}", "}")
        ),
        selection = list(mode = "single", selected = "1"),
        style="bootstrap") %>%
          formatStyle("Type", fontWeight = styleEqual("MT", "bold")) %>%
          formatStyle(c("Peptide Sequence"), "Has ProbPos", border = styleEqual(c(TRUE), c("2px solid red"))) %>%
          formatStyle(c("Problematic Positions"), "Has ProbPos", border = styleEqual(c(TRUE), c("2px solid red"))) %>%
          formatStyle(c("Peptide Sequence"), "Has AnchorResidueFail", border = styleEqual(c(TRUE), c("2px solid red"))) %>%
          formatStyle(c("Anchor Residue Fail"), "Has AnchorResidueFail", border = styleEqual(c(TRUE), c("2px solid red"))) %>%
          formatStyle("Peptide Sequence", backgroundColor = styleEqual(c(best_peptide), c("#98FF98"))) %>%
          formatStyle("Peptide Sequence", fontFamily="monospace")
        dtable$x$data[[1]] <- as.numeric(dtable$x$data[[1]])
        dtable
      }else {
        incProgress(1)
        datatable(data.frame("Peptide Datatable" = character()), selection = list(mode = "single", selected = "1"), style="bootstrap")
      }
    })
  })
  ##update selected peptide data
  selectedPeptideData <- reactive({
    selection <- input$peptideTable_rows_selected
    if (is.null(selection)) {
      selection <- 1
    }
    peptide_names <- names(df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$`peptides`)
    index <- floor((as.numeric(selection) + 1) / 2)
    df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$peptides[[peptide_names[index]]]
  })
  ##Add legend for anchor heatmap
  output$peptideFigureLegend <- renderPlot({
    colors <- colorRampPalette(c("lightblue", "blue"))(99)[seq(1, 99, 7)]
    color_pos <- data.frame(d = as.character(seq(1, 99, 7)), x1 = seq(0.1, 1.5, 0.1), x2 = seq(0.2, 1.6, 0.1), y1 = rep(1, 15), y2 = rep(1.1, 15), colors = colors)
    color_label <- data.frame(x = c(0.1, 0.8, 1.6), y = rep(0.95, 3), score = c(0, 0.5, 1))
    p1 <- ggplot() +
      scale_y_continuous(limits = c(0.90, 1.2), name = "y") + scale_x_continuous(limits = c(0, 1.7), name = "x") +
      geom_rect(data = color_pos, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = colors), color = "black", alpha = 1) +
      scale_fill_identity()
    p1 <- p1 + geom_text(data = color_label, aes(x = x, y = y, label = score), size = 4, fontface = 2) +
      annotate(geom = "text", x = 0.5, y = 1.18, label = "Normalized Anchor Score", size = 4, fontface = 2) +
      coord_fixed() +
      theme_void() + theme(legend.position = "none", panel.border = element_blank(), plot.margin = margin(0, 0, 0, 0, "cm"))
    print(p1)
  })
  ##Anchor Heatmap overlayed on selected peptide sequences
  output$anchorPlot <- renderPlot({
    if (is.null(df$metricsData)) {
      return()
    }
    withProgress(message = "Loading Anchor Heatmap", value = 0, {
      if (type() == 2) {
        p1 <- ggplot() + annotate(geom = "text", x = 10, y = 20, label = "No data available for Class II HLA alleles", size = 6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p1)
      }else if (length(df$metricsData[[selectedID()]]$sets) != 0) {
        peptide_data <- df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$`peptides`
        peptide_names <- names(peptide_data)
        for (i in 1:length(peptide_names)) {
          peptide_data[[peptide_names[[i]]]]$individual_ic50_calls <- NULL
          peptide_data[[peptide_names[[i]]]]$individual_percentile_calls <- NULL
          peptide_data[[peptide_names[[i]]]]$individual_el_calls <- NULL
          peptide_data[[peptide_names[[i]]]]$individual_el_percentile_calls <- NULL
        }
        peptide_data <- as.data.frame(peptide_data)
        p1 <- ggplot() + scale_x_continuous(limits = c(0, 80)) + scale_y_continuous(limits = c(-31, 1))
        all_peptides <- list()
        incProgress(0.1)
        for (i in 1:length(peptide_names)) {
          #set & constrain mutation_pos' to not exceed length of peptide (may happen if mutation range goes off end)
          mutation_pos <- range_str_to_seq(df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$peptides[[peptide_names[i]]]$`mutation_position`)
          mt_peptide_length <- nchar(peptide_names[i])
          mutation_pos <- mutation_pos[mutation_pos <= mt_peptide_length]
          #set associated wt peptide to current mt peptide
          wt_peptide <- as.character(df$metricsData[[selectedID()]]$good_binders[[selectedTranscriptSet()]]$peptides[[peptide_names[i]]]$`wt_peptide`)
          #create dataframes for mt/wt pair
          df_mt_peptide <- data.frame("aa" = unlist(strsplit(peptide_names[i], "", fixed = TRUE)), "x_pos" = c(1:nchar(peptide_names[i])))
          df_mt_peptide$mutation <- "not_mutated"
          df_mt_peptide$type <- "mt"
          df_mt_peptide$y_pos <- (i * 2 - 1) * -1
          df_mt_peptide$length <- mt_peptide_length
          df_mt_peptide[mutation_pos, "mutation"] <- "mutated"
          df_wt_peptide <- data.frame("aa" = unlist(strsplit(wt_peptide, "", fixed = TRUE)), "x_pos" = c(1:nchar(wt_peptide)))
          df_wt_peptide$mutation <- "not_mutated"
          df_wt_peptide$type <- "wt"
          df_wt_peptide$y_pos <- (i * 2) * -1
          df_wt_peptide$length <- nchar(wt_peptide)
          all_peptides[[i]] <- rbind(df_mt_peptide, df_wt_peptide)
        }
        incProgress(0.4)
        all_peptides <- do.call(rbind, all_peptides)
        peptide_table <- do.call("rbind", lapply(peptide_names, table_formatting, peptide_data))
        peptide_table_filtered <- Filter(function(x) length(unique(x)) != 1, peptide_table)
        peptide_table_names <- names(peptide_table_filtered)
        hla_list <- peptide_table_names[grepl("^HLA-*", peptide_table_names)]
        hla_data <- data.frame(hla = hla_list)
        hla_sep <- max(nchar(peptide_table$`Peptide Sequence`))
        hla_data$y_pos <- 1
        hla_data$x_pos <- hla_sep / 2
        pad <- 3
        all_peptides_multiple_hla <- list()
        incProgress(0.1)
        for (i in 1:length(hla_list)) {
          hla_data$x_pos[i] <- hla_data$x_pos[i] + (hla_sep + pad) * (i - 1)
          omit_rows <- which(is.na(peptide_table_filtered[names(peptide_table_filtered) == hla_list[[i]]])) * -1
          all_peptides_multiple_hla[[i]] <- all_peptides[!(all_peptides$y_pos %in% omit_rows), ]
          all_peptides_multiple_hla[[i]]$color_value <- apply(all_peptides_multiple_hla[[i]], 1, function(x) peptide_coloring(hla_list[[i]], x))
          all_peptides_multiple_hla[[i]]$x_pos <- all_peptides_multiple_hla[[i]]$x_pos + (hla_sep + pad) * (i - 1)
        }
        incProgress(0.2)
        all_peptides_multiple_hla <- do.call(rbind, all_peptides_multiple_hla)
        h_line_pos <- data.frame(y_pos = seq(min(all_peptides_multiple_hla["y_pos"]) - 0.5, max(all_peptides_multiple_hla["y_pos"]) - 1.5, 2), x_pos = c(min(all_peptides_multiple_hla["x_pos"]) - 1))
        h_line_pos <- rbind(h_line_pos, data.frame(x_pos = max(all_peptides_multiple_hla["x_pos"]) + 1, y_pos = seq(min(all_peptides_multiple_hla["y_pos"]) - 0.5, max(all_peptides_multiple_hla["y_pos"]) - 1.5, 2)))
        incProgress(0.2)
        p1 <- p1 +
          geom_rect(data = all_peptides_multiple_hla, aes(xmin = x_pos - 0.5, xmax = 1 + x_pos - 0.5, ymin = .5 + y_pos, ymax = -.5 + y_pos), fill = all_peptides_multiple_hla$color_value) +
          geom_text(data = all_peptides_multiple_hla, aes(x = x_pos, y = y_pos, label = aa, color = mutation), size = 4) +
          geom_text(data = hla_data, aes(x = x_pos, y = y_pos, label = hla), size = 4, fontface = "bold") +
          geom_line(data = h_line_pos, (aes(x = x_pos, y = y_pos, group = y_pos)), linetype = "dashed")
        p1 <- p1 + scale_color_manual("mutation", values = c("not_mutated" = "#000000", "mutated" = "#e74c3c"))
        p1 <- p1 + theme_void() + theme(legend.position = "none", panel.border = element_blank())
        print(p1)
      }else {
        p1 <- ggplot() + annotate(geom = "text", x = 10, y = 20, label = "No data available", size = 6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p1)
      }
    })
  }, height = 400, width = 800)
  #anchor score tables for each HLA allele
  output$anchorWeights<- renderDT({
    withProgress(message = "Loading Anchor Weights Table", value = 0, {
      weights <- anchor_weights_for_alleles(df$metricsData$alleles)
      dtable <- datatable(weights, options = list(
        pageLength = 10,
        lengthChange = FALSE,
        rowCallback = JS("function(row, data, index, rowId) {",
                         "if(((rowId+1) % 4) == 3 || ((rowId+1) % 4) == 0) {",
                         'row.style.backgroundColor = "#E0E0E0";', "}", "}")
      ))
      dtable
    })
  })
  ##updating IC50 binding score for selected peptide pair
  bindingScoreDataIC50 <- reactive({
    if (is.null(df$metricsData)) {
      return()
    }
    if (length(df$metricsData[[selectedID()]]$sets) != 0) {
      algorithm_names <- data.frame(algorithms = selectedPeptideData()$individual_ic50_calls$algorithms)
      wt_data <- as.data.frame(selectedPeptideData()$individual_ic50_calls$WT, check.names = FALSE)
      colnames(wt_data) <- paste(colnames(wt_data), "_WT_Score", sep = "")
      mt_data <- as.data.frame(selectedPeptideData()$individual_ic50_calls$MT, check.names = FALSE)
      colnames(mt_data) <- paste(colnames(mt_data), "_MT_Score", sep = "")
      full_data <- cbind(algorithm_names, mt_data, wt_data) %>%
        gather("col", "val", colnames(mt_data)[1]:tail(colnames(wt_data), n = 1)) %>%
        separate(col, c("HLA_allele", "Mutant", "Score"), sep = "\\_") %>%
        spread("Score", val)
      full_data
    }else {
      return()
    }
  })
  ##plotting IC5 binding score violin plot
  output$bindingData_IC50 <- renderPlot({
    withProgress(message = "Loading Binding Score Plot (IC50)", value = 0, {
      if (length(df$metricsData[[selectedID()]]$sets) != 0) {
        line.data <- data.frame(yintercept = c(500, 1000), Cutoffs = c("500nM", "1000nM"), color = c("#28B463", "#EC7063"))
        hla_allele_count <- length(unique(bindingScoreDataIC50()$HLA_allele))
        incProgress(0.5)
        p <- ggplot(data = bindingScoreDataIC50(), aes(x = Mutant, y = Score, color = Mutant), trim = FALSE) + geom_violin() + facet_grid(cols = vars(HLA_allele)) + scale_y_continuous(trans = "log10") + #coord_trans(y = "log10") +
          stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)) +
          geom_jitter(data = bindingScoreDataIC50(), aes(shape = algorithms), sizes = 5, stroke = 1, position = position_jitter(0.3)) +
          scale_shape_manual(values = 0:8) +
          geom_hline(aes(yintercept = yintercept, linetype = Cutoffs), line.data, color = rep(line.data$color, hla_allele_count)) +
          scale_color_manual(values = rep(c("MT" = "#D2B4DE", "WT" = "#F7DC6F"), hla_allele_count)) +
          theme(strip.text = element_text(size = 15), axis.text = element_text(size = 10), axis.title = element_text(size = 15), axis.ticks = element_line(size = 3), legend.text = element_text(size = 15), legend.title = element_text(size = 15))
        incProgress(0.5)
        print(p)
      }else {
        p <- ggplot() + annotate(geom = "text", x = 10, y = 20, label = "No data available", size = 6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p)
      }
    })
  })
  ##updating percentile binding score for selected peptide pair
  bindingScoreDataPercentile <- reactive({
    if (length(df$metricsData[[selectedID()]]$sets) != 0) {
      algorithm_names <- data.frame(algorithms = selectedPeptideData()$individual_percentile_calls$algorithms)
      wt_data <- as.data.frame(selectedPeptideData()$individual_percentile_calls$WT, check.names = FALSE)
      colnames(wt_data) <- paste(colnames(wt_data), "_WT_Score", sep = "")
      mt_data <- as.data.frame(selectedPeptideData()$individual_percentile_calls$MT, check.names = FALSE)
      colnames(mt_data) <- paste(colnames(mt_data), "_MT_Score", sep = "")
      full_data <- cbind(algorithm_names, mt_data, wt_data) %>%
        gather("col", "val", colnames(mt_data)[1]:tail(colnames(wt_data), n = 1)) %>%
        separate(col, c("HLA_allele", "Mutant", "Score"), sep = "\\_") %>%
        spread("Score", val)
      full_data
    }else {
      return()
    }
  })
  ##plotting percentile binding score violin plot
  output$bindingData_percentile <- renderPlot({
    withProgress(message = "Loading Binding Score Plot (Percentile)", value = 0, {
      if (length(df$metricsData[[selectedID()]]$sets) != 0) {
        line.data <- data.frame(yintercept = c(0.5, 2), Cutoffs = c("0.5%", "2%"), color = c("#28B463", "#EC7063"))
        hla_allele_count <- length(unique(bindingScoreDataPercentile()$HLA_allele))
        incProgress(0.5)
        p <- ggplot(data = bindingScoreDataPercentile(), aes(x = Mutant, y = Score, color = Mutant), trim = FALSE) + geom_violin() + facet_grid(cols = vars(HLA_allele)) + scale_y_continuous(trans = "log10") + #coord_trans(y = "log10") +
          stat_summary(fun.y = mean, fun.ymin = mean, fun.ymax = mean, geom = "crossbar", width = 0.25, position = position_dodge(width = .25)) +
          geom_jitter(data = bindingScoreDataPercentile(), aes(shape = algorithms), size = 5, stroke = 1, position = position_jitter(0.3)) +
          scale_shape_manual(values = 0:8) +
          geom_hline(aes(yintercept = yintercept, linetype = Cutoffs), line.data, color = rep(line.data$color, hla_allele_count)) +
          scale_color_manual(values = rep(c("MT" = "#D2B4DE", "WT" = "#F7DC6F"), hla_allele_count)) +
          theme(strip.text = element_text(size = 15), axis.text = element_text(size = 10), axis.title = element_text(size = 15), axis.ticks = element_line(size = 3), legend.text = element_text(size = 15), legend.title = element_text(size = 15))
        incProgress(0.5)
        print(p)
      }else {
        p <- ggplot() + annotate(geom = "text", x = 10, y = 20, label = "No data available", size = 6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p)
      }
    })
  })
  ##plotting binding data table with IC50 and percentile values
  output$bindingDatatable <- renderDT({
    withProgress(message = "Loading binding datatable", value = 0, {
      if (length(df$metricsData[[selectedID()]]$sets) != 0) {
        binding_data <- bindingScoreDataIC50()
        names(binding_data)[names(binding_data) == "Score"] <- "IC50 Score"
        binding_data["% Score"] <- bindingScoreDataPercentile()["Score"]
        binding_data["Score"] <- paste(round(as.numeric(binding_data$`IC50 Score`), 2), " (%: ", round(as.numeric(binding_data$`% Score`), 2), ")", sep = "")
        binding_data["IC50 Score"] <- NULL
        binding_data["% Score"] <- NULL
        binding_reformat <- dcast(binding_data, HLA_allele + Mutant ~ algorithms, value.var = "Score")
        incProgress(1)
        dtable <- datatable(binding_reformat, options = list(
          pageLength = 10,
          lengthChange = FALSE,
          rowCallback = JS("function(row, data, index, rowId) {",
                           "if(((rowId+1) % 4) == 3 || ((rowId+1) % 4) == 0) {",
                           'row.style.backgroundColor = "#E0E0E0";', "}", "}")
        )) %>% formatStyle("Mutant", fontWeight = styleEqual("MT", "bold"), color = styleEqual("MT", "#E74C3C"))
        dtable
      }else {
        incProgress(1)
        datatable(data.frame("Binding Predictions Datatable" = character()))
      }
    })
  })
  ##updating elution score for selected peptide pair
  elutionScoreData <- reactive({
    if (is.null(df$metricsData)) {
      return()
    }
    if (length(df$metricsData[[selectedID()]]$sets) != 0 && length(selectedPeptideData()$individual_el_calls$algorithms) > 0) {
      algorithm_names <- data.frame(algorithms = selectedPeptideData()$individual_el_calls$algorithms)
      wt_data <- as.data.frame(selectedPeptideData()$individual_el_calls$WT, check.names = FALSE)
      colnames(wt_data) <- paste(colnames(wt_data), "_WT_Score", sep = "")
      mt_data <- as.data.frame(selectedPeptideData()$individual_el_calls$MT, check.names = FALSE)
      colnames(mt_data) <- paste(colnames(mt_data), "_MT_Score", sep = "")
      full_data <- cbind(algorithm_names, mt_data, wt_data) %>%
        gather("col", "val", colnames(mt_data)[1]:tail(colnames(wt_data), n = 1)) %>%
        separate(col, c("HLA_allele", "Mutant", "Score"), sep = "\\_") %>%
        spread("Score", val)
      full_data
    }else {
      return()
    }
  })
  ##updating elution percentile for selected peptide pair
  elutionPercentileData <- reactive({
    if (length(df$metricsData[[selectedID()]]$sets) != 0) {
      algorithm_names <- data.frame(algorithms = selectedPeptideData()$individual_el_percentile_calls$algorithms)
      wt_data <- as.data.frame(selectedPeptideData()$individual_el_percentile_calls$WT, check.names = FALSE)
      colnames(wt_data) <- paste(colnames(wt_data), "_WT_Score", sep = "")
      mt_data <- as.data.frame(selectedPeptideData()$individual_el_percentile_calls$MT, check.names = FALSE)
      colnames(mt_data) <- paste(colnames(mt_data), "_MT_Score", sep = "")
      full_data <- cbind(algorithm_names, mt_data, wt_data) %>%
        gather("col", "val", colnames(mt_data)[1]:tail(colnames(wt_data), n = 1)) %>%
        separate(col, c("HLA_allele", "Mutant", "Score"), sep = "\\_") %>%
        spread("Score", val)
      full_data
    }else {
      return()
    }
  })
  ##plotting elution data table
  output$elutionDatatable <- renderDT({
    withProgress(message = "Loading elution datatable", value = 0, {
      if (length(df$metricsData[[selectedID()]]$sets) != 0) {
        elution_data <- elutionScoreData()
        if (!is.null(elution_data)) {
          names(elution_data)[names(elution_data) == "Score"] <- "Elution Score"
          elution_data["% Score"] <- elutionPercentileData()["Score"]
          elution_data["Score"] <- paste(round(as.numeric(elution_data$`Elution Score`), 2), " (%: ", round(as.numeric(elution_data$`% Score`), 2), ")", sep = "")
          elution_data["Elution Score"] <- NULL
          elution_data["% Score"] <- NULL
          elution_reformat <- dcast(elution_data, HLA_allele + Mutant ~ algorithms, value.var = "Score")
          incProgress(1)
          dtable <- datatable(elution_reformat, options = list(
            pageLength = 10,
            lengthChange = FALSE,
            rowCallback = JS("function(row, data, index, rowId) {",
                             "if(((rowId+1) % 4) == 3 || ((rowId+1) % 4) == 0) {",
                             'row.style.backgroundColor = "#E0E0E0";', "}", "}")
          )) %>% formatStyle("Mutant", fontWeight = styleEqual("MT", "bold"), color = styleEqual("MT", "#E74C3C"))
          dtable
        }else {
          incProgress(1)
          datatable(data.frame("Elution Datatable" = character()))
        }
      }else {
        incProgress(1)
        datatable(data.frame("Elution Datatable" = character()))
      }
    })
  })
  
  ##updating reference matches for selected peptide
  output$hasReferenceMatchData <- reactive({
    if (is.null(df$metricsData[[selectedID()]]$reference_matches)) {
      "Reference Similarity not run"
    } else {
      ""
    }
  })
  referenceMatchData <- reactive({
    if (length(df$metricsData[[selectedID()]]$reference_matches$matches) != 0) {
      as.data.frame(df$metricsData[[selectedID()]]$reference_matches$matches, check.names = False)
    }else {
      return()
    }
  })
  output$referenceMatchHitCount <- reactive({
    if (is.null(df$metricsData[[selectedID()]]$reference_matches)) {
      "N/A"
    } else {
      df$metricsData[[selectedID()]]$reference_matches$count
    }
  })
  output$referenceMatchQuerySequence <- reactive({
    if (is.null(df$metricsData[[selectedID()]]$reference_matches)) {
      "N/A"
    } else {
      df$metricsData[[selectedID()]]$reference_matches$query_peptide
    }
  })
  output$referenceMatchDatatable <- renderDT({
    withProgress(message = "Loading reference match datatable", value = 0, {
      reference_match_data <- referenceMatchData()
      if (!is.null(reference_match_data)) {
        incProgress(1)
        dtable <- datatable(reference_match_data, options = list(
          pageLength = 10,
          lengthMenu = c(10)
        ),
        style="bootstrap") %>%
          formatStyle("Matched Peptide", fontFamily="monospace")
        dtable
      } else {
        incProgress(1)
        datatable(data.frame("Reference Matches Datatable" = character()))
      }
    })
  })
  ##Best Peptide with mutated positions marked
  output$referenceMatchPlot <- renderPlot({
    withProgress(message = "Loading Reference Match Best Peptide Plot", value = 0, {
      selectedPosition <- if (is.null(df$selectedRow)) {
        df$mainTable$`Pos`[1]
      }else {
        df$mainTable$`Pos`[df$selectedRow]
      }
      selectedPeptide <- if (is.null(df$selectedRow)) {
        df$mainTable$`Best Peptide`[1]
      }else {
        df$mainTable$`Best Peptide`[df$selectedRow]
      }
      #set & constrain mutation_pos' to not exceed length of peptide (may happen if mutation range goes off end)
      mutation_pos <- range_str_to_seq(selectedPosition)
      peptide_length <- nchar(selectedPeptide)
      mutation_pos <- mutation_pos[mutation_pos <= peptide_length]
      #create dataframes
      df_peptide <- data.frame("aa" = unlist(strsplit(selectedPeptide, "", fixed = TRUE)), "x_pos" = c(1:nchar(selectedPeptide)))
      df_peptide$mutation <- "not_mutated"
      df_peptide$type <- "mt"
      df_peptide$y_pos <- 1.05
      df_peptide$x_pos <- seq(0.05, peptide_length*0.1+0.05, length.out=peptide_length)
      df_peptide$length <- peptide_length
      df_peptide[mutation_pos, "mutation"] <- "mutated"
      ref_match_colors <- rep("white", peptide_length)
      x1_bins = seq(0, peptide_length*0.1, length.out=peptide_length)
      x2_start = peptide_length*0.1/(peptide_length-1)
      x2_bins = seq(x2_start, peptide_length*0.1+x2_start, length.out=peptide_length)
      ref_match_color_pos <- data.frame(d = df_peptide, x1 = x1_bins, x2 = x2_bins, y1 = rep(1, peptide_length), y2 = rep(1.1, peptide_length), colors = ref_match_colors)
      p2 <- ggplot() +
        geom_rect(data = ref_match_color_pos, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = colors), color = "black", alpha = 1) +
        geom_text(data = df_peptide, aes(x = x_pos, y = y_pos, label = aa, color = mutation), size = 5) +
        scale_fill_identity() +
        coord_fixed() +
        theme_void() + theme(legend.position = "none", panel.border = element_blank(), plot.margin = margin(0, 0, 0, 0, "pt"))
      p2 <- p2 + scale_color_manual("mutation", values = c("not_mutated" = "#000000", "mutated" = "#e74c3c"))
      print(p2)
    })
  }, height = 20, width = function(){
    selectedPeptide <- if (is.null(df$selectedRow)) {
      df$mainTable$`Best Peptide`[1]
    }else {
      df$mainTable$`Best Peptide`[df$selectedRow]
    }
    nchar(selectedPeptide) * 20
  } )
  ##Best Peptide with best peptide highlighted and mutated positions marked
  output$referenceMatchQueryPlot <- renderPlot({
    withProgress(message = "Loading Reference Match Query Peptide Plot", value = 0, {
      selectedPosition <- if (is.null(df$selectedRow)) {
        df$mainTable$`Pos`[1]
      }else {
        df$mainTable$`Pos`[df$selectedRow]
      }
      selectedPeptide <- if (is.null(df$selectedRow)) {
        df$mainTable$`Best Peptide`[1]
      }else {
        df$mainTable$`Best Peptide`[df$selectedRow]
      }
      mutation_pos <- range_str_to_seq(selectedPosition)
      #remove leading amino acids from the selectedPeptide that don't occur in
      #the query peptide
      if (mutation_pos[1] > 8) {
        offset <- mutation_pos[1] - 8
        selectedPeptide <- substr(selectedPeptide, offset + 1, nchar(selectedPeptide))
        mutation_pos <- mutation_pos - offset
      }
      if (!is.null(df$metricsData[[selectedID()]]$reference_matches)) {
        queryPeptide <- df$metricsData[[selectedID()]]$reference_matches$query_peptide
        peptide_length <- nchar(queryPeptide)
        ref_match_colors <- rep("white", peptide_length)
        #set & constrain mutation_pos' to not exceed length of peptide (may happen if mutation range goes off end)
        bestPeptidePos <- str_locate(queryPeptide, fixed(selectedPeptide))
        #if the selectedPeptide is not found in the queryPeptide there are
        #trailing amino acids that don't occur in the queryPeptide - remove them
        while (is.na(bestPeptidePos[, 1])) {
          selectedPeptide <- substr(selectedPeptide, 1, nchar(selectedPeptide) - 1)
          bestPeptidePos <- str_locate(queryPeptide, fixed(selectedPeptide))
        }
        best_peptide_positions <- seq(bestPeptidePos[, 1], bestPeptidePos[, 2])
        mutation_pos <- mutation_pos + bestPeptidePos[, 1] - 1
        mutation_pos <- mutation_pos[mutation_pos <= peptide_length]
        ref_match_colors[best_peptide_positions] <- "yellow"
        #create dataframes
        df_peptide <- data.frame("aa" = unlist(strsplit(queryPeptide, "", fixed = TRUE)), "x_pos" = c(1:nchar(queryPeptide)))
        df_peptide$mutation <- "not_mutated"
        df_peptide$type <- "mt"
        df_peptide$y_pos <- 1.05
        df_peptide$x_pos <- seq(0.05, peptide_length*0.1+0.05, length.out=peptide_length)
        df_peptide$length <- peptide_length
        df_peptide[mutation_pos, "mutation"] <- "mutated"
        x1_bins = seq(0, peptide_length*0.1, length.out=peptide_length)
        x2_start = peptide_length*0.1/(peptide_length-1)
        x2_bins = seq(x2_start, peptide_length*0.1+x2_start, length.out=peptide_length)
        ref_match_color_pos <- data.frame(d = df_peptide, x1 = x1_bins, x2 = x2_bins, y1 = rep(1, peptide_length), y2 = rep(1.1, peptide_length), colors = ref_match_colors)
        p3 <- ggplot() +
          geom_rect(data = ref_match_color_pos, aes(xmin = x1, xmax = x2, ymin = y1, ymax = y2, fill = colors), color = "black", alpha = 1) +
          geom_text(data = df_peptide, aes(x = x_pos, y = y_pos, label = aa, color = mutation), size = 5) +
          scale_fill_identity() +
          coord_fixed() +
          theme_void() + theme(legend.position = "none", panel.border = element_blank(), plot.margin = margin(0, 0, 0, 0, "pt"))
        p3 <- p3 + scale_color_manual("mutation", values = c("not_mutated" = "#000000", "mutated" = "#e74c3c"))
        incProgress(1)
        print(p3)
      } else {
        p3 <- ggplot() + annotate(geom = "text", x = 0, y = 0, label = "N/A", size = 5) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p3)
      }
    })
  }, height = 20, width = function(){
    if (is.null(df$metricsData[[selectedID()]]$reference_matches$query_peptide)) {
      return(40)
    } else {
      return(nchar(df$metricsData[[selectedID()]]$reference_matches$query_peptide) * 20)
    }
  } )
  ##############################EXPORT TAB##############################################
  #evalutation overview table
  output$checked <- renderTable({
    if (is.null(df$mainTable)) {
      return()
    }
    Evaluation <- data.frame(selected = shinyValue("selecter_", nrow(df$mainTable), df$mainTable))
    data <- as.data.frame(table(Evaluation))
    data$Count <- data$Freq
    data$Freq <- NULL
    data
  })
  #export table display with options to download
  output$ExportTable <- renderDataTable({
    if (is.null(df$mainTable)) {
      return()
    }
    colsToDrop <- colnames(df$mainTable) %in% c("Evaluation", "Eval", "Select", "Scaled BA", "Scaled percentile", "Tier Count", "Bad TSL",
                                                "Comments", "Gene of Interest", "Bad TSL", "Col RNA Expr", "Col RNA VAF", "Col Allele Expr",
                                                "Col RNA Depth", "Col DNA VAF", "Percentile Fail", "Has Prob Pos")
    data <- df$mainTable[, !(colsToDrop)]
    col_names <- colnames(data)
    data <- data.frame(data, Evaluation = shinyValue("selecter_", nrow(df$mainTable), df$mainTable))
    colnames(data) <- c(col_names, "Evaluation")
    comments <- data.frame("ID" = row.names(df$comments), Comments = df$comments[, 1])
    data <- join(data, comments)
    data[is.na(data)] <- "NA"
    data
  }, escape = FALSE, server = FALSE, rownames = FALSE,
  options = list(dom = "Bfrtip",
                 buttons = list(
                   list(extend = "csvHtml5",
                        filename = input$exportFileName,
                        fieldSeparator = "\t",
                        text = "Download as TSV",
                        extension = ".tsv"),
                   list(extend = "excel",
                        filename = input$exportFileName,
                        text = "Download as excel")
                 ),
                 initComplete = htmlwidgets::JS(
                   "function(settings, json) {",
                   paste0("$(this.api().table().header()).css({'font-size': '", "8pt", "'});"),
                   "}")
  ),
  selection = "none",
  extensions = c("Buttons"))
  


  ### Other Modules ############################################################
  
  
  ############### NeoFox Tab ##########################
  df_neofox <- reactiveValues(
    mainTable_neofox = NULL
  )
  
  # Option 1: User uploaded
  observeEvent(input$neofox_data$datapath, {
    #session$sendCustomMessage("unbind-DT", "mainTable_neofox")
    mainData_neofox <- read.table(input$neofox_data$datapath, sep = "\t",  header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(mainData_neofox) <- mainData_neofox[1, ]
    mainData_neofox <- mainData_neofox[-1, ]
    row.names(mainData_neofox) <- NULL
    
    # Columns that have been reviewed as most interesting
    columns_to_star <- c(
      "dnaVariantAlleleFrequency", "rnaExpression", "imputedGeneExpression",
      "rnaVariantAlleleFrequency", "NetMHCpan_bestRank_rank", "NetMHCpan_bestAffinity_affinity",
      "NetMHCpan_bestAffinity_affinityWT", "NetMHCpan_bestRank_rankWT", "PHBR_I",
      "NetMHCIIpan_bestRank_rank", "NetMHCIIpan_bestRank_rankWT", "PHBR_II", "Amplitude_MHCI_bestAffinity",
      "Pathogensimiliarity_MHCI_bestAffinity9mer", "DAI_MHCI_bestAffinity", "Tcell_predictor",
      "Selfsimilarity_MHCI", "Selfsimilarity_MHCII", "IEDB_Immunogenicity_MHCI", "IEDB_Immunogenicity_MHCII",
      "MixMHCpred_bestScore_score", "MixMHCpred_bestScore_rank", "MixMHC2pred_bestRank_peptide",
      "MixMHC2pred_bestRank_rank", "Dissimilarity_MHCI", "Dissimilarity_MHCII", "Vaxrank_bindingScore",
      "PRIME_bestScore_rank", "PRIME_bestScore_score"
    )
    
    # Check if each column is present in the dataframe and modify the names
    for (col_name in columns_to_star) {
      if (col_name %in% names(mainData_neofox)) {
        new_col_name <- paste0("*", col_name)
        names(mainData_neofox)[names(mainData_neofox) == col_name] <- new_col_name
      }
    }   
    
    df_neofox$mainTable_neofox <- mainData_neofox
  })
  
  # Option 2: Demo Data
  observeEvent(input$loadDefaultneofox, {
    #session$sendCustomMessage("unbind-DT", "mainTable_neofox")
    data_neofox <- "data/neofox_neoantigen_candidates_annotated.tsv"
    mainData_neofox <- read.table(data_neofox, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(mainData_neofox) <- mainData_neofox[1, ]
    mainData_neofox <- mainData_neofox[-1, ]
    row.names(mainData_neofox) <- NULL
    
    # Columns that have been reviewed as most interesting
    columns_to_star <- c(
      "dnaVariantAlleleFrequency", "rnaExpression", "imputedGeneExpression",
      "rnaVariantAlleleFrequency", "NetMHCpan_bestRank_rank", "NetMHCpan_bestAffinity_affinity",
      "NetMHCpan_bestAffinity_affinityWT", "NetMHCpan_bestRank_rankWT", "PHBR_I",
      "NetMHCIIpan_bestRank_rank", "NetMHCIIpan_bestRank_rankWT", "PHBR_II", "Amplitude_MHCI_bestAffinity",
      "Pathogensimiliarity_MHCI_bestAffinity9mer", "DAI_MHCI_bestAffinity", "Tcell_predictor",
      "Selfsimilarity_MHCI", "Selfsimilarity_MHCII", "IEDB_Immunogenicity_MHCI", "IEDB_Immunogenicity_MHCII",
      "MixMHCpred_bestScore_score", "MixMHCpred_bestScore_rank", "MixMHC2pred_bestRank_peptide",
      "MixMHC2pred_bestRank_rank", "Dissimilarity_MHCI", "Dissimilarity_MHCII", "Vaxrank_bindingScore",
      "PRIME_bestScore_rank", "PRIME_bestScore_score"
    )
    
    # Check if each column is present in the dataframe and modify the names
    for (col_name in columns_to_star) {
      if (col_name %in% names(mainData_neofox)) {
        new_col_name <- paste0("*", col_name)
        names(mainData_neofox)[names(mainData_neofox) == col_name] <- new_col_name
      }
    }    
    
    df_neofox$mainTable_neofox <- mainData_neofox
    updateTabItems(session, "neofox_tabs", "neofox_explore")
  })
  
  ##Clear file inputs if demo data load button is clicked
  output$neofox_upload_ui <- renderUI({
    fileInput(inputId = "neofox_data", label = "NeoFox output table (tsv required)",
              accept =  c("text/tsv", "text/tab-separated-values,text/plain", ".tsv"))
  })
  
  observeEvent(input$visualize_neofox, {
    updateTabItems(session, "neofox_tabs", "neofox_explore")
  })
  

  ## NeoFox Explore
  output$neofoxTable <- DT::renderDataTable(
    if (is.null(df_neofox$mainTable_neofox)) {
      return(datatable(data.frame("Annotated Table" = character())))
      
    } else {
      return(datatable(df_neofox$mainTable_neofox,
                escape = FALSE,
                selection = "multiple",
                extensions = c("Buttons")
      ))
    })
  
  output$neofox_selected <- renderText({
    if (is.null(df_neofox$mainTable_neofox)) {
      return()
    }
    input$neofoxTable_rows_selected
  })
  
  ### NeoFox Violin Plots
  
  ## Drop down to select what features to show violin plots for
  output$noefox_features_ui <- renderUI({
    df <- df_neofox$mainTable_neofox
    df <- type.convert(df, as.is = TRUE)
    df[is.na(df)] <- 0

    features <- names(df)[sapply(df, is.numeric)]
    sorted_features <- features[order(!grepl("^\\*", features))]


    default_selection <- c("*IEDB_Immunogenicity_MHCI", "*IEDB_Immunogenicity_MHCII", "*PHBR_I",
      "*MixMHCpred_bestScore_score", "*MixMHCpred_bestScore_rank", "*MixMHC2pred_bestRank_peptide")

    pickerInput(inputId = "neofox_features",
                label = "Plots to Display",
                choices = sorted_features,
                selected = default_selection,
                options = list(`live-search` = TRUE, "max-options" = 6),
                multiple = TRUE
                )
  })

  output$neofox_violin_plots_row1 <- renderPlot({
    withProgress(message = "Loading Violin Plots", value = 0, {
      if (length(input$neofoxTable_rows_selected) != 0 & length(input$neofox_features) != 0) {

        plot_cols_neofox <- c("mutatedXmer", input$neofox_features)
        plot_data_neofox <- df_neofox$mainTable_neofox[, plot_cols_neofox]
        plot_data_neofox <- type.convert(plot_data_neofox, as.is = TRUE)
        plot_data_neofox[is.na(plot_data_neofox)] <- 0

        plot_data_neofox$Selected <- "No"
        plot_data_neofox[input$neofoxTable_rows_selected, "Selected"] <- "Yes"
        reformat_data_neofox <- plot_data_neofox %>%
          gather("Feature", "Value", -c("mutatedXmer", "Selected"))


        p_neofox <- ggplot(reformat_data_neofox, aes(x = "", y = Value)) + geom_violin() +
           geom_jitter(data = reformat_data_neofox[reformat_data_neofox["Selected"] == "No", ], aes(color = Selected), size = 1, alpha = 0.5, stroke = 1, position = position_jitter(0.3)) +
           geom_jitter(data = reformat_data_neofox[reformat_data_neofox["Selected"] == "Yes", ], aes(color = Selected), size = 2, alpha = 1, stroke = 1, position = position_jitter(0.3)) +
           scale_color_manual(values = c("No" = "#939094", "Yes" = "#f42409")) +
          labs(x = NULL) +
          facet_wrap(~Feature, scales="free", ncol=6) +
          theme(strip.text = element_text(size = 15), axis.text = element_text(size = 10), axis.title = element_text(size = 10), axis.ticks = element_line(size = 3), legend.text = element_text(size = 10), legend.title = element_text(size = 10))

        incProgress(0.5)
        print(p_neofox)
      }else {
        p_neofox <- ggplot() + annotate(geom = "text", x = 10, y = 20, label = "No data available", size = 6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
        print(p_neofox)
      }
    })
  })


  ## NeoFox Dynamic scatter plot
  output$xvrbl <- renderUI({
    df <- df_neofox$mainTable_neofox
    df <- type.convert(df, as.is = TRUE)
    df[is.na(df)] <- 0

    features <- names(df)[sapply(df, is.numeric)]
    sorted_features <- features[order(!grepl("^\\*", features))]

    default_selection <- "*NetMHCpan_bestAffinity_affinity"

    pickerInput(inputId = "xvrbl",
                label = "X-Axis Variable",
                choices = sorted_features,
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE
    )
  })

  output$xvrbl_log <- renderUI({
    radioButtons(
      inputId = "LogX",
      choices = c("none", "ln", "log2", "log10", "sqrt"),
      label = "Transform",
      inline = TRUE
      )
  })

  output$yvrbl <- renderUI({
    df <- df_neofox$mainTable_neofox
    df <- type.convert(df, as.is = TRUE)
    df[is.na(df)] <- 0

    features <- names(df)[sapply(df, is.numeric)]
    sorted_features <- features[order(!grepl("^\\*", features))]
    default_selection <- "*NetMHCpan_bestAffinity_affinityWT"

    pickerInput(inputId = "yvrbl",
                label = "Y-Axis Variable",
                choices = sorted_features,
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE
    )
  })

  output$yvrbl_log <- renderUI({
    radioButtons(
      inputId = "LogY",
      choices = c("none", "ln", "log2", "log10", "sqrt"),
      label = "Transform",
      inline = TRUE
    )
  })

  output$xvrbl_scale <- renderUI({
    withProgress(message = "Loading Scale", value = 0, {
      req(input$xvrbl, input$LogX)  # Use req() to check if inputs are not NULL
          df <- df_neofox$mainTable_neofox
          df <- type.convert(df, as.is = TRUE)
          df[is.na(df)] <- 0
          df <- df[is.finite(df[[input$xvrbl]]),]

          # Apply log or sqrt transformation
          if (input$LogX == "ln") {
            df[[input$xvrbl]] <- log(ifelse(df[[input$xvrbl]] == 0, 1e-10, df[[input$xvrbl]]))
          } else if (input$LogX == "log2") {
            df[[input$xvrbl]] <- log2(ifelse(df[[input$xvrbl]] == 0, 1e-10, df[[input$xvrbl]]))
          } else if (input$LogX == "log10") {
            df[[input$xvrbl]] <- log10(ifelse(df[[input$xvrbl]] == 0, 1e-10, df[[input$xvrbl]]))
          } else if (input$LogX == "sqrt") {
            df[[input$xvrbl]] <- sqrt(ifelse(df[[input$xvrbl]] < 0, 1e-10, df[[input$xvrbl]]))
          } else {
            df[[input$xvrbl]] <- df[[input$xvrbl]]
          }

          df <- df[is.finite(df[[input$xvrbl]]),]

          xvrbl_values <- df[[input$xvrbl]]
          range_values <- range(as.numeric(xvrbl_values), na.rm = TRUE)
          min_value <- as.numeric(format(round(range_values[1], 2), nsmall = 2))
          max_value <- as.numeric(format(round(range_values[2], 2), nsmall = 2))


          # Check if min_value and max_value are equal, set default values
          if (min_value == max_value) {
            min_value <- min_value - 1
            max_value <- max_value + 1
          }

          sliderInput(
            inputId = "xvrbl_scale",
            label = "Min/Max",
            min = min_value,
            max = max_value,
            value = c(min_value, max_value),
            step = 0.01,
            dragRange = TRUE  # Allow users to drag the range handles
            )
    })
  })

  output$yvrbl_scale <- renderUI({
    withProgress(message = "Loading Scale", value = 0, {
      req(input$yvrbl)  # Use req() to check if inputs are not NULL
        df <- df_neofox$mainTable_neofox
        df <- type.convert(df, as.is = TRUE)
        df[is.na(df)] <- 0
        df <- df[is.finite(df[[input$yvrbl]]),]

        # Apply log or sqrt transformation
        if (input$LogY == "ln") {
          df[[input$yvrbl]] <- log(ifelse(df[[input$yvrbl]] == 0, 1e-10, df[[input$yvrbl]]))
        } else if (input$LogY == "log2") {
          df[[input$yvrbl]] <- log2(ifelse(df[[input$yvrbl]] == 0, 1e-10, df[[input$yvrbl]]))
        } else if (input$LogY == "log10") {
          df[[input$yvrbl]] <- log10(ifelse(df[[input$yvrbl]] == 0, 1e-10, df[[input$yvrbl]]))
        } else if (input$LogY == "sqrt") {
          df[[input$yvrbl]] <- sqrt(ifelse(df[[input$yvrbl]] < 0, 1e-10, df[[input$yvrbl]]))
        } else {
          df[[input$yvrbl]] <- df[[input$yvrbl]]
        }


        df <- df[is.finite(df[[input$yvrbl]]),]

        yvrbl_values <- df[[input$yvrbl]]
        range_values <- range(as.numeric(yvrbl_values), na.rm = TRUE)
        min_value <- as.numeric(format(round(range_values[1], 2), nsmall = 2))
        max_value <- as.numeric(format(round(range_values[2], 2), nsmall = 2))


        # Check if min_value and max_value are equal, set default values
        if (min_value == max_value) {
          min_value <- min_value - 1
          max_value <- max_value + 1
        }

        sliderInput(
          inputId = "yvrbl_scale",
          label = "Min/Max",
          min = min_value,
          max = max_value,
          value = c(min_value, max_value),
          step = 0.01,
          dragRange = TRUE  # Allow users to drag the range handles
        )
    })
  })

  output$color_neofox <- renderUI({
    df <- df_neofox$mainTable_neofox
    df <- type.convert(df, as.is = TRUE)
    df[is.na(df)] <- 0

    features <- names(df)[sapply(df, is.numeric)]
    sorted_features <- features[order(!grepl("^\\*", features))]

    default_selection <- "*Tcell_predictor"
    pickerInput(inputId = "color_scatter",
                label = "Color",
                choices = sorted_features,
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE
    )
  })

  output$size_neofox <- renderUI({
    df <- df_neofox$mainTable_neofox
    df <- type.convert(df, as.is = TRUE)
    df[is.na(df)] <- 0

    features <- names(df)[sapply(df, is.numeric)]
    sorted_features <- features[order(!grepl("^\\*", features))]

    default_selection <- "*rnaExpression"
    pickerInput(inputId = "size_scatter",
                label = "Size",
                choices = sorted_features,
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE
    )
  })

  output$min_color <- renderUI({
    colourInput("min_col", "Select min color", "#434c4c")
  })

  output$max_color <- renderUI({
    colourInput("max_col", "Select max color", "#14c4c4")
  })

  output$scatter <- renderPlotly({
    withProgress(message = "Loading Scatter Plots", value = 0, {
      incProgress(0.5)
      if (!is.null(input$xvrbl) & !(is.null(input$yvrbl))) {
        df <- df_neofox$mainTable_neofox
        df <- type.convert(df, as.is = TRUE)
        df[is.na(df)] <- 0


        # For input$xvrbl
        if (input$LogX == "ln") {
          df[[input$xvrbl]] <- log(ifelse(df[[input$xvrbl]] == 0, 1e-10, df[[input$xvrbl]]))
        } else if (input$LogX == "log2") {
          df[[input$xvrbl]] <- log2(ifelse(df[[input$xvrbl]] == 0, 1e-10, df[[input$xvrbl]]))
        } else if (input$LogX == "log10") {
          df[[input$xvrbl]] <- log10(ifelse(df[[input$xvrbl]] == 0, 1e-10, df[[input$xvrbl]]))
        } else if (input$LogX == "sqrt") {
          df[[input$xvrbl]] <- sqrt(ifelse(df[[input$xvrbl]] < 0, 1e-10, df[[input$xvrbl]]))
        } else {
          df[[input$xvrbl]] <- df[[input$xvrbl]]
        }

        # For input$yvrbl
        if (input$LogY == "ln") {
          df[[input$yvrbl]] <- log(ifelse(df[[input$yvrbl]] == 0, 1e-10, df[[input$yvrbl]]))
        } else if (input$LogY == "log2") {
          df[[input$yvrbl]] <- log2(ifelse(df[[input$yvrbl]] == 0, 1e-10, df[[input$yvrbl]]))
        } else if (input$LogY == "log10") {
          df[[input$yvrbl]] <- log10(ifelse(df[[input$yvrbl]] == 0, 1e-10, df[[input$yvrbl]]))
        } else if (input$LogY == "sqrt") {
          df[[input$yvrbl]] <- sqrt(ifelse(df[[input$yvrbl]] < 0, 1e-10, df[[input$yvrbl]]))
        } else {
          df[[input$yvrbl]] <- df[[input$yvrbl]]
        }


        df[is.na(df)] <- 0

        # Filter data based on the slider range
        df <- subset(df, df[[input$xvrbl]] >= input$xvrbl_scale[1] & df[[input$xvrbl]] <= input$xvrbl_scale[2])
        df <- subset(df, df[[input$yvrbl]] >= input$yvrbl_scale[1] & df[[input$yvrbl]] <= input$yvrbl_scale[2])


        incProgress(0.5)

        scatter_plot <- ggplot(df , aes(x = .data[[input$xvrbl]], y = .data[[input$yvrbl]],
                                        text = paste("Patient:", .data[["patientIdentifier"]], "<br>",
                                                     "Gene:", .data[["gene"]]))) +
          geom_point(aes(color = .data[[input$color_scatter]], size = .data[[input$size_scatter]])) +  # Correct placement of aes() here
          scale_color_gradient(low = input$min_col, high = input$max_col) +
          theme(strip.text = element_text(size = 15), axis.text = element_text(size = 10), axis.title = element_text(size = 15), axis.ticks = element_line(size = 3), legend.text = element_text(size = 15), legend.title = element_text(size = 15))

        scatter_plot <- ggplotly(scatter_plot)

        print(scatter_plot)
        }
        else {
        p <- ggplot() + annotate(geom = "text", x = 10, y = 20, label = "No data available", size = 6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        scatter_plot <- ggplotly(p)
        incProgress(1)
        print(scatter_plot)
      }
    })
  })

  ############### Custom Tab ##########################
  df_custom <- reactiveValues(
    selectedRow = 1,
    fullData = NULL,
    mainTable_custom = NULL,
    group_inds = NULL,
    metricsData = NULL,
    pageLength = 10,
    groupBy = NULL,
    orderBy = NULL,
    peptide_features = NULL
  )
  
  # Option 1: User uploaded 
  observeEvent(input$custom_data$datapath, {
    mainData_custom <- read.table(input$custom_data$datapath, sep = "\t",  header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(mainData_custom) <- mainData_custom[1, ]
    mainData_custom <- mainData_custom[-1, ]
    row.names(mainData_custom) <- NULL
    mainData_custom <- type.convert(mainData_custom, as.is = TRUE)
    mainData_custom[is.na(mainData_custom)] <- 0
    df_custom$fullData <- mainData_custom
  })

  # Option 2: Demo Data
  observeEvent(input$loadDefault_Vaxrank, {
    data <- "data/vaxrank_output.tsv"
    mainData_custom <- read.table(data, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(mainData_custom) <- mainData_custom[1, ]
    mainData_custom <- mainData_custom[-1, ]
    row.names(mainData_custom) <- NULL
    mainData_custom <- type.convert(mainData_custom, as.is = TRUE)
    mainData_custom[is.na(mainData_custom)] <- 0
    df_custom$fullData <- mainData_custom
  })
  observeEvent(input$loadDefault_Neopredpipe, {
    data <- "data/neopredpipe.neoantigens.txt"
    mainData_custom <- read.table(data, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(mainData_custom) <- mainData_custom[1, ]
    mainData_custom <- mainData_custom[-1, ]
    row.names(mainData_custom) <- NULL
    mainData_custom <- type.convert(mainData_custom, as.is = TRUE)
    mainData_custom[is.na(mainData_custom)] <- 0
    df_custom$fullData <- mainData_custom
  })
  observeEvent(input$loadDefault_antigengarnish, {
    data <- "data/antigengarnish_test_antigen.tsv"
    mainData_custom <- read.table(data, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
    colnames(mainData_custom) <- mainData_custom[1, ]
    mainData_custom <- mainData_custom[-1, ]
    row.names(mainData_custom) <- NULL
    mainData_custom <- type.convert(mainData_custom, as.is = TRUE)
    mainData_custom[is.na(mainData_custom)] <- 0
    df_custom$fullData <- mainData_custom
  })
  
  # group peptides by
  output$custom_group_by_feature_ui <- renderUI({
    
    features <- as.list(names(df_custom$fullData))
    default_selection <- ifelse(length(features) >= 1, features[[1]], "")
    
    pickerInput(inputId = "feature_1",
                label = "Group peptides by",
                choices = features, # a list of strings
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE)
  })
  
  # Sort peptides by
  output$custom_order_by_feature_ui <- renderUI({
    
    features <- names(df_custom$fullData)[sapply(df_custom$fullData, is.numeric)]
    default_selection <- ifelse(length(features) >= 2, features[[2]], "")
    
    pickerInput(inputId = "feature_2",
                label = "Sort peptides by",
                choices = features, 
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE)
  })
  
  ## Clear file inputs if demo data load button is clicked
  output$custom_upload_ui <- renderUI({
    fileInput(inputId = "custom_data", label = "Custom input table (tsv required)",
              accept =  c("text/tsv", "text/tab-separated-values,text/plain", ".tsv"))
  })
  
  observeEvent(input$visualize_custom, {
    updateTabItems(session, "custom_tabs", "custom_explore")
    updateSliderInput(session, "xvrbl_scale_custom", value = c(0, 1))  # Set default values for x-axis slider input
    updateSliderInput(session, "yvrbl_scale_custom", value = c(0, 1))  # Set default values for y-axis slider input
  })

  ## Custom Explore
  output$group_feature <- renderUI({
    if (!is.null(input$feature_1)) {
      HTML(paste("Peptides grouped by: ", input$feature_1))
    } else {
      HTML("No feature selected")
    }
  })
  
  output$sort_feature <- renderUI({
    if (!is.null(input$feature_2)) {
      HTML(paste("Peptides sorted by: ", input$feature_2))
    } else {
      HTML("No feature selected")
    }
  })
  

  # features to display 
  output$custom_peptide_features_ui <- renderUI({
    
    features <- as.list(names(df_custom$fullData))
    default_selection <- features
    
    pickerInput(inputId = "peptide_features",
                label = "Features to display for each group of peptides",
                selected = default_selection,
                options = list(`actions-box` = TRUE,`live-search` = TRUE),
                choices = features, # a list of strings
                multiple = TRUE)
  })
  
  observeEvent(input$visualize_custom, {
    #browser()
    df_custom$groupBy <- input$feature_1
    df_custom$orderBy <- input$feature_2
    
    reformat_data <- df_custom$fullData %>% group_by(across(all_of(df_custom$groupBy))) %>% arrange(across(all_of(df_custom$orderBy)))
    reformat_data <- type.convert(reformat_data, as.is = TRUE)
    reformat_data[is.na(reformat_data)] <- 0
    
    df_custom$fullData <- reformat_data
    row_ind <- reformat_data %>% group_rows()
    row_ind_df <- as.data.frame(row_ind)
    df_custom$group_inds <- row_ind_df
    row_ind_top <- apply(row_ind_df, 1, function(x) {unlist(x[1])[1]})
    df_custom$mainTable_custom <- as.data.frame(reformat_data[row_ind_top, ])
    df_custom$mainTable_custom <- cbind(Select = shinyInputSelect(actionButton, nrow(df_custom$mainTable_custom), "button_", label = "Investigate", onclick = 'Shiny.onInputChange(\"custom_select_button\",  this.id)'), df_custom$mainTable_custom)
    df_custom$metricsData <- get_group_inds(df_custom$fullData, df_custom$group_inds)
    df_custom$peptide_features <- input$peptide_features
    updateTabItems(session, "custom_tabs", "custom_explore")
    
    updateSelectInput(session, "feature_1", selected = df_custom$groupBy)
    updateSelectInput(session, "feature_2", selected = df_custom$orderBy)
  })
  
  output$customTable <- DT::renderDataTable(
    if (is.null(df_custom$mainTable_custom)) {
      return(datatable(data.frame("Annotated Table" = character())))
    }else {
      datatable(df_custom$mainTable_custom,
                escape = FALSE, class = "stripe",
                selection = "single",
                extensions = c("Buttons")
                
      )
    }, server = FALSE)
  

  observeEvent(input$custom_select_button, {
    if (is.null(df_custom$mainTable_custom) | is.null(df_custom$selectedRow)){
      return ()
    }
    #browser()
    df_custom$selectedRow <- as.numeric(strsplit(input$custom_select_button, "_")[[1]][2])
    session$sendCustomMessage('unbind-DT', 'customTable')
    dataTableProxy("customMainTable") %>% 
      selectPage((df_custom$selectedRow-1) %/% df_custom$pageLength + 1)
  })
  
  
  output$customPeptideTable <- renderDT({
    withProgress(message = 'Loading Peptide Table', value = 0, {
      incProgress(0.5)
      #browser()
      if (!is.null(df_custom$selectedRow) & !(is.null(df_custom$mainTable_custom)) & !is.null(df_custom$peptide_features)){
    
        
        display_table <- get_current_group_info(df_custom$peptide_features, df_custom$metricsData, df_custom$fullData, df_custom$selectedRow)
        incProgress(0.5)
        dtable_custom <- datatable(display_table, options =list(
          pageLength = 10,
          rowCallback = JS('function(row, data, index, rowId) {',
                           'console.log(rowId)','if(((rowId+1) % 4) == 3 || ((rowId+1) % 4) == 0) {',
                           'row.style.backgroundColor = "#E0E0E0";','}','}'
                           )
        ), selection = list(mode='single', selected = '1')) 
        dtable_custom
      }
      else{
        incProgress(1)
        datatable(data.frame("Peptide Datatable"=character()), selection = list(mode='single', selected = '1'))
      }})
  })
  
  
  # Dynamic Scatter Plot 
  output$xvrbl_custom <- renderUI({
    #df_custom <- df_custom$fullData
    #df_custom <- type.convert(df_custom, as.is = TRUE)
    #df[is.na(df)] <- 0
    
    features <- as.list(names(df_custom$fullData)[sapply(df_custom$fullData, is.numeric)])
    default_selection <- ifelse(length(features) >= 1, features[[1]], "")
    
    pickerInput(inputId = "xvrbl_custom",
                label = "X-Axis Variable",
                choices = features,
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE
    )
  })
  
  output$xvrbl_log_custom <- renderUI({
    radioButtons(
      inputId = "LogX_custom",
      choices = c("none", "ln", "log2", "log10", "sqrt"),
      label = "X Log Transform",
      inline = TRUE
    )
  })
  
  output$xvrbl_scale_custom <- renderUI({
    if (!is.null(input$xvrbl_custom) && input$xvrbl_custom != "") {
      withProgress(message = "Loading Scale", value = 0, {
          req(input$xvrbl_custom)  # Use req() to check if inputs are not NULL
          df <- df_custom$fullData
          df <- type.convert(df, as.is = TRUE)
          df[is.na(df)] <- 0
          df <- df[is.finite(df[[input$xvrbl_custom]]),]
          
          # Apply log or sqrt transformation
          if (input$LogX_custom == "ln") {
            df[[input$xvrbl_custom]] <- log(ifelse(df[[input$xvrbl_custom]] == 0, 1e-10, df[[input$xvrbl_custom]]))
          } else if (input$LogX_custom == "log2") {
            df[[input$xvrbl_custom]] <- log2(ifelse(df[[input$xvrbl_custom]] == 0, 1e-10, df[[input$xvrbl_custom]]))
          } else if (input$LogX_custom == "log10") {
            df[[input$xvrbl_custom]] <- log10(ifelse(df[[input$xvrbl_custom]] == 0, 1e-10, df[[input$xvrbl_custom]]))
          } else if (input$LogX_custom == "sqrt") {
            df[[input$xvrbl_custom]] <- sqrt(ifelse(df[[input$xvrbl_custom]] < 0, 1e-10, df[[input$xvrbl_custom]]))
          } else {
            df[[input$xvrbl_custom]] <- df[[input$xvrbl_custom]]
          }
          
          df <- df[is.finite(df[[input$xvrbl_custom]]),]
          
          xvrbl_values <- df[[input$xvrbl_custom]]
          range_values <- range(as.numeric(xvrbl_values), na.rm = TRUE)
          min_value <- as.numeric(format(round(range_values[1], 2), nsmall = 2))
          max_value <- as.numeric(format(round(range_values[2], 2), nsmall = 2))
          
          
          # Check if min_value and max_value are equal, set default values
          if (min_value == max_value) {
            min_value <- min_value - 1
            max_value <- max_value + 1
          }
          
          sliderInput(
            inputId = "xvrbl_scale_custom",
            label = "X Min/Max",
            min = min_value,
            max = max_value,
            value = c(min_value, max_value),
            step = 0.01,
            dragRange = TRUE  # Allow users to drag the range handles 
          )
      })
    }
  })
  
  output$yvrbl_custom <- renderUI({
    df <- df_custom$fullData
    df <- type.convert(df, as.is = TRUE)
    df[is.na(df)] <- 0
    
    features <- as.list(names(df)[sapply(df, is.numeric)])
    default_selection <- ifelse(length(features) >= 2, features[[2]], "")
    
    pickerInput(inputId = "yvrbl_custom",
                label = "Y-Axis Variable",
                choices = features,
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE
    )
  })
  
  output$yvrbl_log_custom <- renderUI({
    radioButtons(
      inputId = "LogY_custom",
      choices = c("none", "ln", "log2", "log10", "sqrt"),
      label = "Y Log Transform",
      inline = TRUE
    )
  })
  
  output$yvrbl_scale_custom <- renderUI({
    if (!is.null(input$yvrbl_custom) && input$yvrbl_custom != "") {
          
      withProgress(message = "Loading Scale", value = 0, {
        req(input$yvrbl_custom)  # Use req() to check if inputs are not NULL
        df <- df_custom$fullData
        df <- type.convert(df, as.is = TRUE)
        df[is.na(df)] <- 0
        df <- df[is.finite(df[[input$yvrbl_custom]]),]
        
        # Apply log or sqrt transformation
        if (input$LogY_custom == "ln") {
          df[[input$yvrbl_custom]] <- log(ifelse(df[[input$yvrbl_custom]] == 0, 1e-10, df[[input$yvrbl_custom]]))
        } else if (input$LogY_custom == "log2") {
          df[[input$yvrbl_custom]] <- log2(ifelse(df[[input$yvrbl_custom]] == 0, 1e-10, df[[input$yvrbl_custom]]))
        } else if (input$LogY_custom == "log10") {
          df[[input$yvrbl_custom]] <- log10(ifelse(df[[input$yvrbl_custom]] == 0, 1e-10, df[[input$yvrbl_custom]]))
        } else if (input$LogY_custom == "sqrt") {
          df[[input$yvrbl_custom]] <- sqrt(ifelse(df[[input$yvrbl_custom]] < 0, 1e-10, df[[input$yvrbl_custom]]))
        } else {
          df[[input$yvrbl_custom]] <- df[[input$yvrbl_custom]]
        }
        
        df <- df[is.finite(df[[input$yvrbl_custom]]),]
        yvrbl_values <- df[[input$yvrbl_custom]]
        range_values <- range(as.numeric(yvrbl_values), na.rm = TRUE)
        min_value <- as.numeric(format(round(range_values[1], 2), nsmall = 2))
        max_value <- as.numeric(format(round(range_values[2], 2), nsmall = 2))
        
        # Check if min_value and max_value are equal, set default values
        if (min_value == max_value) {
          min_value <- min_value - 1
          max_value <- max_value + 1
        }
        
        sliderInput(
          inputId = "yvrbl_scale_custom",
          label = "Y Min/Max",
          min = min_value,
          max = max_value,
          value = c(min_value, max_value),
          step = 0.01,
          dragRange = TRUE  # Allow users to drag the range handles 
        )
      })
    }
  })
  
  
  output$color_custom <- renderUI({
    df <- df_custom$fullData
    df <- type.convert(df, as.is = TRUE)
    df[is.na(df)] <- 0
    
    features <- as.list(names(df)[sapply(df, is.numeric)])
    default_selection <- ifelse(length(features) >= 3, features[[3]], "")
    
    pickerInput(inputId = "color_scatter_custom",
                label = "Color",
                choices = features,
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE
    )
  })
  
  output$min_color_custom <- renderUI({
    colourInput("min_col_custom", "Select min color", "#434c4c")
  })
  
  output$max_color_custom <- renderUI({
    colourInput("max_col_custom", "Select max color", "#14c4c4")
  })
  
  output$size_custom<- renderUI({
    df <- df_custom$fullData
    df <- type.convert(df, as.is = TRUE)
    df[is.na(df)] <- 0
    
    features <- as.list(names(df)[sapply(df, is.numeric)])
    default_selection <- ifelse(length(features) >= 4, features[[4]], "")

    pickerInput(inputId = "size_scatter_custom",
                label = "Size",
                choices = features,
                selected = default_selection,
                options = list(`live-search` = TRUE),
                multiple = FALSE
    )
  })
  
  output$scatter_custom <- renderPlotly({
    withProgress(message = "Loading Scatter Plots", value = 0, {
      incProgress(0.5)
      if (!is.null(input$xvrbl_custom) & !is.null(input$yvrbl_custom)
          & !is.null(input$min_col_custom)  & !is.null(input$max_col_custom)
          & !is.null(input$xvrbl_scale_custom)  & !is.null(input$yvrbl_scale_custom)
          & !is.null(input$color_scatter_custom)  & !is.null(input$size_scatter_custom)) {
        
        df <- df_custom$fullData
        df <- type.convert(df, as.is = TRUE)
        df[is.na(df)] <- 0
        
        
        # For input$xvrbl_custom
        if (input$LogX_custom == "ln") {
          df[[input$xvrbl_custom]] <- log(ifelse(df[[input$xvrbl_custom]] == 0, 1e-10, df[[input$xvrbl_custom]]))
        } else if (input$LogX_custom == "log2") {
          df[[input$xvrbl_custom]] <- log2(ifelse(df[[input$xvrbl_custom]] == 0, 1e-10, df[[input$xvrbl_custom]]))
        } else if (input$LogX_custom == "log10") {
          df[[input$xvrbl_custom]] <- log10(ifelse(df[[input$xvrbl_custom]] == 0, 1e-10, df[[input$xvrbl_custom]]))
        } else if (input$LogX_custom == "sqrt") {
          df[[input$xvrbl_custom]] <- sqrt(ifelse(df[[input$xvrbl_custom]] < 0, 1e-10, df[[input$xvrbl_custom]]))
        } else {
          df[[input$xvrbl_custom]] <- df[[input$xvrbl_custom]]
        }
        # For input$yvrbl_custom
        if (input$LogY_custom == "ln") {
          df[[input$yvrbl_custom]] <- log(ifelse(df[[input$yvrbl_custom]] == 0, 1e-10, df[[input$yvrbl_custom]]))
        } else if (input$LogY_custom == "log2") {
          df[[input$yvrbl_custom]] <- log2(ifelse(df[[input$yvrbl_custom]] == 0, 1e-10, df[[input$yvrbl_custom]]))
        } else if (input$LogY_custom == "log10") {
          df[[input$yvrbl_custom]] <- log10(ifelse(df[[input$yvrbl_custom]] == 0, 1e-10, df[[input$yvrbl_custom]]))
        } else if (input$LogY_custom == "sqrt") {
          df[[input$yvrbl_custom]] <- sqrt(ifelse(df[[input$yvrbl_custom]] < 0, 1e-10, df[[input$yvrbl_custom]]))
        } else {
          df[[input$yvrbl_custom]] <- df[[input$yvrbl_custom]]
        }
        
        
        df[is.na(df)] <- 0
        
        # Filter data based on the slider range
        df <- subset(df, df[[input$xvrbl_custom]] >= input$xvrbl_scale_custom[1] & df[[input$xvrbl_custom]] <= input$xvrbl_scale_custom[2])
        df <- subset(df, df[[input$yvrbl_custom]] >= input$yvrbl_scale_custom[1] & df[[input$yvrbl_custom]] <= input$yvrbl_scale_custom[2])
        
        incProgress(0.5)
        
        scatter_plot <- ggplot(df , aes(x = .data[[input$xvrbl_custom]], y = .data[[input$yvrbl_custom]])) +
          geom_point(aes(color = .data[[input$color_scatter_custom]], size = .data[[input$size_scatter_custom]])) +  # Correct placement of aes() here
          scale_color_gradient(low = input$min_col_custom, high = input$max_col_custom) +
          theme(strip.text = element_text(size = 15), axis.text = element_text(size = 10), axis.title = element_text(size = 15), axis.ticks = element_line(size = 3), legend.text = element_text(size = 15), legend.title = element_text(size = 15))
        
        scatter_plot <- ggplotly(scatter_plot)
        
        #print(scatter_plot)
      } 
      else {
          ggplot() + annotate(geom = "text", x = 10, y = 20, label = "No data available", size = 6) +
          theme_void() + theme(legend.position = "none", panel.border = element_blank())
        incProgress(1)
      }
    })
  })
  
})
