### OTHER CODE


############### NeoFox Tab ##########################
df_neofox <- reactiveValues(
  mainTable = NULL
)
observeEvent(input$loadDefaultneofox, {
  #session$sendCustomMessage("unbind-DT", "neofoxTable")
  #data <- getURL("https://raw.githubusercontent.com/griffithlab/pVACtools/f83c52c8b8387beae69be8b200a44dcf199d9af2/pvactools/tools/pvacview/data/H_NJ-HCC1395-HCC1395.Class_I.all_epitopes.aggregated.tsv")
  #mainData <- read.table(text = data, sep = '\t', header = FALSE, stringsAsFactors = FALSE, check.names=FALSE)
  data <- "~/Desktop/Griffith_Lab/R_shiny_visualization/neoantigen_visualization/v11/data/test_pt1_neoantigen_candidates_annotated.tsv"
  mainData <- read.table(data, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(mainData) <- mainData[1, ]
  mainData <- mainData[-1, ]
  row.names(mainData) <- NULL
  df_neofox$mainTable <- mainData
  updateTabItems(session, "neofox_tabs", "neofox_explore")
})
output$neofox_upload_ui <- renderUI({
  fileInput(inputId = "neofox_data", label = "NeoFox output table (tsv required)",
            accept =  c("text/tsv", "text/tab-separated-values,text/plain", ".tsv"))
})
observeEvent(input$neofox_data$datapath, {
  #session$sendCustomMessage("unbind-DT", "neofoxTable")
  mainData <- read.table(input$neofox_data$datapath, sep = "\t",  header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(mainData) <- mainData[1, ]
  mainData <- mainData[-1, ]
  row.names(mainData) <- NULL
  df_neofox$mainTable <- mainData
})
observeEvent(input$visualize_neofox, {
  updateTabItems(session, "neofox_tabs", "neofox_explore")
})
output$neofoxTable <- DT::renderDataTable(
  if (is.null(df_neofox$mainTable)) {
    return(datatable(data.frame("Annotated Table" = character())))
  }else {
    datatable(df_neofox$mainTable,
              escape = FALSE, class = "stripe",
              options = list(lengthChange = FALSE, dom = "Bfrtip", pageLength = input$neofox_page_length,
                             columnDefs = list(list(visible = FALSE, targets = c(-1:-12)),
                                               list(orderable = TRUE, targets = 0)), buttons = list(I("colvis")),
                             initComplete = htmlwidgets::JS(
                               "function(settings, json) {",
                               paste("$(this.api().table().header()).css({'font-size': '", "10pt", "'});"),
                               "}")),
              selection = "multiple",
              extensions = c("Buttons"))
  }, server = FALSE)
output$neofox_selected <- renderText({
  if (is.null(df_neofox$mainTable)) {
    return()
  }
  input$neofoxTable_rows_selected
})
output$neofox_violin_plots_row1 <- renderPlot({
  withProgress(message = "Loading Violin Plots", value = 0, {
    if (length(input$neofoxTable_rows_selected) != 0) {
      plot_cols <- c("mutatedXmer", "imputedGeneExpression", "DAI_MHCI_bestAffinity", "IEDB_Immunogenicity_MHCI")
      plot_data <- df_neofox$mainTable[, plot_cols]
      plot_data[, "imputedGeneExpression"] <- as.numeric(plot_data[, "imputedGeneExpression"])
      plot_data[, "DAI_MHCI_bestAffinity"] <- as.numeric(plot_data[, "DAI_MHCI_bestAffinity"])
      plot_data[, "IEDB_Immunogenicity_MHCI"] <- as.numeric(plot_data[, "IEDB_Immunogenicity_MHCI"])
      plot_data$Selected <- "No"
      plot_data[input$neofoxTable_rows_selected, "Selected"] <- "Yes"
      reformat_data <- plot_data %>%
        gather("Feature", "Value", colnames(plot_data)[2]:tail(colnames(plot_data), n = 2))
      gene_expr_data <- reformat_data[reformat_data["Feature"] == "imputedGeneExpression", ]
      gene_expr_plot <- ggplot(data = gene_expr_data, aes(x = Feature, y = Value)) + geom_violin() +
        geom_jitter(data = gene_expr_data[gene_expr_data["Selected"] == "No", ], aes(color = Selected), size = 1, alpha = 0.5, stroke = 1, position = position_jitter(0.3)) +
        geom_jitter(data = gene_expr_data[gene_expr_data["Selected"] == "Yes", ], aes(color = Selected), size = 2, alpha = 1, stroke = 1, position = position_jitter(0.3)) +
        scale_color_manual(values = c("No" = "#939094", "Yes" = "#f42409")) +
        theme(legend.position = "none")
      DAI_ClassI_data <- reformat_data[reformat_data["Feature"] == "DAI_MHCI_bestAffinity", ]
      DAI_ClassI_plot <- ggplot(DAI_ClassI_data, aes(x = Feature, y = Value)) + geom_violin() +
        geom_jitter(data = DAI_ClassI_data[DAI_ClassI_data["Selected"] == "No", ], aes(color = Selected), size = 1, alpha = 0.5, stroke = 1, position = position_jitter(0.3)) +
        geom_jitter(data = DAI_ClassI_data[DAI_ClassI_data["Selected"] == "Yes", ], aes(color = Selected), size = 2, alpha = 1, stroke = 1, position = position_jitter(0.3)) +
        scale_color_manual(values = c("No" = "#939094", "Yes" = "#f42409")) +
        theme(legend.position = "none")
      IEDB_Immuno_data <- reformat_data[reformat_data["Feature"] == "IEDB_Immunogenicity_MHCI", ]
      IEDB_Immuno_plot <- ggplot(IEDB_Immuno_data, aes(x = Feature, y = Value)) + geom_violin() +
        geom_jitter(data = IEDB_Immuno_data[IEDB_Immuno_data["Selected"] == "No", ], aes(color = Selected), size = 1, alpha = 0.5, stroke = 1, position = position_jitter(0.3)) +
        geom_jitter(data = IEDB_Immuno_data[IEDB_Immuno_data["Selected"] == "Yes", ], aes(color = Selected), size = 2, alpha = 1, stroke = 1, position = position_jitter(0.3)) +
        scale_color_manual(values = c("No" = "#939094", "Yes" = "#f42409"))
      p <- grid.arrange(gene_expr_plot, DAI_ClassI_plot, IEDB_Immuno_plot, ncol = 3)
      incProgress(1)
      print(p)
    }else {
      p <- ggplot() + annotate(geom = "text", x = 10, y = 20, label = "No data available", size = 6) +
        theme_void() + theme(legend.position = "none", panel.border = element_blank())
      incProgress(1)
      print(p)
    }
  })
})
output$neofox_violin_plots_row2 <- renderPlot({
  withProgress(message = "Loading Violin Plots", value = 0, {
    if (length(input$neofoxTable_rows_selected) != 0) {
      plot_cols <- c("mutatedXmer", "MixMHCpred_bestScore_rank", "HexAlignmentScore_MHCI", "PRIME_best_rank")
      plot_data <- df_neofox$mainTable[, plot_cols]
      plot_data[, "MixMHCpred_bestScore_rank"] <- as.numeric(plot_data[, "MixMHCpred_bestScore_rank"])
      plot_data[, "HexAlignmentScore_MHCI"] <- as.numeric(plot_data[, "HexAlignmentScore_MHCI"])
      plot_data[, "PRIME_best_rank"] <- as.numeric(plot_data[, "PRIME_best_rank"])
      plot_data$Selected <- "No"
      plot_data[input$neofoxTable_rows_selected, "Selected"] <- "Yes"
      reformat_data <- plot_data %>%
        gather("Feature", "Value", colnames(plot_data)[2]:tail(colnames(plot_data), n = 2))
      mixpred_data <- reformat_data[reformat_data["Feature"] == "MixMHCpred_bestScore_rank", ]
      mixpred_plot <- ggplot(data = mixpred_data, aes(x = Feature, y = Value)) + geom_violin() +
        geom_jitter(data = mixpred_data[mixpred_data["Selected"] == "No", ], aes(color = Selected), size = 1, alpha = 0.5, stroke = 1, position = position_jitter(0.3)) +
        geom_jitter(data = mixpred_data[mixpred_data["Selected"] == "Yes", ], aes(color = Selected), size = 2, alpha = 1, stroke = 1, position = position_jitter(0.3)) +
        scale_color_manual(values = c("No" = "#939094", "Yes" = "#f42409")) +
        theme(legend.position = "none")
      hex_data <- reformat_data[reformat_data["Feature"] == "HexAlignmentScore_MHCI", ]
      hex_plot <- ggplot(hex_data, aes(x = Feature, y = Value)) + geom_violin() +
        geom_jitter(data = hex_data[hex_data["Selected"] == "No", ], aes(color = Selected), size = 1, alpha = 0.5, stroke = 1, position = position_jitter(0.3)) +
        geom_jitter(data = hex_data[hex_data["Selected"] == "Yes", ], aes(color = Selected), size = 2, alpha = 1, stroke = 1, position = position_jitter(0.3)) +
        scale_color_manual(values = c("No" = "#939094", "Yes" = "#f42409")) +
        theme(legend.position = "none")
      prime_data <- reformat_data[reformat_data["Feature"] == "PRIME_best_rank", ]
      prime_plot <- ggplot(prime_data, aes(x = Feature, y = Value)) + geom_violin() +
        geom_jitter(data = prime_data[prime_data["Selected"] == "No", ], aes(color = Selected), size = 1, alpha = 0.5, stroke = 1, position = position_jitter(0.3)) +
        geom_jitter(data = prime_data[prime_data["Selected"] == "Yes", ], aes(color = Selected), size = 2, alpha = 1, stroke = 1, position = position_jitter(0.3)) +
        scale_color_manual(values = c("No" = "#939094", "Yes" = "#f42409"))
      p <- grid.arrange(mixpred_plot, hex_plot, prime_plot, ncol = 3)
      incProgress(1)
      print(p)
    }else {
      p <- ggplot() + annotate(geom = "text", x = 10, y = 20, label = "No data available", size = 6) +
        theme_void() + theme(legend.position = "none", panel.border = element_blank())
      incProgress(1)
      print(p)
    }
  })
})
############### Custom Tab ##########################
df_custom <- reactiveValues(
  selectedRow = 1,
  fullData = NULL,
  mainTable = NULL,
  group_inds = NULL,
  metricsData = NULL,
  pageLength = 10,
  groupBy = NULL,
  orderBy = NULL,
  peptide_features = NULL
)
observeEvent(input$loadDefault_Vaxrank, {
  data <- "~/Desktop/Griffith_Lab/R_shiny_visualization/neoantigen_visualization/v11/data/vaxrank_output.tsv"
  mainData <- read.table(data, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(mainData) <- mainData[1, ]
  mainData <- mainData[-1, ]
  row.names(mainData) <- NULL
  df_custom$fullData <- mainData
})
observeEvent(input$loadDefault_Neopredpipe, {
  data <- "~/Desktop/Griffith_Lab/R_shiny_visualization/neoantigen_visualization/v11/data/HCC1395Run.neoantigens.txt"
  mainData <- read.table(data, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(mainData) <- mainData[1, ]
  mainData <- mainData[-1, ]
  row.names(mainData) <- NULL
  df_custom$fullData <- mainData
})
observeEvent(input$loadDefault_antigengarnish, {
  data <- "~/Desktop/Griffith_Lab/R_shiny_visualization/neoantigen_visualization/v11/data/ag_test_antigen.tsv"
  mainData <- read.table(data, sep = "\t", header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(mainData) <- mainData[1, ]
  mainData <- mainData[-1, ]
  row.names(mainData) <- NULL
  df_custom$fullData <- mainData
})
output$custom_upload_ui <- renderUI({
  fileInput(inputId = "custom_data", label = "Custom input table (tsv required)",
            accept =  c("text/tsv", "text/tab-separated-values,text/plain", ".tsv"))
})
observeEvent(input$custom_data$datapath, {
  mainData <- read.table(input$custom_data$datapath, sep = "\t",  header = FALSE, stringsAsFactors = FALSE, check.names = FALSE)
  colnames(mainData) <- mainData[1, ]
  mainData <- mainData[-1, ]
  row.names(mainData) <- NULL
  df_custom$fullData <- mainData
})

output$custom_group_by_feature_ui <- renderUI({
  feature <- names(df_custom$fullData)
  pickerInput(inputId = "feature_1",
              label = "Feature to group peptides by",
              choices = feature, # a list of strings
              #options = list(`actions-box` = TRUE, `live-search` = TRUE),
              multiple = FALSE)
})
output$custom_order_by_feature_ui <- renderUI({
  feature <- names(df_custom$fullData)
  pickerInput(inputId = "feature_2",
              label = "Feature to sort peptides by",
              choices = feature, # a list of strings
              #options = list(`actions-box` = TRUE, `live-search` = TRUE),
              multiple = FALSE)
})
output$custom_peptide_features_ui <- renderUI({
  feature <- names(df_custom$fullData)
  pickerInput(inputId = "peptide_features",
              label = "Subset of features to display in peptide subtable",
              choices = feature[((feature != input$feature_2) & (feature != input$feature_1))], # a list of strings
              options = list(`actions-box` = TRUE, `live-search` = TRUE),
              multiple = TRUE)
})
observeEvent(input$visualize_custom, {
  #browser()
  df_custom$groupBy <- input$feature_1
  df_custom$orderBy <- input$feature_2
  reformat_data <- df_custom$fullData %>% group_by(across(all_of(df_custom$groupBy))) %>% arrange(across(all_of(df_custom$orderBy)))
  df_custom$fullData <- reformat_data
  row_ind <- reformat_data %>% group_rows()
  row_ind_df <- as.data.frame(row_ind)
  df_custom$group_inds <- row_ind_df
  row_ind_top <- apply(row_ind_df, 1, function(x) {unlist(x[1])[1]})
  df_custom$mainTable <- as.data.frame(reformat_data[row_ind_top, ])
  #df_custom$mainTable <- cbind("Eval" = shinyInput(df_custom$mainTable, selectInput, nrow(df_custom$mainTable), "custom_selecter_", choices = c("Pending", "Accept", "Reject", "Review"), width = "60px"), df_custom$mainTable)
  df_custom$mainTable <- cbind(Select = shinyInputSelect(actionButton, nrow(df_custom$mainTable), "button_", label = "Investigate", onclick = 'Shiny.onInputChange(\"custom_select_button\",  this.id)'), df_custom$mainTable)
  #if (is.null(df_custom$mainTable$`Evaluation`)) {
  #  df_custom$mainTable$`Evaluation` <- rep("Pending", nrow(df_custom$mainTable))
  #}
  #gene_data <- getURL("https://raw.githubusercontent.com/griffithlab/pVACtools/7c7b8352d81b44ec7743578e7715c65261f5dab7/pvactools/tools/pvacview/data/cancer_census_hotspot_gene_list.tsv")
  #gene_list <- read.table(text = gene_data, sep = '\t',  header = FALSE, stringsAsFactors = FALSE, check.names=FALSE)
  #df_custom$gene_list <- gene_list
  #df_custom$mainTable$`Gene of Interest` <- apply(df_custom$mainTable,1, function(x) {any(x['Gene Name'] == df_custom$gene_list)})
  df_custom$metricsData <- get_group_inds(df_custom$fullData, df_custom$group_inds)
  df_custom$peptide_features <- input$peptide_features
  updateTabItems(session, "custom_tabs", "custom_explore")
})
output$customTable <- DT::renderDataTable(
  if (is.null(df_custom$mainTable)) {
    return(datatable(data.frame("Annotated Table" = character())))
  }else {
    datatable(df_custom$mainTable,
              escape = FALSE, class = "stripe",
              options = list(lengthChange = FALSE, dom = "Bfrtip", pageLength = input$custom_page_length,
                             columnDefs = list(list(visible = FALSE, targets = c(-1:-12)),
                                               list(orderable = TRUE, targets = 0)), buttons = list(I("colvis")),
                             initComplete = htmlwidgets::JS(
                               "function(settings, json) {",
                               paste("$(this.api().table().header()).css({'font-size': '", "10pt", "'});"),
                               "}")),
              selection = "none",
              extensions = c("Buttons"))
  }, server = FALSE)
observeEvent(input$custom_select_button, {
  if (is.null(df_custom$mainTable) | is.null(df_custom$selectedRow)){
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
    if (!is.null(df_custom$selectedRow) & !(is.null(df_custom$mainTable)) & !is.null(df_custom$peptide_features)){
      display_table <- get_current_group_info(df_custom$peptide_features, df_custom$metricsData, df_custom$fullData, df_custom$selectedRow)
      incProgress(0.5)
      dtable <- datatable(display_table, options =list(
        pageLength = 10,
        rowCallback = JS('function(row, data, index, rowId) {',
                         'console.log(rowId)','if(((rowId+1) % 4) == 3 || ((rowId+1) % 4) == 0) {',
                         'row.style.backgroundColor = "#E0E0E0";','}','}')
      ), selection = list(mode='single', selected = '1')) 
      dtable
    }
    else{
      incProgress(1)
      datatable(data.frame("Peptide Datatable"=character()), selection = list(mode='single', selected = '1'))
    }})
})
