## helper function defined for generating shinyInputs in mainTable (Investigate button)
shinyInputSelect <- function(FUN, row_ids, button_label, ...) {
    inputs <- character(nrow(row_ids))
    for (i in 1:nrow(row_ids)) {
        inputs[i] <- as.character(FUN(paste0(button_label, row_ids[i, "ID"]), ...))
    }
    inputs
}

process_main_data <- function(mainData) {
    colnames(mainData) <- mainData[1, ]
    mainData <- mainData[-1, ]
    row.names(mainData) <- NULL
    mainData$Acpt <- shinyInputSelect(actionButton, mainData["ID"], "button-acpt_", icon = icon("thumbs-up"), label = "", onclick = 'Shiny.onInputChange(\"accept_eval\",  this.id, {priority: "event"})', onmousedown = "event.preventDefault(); event.stopPropagation();")
    mainData$Rej <- shinyInputSelect(actionButton, mainData["ID"], "button-rej_", icon = icon("thumbs-down"), label = "", onclick = 'Shiny.onInputChange(\"reject_eval\",  this.id, {priority: "event"})', onmousedown = "event.preventDefault(); event.stopPropagation();")
    mainData$Rev <- shinyInputSelect(actionButton, mainData["ID"], "button-rev_", icon = icon("flag"), label = "", onclick = 'Shiny.onInputChange(\"review_eval\",  this.id, {priority: "event"})', onmousedown = "event.preventDefault(); event.stopPropagation();")
    mainData$`IC50 MT` <- as.numeric(mainData$`IC50 MT`)
    mainData$`%ile MT` <- as.numeric(mainData$`%ile MT`)
    mainData$`RNA Depth` <- as.character(as.integer(mainData$`RNA Depth`))
    mainData$`TSL`[is.na(mainData$`TSL`)] <- "NA"
    return(mainData)
}

process_metrics_data <- function(df) {
    df$binding_threshold <- df$metricsData$`binding_threshold`
    df$use_allele_specific_binding_thresholds <- df$metricsData$`use_allele_specific_binding_thresholds`
    df$allele_specific_binding_thresholds <- df$metricsData$`allele_specific_binding_thresholds`
    df$aggregate_inclusion_binding_threshold <- df$metricsData$`aggregate_inclusion_binding_threshold`
    df$binding_percentile_threshold <- df$metricsData$`binding_percentile_threshold`
    df$immunogenicity_percentile_threshold <- df$metricsData$`immunogenicity_percentile_threshold`
    df$presentation_percentile_threshold <- df$metricsData$`presentation_percentile_threshold`
    df$percentile_threshold_strategy <- df$metricsData$`percentile_threshold_strategy`
    df$scoring_candidate_metric <- df$metricsData$`top_score_metric2`
    df$dna_cutoff <- df$metricsData$vaf_clonal
    df$allele_expr <- df$metricsData$allele_expr_threshold
    df$anchor_mode <- ifelse(df$metricsData$`allele_specific_anchors`, "allele-specific", "default")
    df$allele_specific_anchors <- df$metricsData$`allele_specific_anchors`
    df$anchor_contribution <- df$metricsData$`anchor_contribution_threshold`
    df$maximum_transcript_support_level <- df$metricsData$maximum_transcript_support_level
    df$transcript_prioritization_strategy <- df$metricsData$`transcript_prioritization_strategy`
    hla <- df$metricsData$alleles
    df$converted_hla_names <- unlist(lapply(hla, function(x) {
      if (grepl("HLA-", x)) {
        strsplit(x, "HLA-")[[1]][2]
      } else {
        x
      }
    }))
    return(df)
}

postprocess_inputs <- function(df) {
    if (!("Ref Match" %in% colnames(df$mainTable))) {
        df$mainTable$`Ref Match` <- "Not Run"
    }
    columns_needed <- c("ID", "Index", df$converted_hla_names, "Gene", "AA Change", "Num Passing Transcripts", "Best Peptide", "Best Transcript", "MANE Select", "Canonical", "TSL", "Allele", "Pos", "Prob Pos",
                        "Num Included Peptides", "Num Passing Peptides", "IC50 MT", "IC50 WT", "%ile MT", "%ile WT", "IC50 %ile MT", "IC50 %ile WT", "Pres %ile MT", "Pres %ile WT", "IM %ile MT", "IM %ile WT",
                        "RNA Expr", "RNA VAF", "Allele Expr", "RNA Depth", "DNA VAF", "Tier", "Ref Match", "Acpt", "Rej", "Rev")
    if ("ML Prediction (score)" %in% colnames(df$mainTable)) {
        columns_needed <- c(columns_needed, "ML Prediction (score)")
    }
    df$mainTable <- df$mainTable[, columns_needed]
    df$mainTable$`Gene of Interest` <- apply(df$mainTable, 1, function(x) {any(x["Gene"] == df$gene_list)})
    if ("Comments" %in% colnames(df$mainTable)) {
        df$comments <- data.frame(data = df$mainTable$`Comments`, nrow = nrow(df$mainTable), ncol = 1)
    }else {
        df$comments <- data.frame(matrix("No comments", nrow = nrow(df$mainTable)), ncol = 1)
    }
    rownames(df$comments) <- df$mainTable$ID
    return(df)
}

set_formatting_columns <- function(df) {
    df$mainTable$`Scaled BA` <- apply(df$mainTable, 1, function(x) scale_binding_affinity(df$allele_specific_binding_thresholds, df$use_allele_specific_binding_thresholds, df$binding_threshold, x["Allele"], x["IC50 MT"]))
    df$mainTable$`Scaled binding percentile` <- apply(df$mainTable, 1, function(x) {as.numeric(x["IC50 %ile MT"]) / (df$binding_percentile_threshold)})
    df$mainTable$`Scaled immunogenicity percentile` <- apply(df$mainTable, 1, function(x) {as.numeric(x["IM %ile MT"]) / (df$immunogenicity_percentile_threshold)})
    df$mainTable$`Scaled presentation percentile` <- apply(df$mainTable, 1, function(x) {as.numeric(x["Pres %ile MT"]) / (df$presentation_percentile_threshold)})
    df$mainTable$`Col RNA Expr` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["RNA Expr"]), 0, x["RNA Expr"])})
    df$mainTable$`Col RNA VAF` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["RNA VAF"]), 0, x["RNA VAF"])})
    df$mainTable$`Col Allele Expr` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["Allele Expr"]), 0, x["Allele Expr"])})
    df$mainTable$`Col RNA Depth` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["RNA Depth"]), 0, x["RNA Depth"])})
    df$mainTable$`Col DNA VAF` <- apply(df$mainTable, 1, function(x) {ifelse(is.na(x["DNA VAF"]), 0, x["DNA VAF"])})
    df$mainTable$`IC50 Pass` <- apply(df$mainTable, 1, function(x) {is_ic50_pass(df$use_allele_specific_binding_thresholds, x['Allele'], df$allele_specific_binding_thresholds, as.numeric(x['IC50 MT']), as.numeric(df$binding_threshold))})
    df$mainTable$`Binding Percentile Pass` <- apply(df$mainTable, 1, function(x) {is_percentile_pass(df$binding_percentile_threshold, as.numeric(x["IC50 %ile MT"]))})
    df$mainTable$`Immunogenicity Percentile Pass` <- apply(df$mainTable, 1, function(x) {is_percentile_pass(df$immunogenicity_percentile_threshold, as.numeric(x["IM %ile MT"]))})
    df$mainTable$`Presentation Percentile Pass` <- apply(df$mainTable, 1, function(x) {is_percentile_pass(df$presentation_percentile_threshold, as.numeric(x["Pres %ile MT"]))})
    df$mainTable$`Anchor Pass` <- apply(df$mainTable, 1, function(x) {is_anchor_residue_pass(df$anchor_mode, x['Best Peptide'], x['Allele'], as.numeric(df$anchor_contribution), x['Pos'], x['IC50 WT'], as.numeric(df$binding_threshold))})
    df$mainTable$`VAF Clonal Pass` <- apply(df$mainTable, 1, function(x) {is_vaf_clonal_pass(x["DNA VAF"], as.numeric(df$dna_cutoff))})
    df$mainTable$`Allele Expr Pass` <- apply(df$mainTable, 1, function(x) {is_allele_expr_pass(x["RNA VAF"], x["RNA Expr"], x["Allele Expr"], as.numeric(df$allele_expr))})
    df$mainTable$`RNA Expr Fail` <- apply(df$mainTable, 1, function(x) {!is.na(x['RNA Expr']) && as.numeric(x['RNA Expr']) == 0})
    df$mainTable$`RNA VAF Fail` <- apply(df$mainTable, 1, function(x) {!is.na(x['RNA VAF']) && as.numeric(x['RNA VAF']) <= as.numeric(df$metricsData['trna_vaf'])})
    df$mainTable$`RNA Depth Fail` <- apply(df$mainTable, 1, function(x) {!is.na(x['RNA Depth']) && as.numeric(x['RNA Depth']) <= as.numeric(df$metricsData['trna_cov'])})
    df$mainTable$`Prob Pos Pass` <- apply(df$mainTable, 1, function(x) {is_probaa_pass(x["Prob Pos"])})
    transcript_pass <- apply(df$mainTable, 1, function(x) {
      if ('tsl' %in% df$transcript_prioritization_strategy && is_tsl_pass(x["TSL"], as.numeric(df$maximum_transcript_support_level))) {
        return("True")
      }
      else if ('mane_select' %in% df$transcript_prioritization_strategy && is_mane_select_pass(x["MANE Select"])) {
        return("True")
      }
      else if ('canonical' %in% df$transcript_prioritization_strategy && is_canonical_pass(x["Canonical"])) {
        return("True")
      }
      else {
        return("False")
      }
    })
    df$mainTable <- add_column(df$mainTable, `Transcript Pass` = transcript_pass, .after = "TSL")
    return (df)
}

process_neofox_input <- function(mainData_neofox, df_neofox) {
    colnames(mainData_neofox) <- mainData_neofox[1, ]
    mainData_neofox <- mainData_neofox[-1, ]
    row.names(mainData_neofox) <- NULL

    rename_lookup <- c("PRIME_bestScore_allele" = "PRIME_best_allele", "PRIME_bestScore_peptide" = "PRIME_best_peptide", "PRIME_bestScore_rank" = "PRIME_best_rank", "PRIME_bestScore_score" = "PRIME_best_score")
    mainData_neofox <- mainData_neofox %>% rename(any_of(rename_lookup))
    mainData_neofox <- rename_with(mainData_neofox, ~ gsub("_", " ", .x, fixed = TRUE))

    # Columns that have been reviewed as most interesting
    columns_to_star <- c(
      "dnaVariantAlleleFrequency", "rnaExpression", "imputedGeneExpression",
      "rnaVariantAlleleFrequency", "NetMHCpan bestRank rank", "NetMHCpan bestAffinity affinity",
      "NetMHCpan bestAffinity affinityWT", "NetMHCpan bestRank rankWT", "PHBR I",
      "NetMHCIIpan bestRank rank", "NetMHCIIpan bestRank rankWT", "PHBR II", "Amplitude MHCI bestAffinity",
      "Pathogensimiliarity MHCI bestAffinity9mer", "DAI MHCI bestAffinity", "Tcell predictor",
      "Selfsimilarity MHCI", "Selfsimilarity MHCII", "IEDB Immunogenicity MHCI", "IEDB Immunogenicity MHCII",
      "MixMHCpred bestScore score", "MixMHCpred bestScore rank", "MixMHC2pred bestRank peptide",
      "MixMHC2pred bestRank rank", "Dissimilarity MHCI", "Dissimilarity MHCII", "Vaxrank bindingScore",
      "PRIME bestScore rank", "PRIME bestScore score"
    )

    # Check if each column is present in the dataframe and modify the names
    starred_column_names <- map(names(mainData_neofox), function(x) {
      if (x %in% columns_to_star) {
        paste0("*", x)
      } else {
        x
      }
    })
    names(mainData_neofox) <- starred_column_names
    df_neofox$mainTable_neofox <- mainData_neofox

    # Add scaling columns for coloring and barplots
    # There are no checks if user uploads data without one of these columns
    # Maybe an easy solution would be to just create a dummy column in the
    # for loop above for missing columns?
    df_neofox$mainTable_neofox$`Scaled NetMHCpan_bestAffinity` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(df_neofox$binding_threshold), as.numeric(x["*NetMHCpan bestAffinity affinity"]), as.numeric(x["*NetMHCpan bestAffinity affinity"]) / (df_neofox$binding_threshold))})
    df_neofox$mainTable_neofox$`Scaled NetMHCpan_bestAffinity_WT` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(df_neofox$binding_threshold), as.numeric(x["*NetMHCpan bestAffinity affinityWT"]), as.numeric(x["*NetMHCpan bestAffinity affinityWT"]) / (df_neofox$binding_threshold))})
    df_neofox$mainTable_neofox$`Scaled NetMHCpan_bestRank_rank` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(df_neofox$percentile_threshold), as.numeric(x["*NetMHCpan bestRank rank"]), as.numeric(x["*NetMHCpan bestRank rank"]) / (df_neofox$percentile_threshold))})
    df_neofox$mainTable_neofox$`Scaled NetMHCpan_bestRank_rankWT` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(df_neofox$percentile_threshold), as.numeric(x["*NetMHCpan bestRank rankWT"]), as.numeric(x["*NetMHCpan bestRank rankWT"]) / (df_neofox$percentile_threshold))})
    df_neofox$mainTable_neofox$`Scaled NetMHCIIpan_bestRank_rank` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(df_neofox$percentile_threshold), as.numeric(x["*NetMHCIIpan bestRank rank"]), as.numeric(x["*NetMHCIIpan bestRank rank"]) / (df_neofox$percentile_threshold))})
    df_neofox$mainTable_neofox$`Scaled NetMHCIIpan_bestRank_rankWT` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(df_neofox$percentile_threshold), as.numeric(x["*NetMHCIIpan bestRank rankWT"]), as.numeric(x["*NetMHCIIpan bestRank rankWT"]) / (df_neofox$percentile_threshold))})
    df_neofox$mainTable_neofox$`Scaled MixMHCpred_bestScore_rank` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(df_neofox$percentile_threshold), as.numeric(x["*MixMHCpred bestScore rank"]), as.numeric(x["*MixMHCpred bestScore rank"]) / (df_neofox$percentile_threshold))})
    df_neofox$mainTable_neofox$`Scaled MixMHC2pred_bestRank_rank` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(df_neofox$percentile_threshold), as.numeric(x["*MixMHC2pred bestRank rank"]), as.numeric(x["*MixMHC2pred bestRank rank"]) / (df_neofox$percentile_threshold))})
    df_neofox$mainTable_neofox$`Scaled PRIME_bestScore_rank` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(df_neofox$percentile_threshold), as.numeric(x["*PRIME bestScore rank"]), as.numeric(x["*PRIME bestScore rank"]) / (df_neofox$percentile_threshold))})
    # DAI is a measure of agrotopicity - so we want a a high DAI where the MT BA is low and the WT is BA is high, not sure if this is the correct scale
    df_neofox$mainTable_neofox$`Scaled DAI_MHCI_bestAffinity` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.null(1), as.numeric(x["*DAI MHCI bestAffinity"]), as.numeric(x["*DAI MHCI bestAffinity"]) / 10000)})

    df_neofox$mainTable_neofox$`Col DNA VAF` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.na(x["*dnaVariantAlleleFrequency"]), 0, x["*dnaVariantAlleleFrequency"])})
    df_neofox$mainTable_neofox$`Col RNA Expr` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.na(x["*rnaExpression"]), 0, x["*rnaExpression"])})
    df_neofox$mainTable_neofox$`Col Gene Expr` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.na(x["*imputedGeneExpression"]), 0, x["*imputedGeneExpression"])})
    df_neofox$mainTable_neofox$`Col RNA VAF` <- apply(df_neofox$mainTable_neofox, 1, function(x) {ifelse(is.na(x["*rnaVariantAlleleFrequency"]), 0, x["*rnaVariantAlleleFrequency"])})

    len <- nrow(df_neofox$mainTable_neofox)
    if ('Evaluation' %in% colnames(df_neofox$mainTable_neofox)) {
        setButtonStyling(df_neofox$mainTable_neofox$Evaluation, df_neofox$mainTable_neofox$ID)
    } else {
        df_neofox$mainTable_neofox["Evaluation"] = "Pending"
    }
    df_neofox$mainTable_neofox <- cbind(ID = rownames(df_neofox$mainTable_neofox), df_neofox$mainTable_neofox)
    df_neofox$evaluations <- df_neofox$mainTable_neofox[c("ID", "Evaluation")]
    df_neofox$mainTable_neofox$Evaluation <- NULL
    df_neofox$mainTable_neofox$Acpt <- shinyInputSelect(actionButton, df_neofox$mainTable_neofox["ID"], "button-neofox-acpt_", icon = icon("thumbs-up"), label = "", onclick = 'Shiny.onInputChange(\"accept_neofox_eval\", this.id, {priority: "event"})', onmousedown = "event.preventDefault(); event.stopPropagation();")
    df_neofox$mainTable_neofox$Rej <- shinyInputSelect(actionButton, df_neofox$mainTable_neofox["ID"], "button-neofox-rej_", icon = icon("thumbs-down"), label = "", onclick = 'Shiny.onInputChange(\"reject_neofox_eval\", this.id, {priority: "event"})', onmousedown = "event.preventDefault(); event.stopPropagation();")
    df_neofox$mainTable_neofox$Rev <- shinyInputSelect(actionButton, df_neofox$mainTable_neofox["ID"], "button-neofox-rev_", icon = icon("flag"), label = "", onclick = 'Shiny.onInputChange(\"review_neofox_eval\", this.id, {priority: "event"})', onmousedown = "event.preventDefault(); event.stopPropagation();")

    if ("Comments" %in% colnames(df_neofox$mainTable_neofox)) {
      df_neofox$comments <- data.frame(data = df_neofox$mainTable_neofox$`Comments`, nrow = nrow(df_neofox$mainTable_neofox), ncol = 1)
      df_neofox$mainTable_neofox$Comments <- NULL
    }else {
      df_neofox$comments <- data.frame(matrix("No comments", nrow = nrow(df_neofox$mainTable_neofox)), ncol = 1)
    }
    rownames(df_neofox$comments) <- df_neofox$mainTable_neofox$ID

    df_neofox$default_neofox_columns <- c("patientIdentifier", "gene", "mutatedXmer", "wildTypeXmer", "position", map(columns_to_star, function(x) { paste0("*", x) }), "Acpt", "Rej", "Rev")
    df_neofox$hidden_columns <- setdiff(colnames(df_neofox$mainTable_neofox), df_neofox$default_neofox_columns)

    return (df_neofox)
}
