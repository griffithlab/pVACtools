library(RCurl)
library(curl)
library(data.table)

## Load Anchor data
anchor_data <- list()
anchor_data[[8]] <- read.table(curl("https://raw.githubusercontent.com/griffithlab/pVACtools/ae938113ddbbe6c6eeecebf94459d449facd2c2f/tools/pvacview/data/Normalized_anchor_predictions_8_mer.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anchor_data[[9]] <- read.table(curl("https://raw.githubusercontent.com/griffithlab/pVACtools/ae938113ddbbe6c6eeecebf94459d449facd2c2f/tools/pvacview/data/Normalized_anchor_predictions_9_mer.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anchor_data[[10]] <- read.table(curl("https://raw.githubusercontent.com/griffithlab/pVACtools/ae938113ddbbe6c6eeecebf94459d449facd2c2f/tools/pvacview/data/Normalized_anchor_predictions_10_mer.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)
anchor_data[[11]] <- read.table(curl("https://raw.githubusercontent.com/griffithlab/pVACtools/ae938113ddbbe6c6eeecebf94459d449facd2c2f/tools/pvacview/data/Normalized_anchor_predictions_11_mer.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

#get binding affinity colors cutoffs given HLA

scale_binding_affinity <- function(allele_specific_binding_thresholds, use_allele_specific_binding_thresholds, binding_threshold, hla, current_ba) {
  if (use_allele_specific_binding_thresholds && hla %in% names(allele_specific_binding_thresholds[hla])) {
    threshold <- as.numeric(allele_specific_binding_thresholds[hla])
    return(as.numeric(current_ba) / threshold)
  }else {
    threshold <- as.numeric(binding_threshold)
    return(as.numeric(current_ba) / threshold)
  }
}

#for custom table formatting
get_group_inds <- function(reformat_data, group_inds){
  group_data <- as.data.frame(group_inds[[1]])
  group_data$group_id <- rownames(group_data)
  colnames(group_data)[colnames(group_data) == 'group_inds[[1]]'] <- 'inds'
  rownames(group_data) <- NULL
  return(group_data)
}

get_current_group_info <- function(peptide_features, metricsData, fullData, selectedRow){
  subset_data <- fullData[ ,peptide_features]
  inds <- metricsData[metricsData$group_id == selectedRow, ]$inds
  current_peptide_data <- data.frame(subset_data[unlist(inds), ])
  return(current_peptide_data)
}

#reformat table for display
table_formatting <- function(x, y) {
  y[y == "X"] <- NA
  peptide_ind <- grepl(x, colnames(y))
  peptide_columns <- y[, peptide_ind]
  peptide_columns$Mutant <- x
  colnames(peptide_columns) <- gsub(x, "", colnames(peptide_columns))
  colnames(peptide_columns) <- gsub("\\.", "", colnames(peptide_columns))
  peptide_columns_mt <- peptide_columns
  peptide_columns_mt$wt_peptide <- NULL
  ic50_mt <- dcast(peptide_columns_mt, Mutant ~ hla_types, value.var = "ic50s_MT")
  ic50_mt[, !names(ic50_mt) == "Mutant"] <- round(as.numeric(ic50_mt[, !names(ic50_mt) == "Mutant"]), 2)
  colnames(ic50_mt)[colnames(ic50_mt) == "Mutant"] <- "Peptide Sequence"
  ic50_mt <- add_column(ic50_mt, Type = "MT", .after = "Peptide Sequence")
  ic50_mt <- add_column(ic50_mt, `Problematic Positions` = peptide_columns$problematic_positions[[1]])
  ic50_mt <- add_column(ic50_mt, `Anchor Residue Fail` = peptide_columns$anchor_fails[[1]])
  peptide_columns_wt <- peptide_columns
  peptide_columns_wt$Mutant <- NULL
  ic50_wt <- dcast(peptide_columns_wt, wt_peptide ~ hla_types, value.var = "ic50s_WT")
  ic50_wt[, !names(ic50_wt) == "wt_peptide"] <- round(as.numeric(ic50_wt[, !names(ic50_wt) == "wt_peptide"]), 2)
  colnames(ic50_wt)[colnames(ic50_wt) == "wt_peptide"] <- "Peptide Sequence"
  ic50_wt <- add_column(ic50_wt, Type = "WT", .after = "Peptide Sequence")
  ic50_wt <- add_column(ic50_wt, `Problematic Positions` = "")
  ic50_wt <- add_column(ic50_wt, `Anchor Residue Fail` = "")
  combined_data <- rbind(ic50_mt, ic50_wt)
  combined_data$`Mutation Position` <- peptide_columns$mutation_position[[1]]
  reordered_data <- combined_data %>% select(-one_of("Problematic Positions"), -one_of("Anchor Residue Fail"), one_of("Problematic Positions"), one_of("Anchor Residue Fail"))
  reordered_data$`Has ProbPos` <- apply(reordered_data, 1, function(x) (x["Problematic Positions"] != "") & (x["Problematic Positions"] != "None"))
  reordered_data$`Has AnchorResidueFail` <- apply(reordered_data, 1, function(x) (x["Anchor Residue Fail"] != "") & (x["Anchor Residue Fail"] != "None"))
  reordered_data
}
#generate peptide coloring for hla allele
peptide_coloring <- function(hla_allele, peptide_row) {
  peptide_length <- as.numeric(peptide_row["length"])
  if (peptide_length < 8) {
    return(c("#999999"))
  }
  position <- as.numeric(peptide_row["x_pos"])
  anchor_score <- as.numeric(anchor_data[[peptide_length]][anchor_data[[peptide_length]]["HLA"] == hla_allele][2:(peptide_length + 1)])
  value_bins <- cut(anchor_score, breaks = seq(0, 1, len = 100),
                    include.lowest = TRUE)
  colors <- colorRampPalette(c("lightblue", "blue"))(99)[value_bins]
  return(colors[[position]])
}
#calculate per-length anchor score for HLA allele
anchor_weights_for_alleles <- function(hla_alleles) {
  scores_df <- data.frame()
  for (hla_allele in hla_alleles) {
    if (any(hla_allele == anchor_data[[8]]["HLA"])) {
      eight_mer_scores <- append(anchor_data[[8]][anchor_data[[8]]["HLA"] == hla_allele][1:(8 + 1)], "8mer", 1)
    }
    else {
      eight_mer_scores <- c(hla_allele, "8mer")
    }
    if (any(hla_allele == anchor_data[[9]]["HLA"])) {
      nine_mer_scores <- append(anchor_data[[9]][anchor_data[[9]]["HLA"] == hla_allele][1:(9 + 1)], "9mer", 1)
    }
    else {
      nine_mer_scores <- c(hla_allele, "9mer")
    }
    if (any(hla_allele == anchor_data[[10]]["HLA"])) {
      ten_mer_scores <- append(anchor_data[[10]][anchor_data[[10]]["HLA"] == hla_allele][1:(10 + 1)], "10mer", 1)
    }
    else {
      ten_mer_scores <- c(hla_allele, "10mer")
    }
    if (any(hla_allele == anchor_data[[11]]["HLA"])) {
      eleven_mer_scores <- append(anchor_data[[11]][anchor_data[[11]]["HLA"] == hla_allele][1:(11 + 1)], "11mer", 1)
    }
    else {
      eleven_mer_scores <- c(hla_allele, "11mer")
    }
    scores <- list(eight_mer_scores, nine_mer_scores, ten_mer_scores, eleven_mer_scores)
    scores <- lapply(scores, `length<-`, 13)
    scores <- transpose(data.frame(scores))
    colnames(scores) <- c("HLA Allele", "Peptide Length", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11")
    scores_df <- rbind(scores_df, scores)
  }
  return(scores_df)
}
#calculate anchor list for specific peptide length and HLA allele combo given contribution cutoff
calculate_anchor <- function(hla_allele, peptide_length, anchor_contribution) {
  result <- tryCatch({
    anchor_raw_data <- as.numeric(anchor_data[[peptide_length]][anchor_data[[peptide_length]]["HLA"] == hla_allele][2:(peptide_length + 1)])
    if (any(is.na(anchor_raw_data))) {
      return("NA")
    }
    names(anchor_raw_data) <- as.character(1:length(anchor_raw_data))
    anchor_raw_data <- anchor_raw_data[order(unlist(anchor_raw_data), decreasing = TRUE)]
    count <- 0
    anchor_list <- list()
    for (i in 1:length(anchor_raw_data)) {
      if (count >= anchor_contribution) {
        return(anchor_list)
      }else {
        count <- count + anchor_raw_data[[i]]
        anchor_list <- append(anchor_list, names(anchor_raw_data[i]))
      }
    }
    return(anchor_list)
  }, error = function(e) { return("NA") })
}

#converts string range (e.g. '2-4', '6') to associated list
range_str_to_seq <- function(mutation_position) {
  rnge <- strsplit(mutation_position, "-")[[1]]
  if (length(rnge) == 2) {
    return(seq(rnge[1], rnge[2]))
  }else {
    return(c(strtoi(rnge[1])))
  }
  return(0)
}

#get data from metrics file associated with peptide if available
get_mt_peptide_data <- function(metrics_data_row, mt_peptide) {
  for (trn in metrics_data_row$sets) {
    res <- metrics_data_row$good_binders[[trn]]$peptides[[mt_peptide]]
    if (!is.null(res)) {
      return(res)
    }
  }
  return(c())
}

#calculate the positions different between MT and WT peptide
calculate_mutation_info <- function(metrics_data_row) {
  wt_peptide <- metrics_data_row$best_peptide_wt
  if (is.na(wt_peptide)) {
    return(0)
  }
  mt_peptide <- metrics_data_row$best_peptide_mt
  mt_data <- get_mt_peptide_data(metrics_data_row, mt_peptide)
  # if recorded mutation_position range, use it
  if (length(mt_data) > 0) {
    diff_positions <- range_str_to_seq(mt_data$"mutation_position")
  }else {
    split_positions <- strsplit(c(wt_peptide, mt_peptide), split = "")
    diff_positions <- which(split_positions[[1]] != split_positions[[2]])
  }
  return(diff_positions)
}
##Generate Tiering for given variant with specific cutoffs
tier <- function(variant_info, anchor_contribution, dna_cutoff, allele_expr_cutoff, mutation_pos_list, hla_allele, tsl, meta_data, anchor_mode, use_allele_specific_binding_thresholds, binding_threshold) {
  mt_binding <- as.numeric(variant_info["IC50 MT"])
  wt_binding <- as.numeric(variant_info["IC50 WT"])
  mt_percent <- as.numeric(variant_info["%ile MT"])
  wt_percent <- as.numeric(variant_info["%ile WT"])
  gene_expr <- as.numeric(variant_info["RNA Expr"])
  dna_vaf <- as.numeric(variant_info["DNA VAF"])
  rna_vaf <- as.numeric(variant_info["RNA VAF"])
  rna_depth <- as.numeric(variant_info["RNA Depth"])
  allele_expr <- as.numeric(variant_info["Allele Expr"])
  if (use_allele_specific_binding_thresholds && hla_allele %in% names(meta_data[["allele_specific_binding_thresholds"]][hla_allele])) {
    binding_threshold <- as.numeric(meta_data[["allele_specific_binding_thresholds"]][hla_allele])
  }
  trna_vaf <- as.numeric(meta_data["trna_vaf"])
  trna_cov <- as.numeric(meta_data["trna_cov"])
  percentile_filter <- FALSE
  percentile_threshold <- NULL
  if (!is.null(meta_data[["percentile_threshold"]])) {
    percentile_threshold <- as.numeric(meta_data[["percentile_threshold"]])
    percentile_filter <- TRUE
  }
  tsl_max <- as.numeric(meta_data["maximum_transcript_support_level"])
  mutation_pos_list <- mutation_pos_list[["Pos"]]
  if (anchor_mode == "default") {
    anchor_list <- c(1, 2, nchar(variant_info[["Best Peptide"]]), nchar(variant_info[["Best Peptide"]]) - 1)
  }else {
    anchor_list <- unlist(calculate_anchor(hla_allele, length(unlist(strsplit(variant_info["Best Peptide"][[1]], split = ""))), anchor_contribution))
    if (anchor_list[[1]] == "NA") {
      anchor_list <- c(1, 2, nchar(variant_info[["Best Peptide"]]), nchar(variant_info[["Best Peptide"]]) - 1)
    }
  }
  anchor_residue_pass <- TRUE
  # if all of mutated positions in anchors
  if (grepl("-", mutation_pos_list, fixed = TRUE)) {
    range_start <- as.numeric(strsplit(mutation_pos_list, "-")[[1]][1])
    if (range_start == 0) {
      range_start <- 1
    }
    range_stop <- as.numeric(strsplit(mutation_pos_list, "-")[[1]][2])
    mutation_pos_list <- c(range_start:range_stop)
    if (all(mutation_pos_list %in% anchor_list)) {
      if (is.na(wt_binding)) {
        anchor_residue_pass <- FALSE
      }else if (wt_binding < binding_threshold) {
        anchor_residue_pass <- FALSE
      }
    }
  }else if (!is.na(mutation_pos_list)) {
    if (all(as.numeric(mutation_pos_list) %in% anchor_list)) {
      if (is.na(wt_binding)) {
        anchor_residue_pass <- FALSE
      }else if (wt_binding < binding_threshold) {
        anchor_residue_pass <- FALSE
      }
    }
  }
  tsl_pass <- TRUE
  if ((tsl == "Not Supported")) {
    tsl_pass <- TRUE
  }
  else if ((tsl == "NA") || as.numeric(tsl) > tsl_max) {
    tsl_pass <- FALSE
  }
  allele_expr_pass <- TRUE
  if (!is.na(rna_vaf) && !is.na(gene_expr) && allele_expr <= allele_expr_cutoff) {
    allele_expr_pass <- FALSE
  }
  vaf_clonal_pass <- TRUE
  if (!is.na(dna_vaf) && dna_vaf < dna_cutoff / 2) {
    vaf_clonal_pass <- FALSE
  }
  ## Assign Tiering
  if ((mt_binding < binding_threshold) && allele_expr_pass && vaf_clonal_pass && tsl_pass && anchor_residue_pass) {
    if (percentile_filter) {
      if (mt_percent <= percentile_threshold) {
        return("Pass")
      }
    }else {
      return("Pass")
    }
  }
  
  if ((mt_binding < binding_threshold) && allele_expr_pass && vaf_clonal_pass && tsl_pass && !anchor_residue_pass) {
    if (percentile_filter) {
      if (mt_percent <= percentile_threshold) {
        return("Anchor")
      }
    }else {
      return("Anchor")
    }
  }
  if ((mt_binding < binding_threshold) && allele_expr_pass && !vaf_clonal_pass && tsl_pass && anchor_residue_pass) {
    if (percentile_filter) {
      if (mt_percent <= percentile_threshold) {
        return("Subclonal")
      }
    }else {
      return("Subclonal")
    }
  }
  lowexpr <- FALSE
  if (!is.na(rna_vaf) && !is.na(gene_expr) && !is.na(rna_depth)) {
    if ((allele_expr > 0) || ((gene_expr == 0) && (rna_depth > trna_cov) && (rna_vaf > trna_vaf))) {
      lowexpr <- TRUE
    }
  }
  if ((mt_binding < binding_threshold) && lowexpr && vaf_clonal_pass && tsl_pass && anchor_residue_pass) {
    if (percentile_filter) {
      if (mt_percent <= percentile_threshold) {
        return("LowExpr")
      }
    }else {
      return("LowExpr")
    }
  }
  if (!is.na(allele_expr) && ((gene_expr == 0) || (rna_vaf == 0)) && !lowexpr) {
    return("NoExpr")
  }
  return("Poor")
}
#Determine the Tier Count for given variant with specific cutoffs
tier_numbers <- function(variant_info, anchor_contribution, dna_cutoff, allele_expr_cutoff, mutation_pos_list, hla_allele, tsl, meta_data, anchor_mode, allele_specific_binding_thresholds, use_allele_specific_binding_thresholds, binding_threshold) {
  mt_binding <- as.numeric(variant_info["IC50 MT"])
  wt_binding <- as.numeric(variant_info["IC50 WT"])
  mt_percent <- as.numeric(variant_info["%ile MT"])
  wt_percent <- as.numeric(variant_info["%ile WT"])
  gene_expr <- as.numeric(variant_info["RNA Expr"])
  dna_vaf <- as.numeric(variant_info["DNA VAF"])
  rna_vaf <- as.numeric(variant_info["RNA VAF"])
  rna_depth <- as.numeric(variant_info["RNA Depth"])
  allele_expr <- as.numeric(variant_info["Allele Expr"])
  if (use_allele_specific_binding_thresholds && hla_allele %in% names(meta_data[["allele_specific_binding_thresholds"]][hla_allele])) {
    binding_threshold <- as.numeric(meta_data[["allele_specific_binding_thresholds"]][hla_allele])
  }
  trna_vaf <- as.numeric(meta_data["trna_vaf"])
  trna_cov <- as.numeric(meta_data["trna_cov"])
  percentile_filter <- FALSE
  percentile_threshold <- NULL
  if (!is.null(meta_data[["percentile_threshold"]])) {
    percentile_threshold <- as.numeric(meta_data[["percentile_threshold"]])
    percentile_filter <- TRUE
  }
  tsl_max <- as.numeric(meta_data["maximum_transcript_support_level"])
  mutation_pos_list <- mutation_pos_list[["Pos"]]
  if (anchor_mode == "default") {
    anchor_list <- c(1, 2, nchar(variant_info[["Best Peptide"]]), nchar(variant_info[["Best Peptide"]]) - 1)
  }else {
    anchor_list <- unlist(calculate_anchor(hla_allele, length(unlist(strsplit(variant_info["Best Peptide"][[1]], split = ""))), anchor_contribution))
    if (anchor_list[[1]] == "NA") {
      anchor_list <- c(1, 2, nchar(variant_info[["Best Peptide"]]), nchar(variant_info[["Best Peptide"]]) - 1)
    }
  }
  anchor_residue_pass <- TRUE
  # if all of mutated positions in anchors
  if (grepl("-", mutation_pos_list, fixed = TRUE)) {
    range_start <- as.numeric(strsplit(mutation_pos_list, "-")[[1]][1])
    if (range_start == 0) {
      range_start <- 1
    }
    range_stop <- as.numeric(strsplit(mutation_pos_list, "-")[[1]][2])
    mutation_pos_list <- c(range_start:range_stop)
    if (all(mutation_pos_list %in% anchor_list)) {
      if (is.na(wt_binding)) {
        anchor_residue_pass <- FALSE
      }else if (wt_binding < binding_threshold) {
        anchor_residue_pass <- FALSE
      }
    }
  }else if (!is.na(mutation_pos_list)) {
    if (all(as.numeric(mutation_pos_list) %in% anchor_list)) {
      if (is.na(wt_binding)) {
        anchor_residue_pass <- FALSE
      }else if (wt_binding < binding_threshold) {
        anchor_residue_pass <- FALSE
      }else if (!is.null(percentile_threshold) && (wt_percent) < percentile_threshold) {
        anchor_residue_pass <- FALSE
      }
    }
  }
  tsl_pass <- TRUE
  if ((tsl == "Not Supported")) {
    tsl_pass <- TRUE
  }
  else if ((tsl == "NA") || as.numeric(tsl) > tsl_max) {
    tsl_pass <- FALSE
  }
  allele_expr_pass <- TRUE
  if (!is.na(rna_vaf) && !is.na(gene_expr) && allele_expr <= allele_expr_cutoff) {
    allele_expr_pass <- FALSE
  }
  vaf_clonal_pass <- TRUE
  if (!is.na(dna_vaf) && dna_vaf < dna_cutoff / 2) {
    vaf_clonal_pass <- FALSE
  }
  ## Pass
  if ((mt_binding < binding_threshold) && allele_expr_pass && vaf_clonal_pass && tsl_pass && anchor_residue_pass) {
    if (percentile_filter) {
      if (mt_percent <= percentile_threshold) {
        return(1)
      }
    }else {
      return(1)
    }
  }
  ## Anchor
  if ((mt_binding < binding_threshold) && allele_expr_pass && vaf_clonal_pass && tsl_pass && !anchor_residue_pass) {
    if (percentile_filter) {
      if (mt_percent <= percentile_threshold) {
        return(5)
      }
    }else {
      return(5)
    }
  }
  if ((mt_binding < binding_threshold) && allele_expr_pass && !vaf_clonal_pass && tsl_pass && anchor_residue_pass) {
    if (percentile_filter) {
      if (mt_percent <= percentile_threshold) {
        return(6)
      }
    }else {
      return(6)
    }
  }
  lowexpr <- FALSE
  if (!is.na(rna_vaf) && !is.na(gene_expr) && !is.na(rna_depth)) {
    if ((allele_expr > 0) || ((gene_expr == 0) && (rna_depth > trna_cov) && (rna_vaf > trna_vaf))) {
      lowexpr <- TRUE
    }
  }
  if ((mt_binding < binding_threshold) && (lowexpr) && vaf_clonal_pass && tsl_pass && anchor_residue_pass) {
    if (percentile_filter) {
      if (mt_percent <= percentile_threshold) {
        if (allele_expr > 0) {
          return(7)
        }else if ((gene_expr == 0) && (rna_depth > trna_cov) && (rna_vaf > trna_vaf)) {
          return(8)
        }
      }
    }else {
      if (allele_expr > 0) {
        return(7)
      }else if ((gene_expr == 0) && (rna_depth > trna_cov) && (rna_vaf > trna_vaf)) {
        return(8)
      }
    }
  }
  if (!is.na(allele_expr) && ((gene_expr == 0) || (rna_vaf == 0)) && !lowexpr) {
    if ((gene_expr == 0) && (rna_vaf != 0)) {
      return(9)
    }else if ((gene_expr != 0) && (rna_vaf == 0)) {
      return(10)
    }else {
      return(11)
    }
  }
  count <- 12
  if (!anchor_residue_pass) {
    count <- count + 1
  }
  if (!vaf_clonal_pass) {
    count <- count + 2
  }
  if (!is.na(gene_expr) && !is.na(rna_depth) && !is.na(rna_vaf) && (gene_expr == 0) && (rna_depth > trna_cov) && (rna_vaf > trna_vaf)) {
    count <- count + 4
  }
  if (!is.na(allele_expr) && allele_expr > 0 && allele_expr < allele_expr_cutoff) {
    count <- count + 8
  }
  return(count)
}
