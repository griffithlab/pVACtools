library(RCurl)
library(curl)

## Load Anchor data
anchor_data = list()
anchor_data[[8]] <- read.table(curl("https://raw.githubusercontent.com/griffithlab/pVACtools/ae938113ddbbe6c6eeecebf94459d449facd2c2f/tools/pvacview/data/Normalized_anchor_predictions_8_mer.tsv"), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
anchor_data[[9]] <- read.table(curl("https://raw.githubusercontent.com/griffithlab/pVACtools/ae938113ddbbe6c6eeecebf94459d449facd2c2f/tools/pvacview/data/Normalized_anchor_predictions_9_mer.tsv"), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
anchor_data[[10]] <- read.table(curl("https://raw.githubusercontent.com/griffithlab/pVACtools/ae938113ddbbe6c6eeecebf94459d449facd2c2f/tools/pvacview/data/Normalized_anchor_predictions_10_mer.tsv"), sep = '\t', header = TRUE, stringsAsFactors = FALSE)
anchor_data[[11]] <- read.table(curl("https://raw.githubusercontent.com/griffithlab/pVACtools/ae938113ddbbe6c6eeecebf94459d449facd2c2f/tools/pvacview/data/Normalized_anchor_predictions_11_mer.tsv"), sep = '\t', header = TRUE, stringsAsFactors = FALSE)


#reformat table for display
table_formatting = function(x,y){
  y[y == 'X'] <- NA
  peptide_ind <- grepl(x, colnames(y))
  peptide_columns <- y[,peptide_ind]
  peptide_columns$Mutant <- x
  colnames(peptide_columns) <- gsub(x, "", colnames(peptide_columns))
  colnames(peptide_columns) <- gsub("\\.", "", colnames(peptide_columns))
  peptide_columns_mt <- peptide_columns
  peptide_columns_mt$wt_peptide <- NULL
  ic50_mt <- dcast(peptide_columns_mt, Mutant ~ hla_types, value.var = "ic50s_MT")
  ic50_mt[, !names(ic50_mt) == 'Mutant'] <- round(as.numeric(ic50_mt[, !names(ic50_mt) == 'Mutant']),2)
  colnames(ic50_mt)[colnames(ic50_mt) == "Mutant"] <- "Peptide Sequence"
  ic50_mt <- add_column(ic50_mt, Type = "MT", .after="Peptide Sequence") 
  peptide_columns_wt <- peptide_columns
  peptide_columns_wt$Mutant <- NULL
  ic50_wt <- dcast(peptide_columns_wt, wt_peptide ~ hla_types, value.var = "ic50s_WT")
  ic50_wt[, !names(ic50_wt) == 'wt_peptide'] <- round(as.numeric(ic50_wt[, !names(ic50_wt) == 'wt_peptide']),2)
  colnames(ic50_wt)[colnames(ic50_wt) == "wt_peptide"] <- "Peptide Sequence"
  ic50_wt <- add_column(ic50_wt, Type = "WT", .after="Peptide Sequence") 
  combined_data <- rbind(ic50_mt, ic50_wt)
  combined_data$`Mutation Position` <- peptide_columns$mutation_position[[1]]
  combined_data
}
#generate peptide coloring for hla allele 
peptide_coloring <- function(hla_allele, peptide_row){
  peptide_length <- as.numeric(peptide_row['length'])
  if (peptide_length < 8){
    return(c("#999999"))
  }
  position <- as.numeric(peptide_row['x_pos'])
  anchor_score <- as.numeric(anchor_data[[peptide_length]][anchor_data[[peptide_length]]['HLA'] == hla_allele][2:(peptide_length+1)])
  value_bins <- cut(anchor_score, breaks = seq(0, 1, len = 100), 
                    include.lowest = TRUE)
  colors <- colorRampPalette(c("lightblue", "blue"))(99)[value_bins]
  return (colors[[position]])
}
#calculate anchor list for specific peptide length and HLA allele combo given contribution cutoff
calculate_anchor <- function(hla_allele, peptide_length, anchor_contribution){
  #browser()
  anchor_raw_data <- as.numeric(anchor_data[[peptide_length]][anchor_data[[peptide_length]]['HLA'] == hla_allele][2:(peptide_length+1)])
  if (any(is.na(anchor_raw_data))) {
    return("NA")
  }
  names(anchor_raw_data) <- as.character(1:length(anchor_raw_data))
  anchor_raw_data <- anchor_raw_data[order(unlist(anchor_raw_data), decreasing = TRUE)]
  count <- 0
  anchor_list <- list()
  for (i in 1:length(anchor_raw_data)){
    if (count >= anchor_contribution) {
      return (anchor_list)
    }
    else{
      count <- count + anchor_raw_data[[i]]
      anchor_list <- append(anchor_list,names(anchor_raw_data[i]))
    }
  }
}

#calculate the positions different between MT and WT peptide
calculate_mutation_info <- function(metrics_data_row){
  wt_peptide <- metrics_data_row$best_peptide_wt
  if (is.na(wt_peptide)){
    return (0)
  }
  mt_peptide <- metrics_data_row$best_peptide_mt
  split_positions <- strsplit(c(wt_peptide, mt_peptide), split = "")
  diff_positions <- which(split_positions[[1]] != split_positions[[2]])
  return (diff_positions)
}

##Generate Tiering for given variant with specific cutoffs
tier <- function(variant_info, anchor_contribution, dna_cutoff, allele_expr_high, allele_expr_low, mutation_pos_list, hla_allele, anchor_mode = "allele-specific"){
  mt_binding <- as.numeric(variant_info['IC50 MT'])
  wt_binding <- as.numeric(variant_info['IC50 WT'])
  gene_expr <- as.numeric(variant_info['RNA Expr'])
  dna_vaf <- as.numeric(variant_info['DNA VAF'])
  rna_vaf <- as.numeric(variant_info['RNA VAF'])
  rna_depth <- as.numeric(variant_info['RNA Depth'])
  allele_expr <- as.numeric(variant_info['Allele Expr'])
  if (anchor_mode == "default"){
    anchor_list <- c(1,2,nchar(variant_info['Best Peptide']), nchar(variant_info['Best Peptide'])-1)
  }
  else{
    anchor_list <- unlist(calculate_anchor(hla_allele, length(unlist(strsplit(variant_info['Best Peptide'][[1]], split = ""))), anchor_contribution))
    if (anchor_list[[1]] == 'NA'){
      anchor_list <- c(1,2,nchar(variant_info['Best Peptide']), nchar(variant_info['Best Peptide'])-1)
    }
  }
  anchor_residue_pass <- TRUE
  if (all(as.numeric(mutation_pos_list) %in% anchor_list)){
    if (is.na(wt_binding)){
      anchor_residue_pass <- FALSE
    }
    else if (wt_binding < 1000) {
      anchor_residue_pass <- FALSE
    }
  }
  if ((mt_binding < 500) & (allele_expr > allele_expr_high) & (dna_vaf >= dna_cutoff/2) & anchor_residue_pass){
    return ("Pass")
  }
  if ((mt_binding < 1000) & (allele_expr > allele_expr_low) & (dna_vaf >= dna_cutoff/2) & anchor_residue_pass){
    return ("Relaxed")
  }
  if ((mt_binding < 1000) & (allele_expr > allele_expr_low) & (dna_vaf >= dna_cutoff/2) & !anchor_residue_pass){
    return ("Anchor")
  }
  if ((mt_binding < 1000) & (allele_expr > allele_expr_low) & (dna_vaf < dna_cutoff/2) & anchor_residue_pass){
    return ("Subclonal")
  }
  lowexpr <- FALSE
  if ((allele_expr > 0) | ((gene_expr == 0) & (rna_depth > 50) & (rna_vaf > 0.10))) {
    lowexpr <- TRUE
  }
  if ((mt_binding < 1000) & (lowexpr) & (dna_vaf >= dna_cutoff/2) & anchor_residue_pass){
    return ("LowExpr")
  }
  if (((gene_expr == 0) | (rna_vaf == 0)) & !lowexpr){
    return ("NoExpr")
  }
  
  return ("Poor")
}


tier_numbers <- function(variant_info, anchor_contribution, dna_cutoff, allele_expr_high, allele_expr_low, mutation_pos_list, hla_allele = NULL, anchor_mode = "allele-specific"){
  mt_binding <- as.numeric(variant_info['IC50 MT'])
  wt_binding <- as.numeric(variant_info['IC50 WT'])
  gene_expr <- as.numeric(variant_info['RNA Expr'])
  dna_vaf <- as.numeric(variant_info['DNA VAF'])
  rna_vaf <- as.numeric(variant_info['RNA VAF'])
  rna_depth <- as.numeric(variant_info['RNA Depth'])
  allele_expr <- as.numeric(variant_info['Allele Expr'])
  count = 12
  if (anchor_mode == "default"){
    anchor_list <- c(1,2,nchar(variant_info['Best Peptide']), nchar(variant_info['Best Peptide'])-1)
  }
  else{
    anchor_list <- unlist(calculate_anchor(hla_allele, length(unlist(strsplit(variant_info['Best Peptide'][[1]], split = ""))), anchor_contribution))
    if (anchor_list[[1]] == 'NA'){
      anchor_list <- c(1,2,nchar(variant_info['Best Peptide']), nchar(variant_info['Best Peptide'])-1)
    }
  }
  anchor_residue_pass <- TRUE
  if (all(as.numeric(mutation_pos_list) %in% anchor_list)){
    if (is.na(wt_binding)){
      anchor_residue_pass <- FALSE
    }
    else if (wt_binding < 1000) {
      anchor_residue_pass <- FALSE
    }
  }
  ## Pass
  if ((mt_binding < 500) & (allele_expr > allele_expr_high) & (dna_vaf >= dna_cutoff/2) & anchor_residue_pass){
    return (1)
  }
  ## Relaxed
  if ((mt_binding < 1000) & (allele_expr > allele_expr_low) & (dna_vaf >= dna_cutoff/2) & anchor_residue_pass){
    if ((mt_binding < 500) & (allele_expr < allele_expr_high)){
      return (2) 
    }
    else if ((mt_binding < 1000) & (allele_expr < allele_expr_high)){
      return (3)
    }
    else{
      return(4)
    }
  }
  ## Anchor
  if ((mt_binding < 1000) & (allele_expr > allele_expr_low) & (dna_vaf >= dna_cutoff/2) & !anchor_residue_pass){
    return (5)
  }
  if ((mt_binding < 1000) & (allele_expr > allele_expr_low) & (dna_vaf < dna_cutoff/2) & anchor_residue_pass){
    return (6)
  }
  lowexpr <- FALSE
  if ((allele_expr > 0) | ((gene_expr == 0) & (rna_depth > 50) & (rna_vaf > 0.10))) {
    lowexpr <- TRUE
  }
  if ((mt_binding < 1000) & (lowexpr) & (dna_vaf >= dna_cutoff/2) & anchor_residue_pass){
    if (allele_expr > 0){
      return (7)
    }
    else if ((gene_expr == 0) & (rna_depth > 50) & (rna_vaf > 0.10)){
      return (8)
    }
  }
  if (((gene_expr == 0) | (rna_vaf == 0)) & !lowexpr){
    if ((gene_expr == 0) & (rna_vaf != 0)){
      return (9)
    }
    else if ((gene_expr != 0) & (rna_vaf == 0)){
      return (10)
    }
    else{
      return (11)
    }
  }
  if (!anchor_residue_pass){
    count = count + 1
  }
  if (dna_vaf < dna_cutoff/2){
    count = count + 2
  }
  if ((gene_expr == 0) & (rna_depth > 50) & (rna_vaf > 0.10)){
    count = count + 4
  }
  if (allele_expr > 0 & allele_expr < allele_expr_low){
    count = count + 8
  }
  return (count)
}










