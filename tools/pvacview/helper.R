table_formatting = function(x,y){
  peptide_ind <- grepl(x, colnames(y))
  peptide_columns <- y[,peptide_ind]
  peptide_columns$Mutant <- x
  colnames(peptide_columns) <- gsub(x, "", colnames(peptide_columns))
  colnames(peptide_columns) <- gsub("\\.", "", colnames(peptide_columns))
  peptide_columns_mt <- peptide_columns
  peptide_columns_mt$wt_peptide <- NULL
  ic50_mt <- dcast(peptide_columns_mt, Mutant ~ hla_types, value.var = "ic50s_MT")
  colnames(ic50_mt)[colnames(ic50_mt) == "Mutant"] <- "Peptide Sequence"
  ic50_mt <- add_column(ic50_mt, Type = "MT", .after="Peptide Sequence") 
  peptide_columns_wt <- peptide_columns
  peptide_columns_wt$Mutant <- NULL
  ic50_wt <- dcast(peptide_columns_wt, wt_peptide ~ hla_types, value.var = "ic50s_WT")
  colnames(ic50_wt)[colnames(ic50_wt) == "wt_peptide"] <- "Peptide Sequence"
  ic50_wt <- add_column(ic50_wt, Type = "WT", .after="Peptide Sequence") 
  combined_data <- rbind(ic50_mt, ic50_wt)
  combined_data$`Mutation Position` <- peptide_columns$mutation_position[[1]]
  combined_data
}

