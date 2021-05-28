## Load Anchor data
anchor_data = list()
anchor_data[[8]] <- read.table("~/Desktop/Griffith_Lab/R_shiny_visualization/neoantigen_visualization/data/Anchor_scores/Normalized_anchor_predictions_8_mer.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
anchor_data[[9]] <- read.table("~/Desktop/Griffith_Lab/R_shiny_visualization/neoantigen_visualization/data/Anchor_scores/Normalized_anchor_predictions_9_mer.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
anchor_data[[10]] <- read.table("~/Desktop/Griffith_Lab/R_shiny_visualization/neoantigen_visualization/data/Anchor_scores/Normalized_anchor_predictions_10_mer.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)
anchor_data[[11]] <- read.table("~/Desktop/Griffith_Lab/R_shiny_visualization/neoantigen_visualization/data/Anchor_scores/Normalized_anchor_predictions_11_mer.tsv", sep = '\t', header = TRUE, stringsAsFactors = FALSE)

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
