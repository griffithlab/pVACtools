## server side callback functions
rowcallback <- function(hla_count, row_num) {
  c(
    "function(row, data, displayNum, displayIndex){",
    gsub("0", row_num, "  if (displayIndex == 0){"),
    "  $('td',row).css('border-top','3px solid #0390fc');",
    "  $('td',row).css('border-bottom','3px solid #0390fc');",
    "  }",
    "}")
}

callback <- function(hla_count, score_mode) {
  c(
    "var tips = ['Gene - The Ensembl gene name of the affected gene.',",
    "        'AA Change - The amino acid change for the mutation. Note that FS indicates a frameshift variant.',",
    "        'Num Passing Transcripts - The number of transcripts for this mutation that resulted in at least one well-binding peptide.',",
    "        'Best Peptide - The best mutant epitope sequence, taking into account various criteria such as mutant binding affinity, epitopes arising from a protein_coding transcript, the MANE Select transcript, the Canonical transcript, or a transcript with TSL below the maximum transcript support level (depending on the chosen transcript prioritization strategy), no problematic positions, passes the anchor evaluation.',",
    "        'Best Transcript - Transcript corresponding to the best peptide with the lowest TSL and shortest length.',",
    "        'MANE Select - MANE select status of the best transcript.',",
    "        'Canonical - Canonical status of the best transcript.',",
    "        'TSL - Transcript support level of the best transcript.',",
    "        'Transcript Pass - Reflects whether the transcript giving rise to the Best Peptide passes the transcript evaluation, i.e., the transcript is either the MANE Select transcript, the Canonical transcript or has TSL below the maximum transcript support level. Criteria to be evaluated depend on the selected transcript prioritization strategy.',",
    "        'Allele',",
    "        'Pos - A list of the mutated positions (one-based) in the Best Peptide compared to its matched wild type peptide. NA if there is no matched wild type.',",
    "        'Prob Pos - Problematic positions within the Best Peptide.',",
    "        'Num Included Peptides - The number of top-scoring, unique peptides included for review.',",
    "        'Num Passing Peptides - The number of unique well-binding peptides for this mutation.',",
    paste("      'IC50 MT - ", score_mode, "IC50 binding affinity of the Best Peptide across all binding affinity prediction algorithms used.', "),
    "        'IC50 WT - IC50 binding affinity of the corresponding wild type peptide.',",
    paste("      '%ile MT - ", score_mode, "combined percentile rank of the Best Peptide across all prediction algorithms used (those that provide percentile output).', "),
    "        '%ile WT - Combined percentile rank of the corresponding wild type epitope across all prediction algorithms used (those that provide percentile output).', ",
    paste("      'IC50 %ile MT - ", score_mode, "binding percentile rank of the Best Peptide across all binding prediction algorithms used (those that provide percentile output).', "),
    "        'IC50 %ile WT - binding percentile rank of the corresponding wild type epitope across all binding prediction algorithms used (those that provide percentile output).', ",
    paste("      'Pres %ile MT - ", score_mode, "presentation percentile rank of the Best Peptide across all presentation prediction algorithms used (those that provide percentile output).', "),
    "        'Pres %ile WT - presentation percentile rank of the corresponding wild type epitope across all presentation prediction algorithms used (those that provide percentile output).', ",
    paste("      'IM %ile MT - ", score_mode, "immunogenicity percentile rank of the Best Peptide across all immunogenicity prediction algorithms used (those that provide percentile output).', "),
    "        'IM %ile WT - immunogenicity percentile rank of the corresponding wild type epitope across all immunogenicity prediction algorithms used (those that provide percentile output).', ",
    "        'RNA Expr - Gene expression value for the annotated gene containing the variant.',",
    "        'RNA VAF - Tumor RNA variant allele frequency (VAF) at this position.',",
    "        'Allele Expr - Gene expression value * Tumor RNA VAF. This is used to approximate the expression of the variant allele.',",
    "        'RNA Depth - Tumor RNA depth at this position.',",
    "        'DNA VAF - Tumor DNA variant allele frequency (VAF) at this position.',",
    "        'Tier - A tier suggesting the suitability of variants for use in vaccines.',",
    "        'Ref Match - Indicates if the query sequence has a hit in the reference proteome.',",
    "        'Acpt - Click the thumbs-up button to accept a neoantigen candidate.',",
    "        'Rej - Click the thumbs-down button to reject a neoantigen candidate.',",
    "        'Rev - Click the flag button to mark a neoantigen candidate for review.'],",
    "header = table.columns().header();",
    paste("for (var i = ", hla_count+1, "; i-", hla_count+1, " < tips.length; i++) {"),
    paste("$(header[i]).attr('title', tips[i-", hla_count+1, "]);"),
    "}"
  )
}


#### ui side styling settings
csscode <- HTML("
.sidebar-mini.sidebar-collapse .shiny-bound-input.action-button {
  margin: 6px 6px 6px 3px;
  max-width: 85%;
}
.sidebar-mini.sidebar-collapse .fa {
  font-size: initial;
}
.sidebar-mini.sidebar-collapse #tohide {
  display: none;
}
table.dataTable tbody tr.selected td {
  box-shadow: inset 0 0 0 9999px rgba(55, 55, 55, 0.3) !important;
}
.dataTables_length {
  margin-left: 20px;
  padding-top: 6px;
}
div #neofox_last_selected {
  display: inline;
}
")

# Create the theme
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#4e635c"
  )
)
