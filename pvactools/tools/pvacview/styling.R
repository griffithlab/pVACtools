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
  "        'Best Peptide - The best-binding mutant epitope sequence (lowest mutant binding affinity) prioritizing epitope sequences that resulted from a protein_coding transcript with a TSL below the maximum transcript support level and having no problematic positions.',",
  "        'Best Transcript - Transcript corresponding to the best peptide with the lowest TSL and shortest length.',",
  "        'TSL - Transcript support level of the best peptide.',",
  "        'Allele - HLA allele the best peptide binds well to.',",
  "        'Pos - The one-based position of the start of the mutation within the epitope sequence. 0 if the start of the mutation is before the epitope (as can occur downstream of frameshift mutations).',",
  "        'Prob Pos - Problematic positions within the best peptide.',",
  "        'Num Passing Peptides - The number of unique well-binding peptides for this mutation.',",
  gsub("X", score_mode,"      'IC50 MT - X IC50 binding affinity of the best-binding mutant epitope across all prediction algorithms used.', "),
  "        'IC50 WT - IC50 binding affinity of the corresponding wildtype epitope.',",
  gsub("X", score_mode,"      '%ile MT - X binding affinity percentile rank of the best-binding mutant epitope across all prediction algorithms used (those that provide percentile output).', "),
  "        '%ile WT - Binding affinity percentile rank of the corresponding wildtype epitope across all prediction algorithms used (those that provide percentile output).', ",
  "        'RNA Expr - Gene expression value for the annotated gene containing the variant.',",
  "        'RNA VAF - Tumor RNA variant allele frequency (VAF) at this position.',",
  "        'Allele Expr - Gene expression value * Tumor RNA VAF. This is used to approximate the expression of the variant allele.',",
  "        'RNA Depth - Tumor RNA depth at this position.',",
  "        'DNA VAF - Tumor DNA variant allele frequency (VAF) at this position.',",
  "        'Tier - A tier suggesting the suitability of variants for use in vaccines.',",
  "        'Eval - User-selected evaluation of neoantigen candidate. Options include: Accept, Reject, Review. (Default: Pending)'],",
  "header = table.columns().header();",
  gsub("7", hla_count, "for (var i = 7; i-7 < tips.length; i++) {"),
  gsub("7", hla_count, "$(header[i]).attr('title', tips[i-7]);"),
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
")

# Create the theme
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#4e635c"
  )
)
