
##### SERVER 
rowCallback <- function(hla_count, row_num) {
  gsub("15", hla_count+8,
  gsub("13", hla_count+6,
  c(
  "function(row, data, displayNum, displayIndex){", 
  gsub("0", row_num,"  if (displayIndex == 0){"),
  "  $('td',row).css('border-top','3px solid #0390fc');",
  "  $('td',row).css('border-bottom','3px solid #0390fc');",
  "  }",
  "}"
  ))
)}

callBack <- function(hla_count) {
  c(
  "var tips = ['Gene - The Ensembl gene name of the affected gene.',",
  "        'AA Change - The amino acid change for the mutation. Note that FS indicates a frameshift variant.',",
  "        'Num Passing Transcripts - The number of transcripts for this mutation that resulted in at least one well-binding peptide (median mutant binding affinity < 1000).',",
  "        'Best Peptide - The best-binding mutant epitope sequence (lowest median mutant binding affinity).',",
  "        'Pos - The one-based position of the start of the mutation within the epitope sequence. 0 if the start of the mutation is before the epitope (as can occur downstream of frameshift mutations).',",
  "        'Num Passing Peptides - The number of unique well-binding peptides for this mutation.',",
  "        'IC50 MT - Median ic50 binding affinity of the best-binding mutant epitope across all prediction algorithms used.', ",
  "        'IC50 WT - Median ic50 binding affinity of the corresponding wildtype epitope across all prediction algorithms used.',",
  "        '%ile MT - Median binding affinity percentile rank of the best-binding mutant epitope across all prediction algorithms used (those that provide percentile output).', ",
  "        '%ile WT - Median binding affinity percentile rank of the corresponding wildtype epitope across all prediction algorithms used (those that provide percentile output).', ",
  "        'RNA Expr - Gene expression value for the annotated gene containing the variant.',",
  "        'RNA VAF - Tumor RNA variant allele frequency (VAF) at this position.',",
  "        'Allele Expr - Gene expression value * Tumor RNA VAF. This is used to approximate the expression of the variant allele.',",
  "        'RNA Depth - Tumor RNA depth at this position.',",
  "        'DNA VAF - Tumor DNA variant allele frequency (VAF) at this position.',",
  "        'Tier - A tier suggesting the suitability of variants for use in vaccines.',",
  "        'Eval - User-selected evaluation of neoantigen candidate. Options include: Accept, Reject, Review. (Default: Pending)'],",
  "header = table.columns().header();",
  gsub("7", hla_count,"for (var i = 7; i-7 < tips.length; i++) {"),
  gsub("7", hla_count,"$(header[i]).attr('title', tips[i-7]);"),
  "}"
)
}


#### UI 

css <- "
table.dataTable tbody tr td.zero {background-color: #00FF00 !important}
table.dataTable tbody tr td.fifty {background-color: #00EE00 !important}
table.dataTable tbody tr td.one {background-color: #00D500 !important}
table.dataTable tbody tr td.two {background-color: #00BC00 !important}
table.dataTable tbody tr td.three {background-color: #00A300 !important}
table.dataTable tbody tr td.four {background-color: #008B00 !important}
table.dataTable tbody tr td.five {background-color: #FFFF00 !important}
table.dataTable tbody tr td.six {background-color: #FFEB00 !important}
table.dataTable tbody tr td.seven {background-color: #FFD800 !important}
table.dataTable tbody tr td.eight {background-color: #FFC500 !important}
table.dataTable tbody tr td.nine {background-color: #FFB100 !important}
table.dataTable tbody tr td.ten {background-color: #FF9999 !important}
"

# Create the theme
mytheme <- create_theme(
  adminlte_color(
    light_blue = "#4e635c"
  )
)
