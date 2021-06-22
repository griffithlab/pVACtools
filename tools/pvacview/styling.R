
##### SERVER 
rowCallback <- function(hla_count, row_num) {
  gsub("15", hla_count+8,
  gsub("13", hla_count+6,
  c(
  "function(row, data, displayNum, displayIndex){", 
  "  var x = data[13];", 
  "  var y = data[15];",
  gsub("0", row_num,"  if (displayIndex == 0){"),
  "  $('td',row).css('border-top','3px solid #0390fc');",
  "  $('td',row).css('border-bottom','3px solid #0390fc');",
  "  }",
  "  if (x > 0 & x<= 50){",
  "    $('td:eq(13)', row).addClass('zero');",
  "  } else if(x > 50 & x <= 100){",
  "    $('td:eq(13)', row).addClass('fifty');",
  "  } else if(x > 100 & x <= 200){",
  "    $('td:eq(13)', row).addClass('one');",
  "  } else if(x > 200 & x <= 300){",
  "    $('td:eq(13)', row).addClass('two');",
  "  } else if(x > 300 & x <= 400){",
  "    $('td:eq(13)', row).addClass('three');",
  "  } else if(x > 400 & x <= 500){",
  "    $('td:eq(13)', row).addClass('four');",
  "  } else if(x > 500 & x <= 600){",
  "    $('td:eq(13)', row).addClass('five');",
  "  } else if(x > 600 & x <= 700){",
  "    $('td:eq(13)', row).addClass('six');",
  "  } else if(x > 700 & x <= 800){",
  "    $('td:eq(13)', row).addClass('seven');",
  "  } else if(x > 800 & x <= 900){",
  "    $('td:eq(13)', row).addClass('eight');",
  "  } else if(x > 900 & x <= 1000){",
  "    $('td:eq(13)', row).addClass('nine');",
  "  } else if(x > 1000){",
  "    $('td:eq(13)', row).addClass('ten');",
  "  }",
  "  if (y > 0 & y <= 0.1){",
  "    $('td:eq(15)', row).addClass('fifty');",
  "  } else if(y > 0.1 & y <= 0.2){",
  "    $('td:eq(15)', row).addClass('one');",
  "  } else if(y > 0.2 & y <= 0.3){",
  "    $('td:eq(15)', row).addClass('two');",
  "  } else if(y > 0.3 & y <= 0.4){",
  "    $('td:eq(15)', row).addClass('three');",
  "  } else if(y > 0.4 & y <= 0.5){",
  "    $('td:eq(15)', row).addClass('four');",
  "  } else if(y > 0.5 & y <= 0.75){",
  "    $('td:eq(15)', row).addClass('five');",
  "  } else if(y > 0.75 & y <= 1){",
  "    $('td:eq(15)', row).addClass('six');",
  "  } else if(y > 1 & y <= 1.5){",
  "    $('td:eq(15)', row).addClass('seven');",
  "  } else if(y > 1.5 & y <= 2){",
  "    $('td:eq(15)', row).addClass('eight');",
  "  } else if(y > 2){",
  "    $('td:eq(15)', row).addClass('ten');",
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
  #adminlte_sidebar(
  #  width = "300px",
  #  dark_bg = "#D8DEE9",
  #  #dark_hover_bg = "#81A1C1",
  #  dark_color = "#2E3440"
  #),
  #adminlte_global(
  #  content_bg = "#d8ede2"
  #)
)
