# load shiny library
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(fresh)
library(shinycssloaders)

source("styling.R")

## UPLOAD TAB ##
upload_tab <- tabItem(
    "upload",
    # infoBoxes
    fluidRow(
        column(width = 6,
            box(
                title="Option 1: View demo data", status = "primary", solidHeader = TRUE, width = NULL,
                actionButton("loadDefaultmain", "Load demo data", style = "color: #fff; background-color: #c92424; border-color: #691111"),
                h5("Please wait a couple seconds after clicking and you should be redirected to the Visualize and Explore tab.")
            ),
            box(
                title = "Option 2: Upload your own data Files", status = "primary", solidHeader = TRUE, width = NULL,
                HTML("<h5><b>(Required)</b> Please upload the aggregate report file. Note that this will be the data displayed in the main table in the Explore tab.</h5>"),
                uiOutput("aggregate_report_ui"),
                radioButtons("hla_class", "Does this aggregate report file correspond to Class I or Class II prediction data?",
                            c("Class I data (e.g. HLA-A*02:01) " = "class_i", "Class II data (e.g. DPA1*01:03)" = "class_ii")),
                hr(style = "border-color: white"),
                HTML("<h5><b>(Required)</b> Please upload the corresponding metrics file for the main file that you have chosen.</h5>"),
                uiOutput("metrics_ui"),
                hr(style = "border-color: white"),
                HTML("<h5><b>(Optional)</b> If you would like, you can upload an additional aggregate report file generated with either Class I or Class II results to supplement your main table. (E.g. if you uploaded Class I data as the main table, you can upload your Class II report here as supplemental data)</h5>"),
                uiOutput("add_file_ui"),
                textInput("add_file_label", "Please provide a label for the additional file uploaded (e.g. Class I data or Class II data)"),
                hr(style = "border-color: white"),
                HTML("<h5><b>(Optional)</b> Additionally, you can upload a gene-of-interest list in a tsv format, where each row is a single gene name. These genes (if in your aggregate report) will be highlighted in the Gene Name column.</h5>"),
                fileInput(inputId = "gene_list", label = "4. Gene-of-interest List (tsv required)", accept = c("text/tsv", "text/tab-separated-values,text/plain", ".tsv")),
                actionButton("visualize", "Visualize")
            )
        ),
        column(6,
            box(
                title = "Basic Instructions: How to explore your data using pVACview?", status = "primary", solidHeader = TRUE, width = NULL,
                h4("Step 1: Upload your own data / Load demo data", style = "font-weight: bold"),
                h5("You can either choose to explore a demo dataset that we have prepared from the HCC1395 cell line, or choose to upload your own datasets."),
                HTML("<h5>If you are uploading your own datasets, the two required inputs are output files you obtain after running the pVACseq pipeline. 
                    The <b>aggregated tsv file</b> is a list of all predicted epitopes and their binding affinity scores with additional variant information
                    and the <b>metrics json file</b> contains additional transcript and peptide level information.</h5>"),
                h5("You have the option of uploading an additional file to supplement the data you are exploring. This includes: additional class I or II information and
                    a gene-of-interest tsv file."),
                actionButton("help_doc_upload", "More details", onclick = "window.open('https://pvactools.readthedocs.io/en/latest/pvacview/getting_started.html#upload', '_blank')"),
                h4("Step 2: Exploring your data", style = "font-weight: bold"),
                HTML("<h5>To explore the different aspects of your neoantigen candidates, you will need to navigate to the <b>Aggregate Report of Best Candidate by Variant</b> on the visualize and explore tab.
                    For detailed variant, transcript and peptide information for each candidate listed, you will need to click on the <b>Investigate button</b> for the specific row of interest.
                    This will prompt both the transcript and peptide table to reload with the matching information.</h5>"),
                h5("By hovering over each column header, you will be able to see a brief description of the corresponding column and for more details, you can click on the tooltip located at the top right of the aggregate report table.", br(),
                "After investigating each candidate, you can label the candidate using the dropdown menu located at the second to last column of the table. Choices include:
                    Accept, Reject or Review."),
                actionButton("help_doc_explore", "More details", onclick = "window.open('https://pvactools.readthedocs.io/en/latest/pvacview/getting_started.html#visualize-and-explore', '_blank')"),
                h4("Step 3: Exporting your data", style = "font-weight: bold"),
                h5("When you have either finished ranking your neoantigen candidates or need to pause and would like to save your current evaluations, 
                    you can export the current main aggregate report using the export page."),
                HTML("<h5>Navigate to the export tab, and you will be able to name your file prior to downloading in either tsv or excel format.
                    The excel format is user-friendly for downstream visualization and manipulation. However, if you plan on to continuing editing the aggregate report 
                    and would like to load it back in pVACview with the previous evaluations preloaded, you will <b>need</b> to download the file in a <b>tsv format</b>. 
                    This serves as a way to save your progress as your evaluations are cleared upon closing or refreshing the pVACview app.</h5>"),
                actionButton("help_doc_export", "More details", onclick = "window.open('https://pvactools.readthedocs.io/en/latest/pvacview/getting_started.html#export', '_blank')")
            )
        ),
    )
)

## EXPLORE TAB ##
explore_tab <- tabItem(
    "explore",
    conditionalPanel(
        condition = "output.filesUploaded",
        fluidRow(
            tags$style(
                type = "text/css",
                ".modal-dialog { width: fit-content !important; }"
            ),
            tags$script(HTML("Shiny.addCustomMessageHandler('unbind-DT', function(id) {
            Shiny.unbindAll($('#'+id).find('table').DataTable().table().node());
                                })")),
            box(width = 6,
                title = "Advanced Options: Regenerate Tiering with different parameters",
                status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                "*Please note that the metrics file is required in order to regenerate tiering information with different parameters", br(),
                "Current version of pVACseq results defaults to positions 1, 2, n-1 and n (for a n-mer peptide) when determining anchor positions.
                If you would like to use our allele specific anchor results and regenerate the tiering results for your variants,
                please specify your contribution cutoff and submit for recalculation. ", tags$a(href = "https://www.biorxiv.org/content/10.1101/2020.12.08.416271v1", "More details can be found here.", target = "_blank"), br(),
                checkboxInput("use_anchor", "If you want to use allele-specific anchor calculations, please check this box. Otherwise anchors will be calculated as 1,2 and n-1,n for n-mer peptides.", value = FALSE, width = NULL),
                sliderInput("anchor_contribution", "Contribution cutoff for determining anchor locations", 0.5, 0.9, 0.8, step = 0.1, width = 400),
                uiOutput("binding_threshold_ui"),
                checkboxInput("allele_specific_binding", "If you want to use allele-specific binding thresholds for tiering purposes please check this box.", value = FALSE, width = NULL),
                uiOutput("percentile_threshold_ui"),
                uiOutput("dna_cutoff_ui"),
                uiOutput("allele_expr_ui"),
                h5("For further explanations on these inputs, please refer to the ", tags$a(href = "https://pvactools.readthedocs.io/en/latest/pvacview/getting_started.html#visualize-and-explore", "pVACview documentation.", target = "_blank")),
                actionButton("submit", "Recalculate Tiering with new parameters"),
                style = "overflow-x: scroll;font-size:100%"),
            box(width = 3,
                title = "Original Parameters for Tiering",
                status = "primary", solidHeader = TRUE, collapsible = TRUE,
                column(width = 12,
                h5("These are the original parameters used in the tiering calculations extracted from the metrics data file given as input."),
                tableOutput("paramTable"), style = "height:250px; overflow-y: scroll;overflow-x: scroll;"),
                actionButton("reset_params", "Reset to original parameters"),
                style = "overflow-x: scroll;font-size:100%"),
            box(width = 3,
                title = "Add Comments for selected variant",
                status = "primary", solidHeader = TRUE, collapsible = TRUE,
                textAreaInput("comments", "Please add/update your comments for the variant you are currently examining", value = ""),
                actionButton("comment", "Update Comment Section"),
                h5("Comment:"), verbatimTextOutput("comment_text"),
                style = "overflow-x: scroll;font-size:100%")
        ),
        fluidRow(
            box(width = 12,
                title = "Aggregate Report of Best Candidates by Variant",
                status = "primary", solidHeader = TRUE, collapsible = TRUE,
                enable_sidebar = TRUE, sidebar_width = 25, sidebar_start_open = TRUE,
                dropdownMenu = boxDropdown(boxDropdownItem("Help", id = "help", icon = icon("question-circle"))),
                selectInput("page_length", "Number of variants displayed per page:", selected = "10", c("10", "20", "50", "100"), width = "280px"),
                DTOutput("mainTable") %>% withSpinner(color = "#8FCCFA"),
                span("Currently investigating row: ", verbatimTextOutput("selected")),
                style = "overflow-x: scroll;font-size:100%")
        ),

        fluidRow(
            box(width = 12, title = "Variant Information",  status = "primary", solidHeader = TRUE, collapsible = TRUE,
                tabBox(width = 6, title = " ",
                    tabPanel("Transcript Sets of Selected Variant",
                        DTOutput("transcriptSetsTable") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;font-size:100%"),
                    tabPanel("Additional Data",
                        span("Additional Data Type: ", verbatimTextOutput("type_text")),
                        span("Median MT IC50: ", verbatimTextOutput("addData_IC50")),
                        span("Median MT Percentile: ", verbatimTextOutput("addData_percentile")))
                ),
                box(width = 4, solidHeader = TRUE, title = "Variant & Gene Info",
                    span("DNA VAF", verbatimTextOutput("metricsTextDNA")),
                    span("RNA VAF", verbatimTextOutput("metricsTextRNA")),
                    span("Gene Expression", verbatimTextOutput("metricsTextGene")),
                    span("Genomic Information (chromosome - start - stop - ref - alt)", verbatimTextOutput("metricsTextGenomicCoord")),
                    h5("Additional variant information:"),
                    uiOutput("url"), style = "overflow-x: scroll;font-size:100%"),
                box(width = 2, solidHeader = TRUE, title = "Peptide Evalutation Overview",
                    tableOutput("checked"), style = "overflow-x: scroll;font-size:100%")
            )
        ),
        fluidRow(
            box(width = 12, title = "Transcript Set Detailed Data", solidHeader = TRUE, collapsible = TRUE, status = "primary",
                tabBox(width = 12, title = " ",
                    tabPanel("Peptide Candidates from Selected Transcript Set",
                            DTOutput("peptideTable") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;font-size:100%"),
                    tabPanel("Transcripts in Set",
                            DTOutput("transcriptsTable") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;font-size:100%")
                )
            )
        ),
        fluidRow(
            box(width = 12, title = "Additional Peptide Information",  status = "primary", solidHeader = TRUE, collapsible = TRUE,
                tabBox(title = " ", id = "info",
                    tabPanel("IC50 Plot",
                        h4("Violin Plots showing distribution of MHC IC50 predictions for selected peptide pair (MT and WT)."),
                        plotOutput(outputId = "bindingData_IC50") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;"
                    ),
                    tabPanel("%ile Plot",
                        h4("Violin Plots showing distribution of MHC percentile predictions for selected peptide pair (MT and WT)."),
                        plotOutput(outputId = "bindingData_percentile") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;"
                    ),
                    tabPanel("Binding Data",
                        h4("Prediction score table showing exact MHC binding values for IC50 and percentile calculations."),
                        DTOutput(outputId = "bindingDatatable"), style = "overflow-x: scroll;"
                    ),
                    tabPanel("Elution Table",
                        h4("Prediction score table showing exact MHC binding values for elution and percentile calculations."),
                        DTOutput(outputId = "elutionDatatable"),
                        br(),
                        strong("MHCflurryEL Processing"), span(': An "antigen processing" predictor that attempts to model MHC allele-independent effects such as proteosomal cleavage. ('),
                        a(href = "https://www.sciencedirect.com/science/article/pii/S2405471220302398", "Citation"), span(")"),
                        br(),
                        strong("MHCflurryEL Presentation"), span(': A predictor that integrates processing predictions with binding affinity predictions to give a composite "presentation score." ('),
                        a(href = "https://www.sciencedirect.com/science/article/pii/S2405471220302398", "Citation"), span(")"),
                        br(),
                        strong("NetMHCpanEL / NetMHCIIpanEL"), span(": A predictor trained on eluted ligand data. ("),
                        a(href = "https://academic.oup.com/nar/article/48/W1/W449/5837056", "Citation"), span(")"),
                        style = "overflow-x: scroll;"
                    ),
                    tabPanel("Anchor Heatmap",
                        h4("Allele specific anchor prediction heatmap for top 20 candidates in peptide table."),
                        plotOutput(outputId = "peptideFigureLegend", height = "50px"),
                        plotOutput(outputId = "anchorPlot") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;"
                    )
                ),
                box(
                    column(width = 4,
                        h4("Allele Specific Anchor Prediction Heatmap"),
                        h5(" This tab displays HLA allele specific anchor predictions overlaying good-binding peptide sequences generated from each specific transcript.", br(),
                            " Current version supports the first 15 MT/WT peptide sequence pairs (first 30 rows of the peptide table)."), br(),
                        h4("MHC Binding Prediction Scores"),
                        h5(" This tab contains violin plots that showcase individual binding prediction scores from each algorithm used. A solid line is used to represent the median score.")
                    ),
                    column(width = 8,
                        box(title = "Anchor vs Mutation position Scenario Guide", collapsible = TRUE, collapsed = FALSE, width = 12,
                            img(src = "https://github.com/griffithlab/pVACtools/raw/master/pvactools/tools/pvacview/www/anchor.jpg",
                            align = "center", height = "350px", width = "600px"), style = "overflow-x: scroll;")
                    )
                )
            )
        )
    ),
    conditionalPanel(
        condition = "output.filesUploaded == false",
        h4("Error: Missing required files (both aggregate report and metrics files are required to properly visualize and explore candidates).", style = "font-weight: bold"),
    )
)

## EXPORT TAB ##
export_tab <- tabItem(
    "export",
    fluidRow(
        textInput("exportFileName", "Export filename: ", value = "Annotated.Neoantigen_Candidates", width = NULL, placeholder = NULL)
    ),
    fluidRow(
        column(12,
            DTOutput("ExportTable") %>% withSpinner(color = "#8FCCFA"))
    )
)

## TUTORIAL TAB ##
tutorial_tab <- tabItem("tutorial",
  tabsetPanel(type = "tabs",
    tabPanel("Variant Level",
        ## Aggregate Report Column Descriptions"
        h3("Main table full column descriptions"),
        p("If using pVACview with pVACtools output, the user is required to at least provide two files: ",
        code("all_epitopes.aggregated.tsv"), code("all_epitopes.aggregated.metrics.json")), br(),
        p("The ", code("all_epitopes.aggregated.tsv"),
        "file is an aggregated version of the all_epitopes TSV. 
        It presents the best-scoring (lowest binding affinity) epitope for each variant, and outputs 
        additional binding affinity, expression, and coverage information for that epitope. 
        It also gives information about the total number of well-scoring epitopes for each variant, 
        the number of transcripts covered by those epitopes, as well as the HLA alleles that those 
        epitopes are well-binding to. Lastly, the report will bin variants into tiers that offer 
        suggestions as to the suitability of variants for use in vaccines."), br(),
        p(strong("Column Names : Description")),
        p(code("ID"), " : ", "A unique identifier for the variant"),
        p(code("HLA Alleles"), " : ", "For each HLA allele in the run, the number of this variant’s 
        epitopes that bound well to the HLA allele (with ", code("lowest"), " or ", code("median"),
        " mutant binding affinity < ", code("aggregate_inclusion_binding_threshold"), ")"),
        p(code("Gene"), " : ", "The Ensembl gene name of the affected gene"),
        p(code("AA Change"), " : ", "The amino acid change for the mutation"),
        p(code("Num Passing Transcripts"), " : ", "The number of transcripts 
        for this mutation that resulted in at least one well-binding peptide (", code("lowest"), " or ",
        code("median"), " mutant binding affinity < ", code("aggregate_inclusion_binding_threshold"), ")"),
        p(code("Best Peptide"), " : ", "The best-binding mutant epitope sequence (lowest mutant binding affinity) 
        prioritizing epitope sequences that resulted from a protein_coding transcript with a TSL below the 
        maximum transcript support level and having no problematic positions."),
        p(code("Best Transcript"), " : ", "Transcript corresponding to the best peptide with the lowest TSL and shortest length."),
        p(code("TSL"), " : ", "Transcript support level of the best peptide"),
        p(code("Pos"), " : ", "The one-based position of the start of the mutation within the epitope sequence. ",
        code("0"), " if the start of the mutation is before the epitope (as can occur downstream of frameshift mutations)"),
        p(code("Num Passing Peptides"), " : ", "The number of unique well-binding peptides for this mutation."),
        p(code("IC50 MT"), " : ", code("Lowest"), " or ", code("Median"), " ic50 binding affinity of 
        the best-binding mutant epitope across all prediction algorithms used."),
        p(code("IC50 WT"), " : ", code("Lowest"), " or ", code("Median"), " ic50 binding affinity of 
        the corresponding wildtype epitope across all prediction algorithms used."),
        p(code("%ile MT"), " : ", code("Lowest"), " or ", code("Median"), "binding affinity percentile rank 
        of the best-binding mutant epitope across all prediction algorithms used (those that provide percentile output)"),
        p(code("%ile WT"), " : ", code("Lowest"), " or ", code("Median"), "binding affinity percentile rank of the 
        corresponding wildtype epitope across all prediction algorithms used (those that provide percentile output)"),
        p(code("RNA Expr"), " : ", "Gene expression value for the annotated gene containing the variant."),
        p(code("RNA VAF"), " : ", "Tumor RNA variant allele frequency (VAF) at this position."),
        p(code("Allele Expr"), " : ", "RNA Expr * RNA VAF"),
        p(code("RNA Depth"), " : ", "Tumor RNA depth at this position."),
        p(code("DNA VAF"), " : ", "Tumor DNA variant allele frequency (VAF) at this position."),
        p(code("Tier"), " : ", "A tier suggesting the suitability of variants for use in vaccines."),
        p(code("Evaluation"), " : ", "Column to store the evaluation of each variant when evaluating the run in pVACview.
        Either ", code("Accept"), " ", code("Reject"), "  or ", code("Review"), "."),
        ## Tiering Explained ##
        h3("How is the Tiering column determined / How are the Tiers assigned?"), br(),
        p("Note that if a percentile threshold has been provided, then the ", code("%ile MT"), " column is also required to be lower than
        the given threshold to qualify for tiers, including Pass, Anchor, Subclonal and LowExpr."), br(),
        p(strong("Tier : Criteria")),
        p(code("Pass"), " : ", code(("(MT binding < binding threshold) AND allele expr filter pass AND vaf clonal filter pass 
        AND tsl filter pass AND anchor residue filter pass"))),
        p(code("Anchor"), " : ", code(("(MT binding < binding threshold) AND allele expr filter pass AND vaf clonal filter pass 
        AND tsl filter pass AND anchor residue filter fail"))),
        p(code("Subclonal"), " : ", code(("(MT binding < binding threshold) AND allele expr filter pass AND vaf clonal filter fail 
        AND tsl filter pass AND anchor residue filter pass"))),
        p(code("LowExpr"), " : ", code(("(MT binding < binding threshold) AND low expression criteria met AND allele expr filter pass 
        AND vaf clonal filter pass AND tsl filter pass AND anchor residue filter pass"))),
        p(code("Poor"), " : ", "Best peptide for current variant FAILS in two or more categories"),
        p(code("NoExpr"), " : ", code("((gene expr == 0) OR (RNA VAF == 0)) AND low expression criteria not met")), br(),
        p("Here we list out the exact criteria for passing each respective filter: "),
        p(strong("Allele Expr Filter: "), code("(allele expr >= allele expr cutoff) OR (rna_vaf == 'NA') OR (gene_expr == 'NA')")),
        p(strong("VAF Clonal Filter: "), code("(dna vaf < vaf subclonal) OR (dna_vaf == 'NA')")),
        p(strong("TSL Filter: "), code("(TSL != 'NA') AND (TSL < maximum transcript support level)")),
        p(strong("Anchor Residue Filter: "), br(),
        strong("1. "), code("(Mutation(s) is at anchor(s)) AND
        ((WT binding < binding threshold) OR (WT percentile < percentile threshold))"), br(),
        strong("    OR"), br(), strong("2. "), code("Mutation(s) not or not entirely at anchor(s)")),
        p(strong("Low Expression Criteria: "), code("(allele expr > 0) OR ((gene expr == 0) AND (RNA Depth > RNA Coverage Cutoff) AND (RNA VAF > RNA vaf cutoff))"))
    ),
    tabPanel("Transcript Level",
        h3(" "),
        fluidRow(
            column(width = 6,
                h4("Transcript Set Table", style = "font-weight: bold; text-decoration: underline;"),
                p("Upon selecting a variant for investigation, you may have multiple transcripts covering the region.", br(), br(),
                "These transcripts are first grouped into ", strong("Trancripts Sets"), " , which is based on the good binding peptides
                produced. Transcripts that produce the exact same set of peptides are put into the same group.", br(), br(),
                "The table also lists the number of transcripts and corresponding peptides in each set (each pair of WT and MT peptides are considered 1 when
                counting).", br(), " A sum of the total expression across all transcripts in each set is also shown.", br(), " A light green color is used to
                highlight the ", strong("Transcript Set"), " producing the Best Peptide for the variant in question.")
            ),
            column(width = 6,
                img(src = "Explore - Transcript Set.png",
                align = "center", height = "300px", width = "500px"),
            )
        ),
        fluidRow(
            column(width = 3,
                h4("Transcript Set Detailed Data", style = "font-weight: bold; text-decoration: underline;"),
                p("Upon selecting a specific transcript set, you can now look in more detail which exact transcripts are included.", br(), br(),
                "The ", strong("Transcripts in Set"), "table lists all information regarding each transcript including:", br(), br(),
                "Transcript ID, Gene Name, Amino Acid Change, Mutation Position, individual transcript expression, transcript support level, biotype and transcript length.", br(), br(),
                " A light green color is used to highlight the ", strong("Transcript in Selected Set"), " that produced the Best Peptide for the variant in question.")
            ),
            column(width = 9,
                img(src = "Explore - Transcript in Set.png",
                align = "center", height = "300px", width = "1200px"),
            )
        )
    ),
    tabPanel("Peptide Level",
        h4(" "),
        fluidRow(
            column(width = 12,
                h4("Peptide Table", style = "font-weight: bold; text-decoration: underline;"),
                p("Upon selecting a specific transcript set, you can also now look in detail which good-binding peptides are produced from this set.", br(), br(),
                "Both mutant (", code("MT"), ") and wildtype (", code("WT"), ") sequences are shown, along with ", code("lowest"), " or ", code("median"),
                " binding affinities, depending on how you generated the aggregate report.", br(), br(),
                "An ", code("X"), "is marked for binding affinities higher than the ", code("aggregate_inclusion_binding_threshold"), " set when generating the aggregate report.", br(), br(),
                "We also include two extra columns, one specifying the mutation and position and another providing information on any problematic amino acids.", br(),
                "Note that if users wish to utlitize the ", strong("problematic positions"), " feature, they should run the standalone command ", code("pvacseq identify_problematic_amino_acids"),
                " or run pVACseq with the ", code("--problematic-amino-acids"), " option enabled to generate the needed information."
                )
            )
        ),
        fluidRow(
            column(width = 12,
                img(src = "Explore - Peptide Table.png",
                align = "center", height = "400px", width = "1500px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 4,
                h4("Additional Information",  style = "font-weight: bold; text-decoration: underline;"),
                h5("IC50 Plot", style = "font-weight: bold;"),
                p("By clicking on each MT/WT peptide pair, you can then assess the peptides in more detail by navigating to the ", strong("Additional Peptide Information"), " tab.", br(), br(),
                "There’s five different tabs in this section of the app, providing peptide-level details on the MT/WT peptide pair that you have selected", br(),
                "The ", strong("IC50 Plot"), "tab shows violin plots of the individual IC50-based binding affinity predictions of the MT and WT peptides for HLA
                alleles were the MT binds well to. These peptides each have up to 8 binding algorithm scores (for Class I alleles with pVACseq version 3.0) or up
                to 4 algorithm scores (for Class II alleles with pvacseq version 3.0).", br())
            ),
            column(width = 8,
                img(src = "Explore - IC50 Plots.png",
                align = "center", height = "350px", width = "700px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 4,
                h5("%ile Plot", style = "font-weight: bold;"),
                p("The ", strong("%ile Plot"), "tab shows violin plots of the individual percentile-based binding affinity predictions of the MT and WT peptides
                for HLA alleles were the MT binds well to.", br())
            ),
            column(width = 8,
                img(src = "Explore - Percentile Plots.png",
                align = "center", height = "350px", width = "700px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 4,
                h5("Binding Data", style = "font-weight: bold;"),
                p("The ", strong("Binding Data"), "tab shows the specific IC50 and percentile binding affinity predictions generated from each individual algorithm.
                Each cell shows the IC50 prediction followed by the percentile predictions in parenthesis.", br())
            ),
            column(width = 8,
                img(src = "Explore - Binding Data.png",
                align = "center", height = "350px", width = "720px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 4,
                h5("Elution Table", style = "font-weight: bold;"),
                p("The ", strong("Elution Table"), "tab shows prediction results based on algorithms trained from peptide elution data. This includes algorithms
                such as NetMHCpanEL/NetMHCIIpanEL, MHCflurryELProcessing and MHCflurryELPresentation.", br())
            ),
            column(width = 8,
                img(src = "Explore - Elution Data.png",
                align = "center", height = "350px", width = "720px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 4,
                h5("Anchor Heatmap", style = "font-weight: bold;"),
                p("The ", strong("Anchor Heatmap"), "tab shows the top 30 MT/WT peptide pairs from the peptide table with anchor probabilities overlaying as a heatmap.
                The anchor probabilities shown are both allele and peptide length specific. The mutated amino acid is marked in red (for missense mutations) and each 
                MT/WT pair are separated from others using a dotted line. ", br(),
                "For peptide sequences with no overlaying heatmap, we currently do not have allele-specific predictions for them in our database.", br(), br(),
                "For more details and explanations regarding anchor positions and its influence on neoantigen prediction and prioritization, please refer to the next section: ",
                strong("Advanced Options: Anchor Contribution"))
            ),
            column(width = 8,
                img(src = "Explore - Anchor Heatmap.png",
                align = "center", height = "350px", width = "720px"), br(), br()
            )
        )
    )
  )
)

## CONTACT TAB ##
contact_tab <- tabItem("contact",
 p("Bug reports or feature requests can be submitted on the ", tags$a(href = "https://github.com/griffithlab/pVACtools", "pVACtools Github page."),
 "You may also contact us by email at ", code("help@pvactools.org", "."))

)

ui <- dashboardPage(
  ## HEADER ##
  header = dashboardHeader(
    title = tagList(tags$a(class = "logo",
                         span(class = "logo-mini", tags$img(src = "https://github.com/griffithlab/pVACtools/raw/master/pvactools/tools/pvacview/www/pVACview_logo_mini.png")),
                         span(class = "logo-lg", tags$img(src = "https://github.com/griffithlab/pVACtools/raw/master/pvactools/tools/pvacview/www/pVACview_logo.png"))
                  )),
    tags$li(class = "dropdown", tags$a(href = "https://pvactools.readthedocs.io/en/latest/", class = "my_class", "Help", target = "_blank"))
    ),
  ## SIDEBAR ##
  sidebar = dashboardSidebar(
    sidebarMenu(
      tags$head(tags$style(csscode)),
      id = "tabs",
      menuItem("pVACtools Output", tabName = "pvactools", startExpanded = TRUE, icon = icon("far fa-chart-bar"),
               br(),
               menuSubItem("Upload", tabName = "upload", icon = icon("upload")),
               br(),
               menuSubItem("Visualize and Explore", tabName = "explore", icon = icon("digital-tachograph")),
               br(),
               menuSubItem("Export", tabName = "export", icon = icon("file-export")),
               br()
      ),
      menuItem("Tutorials", tabName = "tutorial", startExpanded = TRUE, icon = icon("fas fa-book-open")),
      menuItem("pVACview Documentation", icon = icon("fas fa-file-invoice"), href = "https://pvactools.readthedocs.io/en/latest/pvacview.html"),
      menuItem("Submit Github Issue", tabName = "contact", icon = icon("far fa-question-circle"))
    )
  ),
  body = dashboardBody(
    use_theme(mytheme),
    tags$head(
      tags$style(HTML(css)),
      tags$style(HTML("table.dataTable tr.selected td, table.dataTable td.hover {background-color: #EAF2F8 !important;}")),
      tags$style(HTML("table.dataTable { border-collapse: collapse;}")),
      tags$style(HTML("table.dataTable.hover tbody tr:hover, table.dataTable.display tbody tr:hover {
                              background-color: #92c8f0 !important; } ")),
      tags$style(HTML(".skin-blue .main-header .logo {background-color: #dff5ee;}")),
      tags$style(HTML(".skin-blue .main-header .navbar { background-color: #739187;}")),
      tags$style(HTML("element.style {}.skin-blue .wrapper, .skin-blue .main-sidebar, .skin-blue .left-side {background-color: #739187;}")),
      tags$style(HTML(".main-header .sidebar-toggle {background-color: #b6d1c8}")),
      tags$style(HTML(".box-header.with-border {border-bottom: 1px solid #f4f4f4;}")),
      tags$style(HTML(".skin-blue .main-header .navbar .sidebar-toggle {color: #4e635c;}")),
      tags$style(HTML(".content-wrapper {background-color: #ecf0f5;}")),
      tags$style(HTML(".main-header .logo {padding-right : 5px; padding-left : 5px;}")),
      tags$style(HTML(".box.box-solid.box-primary {border-radius: 12px}")),
      tags$style(HTML(".box-header.with-border {border-radius: 10px}"))
    ),

    tabItems(
      ## UPLOAD TAB ##
      upload_tab,
      ## EXPLORE TAB ##
      explore_tab,
      ## EXPORT TAB ##
      export_tab,
      ## TUTORIAL TAB ##
      tutorial_tab,
      ## CONTACT TAB ##
      contact_tab
    )
  )
)