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
            column(width = 6,
                box(width = 12,
                    title = "Advanced Options: Regenerate Tiering with different parameters",
                    status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                    "*Please note that the metrics file is required in order to regenerate tiering information with different parameters", br(),
                    "Current version of pVACseq results defaults to positions 1, 2, n-1 and n (for a n-mer peptide) when determining anchor positions.
                    If you would like to use our allele specific anchor results and regenerate the tiering results for your variants,
                    please specify your contribution cutoff and submit for recalculation. ", tags$a(href = "https://www.biorxiv.org/content/10.1101/2020.12.08.416271v1", "More details can be found here.", target = "_blank"), br(),
                    uiOutput("allele_specific_anchors_ui"),
                    uiOutput("anchor_contribution_ui"),
                    uiOutput("binding_threshold_ui"),
                    uiOutput("allele_specific_binding_ui"),
                    uiOutput("percentile_threshold_ui"),
                    uiOutput("dna_cutoff_ui"),
                    uiOutput("allele_expr_ui"),
                    h5("For further explanations on these inputs, please refer to the ", tags$a(href = "https://pvactools.readthedocs.io/en/latest/pvacview/getting_started.html#visualize-and-explore", "pVACview documentation.", target = "_blank")),
                    actionButton("submit", "Recalculate Tiering with new parameters"),
                    style = "overflow-x: scroll;font-size:100%"),
                style = "padding:0px;"
            ),
            column(width = 3,
                fluidRow(
                    box(width = 12,
                        title = "Original Parameters for Tiering",
                        status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
                        column(width = 12,
                            h5("These are the original parameters used in the tiering calculations extracted from the metrics data file given as input."),
                            tableOutput("paramTable"),
                            tableOutput("bindingParamTable"),
                            style = "height:250px; overflow-y: scroll;overflow-x: scroll;"
                        ),
                        style = "font-size:100%"
                    )
                ),
                fluidRow(
                    box(width = 12,
                        title = "Current Parameters for Tiering",
                        status = "primary", solidHeader = TRUE, collapsible = TRUE, collapsed = FALSE,
                        column(width = 12,
                            h5("These are current parameters used in the tiering calculaions which may be different from the original parameters if candidates were re-tiered."),
                            tableOutput("currentParamTable"),
                            tableOutput("currentBindingParamTable"),
                            style = "height:250px; overflow-y: scroll;overflow-x: scroll;"
                        ),
                        style = "font-size:100%"
                    )
                ),
                fluidRow(
                    column(
                        width = 12,
                        actionButton("reset_params", "Reset to original parameters", style = "width: 100%"),
                        align = "center",
                        style = "padding-bottom: 20px"
                    )
                ),
                style = "padding:0px;"
            ),
            column(width = 3,
                box(width = 12,
                    title = "Add Comments for selected variant",
                    status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    textAreaInput("comments", "Please add/update your comments for the variant you are currently examining", value = ""),
                    actionButton("comment", "Update Comment Section"),
                    h5("Comment:"), htmlOutput("comment_text"),
                    style = "font-size:100%"),
                style = "padding:0px;"
            )
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
                    tabPanel("Reference Matches",
                        h4("Best Peptide Data"),
                        column(6,
                            span("Best Peptide: "),
                            plotOutput(outputId = "referenceMatchPlot", height="20px")
                        ),
                        column(2,
                            span("AA Change: ", verbatimTextOutput("selectedAAChange"))
                        ),
                        column(2,
                            span("Pos: ", verbatimTextOutput("selectedPos"))
                        ),
                        column(2,
                            span("Gene: ", verbatimTextOutput("selectedGene"))
                        ),
                        h4("Query Data"),
                        h5(uiOutput("hasReferenceMatchData")),
                        column(10,
                            span("Query Sequence: "),
                            plotOutput(outputId = "referenceMatchQueryPlot", height="20px")
                        ),
                        column(2,
                            span("Hits: ", verbatimTextOutput("referenceMatchHitCount"))
                        ),
                        h4("Hits"),
                        DTOutput(outputId = "referenceMatchDatatable") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;"
                    ),
                    tabPanel("Additional Data",
                        span("Additional Data Type: ", verbatimTextOutput("type_text")),
                        span("Median MT IC50: ", verbatimTextOutput("addData_IC50")),
                        span("Median MT Percentile: ", verbatimTextOutput("addData_percentile")),
                        span("Best Peptide: ", verbatimTextOutput("addData_peptide")),
                        span("Corresponding HLA allele: ", verbatimTextOutput("addData_allele")),
                        span("Best Transcript: ", verbatimTextOutput("addData_transcript")))
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
            box(width = 12, title = "Transcript and Peptide Set Data", solidHeader = TRUE, collapsible = TRUE, status = "primary",
                tabBox(width = 12, title = " ",
                    tabPanel("Peptide Candidates from Selected Transcript Set",
                            DTOutput("peptideTable") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;font-size:100%"),
                    tabPanel("Anchor Heatmap",
                        fluidRow(
                            column(width = 6,
                                div(style = 'overflow-y:scroll;height: 500px;',
                                    h4("Allele specific anchor prediction heatmap for top 20 candidates in peptide table."),
                                    h5("HLA allele specific anchor predictions overlaying good-binding peptide sequences generated from each specific transcript.", br(),
                                        " Current version supports the first 15 MT/WT peptide sequence pairs (first 30 rows of the peptide table)."), br(),
                                    plotOutput(outputId = "peptideFigureLegend", height = "50px"),
                                    plotOutput(outputId = "anchorPlot")
                                )
                            ) %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll; overflow-y: scroll",
                            column(width = 6,
                                h4("Anchor vs Mutation position Scenario Guide",
                                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/anchor.jpg",
                                    align = "center", width = "100%")
                                )
                            )
                        ),
                        fluidRow(
                            box(width = 12, title = "Anchor Weights", solidHeader = TRUE, collapsible = TRUE, status = "primary", collapsed = TRUE,
                                DTOutput("anchorWeights") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll"
                            )
                        )
                    ),
                    tabPanel(
                      "Transcripts in Set",
                      DTOutput("transcriptsTable") %>% withSpinner(color = "#8FCCFA"),
                      style = "overflow-x: scroll;font-size:100%"
                    )
                )
            )
        ),
        fluidRow(
            box(width = 12, title = "Additional Peptide Information",  status = "primary", solidHeader = TRUE, collapsible = TRUE,
                tabBox(width = 12, title = " ", id = "info",
                    tabPanel("IC50 Plot",
                        h4("Violin Plots showing distribution of MHC IC50 predictions for selected peptide pair (MT and WT)."),
                        h5("Showcases individual binding prediction scores from each algorithm used. A solid line is used to represent the median score."),
                        plotOutput(outputId = "bindingData_IC50") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;"
                    ),
                    tabPanel("%ile Plot",
                        h4("Violin Plots showing distribution of MHC percentile predictions for selected peptide pair (MT and WT)."),
                        h5("Showcases individual percentile scores from each algorithm used. A solid line is used to represent the median percentile score."),
                        plotOutput(outputId = "bindingData_percentile") %>% withSpinner(color = "#8FCCFA"), style = "overflow-x: scroll;"
                    ),
                    tabPanel("Binding Data",
                        h4("Prediction score table showing exact MHC binding values for IC50 and percentile calculations."),
                        DTOutput(outputId = "bindingDatatable"), style = "overflow-x: scroll;"
                    ),
                    tabPanel("Elution and Immunogenicity Data",
                        h4("Prediction score table showing exact MHC scpres for elution, immunogenicity, and percentile calculations."),
                        DTOutput(outputId = "elutionDatatable"),
                        br(),
                        strong("BigMHC_EL / BigMHC_IM"), span(": A deep learning tool for predicting MHC-I (neo)epitope presentation and immunogenicity. ("),
                        a(href = "https://www.nature.com/articles/s42256-023-00694-6", "Citation", target = "_blank"), span(")"),
                        br(),
                        strong("DeepImmuno"), span(": Deep-learning empowered prediction of immunogenic epitopes for T cell immunity. ("),
                        a(href = "https://academic.oup.com/bib/article/22/6/bbab160/6261914?login=false", "Citation", target = "_blank"), span(")"),
                        br(),
                        strong("MHCflurryEL Processing"), span(': An "antigen processing" predictor that attempts to model MHC allele-independent effects such as proteosomal cleavage. ('),
                        a(href = "https://www.sciencedirect.com/science/article/pii/S2405471220302398", "Citation", target = "_blank"), span(")"),
                        br(),
                        strong("MHCflurryEL Presentation"), span(': A predictor that integrates processing predictions with binding affinity predictions to give a composite "presentation score." ('),
                        a(href = "https://www.sciencedirect.com/science/article/pii/S2405471220302398", "Citation", target = "_blank"), span(")"),
                        br(),
                        strong("NetMHCpanEL / NetMHCIIpanEL"), span(": A predictor trained on eluted ligand data. ("),
                        a(href = "https://academic.oup.com/nar/article/48/W1/W449/5837056", "Citation", target = "_blank"), span(")"),
                        style = "overflow-x: scroll;"
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
        p("If using pVACview with pVACtools output, the user is required to provide at least the following two files: ",
        code("all_epitopes.aggregated.tsv"), code("all_epitopes.aggregated.metrics.json")), br(),
        p("The ", code("all_epitopes.aggregated.tsv"),
        "file is an aggregated version of the all_epitopes TSV. 
        It presents the best-scoring (lowest binding affinity) epitope for each variant, along with 
        additional binding affinity, expression, and coverage information for that epitope. 
        It also gives information about the total number of well-scoring epitopes for each variant, 
        the number of transcripts covered by those epitopes, and the HLA alleles that those 
        epitopes are well-binding to. Here, a well-binding or well-scoring epitope is any epitope that has a stronger 
        binding affinity than the ", code("aggregate_inclusion_binding_threshold"), "described below. The report then bins variants into 
        tiers that offer suggestions about the suitability of variants for use in vaccines."), br(),
        p("The ", code("all_epitopes.aggregated.metrics.json"),
        "complements the ", code("all_epitopes_aggregated.tsv"), "and is required for the tool's proper functioning."), br(),
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
        p(code("Best Peptide"), " : ", "The best-binding mutant epitope sequence (lowest binding affinity) 
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
        Can be ", code("Accept,"), " ", code("Reject"), "  or ", code("Review"), "."),
        ## Tiering Explained ##
        h3("How is the Tiering column determined / How are the Tiers assigned?"), br(),
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
        p(strong("Low Expression Criteria: "), code("(allele expr > 0) OR ((gene expr == 0) AND (RNA Depth > RNA Coverage Cutoff) AND (RNA VAF > RNA vaf cutoff))")),br(),
        p("Note that if a percentile threshold has been provided, then the ", code("%ile MT"), " column is also required to be lower than
        the given threshold to qualify for tiers, including Pass, Anchor, Subclonal and LowExpr.") 
    ),
    tabPanel("Transcript Level",
        h3(" "),
        fluidRow(
            column(width = 6,
                h4("Transcript Set Table", style = "font-weight: bold; text-decoration: underline;"),
                p("Upon selecting a variant for investigation, you may have multiple transcripts covering the region.", br(), br(),
                "These transcripts are grouped into ", strong("Trancripts Sets"), " , based on the good-binding peptides
                produced. (Transcripts that produce the exact same set of peptides are grouped together.)", br(), br(),
                "The table also lists the number of transcripts and corresponding peptides in each set (each pair of WT and MT peptides are considered 1 when
                counting).", br(), " A sum of the total expression across all transcripts in each set is also shown.", br(), " A light green color is used to
                highlight the ", strong("Transcript Set"), " producing the Best Peptide for the variant in question.")
            ),
            column(width = 6,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_Transcript_Set.png",
                align = "center", height = "300px", width = "500px"),
            )
        ),
        fluidRow(
            column(width = 3,
                h4("Transcript Set Detailed Data", style = "font-weight: bold; text-decoration: underline;"),
                p("Upon selecting a specific transcript set, you can see more details about the exact transcripts that are included.", br(), br(),
                "The ", strong("Transcripts in Set"), "table lists all information regarding each transcript including:", br(), br(),
                "Transcript ID, Gene Name, Amino Acid Change, Mutation Position, individual transcript expression, transcript support level, biotype and transcript length.", br(), br(),
                " A light green color is used to highlight the specific", strong("Transcript in Selected Set"), " that produced the Best Peptide for the variant in question.")
            ),
            column(width = 9,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_Transcript_in_Set.png",
                align = "center", height = "300px", width = "1200px"),
            )
        )
    ),
    tabPanel("Peptide Level",
        h4(" "),
        fluidRow(
            column(width = 12,
                h4("Peptide Table", style = "font-weight: bold; text-decoration: underline;"),
                p("Upon selecting a specific transcript set, you can also visualize which well-binding peptides are produced from this set. The best peptide is highlighted in light green.", br(), br(),
                "Both, mutant (", code("MT"), ") and wildtype (", code("WT"), ") sequences are shown, along with either the", code("lowest"), " or ", code("median"),
                " binding affinities, depending on how you generated the aggregate report.", br(), br(),
                "An ", code("X"), "is marked for binding affinities higher than the ", code("aggregate_inclusion_binding_threshold"), " set when generating the aggregate report.", br(), br(),
                "We also include three extra columns, one specifying the mutated position(s) in the peptide, one providing information on any problematic amino acids in the mutant sequence, and one identifying whether the peptide failed the anchor criteria for any of the HLA alleles.", br(),
                "Note that if users wish to utlitize the ", strong("problematic positions"), " feature, they should run the standalone command ", code("pvacseq identify_problematic_amino_acids"),
                " or run pVACseq with the ", code("--problematic-amino-acids"), " option enabled to generate the needed information."
                )
            )
        ),
        fluidRow(
            column(width = 12,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_Peptide_Table.png",
                align = "center", width = "1500px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 12,
                h4("Anchor Heatmap", style = "font-weight: bold; text-decoration: underline;"),
                p("The ", strong("Anchor Heatmap"), "tab shows the top 30 MT/WT peptide pairs from the peptide table with anchor probabilities overlayed as a heatmap.
                The anchor probabilities shown are both allele and peptide length specific. The mutated amino acid is marked in red (for missense mutations) and each 
                MT/WT pair are separated from others using a dotted line. ", br(),
                "For peptide sequences with no overlaying heatmap, we currently do not have allele-specific predictions in our database.", br(), br(),
                "The Anchor Weights section shows a table of the per-allele per-length anchor weights for each peptide position.", br(), br(),
                "For more details and explanations regarding anchor positions and its influence on neoantigen prediction and prioritization, please refer to the next section: ",
                strong("Advanced Options: Anchor Contribution"))
            ),
            column(width = 12,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_Anchor_Heatmap.png",
                align = "center", width = "1500px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 12,
                h4("Additional Information",  style = "font-weight: bold; text-decoration: underline;"),
                h5("IC50 Plot", style = "font-weight: bold;"),
                p("By clicking on each MT/WT peptide pair, you can then assess the peptides in more detail by navigating to the ", strong("Additional Peptide Information"), " tab.", br(), br(),
                "There are five different tabs in this section of the app, providing peptide-level details on the MT/WT peptide pair that you have selected.", br(),
                "The ", strong("IC50 Plot"), "tab shows violin plots of the individual IC50-based binding affinity predictions of the MT and WT peptides for HLA
                alleles that the MT binds well to. These peptides each have up to 8 binding algorithm scores for Class I alleles or up
                to 4 algorithm scores for Class II alleles.", br())
            ),
            column(width = 12,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_IC50_Plots.png",
                align = "center", width = "1500px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 12,
                h5("%ile Plot", style = "font-weight: bold;"),
                p("The ", strong("%ile Plot"), "tab shows violin plots of the individual percentile-based binding affinity predictions of the MT and WT peptides
                for HLA alleles that the MT binds well to.", br())
            ),
            column(width = 12,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_Percentile_Plots.png",
                align = "center", width = "1500px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 12,
                h5("Binding Data", style = "font-weight: bold;"),
                p("The ", strong("Binding Data"), "tab shows the specific IC50 and percentile binding affinity predictions generated from each individual algorithm.
                Each cell shows the IC50 prediction followed by the percentile predictions in parenthesis.", br())
            ),
            column(width = 12,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_Binding_Data.png",
                align = "center", width = "1500px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 12,
                h5("Elution and Immunogenicity Table", style = "font-weight: bold;"),
                p("The ", strong("Elution and Immunigenicity Table"), "tab shows prediction results based on algorithms trained from peptide elution data. This includes algorithms
                such as NetMHCpanEL/NetMHCIIpanEL, MHCflurryELProcessing and MHCflurryELPresentation.", br())
            ),
            column(width = 12,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_Elution_Data.png",
                align = "center", width = "1500px"), br(), br()
            )
        )
    ),
    tabPanel("Advanced Options: Anchor Contribution",
        h4(" "),
        fluidRow(
            column(width = 6,
                h4("Anchor vs Mutation Positions", style = "font-weight: bold; text-decoration: underline;"),
                p("Neoantigen identification and prioritization relies on correctly predicting whether the presented 
                peptide sequence can successfully induce an immune response. As the majority of somatic mutations are single nucleotide variants,
                changes between wildtype and mutated peptides are typically subtle and require cautious interpretation. ", br(), br(),
                "In the context of neoantigen presentation by specific MHC alleles, researchers have noted that a subset of 
                peptide positions are presented to the T-cell receptor for recognition, while others are responsible for anchoring 
                to the MHC, making these positional considerations critical for predicting T-cell responses.", br(), br(),
                "Multiple factors should be considered when prioritizing neoantigens, including mutation location, anchor position, predicted MT
                and WT binding affinities, and WT/MT fold change, also known as agretopicity.", br(), br(),
                "Examples of the four distinct possible scenarios for a predicted strong MHC binding peptide involving these factors are illustrated
                in the figure on the right. There are other possible scenarios where the MT is a poor binder, however those are not listed as 
                they would not pertain to our goal of neoantigen identification.", br(), br(),
                strong("Scenario 1"), "shows the case where the WT is a poor binder and the MT peptide is a strong binder,
                containing a mutation at an anchor location. Here, the mutation results in a tighter binding of the MHC and allows for
                better presentation and potential for recognition by the TCR. As the WT does not bind (or is a poor binder), this neoantigen 
                remains a good candidate since the sequence presented to the TCR is novel.", br(), br(),
                strong("Scenario 2"), " and ", strong("Scenario 3"), " both have strong binding WT and MT peptides. In ", strong("Scenario 2"),
                ", the mutation of the peptide is located at a non-anchor location, creating a difference in the sequence participating in TCR
                recognition compared to the WT sequence. In this case, although the WT is a strong binder, the neoantigen remains a good candidate
                that should not be subject to central tolerance.", br(), br(),
                "However, as shown in ", strong("Scenario 3"), ", there are neoantigen candidates where the mutation is located at the anchor position
                and both peptides are strong binders. Although anchor positions can themselves influence TCR recognition, a mutation at a strong
                anchor location generally implies that both WT and MT peptides will present the same residues for TCR recognition. As the WT peptide
                is a strong binder, the MT neoantigen, while also a strong binder, will likely be subject to central tolerance and should not be 
                considered for prioritization.", br(), br(),
                strong("Scenario 4"), " is similar to the first scenario where the WT is a poor binder. However, in this case, the mutation is
                located at a non-anchor position, likely resulting in a different set of residues presented to the TCR and thus making the neoantigen a good candidate."
                )
            ),
            column(width = 6,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Anchor_Scenarios.png",
                    align = "center", height = "800px", width = "400px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 6,
                h4("Anchor Guidance", style = "font-weight: bold; text-decoration: underline;"),
                p("To summarize, here are the specific criteria for prioritizing (accept) and not prioritizing (reject) a neoantigen candidate:", br(),
                "Note that in all four cases, we are assuming a strong MT binder which means ",
                code("(MT IC50 < binding threshold) OR (MT percentile < percentile threshold)"), br(), br()),
                p(code("I: WT Weak binder"), " : ", code("(WT IC50 < binding threshold) OR (WT percentile < percentile threshold)")),
                p(code("II: WT Strong binder"), " : ", code("(WT IC50 > binding threshold) AND (WT percentile > percentile threshold)")),
                p(code("III: Mutation at Anchor"), " : ", code("set(All mutated positions) is a subset of set(Anchor Positions of corresponding HLA allele)")),
                p(code("IV: Mutation not at Anchor"), " : ", code("There is at least one mutated position between the WT and MT that is not at an anchor position")),
                p(strong("Scenario 1 : "), code(" I + IV"), strong(" -> Accept")),
                p(strong("Scenario 2 : "), code(" II + IV"), strong(" -> Accept")),
                p(strong("Scenario 3 : "), code(" II + III"), strong(" -> Reject")),
                p(strong("Scenario 4 : "), code(" I + III"), strong(" -> Accept"))
            ),
            column(width = 6,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/anchor.jpg",
                    align = "center", height = "350px", width = "600px"), br(), br()
            )
        )
    ),
    tabPanel("Advanced Options: Regenerate Tiering",
        h4(" "),
        fluidRow(
            column(width = 6,
                h4("Reassigning Tiers for all variants after adjusting parameters", style = "font-weight: bold; text-decoration: underline;"),
                p("The Tier column generated by pVACtools is aimed at helping users group and prioritize neoantigens in a more efficient manner.
                For details on how Tiering is done, please refer to the Variant Level tutorial tab where we break down each
                specific Tier and its criteria.", br(), br(),
                "While we try to provide a set of reasonable default parameters, we fully understand the need for flexible changes to the
                parameters used in the underlying Tiering algorithm. Thus, we provide an Advanced Options tab in our app where users can change the following
                cutoffs custom to their individual analysis: ", br(), br(),
                code("Binding Threshold"), p("IC50 cutoff for a peptide to be considered a strong binder. Note that if allele-specific binding thresholds are
                in place, those will stay the same and not overwritten by this parameter value change."), br(),
                code("Percentile Threshold"), p("Percentile cutoff for a peptide to be considered a strong binder."), br(),
                code("Clonal DNA VAF"), p("VAF cutoff that is taken into account when deciding subclonal variants. Note that variants with a DNA VAF lower
                than half of the clonal VAF cutoff will be considered subclonal (e.g. setting a 0.6 clonal VAF cutoff means anything under 0.3 VAF is subclonal)."), br(),
                code("Allele Expr"), p("Allele expression cutoff for a peptide to be considered expressed. Note for each variant, the allele expression
                is calculated by multiplying gene expression and RNA VAF."), br(),
                code("Default Anchors vs Allele-specific Anchors"), br(),
                "By default, pVACtools considers positions 1, 2, n-1, and n to be anchors for an n-mer allele. However, a recent study has shown that anchors should be
                considered on an allele-specific basis and different anchor patterns exist among HLA alleles.",
                "Here, we provide users with the option to utilize allele-specific anchors when generating the Anchor Tier. However, to objectively determine
                which positions are anchors for each individual allele, the users need to set a contribution percentage threshold (X).",
                "Per anchor calculation results from the described computational workflow in the cited paper, each position of the n-mer peptide is assigned a
                score based on how its binding to a certain HLA allele was influenced by mutations. These scores can then be used to calculate the relative
                contribution of each position to the overall binding affinity of the peptide. Given the contribution threshold X, we rank the normalized score
                across the peptide in descending order (e.g. [2,9,1,3,2,8,7,6,5] for a 9-mer peptide) and start summing the scores from top to bottom.
                Positions that together account for X% of the overall binding affinity change (e.g. 2,9,1) will be assigned as anchor locations for tiering purposes.", br(), br(),
                "However, we recommend users also navigating to the Anchor Heatmap Tab in the peptide level description for a less binary approach."
                )
            ),
            column(width = 6,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_Regenerate_Tiering.png",
                    align = "center", width = "700px"), br(), br()
            )
        ),
        fluidRow(
            column(width = 6,
                h4("Original Parameters", style = "font-weight: bold; text-decoration: underline;"),
                p(" In this box, we provide users with the original parameters they had used to generate the currently loaded aggregate report and metrics file."),
                p("Note that the app will keep track of your peptide evaluations and comments accordingly even when changing or reseting the parameters."),
                p("If you see a parameter in the original parameter box but did not see an option to change it in the advanced options section, it is likely that you
                will be required to rerun the", code("pvacseq generate-aggregate-report"), " command. This is likely due to the current metrics file not 
                having the necessary peptide information to perform this request."
                ),
                h4("Current Parameters", style = "font-weight: bold; text-decoration: underline;"),
                p(" In this box, we provide users with the tiering parameters that currently applied to the aggregate report.",
                "This not only allows users to compare their current parameters (if changed) with the original parameters."
                ),
                h4("Resetting Parameters", style = "font-weight: bold; text-decoration: underline;"),
                p("The ", strong("reset"),
                " button allows the user to restore the original tiering when desired."
                )
            ),
            column(width = 6,
                img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/Explore_Original_Parameters.png",
                    align = "center", width = "300px"), br(), br()
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
                         span(class = "logo-mini", tags$img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/pVACview_logo_mini.png")),
                         span(class = "logo-lg", tags$img(src = "https://github.com/griffithlab/pVACtools/raw/5834def4/pvactools/tools/pvacview/www/pVACview_logo.png"))
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
    ),
    div(textOutput("version"), style="margin-left:20px;position:fixed;bottom:20px;color:#b8c7ce")
  ),
  body = dashboardBody(
    use_theme(mytheme),
    tags$head(
      tags$style(HTML("table.dataTable tr.active td, table.dataTable td.active {color: black !important}")),
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
