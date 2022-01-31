# load shiny library
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(fresh)
library(shinycssloaders)

source("styling.R")


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


ui <- dashboardPage(
  ## HEADER ##
  header = dashboardHeader(
    title=tagList(tags$a(class="logo",
                         span(class= "logo-mini", tags$img(src='pVACview_logo_mini.png')),
                         span(class= "logo-lg", tags$img(src='pVACview_logo.png'))
                  )),
    tags$li(class = "dropdown", tags$a(href = "http://pvactools.org", class = "my_class", "Help", target="_blank"))
    ),
  ## SIDEBAR ##
  sidebar = dashboardSidebar(
    sidebarMenu(
      tags$head(tags$style(csscode)),
      id="tabs",
      br(),
      menuItem("Upload", tabName = "upload", icon = icon("upload")),
      br(),
      menuItem("Visualize and Explore", tabName = "explore", icon = icon("digital-tachograph")),
      br(),
      menuItem("Export", tabName = "export", icon = icon("file-export"))
    )
  ),
  body = dashboardBody(
    
    use_theme(mytheme), 
    tags$head(
      tags$style(HTML(css)),
      tags$style(HTML('table.dataTable tr.selected td, table.dataTable td.hover {background-color: #EAF2F8 !important;}')),
      tags$style(HTML('table.dataTable { border-collapse: collapse;}')),
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
      tabItem(
        "upload",
        
        # infoBoxes
        fluidRow(
          column(width = 6,
            box(
              title="Option 1: View demo data", status='primary', solidHeader = TRUE, width = NULL,
              actionButton('loadDefaultmain', "Load demo data", style="color: #fff; background-color: #c92424; border-color: #691111"),
              h5("Please wait a couple seconds after clicking and you should be redirected to the Visualize and Explore tab.")
            ),
            box(
              title="Option 2: Upload your own data Files", status='primary', solidHeader = TRUE, width = NULL,
              h5("Please upload the aggregate report file. Note that this will be the data displayed in the main table in the Explore tab."),
              fileInput(inputId="mainDataInput", label="1. Neoantigen Candidate Aggregate Report (tsv required)", accept =  c("text/tsv",
                                                                                                                              "text/tab-separated-values,text/plain",
                                                                                                                              ".tsv")),
              radioButtons("hla_class", "Does this aggregate report file correspond to Class I or Class II prediction data?",  
                           c("Class I data (e.g. HLA-A*02:01) " = "class_i", "Class II data (e.g. DPA1*01:03)" = "class_ii")),
              
              hr(style="border-color: white"),
              h5("Please upload the corresponding metrics file for the main file that you have chosen."),
              fileInput(inputId="metricsDataInput", label="2. Neoantigen Candidate Metrics file (json required)", accept = c("application/json",
                                                                                                                             ".json")),
              hr(style="border-color: white"),
              h5("If you would like, you can upload an additional aggregate report file generated with either Class I or Class II results to supplement your main table. (E.g. if you uploaded Class I data as the main table, you can upload your Class II report here as supplemental data)"),
              fileInput(inputId="additionalDataInput", label="3. Additional Neoantigen Candidate Aggregate Report (tsv required)", accept =  c("text/tsv",
                                                                                                                                               "text/tab-separated-values,text/plain",
                                                                                                                                               ".tsv")),
              textInput("add_file_label", "Please provide a label for the additional file uploaded (e.g. Class I data or Class II data)"),
              hr(style="border-color: white"),
              h5("Additionally, you can upload a gene-of-interest list in a tsv format, where each row is a single gene name. These genes (if in your aggregate report) will be highlighted in the Gene Name column."),
              fileInput(inputId="gene_list", label="4. Gene-of-interest List (tsv required)", accept =  c("text/tsv","text/tab-separated-values,text/plain",".tsv")),
              actionButton('visualize', "Visualize")
            )
          ),
          column(6,
            box(
              title="Basic Instructions: How to explore your data using pVACview?", status='primary', solidHeader = TRUE, width = NULL,
              h4("Step 1: Upload your own data / Load demo data", style="font-weight: bold"),
              h5("You can either choose to explore a demo dataset that we have prepared from the HCC1395 cell line, or choose to upload your own datasets."),
              HTML("<h5>If you are uploading your own datasets, the two required inputs are output files you obtain after running the pVACseq pipeline. 
                 The <b>aggregated tsv file</b> is a list of all predicted epitopes and their binding affinity scores with additional variant information
                 and the <b>metrics json file</b> contains additional transcript and peptide level information.</h5>"),
              h5("You have the option of uploading an additional file to supplement the data you are exploring. This includes: additional class I or II information and
                 a gene-of-interest tsv file."),
              actionButton('help_doc_upload', "More details", onclick ="window.open('https://pvactools.readthedocs.io/en/latest/pvacview/getting_started.html#upload', '_blank')"),
              h4("Step 2: Exploring your data", style="font-weight: bold"),
              HTML("<h5>To explore the different aspects of your neoantigen candidates, you will need to navigate to the <b>Aggregate Report of Best Candidate by Variant</b> on the visualize and explore tab.
                 For detailed variant, transcript and peptide information for each candidate listed, you will need to click on the <b>Investigate button</b> for the specific row of interest.
                 This will prompt both the transcript and peptide table to reload with the matching information.</h5>"),
              h5("By hovering over each column header, you will be able to see a brief description of the corresponding column and for more details, you can click on the tooltip located at the top right of the aggregate report table.", br(),
              "After investigating each candidate, you can label the candidate using the dropdown menu located at the second to last column of the table. Choices include:
                 Accept, Reject or Review."),
              actionButton('help_doc_explore', "More details", onclick ="window.open('https://pvactools.readthedocs.io/en/latest/pvacview/getting_started.html#visualize-and-explore', '_blank')"),
              h4("Step 3: Exporting your data", style="font-weight: bold"),
              h5("When you have either finished ranking your neoantigen candidates or need to pause and would like to save your current evaluations, 
                 you can export the current main aggregate report using the export page. "),
              HTML("<h5>Navigate to the export tab, and you will be able to name your file prior to downloading in either tsv or excel format.
                 The excel format is user-friendly for downstream visualization and manipulation. However, if you plan on to continuing editing the aggregate report 
                 and would like to load it back in pVACview with the previous evaluations preloaded, you will <b>need</b> to download the file in a <b>tsv format</b>. 
                 This serves as a way to save your progress as your evaluations are cleared upon closing or refreshing the pVACview app.</h5>"),
              actionButton('help_doc_export', "More details", onclick ="window.open('https://pvactools.readthedocs.io/en/latest/pvacview/getting_started.html#export', '_blank')")
            )
          ),
        )
      ),
      ## EXPLORE TAB ##
      tabItem(
        "explore",
        fluidRow(
          tags$style(
            type = 'text/css',
            '.modal-dialog { width: fit-content !important; }'
          ),
          tags$script(HTML("Shiny.addCustomMessageHandler('unbind-DT', function(id) {
          Shiny.unbindAll($('#'+id).find('table').DataTable().table().node());
                           })")),
          
          box(width= 6, 
              title="Advanced Options: Regenerate Tiering with different parameters",
              status='primary', solidHeader = TRUE, collapsible = TRUE, collapsed = TRUE,
              "*Please note that the metrics file is required in order to regenerate tiering information with different parameters", br(),
              "Current version of pVACseq results defaults to positions 1, 2, n-1 and n (for a n-mer peptide) when determining anchor positions.
              If you would like to use our allele specific anchor results and regenerate the tiering results for your variants,
              please specify your contribution cutoff and submit for recalculation. ", tags$a(href="https://www.biorxiv.org/content/10.1101/2020.12.08.416271v1", "More details can be found here."), br(),
              checkboxInput("use_anchor", "If you want to use allele-specific anchor calculations, please check this box. Otherwise anchors will be calculated as 1,2 and n-1,n for n-mer peptides.", value = FALSE, width = NULL),
              sliderInput("anchor_contribution", "Contribution cutoff for determining anchor locations", 0.5, 0.9, 0.8, step = 0.1, width = 400),
              numericInput("dna_cutoff", "Clonal DNA VAF (Anything lower than 1/2 of chosen VAF level will be considered subclonal)", 0.5, min = 0, max = 1, step = 0.01, width = 500),
              #numericInput("rna_cutoff", "RNA low gene expression cutoff (Anything lower than chosen expression level will be considered low expression)", 1, min = 0, max = 100, step = 1, width = 500),
              h5("For your reference, the max DNA VAF under 0.6 in the current main table is: "), verbatimTextOutput("max_dna"), br(),
              numericInput("allele_expr_high", "Allele Expression cutoff to be considered a Pass variant (default: 3)", 3, min = 0, max = 100, step = 0.1, width = 500),
              numericInput("allele_expr_low", "Allele Expression cutoff to be considered a Relaxed variant (default: 1). Note that this criteria is also used in determining Anchor and Subclonal variants.", 1, min = 0, max = 100, step = 0.1, width = 500),
              h5("For further explanations on these inputs, please refer to the ", tags$a(href="https://pvactools.readthedocs.io/en/latest/pvacview/getting_started.html#visualize-and-explore", "pVACview documentation.")),
              actionButton('submit','Recalculate Tiering with new parameters'),
              style = "overflow-x: scroll;font-size:100%"),
          
          box(width= 6, 
              title="Add Comments for selected variant",
              status='primary', solidHeader = TRUE, collapsible = TRUE,
              textAreaInput("comments", "Please add/update your comments for the variant you are currently examining", value=""),
              actionButton('comment','Update Comment Section'),
              h5("Comment:"), verbatimTextOutput("comment_text"),
              style = "overflow-x: scroll;font-size:100%")
        ),
        
        fluidRow(
          box(width= 12, 
                  title="Aggregate Report of Best Candidates by Variant", 
                  status='primary', solidHeader = TRUE, collapsible = TRUE,
                  enable_sidebar = TRUE, sidebar_width = 25, sidebar_start_open = TRUE,
                  dropdownMenu = boxDropdown(boxDropdownItem("Help", id = "help", icon = icon("question-circle"))),
              DTOutput('mainTable')%>% withSpinner(color="#8FCCFA"),
              span("Currently investigating row: ", verbatimTextOutput('selected')),
              style = "overflow-x: scroll;font-size:100%")
        ),

        fluidRow(
          tabBox(width = 6, title="Variant Information", 
                 tabPanel("Transcripts of Selected Variant",
                      DTOutput('transcriptsTable')%>% withSpinner(color="#8FCCFA"), style = "overflow-x: scroll;font-size:100%"),
                 tabPanel("Additional Data",
                         span("Additional Data Type: ", verbatimTextOutput('type_text')),
                         span("Median MT IC50: ", verbatimTextOutput('addData_IC50')),
                         span("Median MT Percentile: ", verbatimTextOutput('addData_percentile')) )
          ),
          box(width = 4, solidHeader = TRUE, title="Variant & Gene Info",
                 span("DNA VAF", verbatimTextOutput('metricsTextDNA')),
                 span("RNA VAF", verbatimTextOutput('metricsTextRNA')),
                 span("Gene Expression", verbatimTextOutput('metricsTextGene')),
                 span("Genomic Information (chromosome - start - stop - ref - alt)", verbatimTextOutput('metricsTextGenomicCoord')),
                 h5("Additional variant information:"),
                 uiOutput("url"),style = "overflow-x: scroll;font-size:100%"),
          box(width = 2, solidHeader = TRUE, title="Peptide Evalutation Overview",
                 tableOutput("checked"), style = "overflow-x: scroll;font-size:100%")
        ),
        fluidRow(
          box(width = 12, title="Peptide Candidates from Selected Transcript", status='primary', solidHeader = TRUE, collapsible = TRUE,
                 DTOutput('peptideTable')%>% withSpinner(color="#8FCCFA"), style = "overflow-x: scroll;font-size:100%")
        ),
        fluidRow(
          tabBox(
            title = "Additional Info",
            id = "info",

            tabPanel("MHC Binding Prediction Scores (IC50)", 
              plotOutput(outputId = "bindingData_IC50")%>% withSpinner(color="#8FCCFA"), style = "overflow-x: scroll;"
            ),
            tabPanel("MHC Binding Prediction Scores (%ile)", 
                     plotOutput(outputId = "bindingData_percentile")%>% withSpinner(color="#8FCCFA"), style = "overflow-x: scroll;"
            ),
            tabPanel("MHC Binding Predictions Table", 
                     DTOutput(outputId = "bindingDatatable"), style = "overflow-x: scroll;"
            ),
            tabPanel("Allele Specific Anchor Prediction Heatmap",
                     plotOutput(
                       outputId = "peptideFigureLegend", height = "50px"),
                     plotOutput(
                       outputId = "anchorPlot"
                     )%>% withSpinner(color="#8FCCFA"), style = "overflow-x: scroll;"
            )
          ),
          box(
            column(width=4,
            h4("Allele Specific Anchor Prediction Heatmap"),
            h5(" This tab displays HLA allele specific anchor predictions overlaying good-binding peptide sequences generated from each specific transcript.", br(),
               " Current version supports the first 15 MT/WT peptide sequence pairs (first 30 rows of the peptide table)."), br(),
            h4("MHC Binding Prediction Scores"), 
            h5(" This tab contains violin plots that showcase individual binding prediction scores from each algorithm used. A solid line is used to represent the median score.
               ")
            ),
            column(width=8, 
                   box(title = 'Anchor vs Mutation position Scenario Guide', collapsible = TRUE, collapsed = FALSE, width = 12,
                       img(src='anchor.jpg', align = "center", height='350px', width='600px'), style = "overflow-x: scroll;"))
          )
          )
      ),
      ## EXPORT TAB ##
      tabItem(
        "export",
        fluidRow(
          textInput("exportFileName", "Export filename: ", value = "Annotated.Neoantigen_Candidates", width = NULL, placeholder = NULL)
        ),
        fluidRow(
          column(12,
                 DTOutput('ExportTable')%>% withSpinner(color="#8FCCFA"))
        )
      )
    )
  )
)








