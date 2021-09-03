# load shiny library
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(fresh)
library(shinycssloaders)

source("styling.R")

ui <- dashboardPage(
  #skin = "black",
  #, div(class= "logo-lg", tags$img(src='pVACview_logo.png'))
  #tags$a(href='http://pvactools.org',
  header = dashboardHeader(
    title=tagList(tags$a(href='http://pvactools.org', class="logo",
                         span(class= "logo-mini", tags$img(src='pVACview_logo_mini.png')),
                         span(class= "logo-lg", tags$img(src='pVACview_logo.png'))
                  ))
    ),
  sidebar = dashboardSidebar(
    sidebarMenu(
      id="tabs",
      menuItem("Upload", tabName = "upload", icon = icon("upload")),
      menuItem("Visualize and Explore", tabName = "explore", icon = icon("digital-tachograph")),
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
      tabItem(
        "upload",
        
        # infoBoxes
        fluidRow(
          box(
            title="Load Default Data", status='primary', solidHeader = TRUE,
            actionButton('loadDefaultmain', "Load Default Data", style="color: #fff; background-color: #c92424; border-color: #691111")
          )
        ),
        fluidRow(
          box(
            title="Upload Data Files", status='primary', solidHeader = TRUE,
            h5("Please upload the aggregate report file. Note that this will be the data displayed in the main table in the Explore tab."),
            fileInput(inputId="mainDataInput", label="Neoantigen Candidate Aggregate Report (tsv required)", accept =  c("text/tsv",
                                                                                                                         "text/tab-separated-values,text/plain",
                                                                                                                         ".tsv")),
            radioButtons("hla_class", "Does this aggregate report file correspond to Class I or Class II prediction data?",  
                         c("Class I data (e.g. HLA-A*02:01) " = "class_i", "Class II data (e.g. DPA1*01:03)" = "class_ii")),
            
            hr(style="border-color: white"),
            h5("Please upload the corresponding metrics file for the main file that you have chosen."),
            fileInput(inputId="metricsDataInput", label="Neoantigen Candidate Metrics file (json required)", accept = c("application/json",
                                                                                                                        ".json")),
            hr(style="border-color: white"),
            h5("If you would like, you can upload an additional aggregate report file generated with either Class I or Class II results to supplement your main table. (E.g. if you uploaded Class I data as the main table, you can upload your Class II report here as supplemental data)"),
            fileInput(inputId="additionalDataInput", label="Additional Neoantigen Candidate Aggregate Report (tsv required)", accept =  c("text/tsv",
                                                                                                                               "text/tab-separated-values,text/plain",
                                                                                                                               ".tsv")),
            textInput("add_file_label", "Please provide a label for the additional file uploaded (e.g. Class I data or Class II data)"),
            actionButton('visualize', "Visualize")
          ),
        )
      ),
      
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
              title="Regenerate Tiering with different parameters",
              status='primary', solidHeader = TRUE, collapsible = TRUE,
              "*Please note that the metrics file is required in order to regenerate tiering information with different parameters", br(),
              sliderInput("anchor_contribution", "Contribution cutoff for determining anchor locations", 0.5, 0.9, 0.8, step = 0.1, width = 400),
              numericInput("dna_cutoff", "DNA VAF clonal VAF (Anything lower than 1/2 of chosen VAF level will be considered subclonal)", 0.5, min = 0, max = 1, step = 0.01, width = 500),
              #numericInput("rna_cutoff", "RNA low gene expression cutoff (Anything lower than chosen expression level will be considered low expression)", 1, min = 0, max = 100, step = 1, width = 500),
              actionButton('submit','Recalculate Tiering with new parameters'),
              style = "overflow-x: scroll;font-size:100%")
        ),
        
        fluidRow(
          box(width= 12, 
                  title="Aggregate Report of Best Candidates by Mutation", 
                  status='primary', solidHeader = TRUE, collapsible = TRUE,
                  enable_sidebar = TRUE, sidebar_width = 25, sidebar_start_open = TRUE,
                  dropdownMenu = boxDropdown(boxDropdownItem("Help", id = "help", icon = icon("question-circle"))),
              DTOutput('mainTable')%>% withSpinner(color="#8FCCFA"),
              span("Currently investigating row: ", verbatimTextOutput('selected')),
              style = "overflow-x: scroll;font-size:100%")
        ),

        fluidRow(
          tabBox(width = 6, title="Variant Information", 
                 tabPanel("Transcripts of Selected Mutation",
                      DTOutput('transcriptsTable')%>% withSpinner(color="#8FCCFA"), style = "overflow-x: scroll;font-size:100%"),
                 tabPanel("Additional Data",
                         span("Additional Data Type: ", verbatimTextOutput('type_text')),
                         span("Median MT IC50: ", verbatimTextOutput('addData_IC50')),
                         span("Median MT Percentile: ", verbatimTextOutput('addData_percentile')) )
          ),
          box(width = 4, solidHeader = TRUE, title="Mutation & Gene Info",
                 span("DNA VAF", verbatimTextOutput('metricsTextDNA')),
                 span("RNA VAF", verbatimTextOutput('metricsTextRNA')),
                 span("Gene Expression", verbatimTextOutput('metricsTextGene')),
                 span("Genomic Information (chromosome - start - stop - ref - alt)", verbatimTextOutput('metricsTextGenomicCoord')), style = "overflow-x: scroll;font-size:100%"),
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








