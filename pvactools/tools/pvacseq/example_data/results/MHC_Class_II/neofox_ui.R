# ui.R
library(shiny)
library(plotly)

neofox_tab <- tabItem("neofox",
    tabsetPanel(type = "tabs", id = "neofox_tabs",
        tabPanel(title = "Upload Data", value = "neofox_upload",
            fluidRow(
                column(width = 6,
                    box(
                        title = "Option 1: View NeoFox demo data", status = "primary", solidHeader = TRUE, width = NULL,
                        actionButton("loadDefaultneofox", "Load neofox output data for HCC1395", style = "color: #fff; background-color: #c92424; border-color: #691111"),
                        h5("Please wait a couple seconds after clicking for the data to load.")
                    ),
                    box(
                        title = "Option 2: Upload your own neofox data files", status = "primary", solidHeader = TRUE, width = NULL,
                        HTML("<h5><b>(Required)</b> Please upload your neofox output file. This file should be a table generated by NeoFox with the suffix “_neoantigen_candidates_annotated.tsv“"),
                        br(), br(),
                        uiOutput("neofox_upload_ui"),
                        actionButton("visualize_neofox", "Visualize")
                    )
                ),
                column(6,
                    box(
                        title = "NeoFox (NEOantigen Feature toolbOX)", status = "primary", solidHeader = TRUE, width = NULL,
                        h5("NeoFox (NEOantigen Feature toolbOX) is a python package that annotates a given set of neoantigen candidate sequences with relevant neoantigen features."),
                        h5("The tool covers neoepitope prediction by MHC binding and ligand prediction, similarity/foreignness of a neoepitope candidate sequence, combinatorial features and machine learning approaches by
                        running a wide range of published toolsets on the given input data. For more detailed information on the specific neoantigen-related algorithms and how to generate your own NeoFox results, please
                        refer to the link below: "), br(),
                        actionButton("neofox_help_doc", "NeoFox Website", onclick = "window.open('https://neofox.readthedocs.io/en/latest/index.html', '_blank')")
                    )
                )
            )
        ),
        tabPanel(title = "Explore Data", value = "neofox_explore",
            fluidRow(
                box(width = 4, solidHeader = TRUE, title = "Peptide Evaluation Overview", status = "primary",
                    tableOutput("neofox_checked"), style = "overflow-x: scroll;font-size:100%"),
                box(width = 8,
                    title = "Add Comments for last selected variant(s)",
                    status = "primary", solidHeader = TRUE, collapsible = TRUE,
                    textAreaInput("neofox_comments", label = "Please add/update your comments for the selected variant(s)", value = ""),
                    actionButton("neofox_comment", "Update Comment Section"),
                    h5("Comment:"), tableOutput("neofox_comment_text"),
                    style = "font-size:100%")
            ),
            fluidRow(
                    box(width = 12,
                        title = "Annotated Neoantigen Candidates using NeoFox",
                        status = "primary", solidHeader = TRUE, collapsible = TRUE,
                        enable_sidebar = TRUE, sidebar_width = 25, sidebar_start_open = TRUE,
                        #selectInput("neofox_page_length", "Number of entries displayed per page:", selected = "10", c("10", "20", "50", "100"), width = "280px"),
                        DTOutput("neofoxTable") %>% withSpinner(color = "#8FCCFA"),
                        span("Currently investigating row(s): ", verbatimTextOutput("neofox_selected")),
                        style = "overflow-x: scroll;font-size:100%",
                        "* indicates variable of interest designated by authors"
                        )
            ),
            fluidRow(
              box(width = 12, 
                  title = "Comparative Violin Plots",  status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  h4("Violin Plots showing distribution of various neoantigen features for selected variants."),
                  uiOutput("noefox_features_ui"),
                  plotOutput(outputId = "neofox_violin_plots_row1") %>% withSpinner(color = "#8FCCFA"),
                  "* indicates variable of interest designated by authors"
                  )
            ),
            fluidRow(
              box(width = 12,
                  title = "Dynamic Scatter Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                  h4("Scatter plot to explore characteristics of data"),
                  sidebarPanel(
                      # variable selection for x-axis
                      uiOutput("xvrbl"),
                      uiOutput("xvrbl_log"),
                      uiOutput("xvrbl_scale"),
                      # variable selection for y-axis
                      uiOutput("yvrbl"),
                      uiOutput("yvrbl_log"),
                      uiOutput("yvrbl_scale"),
                      # color
                      uiOutput("color_neofox"),
                      uiOutput("min_color"),
                      uiOutput("max_color"),
                      # size
                      uiOutput("size_neofox"),
                      "* indicates variable of interest designated by authors"
                    ),
                    mainPanel(
                      align = "center",
                      plotlyOutput(outputId = "scatter", height = "800px") %>% withSpinner(color = "#8FCCFA"),
                    )
              )
            )
        ),
        ## EXPORT TAB ##
        tabPanel(title = "Export Data", value = "neofox_export",
          fluidRow(
            textInput("exportNeofoxFileName", "Export filename: ", value = "Annotated.Neoantigen_Candidates", width = NULL, placeholder = NULL)
          ),
          fluidRow(
            column(12,
                   DTOutput("NeofoxExportTable") %>% withSpinner(color = "#8FCCFA"))
          )
        )
    )
)
