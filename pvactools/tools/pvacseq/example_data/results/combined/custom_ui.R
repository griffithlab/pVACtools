custom_tab <- tabItem("custom",
    tabsetPanel(type = "tabs", id = "custom_tabs",
        tabPanel(title = "Upload Data", value = "custom_upload",
            fluidRow(
                column(width = 6,
                    box(
                        title = "Option 1: View Demo data", status = "primary", solidHeader = TRUE, width = NULL,
                        actionButton("loadDefault_Vaxrank", "Load demo data (VaxRank output)", style = "color: #fff; background-color: #c92424; border-color: #691111"),
                        h5("After clicking the \"Load demo data\" button, select your desired grouping and sorting parameters in the \"Choose How to Visualize Data\" panel
                           and click \"Visualize\". Please wait a couple seconds for the data to load.")
                    ),
                    box(
                        title = "Option 2: View NeoPredPipe demo data", status = "primary", solidHeader = TRUE, width = NULL,
                        actionButton("loadDefault_Neopredpipe", "Load demo data (NeoPredPipe output)", style = "color: #fff; background-color: #c92424; border-color: #691111"),
                        h5("After clicking the \"Load demo data\" button, select your desired grouping and sorting parameters in the \"Choose How to Visualize Data\" panel
                           and click \"Visualize\". Please wait a couple seconds for the data to load.")
                    ),
                    box(
                        title = "Option 3: View antigen.garnish demo data", status = "primary", solidHeader = TRUE, width = NULL,
                        actionButton("loadDefault_antigengarnish", "Load demo data (antigen.garnish output)", style = "color: #fff; background-color: #c92424; border-color: #691111"),
                        h5("After clicking the \"Load demo data\" button, select your desired grouping and sorting parameters in the \"Choose How to Visualize Data\" panel
                           and click \"Visualize\". Please wait a couple seconds for the data to load.")
                    ),
                    box(
                        title = "Option 4: Upload your own custom data files", status = "primary", solidHeader = TRUE, width = NULL,
                        HTML("<h5><b>(Required)</b> Please upload your TSV file."),
                        br(), br(),
                        uiOutput("custom_upload_ui")
                    ),
                    box(
                      title = "Choose How to Visualize Data", status = "primary", solidHeader = TRUE, width = NULL,
                      uiOutput("custom_group_by_feature_ui"),
                      h5("Group peptides together by a certain feature. For example, grouping by variant would allow user to explore all proposed peptides for one variant at a time."),
                      uiOutput("custom_order_by_feature_ui"),
                      h5("Order peptides by a certain feature. For example, ordering peptides by binding scores to find the best binders."),
                      uiOutput("custom_peptide_features_ui"),
                      h5("Choose what features you would like to consider for each group of peptides."),
                      actionButton("visualize_custom", "Visualize")
                    ),
                ),
                column(6,
                    box(
                        title = "Example Neoantigen Prediction Pipelines", status = "primary", solidHeader = TRUE, width = NULL,
                        h4("Vaxrank: A computational tool for designing personalized cancer vaccines", style = "font-weight: bold; text-decoration: underline;"),
                        h5("Therapeutic vaccines targeting mutant tumor antigens (“neoantigens”) are an increasingly popular form of personalized cancer immunotherapy.
                        Vaxrank is a computational tool for selecting neoantigen vaccine peptides from tumor mutations, tumor RNA data, and patient HLA type.
                        Vaxrank is freely available at www.github.com/openvax/vaxrank under the Apache 2.0 open source license and can also be installed from the Python Package Index."),
                        actionButton("vaxrank_help_doc", "Vaxrank Gtihub", onclick = "window.open('https://github.com/openvax/vaxrank', '_blank')"),
                        hr(style = "border-color: white"),
                        h4("NeoPredPipe: high-throughput neoantigen prediction and recognition potential pipeline", style = "font-weight: bold; text-decoration: underline;"),
                        h5("NeoPredPipe (Neoantigen Prediction Pipeline) is offered as a contiguous means of predicting putative neoantigens and their corresponding recognition potentials for
                        both single and multi-region tumor samples. This tool allows a user to process neoantigens predicted from single- or multi-region vcf files using ANNOVAR and netMHCpan."),
                        actionButton("neopredpipe_help_doc", "NeoPredPipe Gtihub", onclick = "window.open('https://github.com/MathOnco/NeoPredPipe', '_blank')"),
                        hr(style = "border-color: white"),
                        h4("antigen.garnish.2: Tumor neoantigen prediction", style = "font-weight: bold; text-decoration: underline;"),
                        h5("Human and mouse ensemble tumor neoantigen prediction from SNVs and complex variants. Immunogenicity filtering based on the Tumor Neoantigen Selection Alliance (TESLA)."),
                        actionButton("antigen_garnish_help_doc", "antigen.garnish Gtihub", onclick = "window.open('https://github.com/andrewrech/antigen.garnish', '_blank')"),
                        hr(style = "border-color: white")
                    )
                )
            )
        ),
        tabPanel(title = "Explore Data", value = "custom_explore",
            fluidRow(
                box(width = 12,
                    title="Overview of Neoantigen Features", 
                    status='primary', solidHeader = TRUE, collapsible = TRUE,
                    enable_sidebar = TRUE, sidebar_width = 25, sidebar_start_open = TRUE,
                    uiOutput("group_feature"),
                    DTOutput('customTable')%>% withSpinner(color="#8FCCFA"),
                    span("Currently investigating row: ", verbatimTextOutput("customSelected")),
                    style = "overflow-x: scroll;font-size:100%"),
          ),
          fluidRow(
            box(width = 12, title = "Detailed Data", solidHeader = TRUE, collapsible = TRUE, status = "primary",
                uiOutput("sort_feature"),
                DTOutput('customPeptideTable')%>% withSpinner(color="#8FCCFA"), style = "overflow-x: scroll;font-size:100%"
            )
          ),
          fluidRow(
            box(width = 12,
                title = "Dynamic Scatter Plot", status = "primary", solidHeader = TRUE, collapsible = TRUE,
                h4("Scatter plot to explore characteristics of data"),
                sidebarPanel(
                  #variable selection for x-axis
                  uiOutput("xvrbl_custom"),
                  uiOutput("xvrbl_log_custom"),
                  uiOutput("xvrbl_scale_custom"),
                  # variable selection for y-axis
                  uiOutput("yvrbl_custom"),
                  uiOutput("yvrbl_log_custom"),
                  uiOutput("yvrbl_scale_custom"),
                  # color
                  uiOutput("color_custom"),
                  uiOutput("min_color_custom"),
                  uiOutput("max_color_custom"),
                  # size
                  uiOutput("size_custom")
                ),
                mainPanel(
                  align = "center",
                  plotlyOutput(outputId = "scatter_custom", height = "800px") %>% withSpinner(color = "#8FCCFA"),
                )
                
            )
          )
        )
    )
)
