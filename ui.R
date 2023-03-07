

source("./mod/mod_plot_output.R")
source("./mod/mod_dnCNV.R")

ui <- dashboardPage(
  
  dashboardHeader(title = "VizCNV-dev"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem(tabName = "input", text = "Input/Filtering", icon = icon("upload")),
      menuItem(tabName = "plot", text = "Plots", icon = icon("search")),
      menuItem(tabName = "table", text = "Tables", icon = icon("table")),
      menuItem(tabName = "help", text = "Help", icon = icon("question-circle")),
      menuItem(tabName = "about", text = "About", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
      tabItems(
        tabItem(tabName = "input",
          fluidRow(box(title=" ",status="primary",width = 12,solidHeader = T,collapsible = T,
                       shinyWidgets::switchInput(inputId = "file_source",onLabel = "Local",offLabel = "Server", labelWidth = 20, handleWidth = 300),
                       radioButtons(inputId = "ref",label = h3("Reference Genome Build"),choices = list("GRCh37"="GRCh37","GRCh38"="GRCh38"),inline = T,selected = "GRCh37")
                    ),
          box(title="File upload",status="primary",width = 12,solidHeader = T,collapsible = T,
              fluidRow(
                box(width = 6,uiOutput("file_source_ui1")),
                box(width = 6,uiOutput("file_source_ui2"))
              )),
          box(title="Filter parameters",status="primary",width = 12,solidHeader = T,collapsible = T,
                column(2, 
                       p(HTML("<b>Chromosome</b>"),span(shiny::icon("info-circle"), id = "info_chr"),selectInput('chr', choice=c(paste0("chr", c(seq(1,22), "X"))),label = NULL, multiple = F),
                         tippy::tippy_this(elementId = "info_chr",tooltip = "Selected chromosome",placement = "right")
                       )),

              fluidRow(
                column(2,
                       p(HTML("<b>Segment option</b>"),span(shiny::icon("info-circle"), id = "info_seg"),radioButtons('seg_option', 
                                                                                                                      label = NULL,
                                                                                                                      choiceNames = list("SLM","CBS"),
                                                                                                                      choiceValues = list("slm","cbs")),
                         tippy::tippy_this(elementId = "info_seg",tooltip = "Select options for segment. 
                                                             SLM is fast and tend to give more segment,can be used for high-quality data; 
                                                             CBS is slow and give less segment, can be used for noisy data ",placement = "right")
                       )),
                column(3, 
                       p(HTML("<b>Segments to be included</b>"),span(shiny::icon("info-circle"), id = "info_include"),checkboxGroupInput(inputId="include_seg",
                                                                                                                                               label = NULL,
                                                                                                                                               c("index case"="Proband",
                                                                                                                                                 "Mom"="Mother",
                                                                                                                                                 "Dad"="Father"),selected = "Proband"),
                         tippy::tippy_this(elementId = "info_include",tooltip = "Choose to show segment from either or both parents",placement = "right")
                       )),
              ),
              fluidRow(
                use_waiter(),
                column(1, actionButton("btn_filter", "Filter"))
              )
            )
          )
        ),
        tabItem(tabName = "plot",
                shinyjs::useShinyjs(),
                fluidRow(
                  box(title = "Customizations",width = 12,solidHeader = T, status = "success",collapsible = T,
                    use_waiter(),
                    column(4, 
                           HTML("<b>Basic Plot Options: </b>"),
                           mod_checkbox_UI("RD-static"),
                           mod_checkbox_UI("RD-dynamic"),
                           mod_checkbox_UI("Baf-B_allele"),
                           mod_checkbox_UI("Baf-A_allele", value = F),
                           uiOutput("blt_dnSNV_ui"),
                           actionButton("btn_plot", "Plot")),
                    column(6,
                           HTML("<b>Annotation Track Options:</b>"),
                           mod_checkbox_UI("RefSeq"),
                           mod_checkbox_UI("IDR", value = F),
                           mod_checkbox_UI("SegDup"),
                           mod_checkbox_UI("OMIM"),
                           mod_checkbox_UI("gnomAD"),
                           mod_checkbox_UI("RMSK", value = F),
                           actionButton("btn_anno", "Annotate")),
                    actionButton("btn_dnCNV", "Show potential dnCNVs")
                    )
                  ),
                  fluidRow(box(title = "Controls",width = 12,solidHeader = T, status = "success",collapsible = T,

                             column(4,shiny::textInput("goto_reg",label = NULL,placeholder = "gene, chromosome range")),
                             column(2,shiny::actionButton("btn_go","go")),
                             verbatimTextOutput("cur_range")
                          
                        )
                      ),
                  fluidRow(box(title = "Plots",width = 12,solidHeader = T, status = "success",collapsible = T,
                    mod_plot_switch_UI("RD-static", height = 200),
                    mod_plot_switch_UI("RD-dynamic", height = 200),
                    mod_plot_switch_UI("Baf-B_allele", height = 200),
                    mod_plot_switch_UI("Baf-A_allele", height = 200),
                    # uiOutput("ui_plot_anno"),
                    mod_plot_switch_UI("RefSeq", height = 120),
                    mod_plot_switch_UI("IDR"),
                    mod_plot_switch_UI("Segdup"),
                    mod_plot_switch_UI("OMIM"),
                    mod_plot_switch_UI("gnomAD"),
                    mod_plot_switch_UI("RMSK"),
                    fluidRow(column(1,uiOutput("ui_dlbtn_plt")),
                             column(1,uiOutput("ui_clbtn_plt")),
                             column(1,uiOutput("ui_dlbtn_dnsnv")))
                  )
                 )
               ),
        tabItem(tabName = "table",
                fluidRow(box(title = "dnCNV Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                             mod_dnCNV_UI("dnCNV")),
                box(title = "SV Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                      DT::dataTableOutput("filter_sv_table"),
                      fluidRow(
                        column(1,uiOutput("ui_dlbtn_tbl"))
                      ),
                      DT::dataTableOutput("Select_table")),
                box(title = "Annotation Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                tabsetPanel(
                  tabPanel("RefSeq",anno_table_UI("RefSeq")),
                  tabPanel("IDR",anno_table_UI("IDR")),
                  tabPanel("SegDup",anno_table_UI("SegDup")),
                  tabPanel("OMIM",anno_table_UI("OMIM")),
                  tabPanel("gnomAD",anno_table_UI("gnomAD")),
                  tabPanel("rmsk",anno_table_UI("rmsk"))
                )
              )
            )
        ),
        tabItem(tabName = "help",
                "Lorem ipsum dolor sit amet, consectetur adipiscing elit, sed do eiusmod tempor incididunt ut labore et dolore magna aliqua. Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat. Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur. Excepteur sint occaecat cupidatat non proident, sunt in culpa qui officia deserunt mollit anim id est laborum"),
        tabItem(tabName = "about",
                "Sed ut perspiciatis unde omnis iste natus error sit voluptatem accusantium doloremque laudantium, totam rem aperiam, eaque ipsa quae ab illo inventore veritatis et quasi architecto beatae vitae dicta sunt explicabo. Nemo enim ipsam voluptatem quia voluptas sit aspernatur aut odit aut fugit, sed quia consequuntur magni dolores eos qui ratione voluptatem sequi nesciunt. Neque porro quisquam est, qui dolorem ipsum quia dolor sit amet, consectetur, adipisci velit, sed quia non numquam eius modi tempora incidunt ut labore et dolore magnam aliquam quaerat voluptatem. Ut enim ad minima veniam, quis nostrum exercitationem ullam corporis suscipit laboriosam, nisi ut aliquid ex ea commodi consequatur? Quis autem vel eum iure reprehenderit qui in ea voluptate velit esse quam nihil molestiae consequatur, vel illum qui dolorem eum fugiat quo voluptas nulla pariatur?")
    )
  )

)

