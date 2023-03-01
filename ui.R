
source("./mod/mod_anno.R")

ui <- dashboardPage(
  
  dashboardHeader(title = "GREGoR VizCNV"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem(tabName = "input", text = "Input/Filtering", icon = icon("upload")),
      menuItem(tabName = "plot", text = "Plots", icon = icon("search")),
      menuItem(tabName = "table", text = "Tables", icon = icon("table"))
    )
  ),
  
  dashboardBody(
      tabItems(
        tabItem(tabName = "input",
          fluidRow(
            box(title="File upload or select",status="primary",width = 12,solidHeader = T,collapsible = T,
                fluidRow(
                  shinyWidgets::switchInput(inputId = "file_source",onLabel = "Local",offLabel = "Server"),
                  box(width = 6,uiOutput("file_source_ui1")),
                  box(width = 6,uiOutput("file_source_ui2"))
                ))
          ),
          fluidRow(
            box(title="Filter parameters",status="primary",width = 12,solidHeader = T,collapsible = T,
                
                fluidRow(
                  column(4, 
                         p(HTML("<b>Chromosome</b>"),span(shiny::icon("info-circle"), id = "info_chr"),selectizeInput('chr', choice=NULL,label = NULL,options = list(maxItems = 1,placeholder="chr6")),
                           tippy::tippy_this(elementId = "info_chr",tooltip = "Selected chromosome",placement = "right")
                         )),
                  column(4,uiOutput("blt_dnSNV_ui"))
                ),
                fluidRow(
                  column(4, 
                         p(HTML("<b>Segment option</b>"),span(shiny::icon("info-circle"), id = "info_seg"),radioButtons('seg_option', 
                                                                                                                        label = NULL,
                                                                                                                        choiceNames = list("SLM","CBS"),
                                                                                                                        choiceValues = list("slm","cbs")),
                           tippy::tippy_this(elementId = "info_seg",tooltip = "Select options for segment. 
                                                               SLM is fast and tend to give more segment,can be used for high-quality data; 
                                                               CBS is slow and give less segment, can be used for noisy data ",placement = "right")
                         )),
                  column(4, 
                         p(HTML("<b>Choose segment to be included</b>"),span(shiny::icon("info-circle"), id = "info_include"),checkboxGroupInput(inputId="include_seg",
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
                box(title = "Plot Options",width = 12,solidHeader = T, status = "success",collapsible = T,
                  use_waiter(),
                  column(2, 
                         actionButton("btn_plot", "Plot"),
                         actionButton("btn_anno", "Annotate")),
                  column(6,
                         mod_anno_check_UI("IDR"),
                         mod_anno_check_UI("SegDup"),
                         mod_anno_check_UI("OMIM"),
                         mod_anno_check_UI("gnomAD"),
                         mod_anno_check_UI("RMSK"))

                  ),
                
                   
                fluidRow(
                      column(width = 4,uiOutput("ui_btn_goto")),
                      column(width = 4,uiOutput("ui_btn_add")),
                      column(width = 4,verbatimTextOutput("brush_info"))),
                    plotOutput(
                        "plot1",
                        height = 400,
                        dblclick = "plot1_dblclick",
                        brush = brushOpts(id = "plot1_brush",direction = "x",
                                          resetOnNew = TRUE)),
                    uiOutput("ui_plot_snp"),
                    uiOutput("ui_plot_anno"),
                    mod_anno_plot_switch_UI("IDR"),
                    mod_anno_plot_switch_UI("Segdup"),
                    mod_anno_plot_switch_UI("OMIM"),
                    mod_anno_plot_switch_UI("gnomAD"),
                    mod_anno_plot_switch_UI("RMSK"),
                    fluidRow(column(1,uiOutput("ui_dlbtn_plt")),
                             column(1,uiOutput("ui_clbtn_plt")),
                             column(1,uiOutput("ui_dlbtn_dnsnv")))
                 ),
        tabItem(tabName = "table",
                box(title = "SV Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                      DT::dataTableOutput("filter_sv_table"),
                      fluidRow(
                        column(1,uiOutput("ui_dlbtn_tbl"))
                      ),
                      DT::dataTableOutput("Select_table")),
                box(title = "Annotation Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                tabsetPanel(
                  tabPanel("IDR",anno_table_UI("IDR")),
                  tabPanel("SegDup",anno_table_UI("SegDup")),
                  tabPanel("OMIM",anno_table_UI("OMIM")),
                  tabPanel("gnomAD",anno_table_UI("gnomAD")),
                  tabPanel("rmsk",anno_table_UI("rmsk"))
                )
              )
        )
      )
    )
  )


