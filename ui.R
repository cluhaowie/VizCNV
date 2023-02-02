ui <- dashboardPage(
  dashboardHeader(title = "GREGoR VizCNV"),
  dashboardSidebar(sidebarMenu(
    #        menuItem(tabName = "input", text = "Input", icon = icon("upload")),
    menuItem(tabName = "zoom", text = "Zoom", icon = icon("search"))
    #        menuItem(tabName = "qc", text = "QC", icon = icon("check")),
    #        menuItem(tabName = "wgv", text = "Genome", icon = icon("info")),
    #        menuItem(tabName = "pheno", text = "HPO", icon = icon("address-card"))
  )),
  dashboardBody(
    fluidRow(
      tabItems(
        #Tab1
        #------
        tabItem(
          tabName = "input",
          h1("Upload")
        ),
        #-------
        #Tab2
        #--------
        tabItem(tabName = "qc",
                h1("Quality Control")
        ),
        #-------
        #Tab3_wgv
        #-------
        tabItem(tabName = "wgv",
                h1("Genome wide view")
        ),
        #-------
        #Tab4_zoom
        #-------
        tabItem(tabName = "zoom",
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
                        column(4,uiOutput("dnSNV_ui"))
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
                               ))
                      ),
                      fluidRow(
                        useWaiter(),
                        uiOutput("btn_filter_plot")
                      ))),
                # fluidRow(
                #     box(title = "plot parameter",width = 10,solidHeader = T, status = "primary",collapsible = T,
                #         fluidRow(
                #            # column(1,actionButton("btn_filter","Filter")),
                #             column(1,actionButton("btn_plot","Plot")))
                #         )),
                
                fluidRow(
                  box(title = "Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                      DT::dataTableOutput("filter_sv_table"),
                      fluidRow(
                        column(1,uiOutput("ui_dlbtn_tbl"))
                      ),
                      DT::dataTableOutput("Select_table"))),
                fluidRow(
                  box(title ="Plot",width = 12,solidHeader = T, status = "success",collapsible = T,
                      fluidRow(
                        column(1,uiOutput("ui_dlbtn_goto"))
                      ),
                      fluidRow(
                        plotOutput(
                          "plot1",
                          height = 400,
                          dblclick = "plot1_dblclick",
                          brush = brushOpts(id = "plot1_brush",direction = "x",
                                            resetOnNew = TRUE))),
                      fluidRow(                                     
                        plotOutput(
                          "plot2",
                          height = 400,
                          dblclick = "plot2_dblclick",
                          brush = brushOpts(id = "plot2_brush",direction = "x",
                                            resetOnNew = TRUE))),
                      fluidRow(
                        column(1,uiOutput("ui_dlbtn_plt")),
                        column(1,uiOutput("ui_clbtn_plt")),
                        column(1,uiOutput("ui_dlbtn_dnsnv"))
                      )
                  )),
        ),
        
        #--------
        #Tab5_hpo
        #-------
        tabItem(tabName = "pheno",
                h1("Phenotype based priotization of SV")
        )
      )
    )
  )
)




#shiny::shinyApp(ui, server)





