ui <- dashboardPage(
  dashboardHeader(title = "VizCNV",
                  leftUi = tagList(
                    dropdownBlock(
                      id = "input_dropdown_upload",
                      title = "Upload",
                      icon = icon("upload"),
                      uiOutput("file_source_ui2"))
                    )),
  dashboardSidebar(sidebarMenu(
    menuItem(tabName = "zoom", text = "Zoom", icon = icon("search"))
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
                        box(width = 6,uiOutput("file_source_ui1"))
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
                        column(4, 
                               p(HTML("<b>Annotation option</b>"),span(shiny::icon("info-circle"), id = "info_anno"),checkboxGroupInput(inputId="include_anno",
                                                                                                                                        label = NULL,
                                                                                                                                        c("GRCh38 v1.0 MANE refseq"="hg38_gene"),selected = "hg38_gene"),
                                 tippy::tippy_this(elementId = "include_anno",tooltip = "Include MANE 1.0 transcripts
                                                   in the annotation track",placement = "right")
                               ))
                      ),
                      fluidRow(
                        use_waiter(),
                        uiOutput("btn_filter_plot")
                      ))),
                fluidRow(
                  box(title = "Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                      DT::dataTableOutput("filter_sv_table"),
                      fluidRow(
                        column(1,uiOutput("ui_dlbtn_tbl"))
                      ),
                      DT::dataTableOutput("Select_table"))),
                fluidRow(
                  box(title ="Plot",width = 12,solidHeader = T, status = "success",collapsible = T,
                      fluidRow(uiOutput("ui_plot_anno")),
                      fluidRow(
                        column(width = 4,uiOutput("ui_btn_goto")),
                        column(width = 4,uiOutput("ui_btn_add")),
                        column(width = 4,verbatimTextOutput("brush_info"))),
                      fluidRow(
                        plotOutput(
                          "plot1",
                          height = 400,
                          dblclick = "plot1_dblclick",
                          brush = brushOpts(id = "plot1_brush",direction = "x",
                                            resetOnNew = TRUE))),
                      fluidRow(uiOutput("ui_plot_snp")),
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





