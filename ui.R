source("./mod/mod_plot_output.R")
source("./mod/mod_dnCNV.R")
source("./mod/mod_upload.R")
source("./mod/mod_UCSC.R")

ui <- dashboardPage(
  
  dashboardHeader(title = HTML("VizCNV-dev"),
                  tags$span(style = "width: 90%"),
                  tagList(
                    tags$a(href="https://github.com/cluhaowie/VizCNV",icon("github"))
                    )
                  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem(tabName = "input", text = "Input/Filtering", icon = icon("upload")),
      menuItem(tabName = "wg_plot", text = "WG Plots", icon = icon("chart-column")),
      menuItem(tabName = "chr_plot", text = "Chromosomal Plots", icon = icon("search")),
      menuItem(tabName = "table", text = "Tables", icon = icon("table")),
      menuItem(tabName = "help", text = "Help", icon = icon("question-circle")),
      menuItem(tabName = "about", text = "About", icon = icon("info-circle"))
    )
  ),
  
  dashboardBody(
      tabItems(
        tabItem(tabName = "input",
                tags$head(
                  tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css")
                ),
          fluidRow(
            column(width=8,
                   box(title="Step 0: File import/upload",status="primary",width = 12,solidHeader = T,collapsible = T,
                       radioButtons(inputId = "ref",label = h4("Reference Genome Build"),choices = list("hg19"="hg19","hg38"="hg38"),inline = T,selected = "hg19"),
                       h5("Upload or select proband read depth file",dashboardBadge("required", color = "primary")),
                       mod_rd_upload_UI("pr_rd"),
                       h5("Upload or select mother read depth file",dashboardBadge("optional", color = "secondary")),
                       mod_rd_upload_UI("m_rd"),
                       h5("Upload or select father read depth file",dashboardBadge("optional", color = "secondary")),
                       mod_rd_upload_UI("f_rd")
                       )
                   ),
            column(width=4,
                   box(title = "Help",icon=icon("question-circle"),status = "secondary",width = 12,solidHeader = T,collapsible = T,
                       bs4Dash::accordion(
                         id = "help_panel1",
                         bs4Dash::accordionItem(
                           title = "Sample read depth file (.bed or .bed.gz)",
                           collapsed = F,
                           tags$table(style = "border: 1px gray; padding: 1%; width: 100%;",
                             tags$tr(
                             tags$th("chr"),
                             tags$th("start"),
                             tags$th("end"),
                             tags$th("depth")
                           ),
                           tags$tr(
                             tags$td("chr1"),
                             tags$td("10000"),
                             tags$td("11000"),
                             tags$td("30.2")
                           ),
                           tags$tr(
                             tags$td("chr1"),
                             tags$td("11000"),
                             tags$td("12000"),
                             tags$td("32.2")
                           )
                           )
                         )),
                       bs4Dash::accordion(
                         id = "help_panel2",
                         bs4Dash::accordionItem(
                           title = "How to make .bed files",
                           collapsed = F,
                           "A output from",
                           tags$a(href="https://github.com/brentp/mosdepth","mosedepth"), 
                           "can be used as the input. An example of generating the read depth file for 1Kb window size would be:",
                           br(),
                           code("mosdepth -n --fast-mode --by 1000 sample.wgs $sample.wgs.bam")
                         )
                       )

                       )
                       )
                   ),
          fluidRow(
            column(width=8,
                   box(title="Advanced",status="primary",width = 12,solidHeader = T,collapsible = T,
                       h5("Upload or select proband SV vcf",dashboardBadge("optional", color = "secondary")),
                       mod_sv_upload_UI("pr_sv"),
                       h5("Upload or select mother SV vcf",dashboardBadge("optional", color = "secondary")),
                       mod_sv_upload_UI("m_sv"),
                       h5("Upload or select father SV vcf",dashboardBadge("optional", color = "secondary")),
                       mod_sv_upload_UI("f_sv"),
                       h5("Upload or select joint snp vcf (maximun file size 3Gb)",dashboardBadge("optional", color = "secondary")),
                       mod_snp_upload_UI("snp_file")
                       )
                   )
          )
          ),
        tabItem( tabName = "wg_plot",
                 fluidRow(
                   box(title = "Options",width = 12,solidHeader = T, status = "success",collapsible = T,
                       fluidRow(column(3, 
                                       actionButton("btn_wg_rd", "Show whole genome read depth"),
                                       actionButton("btn_wg_dnCNV", "Show potential dnCNVs")),
                       column(3, 
                              radioButtons("wg_norm_options",
                                           label = "Normalization options", 
                                           choiceNames = list("chromosomal median", "whole genome median"),
                                           choiceValues = list("chr_med", "wg_med")))
                   )),
                   box(title = "WG Plots",width = 12,solidHeader = T, status = "success",collapsible = T,
                       mod_plot_wg_UI("wg_pr_rd", height = 270),
                       mod_plot_wg_UI("wg_m_rd", height = 270),
                       mod_plot_wg_UI("wg_f_rd", height = 270)),
                   box(title = "WG Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                       mod_dnCNV_UI("wg_dnCNV"))
                 )
        ),
        tabItem(tabName = "chr_plot",
                shinyjs::useShinyjs(),
                fluidRow(
                  box(title = "Customizations",width = 12,solidHeader = T, status = "success",collapsible = T,
                    use_waiter(),
                    fluidRow(
                      use_waiter(),
                      column(1, 
                             p(HTML("<b>Chromosome</b>"),span(shiny::icon("info-circle"), id = "info_chr"),selectInput('chr', choice=c(paste0("chr", c(seq(1,22), "X"))),label = NULL, multiple = F),
                               tippy::tippy_this(elementId = "info_chr",tooltip = "Selected chromosome",placement = "right")
                             )),
                      column(2,
                             offset = 1, 
                             p(HTML("<b>Segment option</b>"),span(shiny::icon("info-circle"), id = "info_seg"),radioButtons('seg_option', 
                                                                                                                            label = NULL,
                                                                                                                            choiceNames = list("SLM","CBS"),
                                                                                                                            choiceValues = list("slm","cbs")),
                               tippy::tippy_this(elementId = "info_seg",tooltip = "Select options for segment. 
                                                             SLM is fast and tend to give more segment,can be used for high-quality data; 
                                                             CBS is slow and give less segment, can be used for noisy data ",placement = "right")
                             )),
                      column(3, 
                             radioButtons("norm_options",
                                          label = "Normalization options", 
                                          choiceNames = list("chromosomal median", "whole genome median"),
                                          choiceValues = list("chr_med", "wg_med"))),
                      column(2, 
                             p(HTML("<b>Segments to be included</b>"),span(shiny::icon("info-circle"), id = "info_include"),
                               checkboxGroupInput(inputId="include_seg",
                                                  label = NULL,
                                                  c("Proband"="Proband",
                                                    "Mom"="Mother",
                                                    "Dad"="Father"),selected = "Proband"),
                               tippy::tippy_this(elementId = "info_include",tooltip = "Choose to show segment from either or both parents",placement = "right")
                             ),
                             uiOutput("blt_dnSNV_ui"))
                      ),
                    fluidRow(
                      column(1, actionButton("btn_filter", "Filter")),
                      column(1, actionButton("btn_plot", "Plot")),
                      column(1, actionButton("btn_anno", "Annotate")),
                      column(2, actionButton("btn_dnCNV", "Show potential dnCNVs"))
                      )
                    )
                  ),
                  fluidRow(box(title = "Plots",closable = TRUE,width = 12,solidHeader = T, status = "success",collapsible = T,
                               sidebar = boxSidebar(
                                 width = 25,
                                 startOpen = TRUE,
                                 id = "mycardsidebar",
                                 HTML("<b>Basic Plot Options: </b>"),
                                 mod_checkbox_UI("RD-static"),
                                 mod_checkbox_UI("RD-dynamic"),
                                 mod_checkbox_UI("Baf-B_allele"),
                                 mod_checkbox_UI("Baf-A_allele", value = F),
                                 HTML("<b>Annotation Track Options:</b>"),
                                 mod_checkbox_UI("RefSeq"),
                                 mod_checkbox_UI("SegDup"),
                                 mod_checkbox_UI("OMIM"),
                                 mod_checkbox_UI("gnomAD"),
                                 mod_checkbox_UI("IDR", value = F),
                                 mod_checkbox_UI("RMSK", value = F),
                                 mod_checkbox_UI("pr_sv", value = F)
                               ),
                               fluidRow(
                                 column(4,shiny::textInput("goto_reg",label = NULL,placeholder = "gene symbol, chromosome location/range")),
                                 column(1,shiny::actionButton("btn_go","go")),
                                 column(2,verbatimTextOutput("cur_loc")),
                                 column(4,verbatimTextOutput("cur_range"))
                               ),
                               fluidRow(
                                 column(4,mod_col_pick_UI("highlight")),
                                 column(4,mod_UCSC_UI("UCSC"))
                               ),
                    mod_plot_switch_UI("RD-static", height = 200),
                    mod_plot_switch_UI("RD-dynamic", height = 200),
                    mod_plot_switch_UI("Baf-B_allele", height = 200),
                    mod_plot_switch_UI("Baf-A_allele", height = 200),
                    
                    mod_plot_switch_UI("pr_sv"),
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
                  tabPanel("pr_sv",anno_table_UI("pr_sv")),
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

