source("./mod/mod_plot_output.R")
source("./mod/mod_dnCNV.R")
source("./mod/mod_findCNV.R")
source("./mod/mod_hmzcnv.R")
source("./mod/mod_upload.R")
source("./mod/mod_UCSC.R")

ui <- dashboardPage(
  
  dashboardHeader(title = "VizCNV",
                  tags$span(style = "width: 90%"),
                  tagList(
                    tags$a(href="https://github.com/cluhaowie/VizCNV",icon("github"))
                    )
                  ),
  
  dashboardSidebar(
    sidebarMenu(
      id = "tabs",
      menuItem(tabName = "home", text = "Home", icon = icon("home")),
      menuItem(tabName = "input", text = "Input", icon = icon("upload")),
      menuItem(tabName = "wg_plot", text = "Genome-wide view", icon = icon("chart-column")),
      menuItem(tabName = "chr_plot", text = "Chromosomal view", icon = icon("search")),
      menuItem(tabName = "table", text = "Tables", icon = icon("table")),
      menuItem(tabName = "help", text = "Help", icon = icon("question-circle"))
      
    )
  ),
  
  dashboardBody(
      tabItems(
        tabItem(tabName = "home", 
                h1("VizCNV"),
                "Welcome to ",strong("VizCNV!"),
                p("VizCNV serves as a dynamic interface for analyzing and representing copy number variations (CNVs) derived from short-read whole genome sequencing (WGS) data, specifically in the realm of rare disease investigation.
                  Constructed using the interactive web application framework R Studio Shiny, VizCNV simplifies and accelerates the process of detecting intricate alterations in the genome structure. By offering a comprehensive, user-friendly platform for both computational analysis and visual representation of genomic data, VizCNV aids in enhancing our understanding of rare genetic disorders, thereby contributing to the improvement of diagnostic strategies.")),
        tabItem(tabName = "input",
                tags$head(
                  tags$link(rel = "stylesheet", href = "https://cdnjs.cloudflare.com/ajax/libs/font-awesome/4.7.0/css/font-awesome.min.css")
                ),
          fluidRow(
            column(width=8,
                   box(title="File(s) upload",status="primary",width = 12,solidHeader = T,collapsible = T,
                       radioButtons(inputId = "ref",label = "Reference Genome Build",choices = list("hg19"="hg19","hg38"="hg38","T2T"="T2T"),inline = T,selected = "hg19"),
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
                       h5("Upload or select joint snp vcf (maximum file size 3Gb)",dashboardBadge("optional", color = "secondary")),
                       mod_snp_upload_UI("snp_file")
                       )
                   )
          )
          ),
        tabItem( tabName = "wg_plot",
                 fluidRow(
                   box(title = "Options",width = 12,solidHeader = T, status = "success",collapsible = T,
                       fluidRow(column(3,
                                       radioButtons("wg_norm_options",
                                                    label = "Normalization options:", 
                                                    choiceNames = list("chromosomal median", "whole genome median"),
                                                    choiceValues = list("chr_med", "wg_med"),
                                                    selected = "wg_med"))),
                       fluidRow(column(3, actionButton("btn_wg_rd", "Plot"))
                                # column(3, actionButton("btn_wg_dnCNV", "Filter dnCNVs"))
                                )
                   ),
                   box(title = "WG Plots",width = 12,solidHeader = T, status = "success",collapsible = T,
                       fluidRow(
                       verbatimTextOutput("pr_name"),
                       verbatimTextOutput("m_name"),
                       verbatimTextOutput("f_name")),
                       mod_plot_wg_UI("wg_pr_rd", height = 270),
                       mod_plot_wg_UI("wg_m_rd", height = 270),
                       mod_plot_wg_UI("wg_f_rd", height = 270))
                   # box(title = "WG Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                   #     mod_dnCNV_UI("wg_dnCNV"))
                 )
        ),
        tabItem(tabName = "chr_plot",
                shinyjs::useShinyjs(),
                fluidRow(
                  box(title = "Customizations",width = 12,solidHeader = T, status = "success",collapsible = T,
                    use_waiter(),
                    fluidRow(
                      use_waiter(),
                      column(2, 
                             p(HTML("<b>Chromosome</b>"),span(shiny::icon("info-circle"), id = "info_chr"),selectInput('chr', choice=c(paste0("chr", c(seq(1,22), "X","Y"))),label = NULL, multiple = F, selected = "chr1"),
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
                      column(2, 
                             radioButtons("norm_options",
                                          label = "Normalization options", 
                                          choiceNames = list("chromosomal median", "whole genome median"),
                                          choiceValues = list("chr_med", "wg_med")),
                             p(HTML("<b>Is this a trio?</b>")),
                             checkboxInput("is_trio", "Yes", value = TRUE)),
                      column(2, 
                             p(HTML("<b>Segments to be included</b>"),span(shiny::icon("info-circle"), id = "info_include"),
                               checkboxGroupInput(inputId="include_seg",
                                                  label = NULL,
                                                  c("Proband"="Proband",
                                                    "Mother"="Mother",
                                                    "Father"="Father"),selected = c("Proband", "Mother", "Father")),
                               tippy::tippy_this(elementId = "info_include",tooltip = "Choose to show segment from either or both parents",placement = "right")
                             ),
                             uiOutput("blt_dnSNV_ui")),
                      column(2, p(HTML("<b>Other options</b>")),
                             checkboxInput(inputId = "mask_option",
                                           label = "remove gaps in ref genome",
                                           value = T),
                             checkboxInput(inputId = "baf_seg",
                                           label = "show Baf plot segments", 
                                           value = T),
                             uiOutput("blt_baf_seg")
                             
                      )
                      ),
                    fluidRow(
                      p(HTML("<b>Step1:</b>")),
                      column(2, actionButton("btn_filter", "Filter")),
                      p(HTML("<b>Step2:</b>")),
                      column(2, actionButton("btn_plot", "Plot")),
                      p(HTML("<b>Step3:</b>")),
                      column(2, actionButton("btn_anno", "Annotate")),
                      # column(2, actionButton("btn_dnCNV", "Filter dnCNVs")),
                      # column(2, actionButton("btn_hmzCNV", "Filter hmzCNV")),
                      p(HTML("<b>Optional:</b>")),
                      column(2, actionButton("btn_findCNV", "Find CNVs")),
                      
                      )
                    )
                  ),
                  fluidRow(box(title = "Plots",closable = F, width = 12,solidHeader = T, status = "success",collapsible = T, 
                               sidebar = boxSidebar(
                                 width = 25,
                                 startOpen = F,
                                 id = "mycardsidebar",
                                 HTML("<b>Basic Plot Options: </b>"),
                                 mod_checkbox_UI("RD-static"),
                                 mod_checkbox_UI("RD-dynamic"),
                                 mod_checkbox_UI("Baf-B_allele"),
                                 mod_checkbox_UI("Baf-A_allele", value = F),
                                 HTML("<b>Annotation Track Options:</b>"),
                                 mod_checkbox_UI("pr_sv", value = F),
                                 mod_checkbox_UI("m_sv", value = F),
                                 mod_checkbox_UI("f_sv", value = F),
                                 uiOutput("ui_chkbox_RefSeq"),
                                 uiOutput("ui_chkbox_SegDup"),
                                 uiOutput("ui_chkbox_OMIM"),
                                 uiOutput("ui_chkbox_gnomAD"),
                                 uiOutput("ui_chkbox_IDR"),
                                 uiOutput("ui_chkbox_RMSK")
                               ),
                               fluidRow(
                                 column(4,shiny::textInput("goto_reg",label = "Search: gene symbol, location/range (eg. 200000/200000-300000)")),
                                 column(1,shiny::actionButton("btn_go","go")),
                                 column(2,verbatimTextOutput("cur_loc")),
                                 column(4,verbatimTextOutput("cur_range"))
                               ),
                               fluidRow(
                                 column(6,mod_col_pick_UI("highlight")),
                                 column(6,mod_UCSC_UI("UCSC"))
                               ),
                               fluidRow(
                                 verbatimTextOutput("pr_name2"),
                                 verbatimTextOutput("m_name2"),
                                 verbatimTextOutput("f_name2"),
                               ),
                    mod_plot_switch_UI("RD-static", height = 200),
                    mod_plot_switch_UI("RD-dynamic", height = 200),
                    mod_plot_switch_UI("Baf-B_allele", height = 200),
                    mod_plot_switch_UI("Baf-A_allele", height = 200),
                    
                    mod_plot_switch_UI("pr_sv"),
                    mod_plot_switch_UI("m_sv"),
                    mod_plot_switch_UI("f_sv"),
                    mod_plot_switch_UI("RefSeq", height = 120),
                    mod_plot_switch_UI("Segdup"),
                    mod_plot_switch_UI("OMIM"),
                    mod_plot_switch_UI("gnomAD"),
                    mod_plot_switch_UI("IDR"),
                    mod_plot_switch_UI("RMSK"),
                    fluidRow(column(1,uiOutput("ui_dlbtn_plt")),
                             column(1,uiOutput("ui_clbtn_plt")),
                             column(1,uiOutput("ui_dlbtn_dnsnv")))
                  )
                 )
               ),
        tabItem(tabName = "table",
                fluidRow(box(title = "find CNV Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                             mod_findCNV_UI("findCNV")),
                         # box(title = "dnCNV Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                         #     mod_dnCNV_UI("dnCNV")),
                         # box(title = "Homozygous CNVs Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                         #     mod_hmzcnv_UI("hmzCNV")),
                         box(title = "GATK Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                             anno_table_UI("gatk_table")),
                         box(title = "SV Table",width = 12,solidHeader = T, status = "success",collapsible = T,
                             tabsetPanel(
                               tabPanel("Proband", anno_table_UI("pr_sv")),
                               tabPanel("Mother", anno_table_UI("m_sv")),
                               tabPanel("Father", anno_table_UI("f_sv"))
                              )
                             ),
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
        tabItem(tabName = "help")
    )
  )

)


