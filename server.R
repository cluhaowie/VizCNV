

server <- function(input, output,session) {
  # Reavtive Values --------------------------
  values <- reactiveValues()
  
  values$data <- data.frame(stringsAsFactors = F)
  values$work_data <- data.frame(stringsAsFactors = F)
  values$pr_rd <- data.frame(stringsAsFactors = F)
  values$m_rd <- data.frame(stringsAsFactors = F)
  values$f_rd <- data.frame(stringsAsFactors = F)
  values$snp_gvcf_file <- data.frame(stringsAsFactors = F)
  values$local_file_paths <- data.frame(file=c("CNV call file",
                                               "Proband read depth file",
                                               "Mom's read depth file",
                                               "Dad's read depth file",
                                               "Joint SNP vcf file"),datapath=rep("None",5),stringsAsFactors = F)
  values$selected_record <- data.frame(stringsAsFactors = F)
  
  plots <- reactiveValues()
  plots$pr_rd <- data.frame(stringsAsFactors = F)
  plots$m_rd <- data.frame(stringsAsFactors = F)
  plots$f_rd <- data.frame(stringsAsFactors = F)
  plots$pr_seg <- data.frame(stringsAsFactors = F)
  plots$m_seg <- data.frame(stringsAsFactors = F)
  plots$f_seg <- data.frame(stringsAsFactors = F)
  plots$snp_chr <- data.frame(stringsAsFactors = F)
  plots$xlabel <- character()
  plots$SNPcols <- vector(length = 3) ## placeholder for color in SNP plot
  
  plots$plot1 <- list()
  plots$plot2 <- list()
  toListen <- reactive({
    list(values$data,values$pr_rd)
  })
  
  volumes <- c(Home="~/Downloads/test","R installation" = R.home(),shinyFiles::getVolumes()())
  shinyFileChoose(input, "local_sv_file", roots = volumes, session = session)
  shinyFileChoose(input, "local_pr_rd_file", roots = volumes, session = session)
  shinyFileChoose(input, "local_m_rd_file", roots = volumes, session = session)
  shinyFileChoose(input, "local_f_rd_file", roots = volumes, session = session)
  shinyFileChoose(input, "local_pr_snv_file", roots = volumes, session = session)
  output$file_source_ui1 <- renderUI({
    if(input$file_source=="FALSE"){
      tagList(
        radioButtons(inputId = "ref",label = h3("Choose one reference genome"),choices = list("GRCh37"="GRCh37","GRCh38"="GRCh38"),inline = T,selected = "GRCh38"),
        fileInput("sv_vcf_file",label = "Structual variant vcf files",accept=c("*.vcf","*.vcf.gz"),multiple = F,buttonLabel = "Browse..."),
        fileInput("snp_gvcf_file",label = "SNV gVCF file",accept=c("*.vcf","*.vcf.gz"),multiple = F,buttonLabel = "Browse...")
      )
    }else{
      tagList(
        radioButtons(inputId = "ref",label = h3("Choose one reference genome"),choices = list("GRCh37"="GRCh37","GRCh38"="GRCh38"),inline = T,selected = "GRCh38"),
        h5(strong("Select local sv file (optional):")),
        shinyFilesButton(id = "local_sv_file", label = "Browse...", title = "Please select a file", multiple = F, viewtype = "detail"),
        h5(strong("Select local proband read depth file (required)")),
        shinyFilesButton(id = "local_pr_rd_file", label = "Browse...", title = "Please select a file", multiple = F, viewtype = "detail"),
        h5(strong("Select local mom's read depth file (optional)")),
        shinyFilesButton(id = "local_m_rd_file", label = "Browse...", title = "Please select a file", multiple = F, viewtype = "detail"),
        h5(strong("Select local dad's read depth file (optional)")),
        shinyFilesButton(id = "local_f_rd_file", label = "Browse...", title = "Please select a file", multiple = F, viewtype = "detail"),
        h5(strong("Select local proband snv file (optional)")),
        shinyFilesButton(id = "local_pr_snv_file", label = "Browse...", title = "Please select a file", multiple = F, viewtype = "detail")
      )
      
    }
  })
  output$file_source_ui2 <- renderUI({
    if(input$file_source=="FALSE"){
      tagList(
        fileInput("pr_rd_file",label = "Proband's read depth files (required)",accept=c("*.bed","*.bed.gz"),multiple = F,buttonLabel = "Browse..."),
        fileInput("m_rd_file",label = "Mom's read depth files (optional)",accept=c("*.bed","*.bed.gz"),multiple = F,buttonLabel = "Browse..."),
        fileInput("f_rd_file",label = "Dad's read depth files (optional)",accept=c("*.bed","*.bed.gz"),multiple = F,buttonLabel = "Browse...")
      )
    }else{
      verbatimTextOutput("filepaths")
    }
  })
  output$btn_filter_plot <- renderUI({
    tagList(
      column(1,actionButton("btn_filter","Filter")),
      column(1,actionButton("btn_plot","Plot"))
    )
  })
  output$filepaths <- renderPrint({
    ifelse(is.integer(input$local_sv_file),
           values$local_file_paths$datapath[1] <- "None",
           values$local_file_paths$datapath[1] <- parseFilePaths(volumes, input$local_sv_file)$name)
    ifelse(is.integer(input$local_pr_rd_file),
           values$local_file_paths$datapath[2] <- "None",
           values$local_file_paths$datapath[2] <- parseFilePaths(volumes, input$local_pr_rd_file)$name)
    ifelse(is.integer(input$local_m_rd_file),
           values$local_file_paths$datapath[3] <- "None",
           values$local_file_paths$datapath[3] <- parseFilePaths(volumes, input$local_m_rd_file)$name)
    ifelse(is.integer(input$local_f_rd_file),
           values$local_file_paths$datapath[4] <- "None",
           values$local_file_paths$datapath[4] <- parseFilePaths(volumes, input$local_f_rd_file)$name)
    ifelse(is.integer(input$local_pr_snv_file),
           values$local_file_paths$datapath[5] <- "None",
           values$local_file_paths$datapath[5] <- parseFilePaths(volumes, input$local_pr_snv_file)$name)
    values$local_file_paths
  })
  
  # observe file uploaded and save in SQLdatabase---------
  # local option
  observeEvent(input$local_sv_file,{ 
    if(is.integer(input$local_sv_file)){cat("no file\n")}else{
      local_sv_file <- parseFilePaths(volumes, input$local_sv_file)
      values$data <- bedr::read.vcf(local_sv_file$datapath,split.info = T)$vcf
      showModal(modalDialog(title = "File upload",
                            "The sv file has been uploaded"))
    }
    
  },ignoreInit = T)
  observeEvent(input$local_pr_rd_file,{ 
    if(is.integer(input$local_pr_rd_file)){cat("no file\n")}else{
      local_pr_rd_file <- parseFilePaths(volumes, input$local_pr_rd_file)
      values$pr_rd <- data.table::fread(local_pr_rd_file$datapath,header = F)%>%mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))
      showModal(modalDialog(title = "File upload",
                            "The proband read depth file has been uploaded"))
    }
  },ignoreInit = T)
  observeEvent(input$local_m_rd_file,{ 
    if(is.integer(input$local_m_rd_file)){cat("no file\n")}else{
      local_m_rd_file <- parseFilePaths(volumes, input$local_m_rd_file)
      values$m_rd <- data.table::fread(local_m_rd_file$datapath,header = F)%>%mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))
      showModal(modalDialog(title = "File upload",
                            "The mom's read depth file has been uploaded"))
    }
  },ignoreInit = T)
  observeEvent(input$local_f_rd_file,{ 
    if(is.integer(input$local_f_rd_file)){cat("no file\n")}else{
      local_f_rd_file <- parseFilePaths(volumes, input$local_f_rd_file)
      values$f_rd <- data.table::fread(local_f_rd_file$datapath,header = F)%>%mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))
      showModal(modalDialog(title = "File upload",
                            "The dad's read depth file has been uploaded"))
    }
  },ignoreInit = T)
  observeEvent(input$local_pr_snv_file,{ 
    if(is.integer(input$local_pr_snv_file)){cat("no file\n")}else{
      values$snp_gvcf_file <- parseFilePaths(volumes, input$local_pr_snv_file)
      index.file <- paste0(values$snp_gvcf_file$datapath,".tbi")
      if(!file.exists(index.file)){
        Rsamtools::indexTabix(values$snp_gvcf_file$datapath,format = "vcf")
      }
      showModal(modalDialog(
        title = "File upload",
        "The joint called SNP file has been indexed"
      ))
    }
  },ignoreInit = T)
  # cloud option
  observe({
    sv_vcf_file=input$sv_vcf_file
    if(is.null(sv_vcf_file)){return(NULL)}
    req(sv_vcf_file)
    values$data <- bedr::read.vcf(sv_vcf_file$datapath,split.info = T)$vcf
    #saveData(values$data,"sv_vcf_file")
    showModal(modalDialog(
      title = "File upload",
      "The sv file has been uploaded"
    ))
  })
  observe({
    pr_rd_file=input$pr_rd_file
    if(is.null(pr_rd_file)){return(NULL)}
    req(pr_rd_file)
    values$pr_rd <- data.table::fread(pr_rd_file$datapath,header = F)%>%mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))
    #saveData(values$pr_rd,"pr_rd_file")
    showModal(modalDialog(
      title = "File upload",
      "The proband bed file has been uploaded"
    ))
  }) 
  observe({
    m_rd_file=input$m_rd_file
    if(is.null(m_rd_file)){return(NULL)}
    req(m_rd_file)
    values$m_rd <- data.table::fread(m_rd_file$datapath,header = F)%>%mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))
    #saveData(values$m_rd,"m_rd_file")
    showModal(modalDialog(
      title = "File upload",
      "The mom's bed file has been uploaded"
    ))
  }) 
  observe({
    f_rd_file=input$f_rd_file
    if(is.null(f_rd_file)){return(NULL)}
    req(f_rd_file)
    values$f_rd <- data.table::fread(f_rd_file$datapath,header = F)%>%mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))
    #saveData(values$f_rd,"f_rd_file")
    showModal(modalDialog(
      title = "File upload",
      "The dad's bed file has been uploaded"
    ))
  }) 
  observeEvent(input$snp_gvcf_file,{
    values$snp_gvcf_file <- input$snp_gvcf_file
    if(is.null(values$snp_gvcf_file)){return(NULL)}
    req(values$snp_gvcf_file)
    Rsamtools::indexTabix(values$snp_gvcf_file$datapath,format = "vcf")
    showModal(modalDialog(
      title = "File upload",
      "The joint called SNP file has been uploaded and indexed"
    ))
  }) 
  
  observeEvent(values$data,{
    data <- values$data
    ALTs <- unique(data$ALT)
    updateSelectizeInput(session, "type", choices = ALTs, server = TRUE)
  })
  
  observeEvent(toListen(),{
    if(nrow(values$data)!=0){
      data <- values$data
      chroms <- unique(data$CHROM)
      updateSelectizeInput(session, "chr", choices = chroms, server = TRUE)
    }else{
      data <- values$pr_rd
      chroms <- unique(data$V1)
      updateSelectizeInput(session, "chr", choices = chroms, server = TRUE)
    }
    
  })
  
  observeEvent(input$chr,{
    chr <- input$chr
    if(nrow(values$data) == 0){
      return(NULL)
    }else{
      values$work_data <- values$data%>%filter(CHROM==chr)
      return(values$work_data)
    }
    
  })
  
  # button to filter range---------
  
  
  output$filter_sv_table <- DT::renderDataTable({ 
    values$work_data
  },
  extensions=c("Responsive","Buttons"),
  server = T,editable = TRUE,filter = list(position = 'top', clear = T),options = list(dom = 'Bfrtip',buttons = c('txt','csv', 'excel')))
  w <- waiter::Waiter$new(html = spin_3(), 
                  color = transparent(.5))
  output$Select_table <- DT::renderDataTable({
    values$selected_record
  })
  ## keep the selected record when click the btl_select
  observeEvent(input$btl_select,{
    rows_selected <- input$filter_sv_table_rows_selected
    if(length(rows_selected)){
      values$selected_record <- rbind(values$selected_record,values$work_data[rows_selected,])%>%
        distinct(across(everything()))
    }
  })
  
  observeEvent(input$btn_filter,{
    seg_option <- input$seg_option
    chr <- input$chr
    if(nrow(values$pr_rd)==0){return(NULL)
    }else{
      w$show()
      plots$pr_rd <- values$pr_rd%>%
        filter(V1==chr)%>%
        mutate(ratio=V4/median(V4+0.00001))
      plots$pr_seg <- SegNormRD(plots$pr_rd,id="index",seg.method = seg_option)
      w$hide()
    }
    if(nrow(values$m_rd)==0){return(NULL)
    }else{
      plots$m_rd <- values$m_rd%>%
        filter(V1==chr)%>%
        mutate(ratio=V4/median(V4+0.00001))
      plots$m_seg <- SegNormRD(plots$m_rd,id="m",seg.method = seg_option)
    }
    if(nrow(values$f_rd)==0){return(NULL)
    }else{
      plots$f_rd <- values$f_rd%>%
        filter(V1==chr)%>%
        mutate(ratio=V4/median(V4+0.00001))
      plots$f_seg <- SegNormRD(plots$f_rd,id="f",seg.method = seg_option)
    }
  })
  observeEvent(input$btn_filter,{
    snp_gvcf_file=values$snp_gvcf_file
    ref_genome <- input$ref
    chr <- input$chr
    if(is.null(snp_gvcf_file$datapath)){
      return(NULL)
    }else{
      w$show()
      loc.start <- 0
      loc.end <- hg38.info%>%filter(chrom==chr)%>%dplyr::select(seqlengths)%>%unlist
      range.gr <- GenomicRanges::GRanges(chr,ranges = IRanges(loc.start,loc.end))
      plots$snp_chr <- ReadGVCF(snp_gvcf_file$datapath,ref_genome=input$ref,param = range.gr)%>%as.data.frame()
      InhFrom <- unique(plots$snp_chr$InhFrom)
      if(length(InhFrom)==3){
        names(plots$SNPcols) <- InhFrom
        plots$SNPcols[names(plots$SNPcols)!="Notphased"] <- SNPCOLOR2
        plots$SNPcols["Notphased"] <- c("#999999")
      }
      w$hide()
    }
  })
  
  
  ext1 <- eventReactive(input$btn_plot,{
    # from input
    if(nrow(plots$pr_rd) == 0){
      return(NULL)
    }
    include_seg <- input$include_seg
    df <- rbindlist(list(plots$pr_seg,plots$m_seg,plots$f_seg))%>%filter(ID%in%include_seg)
    plots$xlabel=unique(df$chrom)[1]
    ggplot(plots$pr_rd, aes(V2, log2(ratio+0.00001))) +
      geom_point(shape=".")+
      geom_segment(data = df,aes(x=loc.start,y=seg.mean,xend=loc.end,yend=seg.mean,color=ID),size=1)+
      #geom_segment(data=work_data,aes(x=POS,xend=END,y=value,yend=value,color=ALT),size=2,alpha=0.5)+
      ylim(-4,4)+xlab(plots$xlabel)+
      scale_rd+style_rd+scale_x_continuous(labels = scales::label_number())
  },ignoreInit = T)
  ext2 <- eventReactive(input$btn_plot,{
    if(nrow(plots$snp_chr) == 0){
      return(NULL)
    }
    df <- plots$snp_chr%>%filter(likelyDN%in%c(input$include_dnSNV,"FALSE"))
    cols <- plots$SNPcols
    xlabel=unique(df$chrom)[1]
    df %>% ggplot(aes(x=start,y=pr_ALT_Freq,col=InhFrom))+
      geom_point(shape=".")+
      geom_point(data = subset(df, likelyDN %in%c("TRUE")),size = 2,shape=8,color="red")+
      scale_fill_manual("LikelyDN",limits=c("dnSNV"),values = "red")+
      xlab(xlabel)+
      scale_snp+
      style_snp+
      scale_colour_manual(values = cols)+
      guides(color = guide_legend(override.aes = list(size = 4)))+scale_x_continuous(labels = scales::label_number())
  })
  
  # interactive plot regions-------
  ranges <- reactiveValues(x = NULL, y = NULL)
  output$plot1 <- renderPlot({
    if(nrow(plots$pr_rd) == 0){
      return(NULL)
    }
    plots$plot1 <- ext1()+coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
      ggtitle(paste0(plots$xlabel,":",paste0(round(as.numeric(ranges$x)),collapse = "-")))
    plots$plot1
  })
  output$plot2 <- renderPlot({
    if(nrow(plots$snp_chr) == 0){
      return(NULL)
    }
    plots$plot2 <- ext2()+coord_cartesian(xlim = ranges$x, expand = FALSE)
    plots$plot2
  })
  #When a double-click happens, check if there's a brush on the plot.
  #If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    chr <- input$chr
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      if(nrow(values$data)!=0){
        values$work_data <- values$data%>%filter(CHROM==chr)%>%filter(POS>=brush$xmin,POS < brush$xmax)
      }
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
      if(nrow(values$data)!=0){
        values$work_data <- values$data%>%filter(CHROM==chr)
      }
    }
  })
  observeEvent(input$plot2_dblclick, {
    brush <- input$plot2_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  ## buttons 
  output$ui_dlbtn_tbl <- renderUI({
    if(nrow(values$data) > 0){
      actionButton("btl_select", "Select Record")
    }
  })
  output$ui_dlbtn_plt <- renderUI({
    if(length(plots$plot1) > 0){
      downloadButton("dl_plt", "Download")
    }
  })
  output$ui_clbtn_plt <- renderUI({
    if(length(plots$plot1) > 0){
      shiny::actionButton("cl_btn","Clear plot",icon("trash"))
    }
  })
  observeEvent(input$cl_btn,{
    plots$snp_chr <- data.frame(stringsAsFactors = F)
    plots$pr_rd <- data.frame(stringsAsFactors = F)
  })
  ## Download handler
  output$dl_plt <- downloadHandler(
    filename = function(){
      paste("index",input$chr,".pdf")
    },
    content = function(file){
      p <- cowplot::plot_grid(plots$plot1,plots$plot2,ncol = 1)
      ggplot2::ggsave(filename =file, plot = p,device = "pdf",width =12 ,height = 8,units = "in")
    }
  )
  
}