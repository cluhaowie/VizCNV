

source("./mod/mod_plot_output.R")
source("./mod/mod_dnCNV.R")

server <- function(input, output,session) {
  # Reavtive Values --------------------------
  values <- reactiveValues()
  
  values$pr_sv <- data.frame(stringsAsFactors = F)
  values$pr_rd <- data.frame(stringsAsFactors = F)
  values$m_rd <- data.frame(stringsAsFactors = F)
  values$f_rd <- data.frame(stringsAsFactors = F)
  values$snp_gvcf_file <- data.frame(stringsAsFactors = F)
  values$snp_gvcf_file_ref <- vector()
  values$local_file_paths <- data.frame(file=c("CNV call file",
                                               "Proband read depth file",
                                               "Mom's read depth file",
                                               "Dad's read depth file",
                                               "Joint SNP vcf file"),datapath=rep("None",5),stringsAsFactors = F)
  values$selected_record <- data.frame(stringsAsFactors = F)
  values$ref_info <- data.frame(stringsAsFactors = F)
  values$anno_rect <- data.frame(stringsAsFactors = F)
  
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
  plots$genelabel <- data.frame(stringsAsFactors = F)
  
  plots$plot1 <- list()
  plots$plot2 <- list()
  plots$plot3 <- list()
  plots$plot3_dl <- list()
  
  toListen <- reactive({
    list(values$pr_sv,values$pr_rd)
  })
  
  volumes <- c(Home="~/Downloads/","R installation" = R.home(),shinyFiles::getVolumes()())
  shinyFileChoose(input, "local_sv_file", roots = volumes, session = session)
  shinyFileChoose(input, "local_pr_rd_file", roots = volumes, session = session)
  shinyFileChoose(input, "local_m_rd_file", roots = volumes, session = session)
  shinyFileChoose(input, "local_f_rd_file", roots = volumes, session = session)
  shinyFileChoose(input, "local_pr_snv_file", roots = volumes, session = session)
  output$file_source_ui1 <- renderUI({
    if(input$file_source=="TRUE"){
      tagList(
         fileInput("sv_vcf_file",label = "Structual variant vcf files",accept=c("*.vcf","*.vcf.gz"),multiple = F,buttonLabel = "Browse..."),
        fileInput("snp_gvcf_file",label = "SNV gVCF file",accept=c("*.vcf","*.vcf.gz"),multiple = F,buttonLabel = "Browse...")
      )
    }else{
      tagList(
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
    if(input$file_source=="TRUE"){
      tagList(
        fileInput("pr_rd_file",label = "Proband's read depth files (required)",accept=c("*.bed","*.bed.gz"),multiple = F,buttonLabel = "Browse..."),
        fileInput("m_rd_file",label = "Mom's read depth files (optional)",accept=c("*.bed","*.bed.gz"),multiple = F,buttonLabel = "Browse..."),
        fileInput("f_rd_file",label = "Dad's read depth files (optional)",accept=c("*.bed","*.bed.gz"),multiple = F,buttonLabel = "Browse...")
      )
    }else{
      verbatimTextOutput("filepaths")
    }
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
  output$blt_dnSNV_ui <- shiny::renderUI({
    if(is.null(input$snp_gvcf_file)&is.integer(input$local_pr_snv_file)){
      return(NULL)
    }
    else{
      shiny::tagList(
        p(HTML("<b>Show de novo SNV ?</b>"),
          span(shiny::icon("info-circle"), id = "info_nor"),
          checkboxGroupInput(inputId="include_dnSNV",label = NULL,c("Show"="TRUE"))))
    }
  })

  
  # observe file uploaded and save in SQLdatabase---------
  # local option
  observeEvent(input$local_sv_file,{ 
    if(is.integer(input$local_sv_file)){cat("no file\n")}else{
      local_sv_file <- parseFilePaths(volumes, input$local_sv_file)
      values$pr_sv <- bedr::read.vcf(local_sv_file$datapath,split.info = T)$vcf
      showModal(modalDialog(title = "File upload",
                            "The sv file has been uploaded"))
    }
    
  },ignoreInit = T)
  observeEvent(input$local_pr_rd_file,{ 
    if(is.integer(input$local_pr_rd_file)){cat("no file\n")}else{
      local_pr_rd_file <- parseFilePaths(volumes, input$local_pr_rd_file)
      values$pr_rd <- data.table::fread(local_pr_rd_file$datapath,header = F)%>%
        mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))%>%regioneR::toGRanges()
      values$pr_rd <- values$pr_rd[!values$pr_rd%over%blacklist,]%>%
        as.data.table()%>%
        dplyr::select(-c("strand","width"))
      setnames(values$pr_rd,c("seqnames","start","end"),c("V1","V2","V3"))
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
      #read the header of p.VCF file to checking the chr id
      values$snp_gvcf_file_ref <- VariantAnnotation::scanVcfHeader(values$snp_gvcf_file$datapath)@reference
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
    values$pr_sv <- bedr::read.vcf(sv_vcf_file$datapath,split.info = T)$vcf
    #saveData(values$pr_sv,"sv_vcf_file")
    showModal(modalDialog(
      title = "File upload",
      "The sv file has been uploaded"
    ))
  })
  observe({
    pr_rd_file=input$pr_rd_file
    if(is.null(pr_rd_file)){return(NULL)}
    req(pr_rd_file)
    values$pr_rd <- data.table::fread(pr_rd_file$datapath,header = F)%>%
      mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))
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
    #read the header of p.VCF file to checking the chr id
    values$snp_gvcf_file_ref <- VariantAnnotation::scanVcfHeader(values$snp_gvcf_file$datapath)@reference
    showModal(modalDialog(
      title = "File upload",
      "The joint called SNP file has been uploaded and indexed"
    ))
  }) 
  
  observeEvent(values$pr_sv,{
    data <- values$pr_sv
    ALTs <- unique(data$ALT)
    updateSelectizeInput(session, "type", choices = ALTs, server = TRUE)
  })
  
  observeEvent(toListen(),{
    if(nrow(values$pr_sv)!=0){
      data <- values$pr_sv
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
    if(nrow(values$pr_sv) == 0){
      return(NULL)
    }else{
      values$pr_sv <- values$pr_sv%>%filter(CHROM==chr)
      return(values$pr_sv)
    }
  })
  observeEvent(input$ref,{
    if(input$ref=="GRCh38"){
      blacklist <- data.table::fread("GRCh38_unified_blacklist.bed.gz")%>%
        regioneR::toGRanges()
      values$ref_info <- data.table::fread("hg38.info.txt")
    }else if(input$ref=="GRCh37"){
      blacklist <- data.table::fread("ENCFF001TDO.bed.gz")%>%
        regioneR::toGRanges()
      values$ref_info <- data.table::fread("hg19.info.txt")
    }
  })
  
  # button to filter range---------
  
  
  output$filter_sv_table <- DT::renderDataTable({ 
    values$pr_sv
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
      values$selected_record <- rbind(values$selected_record,values$pr_sv[rows_selected,])%>%
        distinct(across(everything()))
    }
  })
  observeEvent(input$btl_select2,{
    values$selected_record <- data.frame(stringsAsFactors = F)
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
      plots$pr_seg <- SegNormRD(plots$pr_rd,id="Proband",seg.method = seg_option)
      # print(class(plots$pr_seg))
      # print(head(plots$pr_seg))
      w$hide()
    }
    if(nrow(values$m_rd)==0){return(NULL)
    }else{
      w$show()
      plots$m_rd <- values$m_rd%>%
        filter(V1==chr)%>%
        mutate(ratio=V4/median(V4+0.00001))
      plots$m_seg <- SegNormRD(plots$m_rd,id="Mother",seg.method = seg_option)
      w$hide()
    }
    if(nrow(values$f_rd)==0){return(NULL)
    }else{
      w$show()
      plots$f_rd <- values$f_rd%>%
        filter(V1==chr)%>%
        mutate(ratio=V4/median(V4+0.00001))
      plots$f_seg <- SegNormRD(plots$f_rd,id="Father",seg.method = seg_option)
      w$hide()
    }
  })
  observeEvent(input$btn_filter,{
    snp_gvcf_file=values$snp_gvcf_file
    ref_genome <- input$ref
    chr <- input$chr
    if(!chr%in%values$snp_gvcf_file_ref){
      if(chr%in%chrom_id){
        chr <- names(chrom_id)[which(chrom_id==chr)]
      }else{
        chr <- chrom_id[which(names(chrom_id)==chr)]
      }
    }
    if(is.null(snp_gvcf_file$datapath)){
      return(NULL)
    }else{
      w$show()
      loc.start <- 0
      loc.end <- values$ref_info%>%
        filter(chrom%in%c(chr,paste0("chr",chr)))%>%
        dplyr::select(seqlengths)%>%unlist
      range.gr <- GenomicRanges::GRanges(chr,ranges = IRanges(loc.start,loc.end))
      range.gr <- GenomicRanges::setdiff(range.gr, blacklist)
      plots$snp_chr <- ReadGVCF(snp_gvcf_file$datapath,ref_genome=input$ref,param = range.gr)%>%
        as.data.frame()
      InhFrom <- unique(plots$snp_chr$B_InhFrom)
      if(length(InhFrom)==3){
        names(plots$SNPcols) <- InhFrom
        plots$SNPcols[names(plots$SNPcols)!="Notphased"] <- SNPCOLOR2
        plots$SNPcols["Notphased"] <- c("#999999")
      }
      w$hide()
    }
  })
  


  


  observeEvent(input$btl_add,{
    brush <- input$plot1_brush
    if(!is.null(brush)){
      xmin <- round(brush$xmin,0)
      xmax <- round(brush$xmax,0)
      fill <- input$add_col
      values$anno_rect <- rbindlist(list(values$anno_rect,data.frame(xmin,xmax,fill)))
    }
  })
  observeEvent(input$btl_reset,{
    values$anno_rect <- data.frame(stringsAsFactors = F)
  })
  
  ## buttons 
  output$ui_dlbtn_tbl <- renderUI({
    if(nrow(values$pr_sv) > 0){
      tagList(shiny::actionButton("btl_select", "Select",icon("check")))
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
  output$ui_dlbtn_dnsnv <- renderUI({
    if(length(plots$plot2) > 0){
      shiny::downloadButton("dl_btn_dnsnv","Download dnSNV")
    }
  })
  output$ui_btn_goto <- renderUI({
    if(length(plots$plot1)>0){
      tagList(
        fluidRow(
          column(1,shiny::actionButton("btl_goto","goto")),
          column(8,shiny::textInput("goto_reg",label = NULL,placeholder = "e.g, 10000-20000 or GENE"))))
    }
  })
  output$ui_btn_add <- renderUI({
    if(length(plots$plot1)>0){
      tagList(
        fluidRow(column(1,shiny::actionButton("btl_add","Add")),
                 column(4,colourpicker::colourInput("add_col",NULL,"yellow",palette = "limited")),
                 column(1,shiny::actionButton("btl_reset","Reset"))))
    }
  })
  observeEvent(input$cl_btn,{
    plots$snp_chr <- data.frame(stringsAsFactors = F)
    plots$pr_rd <- data.frame(stringsAsFactors = F)
    input$filter_sv_table_rows_selected <- NULL
  })
  
  ## Download handler
  output$dl_plt <- downloadHandler(
    filename = function(){
      paste0(input$chr,".pdf")
    },
    content = function(file){
      
      mylist <- list(plots$plot1,plots$plot3_dl,plots$plot2)
      mylist <- mylist[lengths(mylist)!= 0]
      n <- length(mylist)
      p <- cowplot::plot_grid(plotlist=mylist,ncol = 1,align = 'v',axis = 'lr')
      ggplot2::ggsave(filename =file, plot = p,device = "pdf",width =12 ,height = n*4,units = "in")
    }
  )
  output$dl_btn_dnsnv <- downloadHandler(
    filename = function(){paste("dnSNV_",input$chr,".csv")},
    content = function(file){
      df <- plots$snp_chr%>%filter(likelyDN%in%c(input$include_dnSNV))
      write.csv(df,file,row.names = F)
    }
  )
  
  ### "Global" reactive values
  ranges <- reactiveValues(x = NULL, y = NULL)
  dnCNV_table <- reactiveValues(t = data.frame(start = c(0), end = c(0), stringsAsFactors = F))
  
  ## WG Plot section
  observeEvent(input$btn_wg_rd, {
    req(!is.null(values$pr_rd))
    
    w$show()
    wg_pr_rd <- values$pr_rd
    names(wg_pr_rd) <- c("chr", "start", "end", "coverage")
    wg_pr_rd <- wg_pr_rd %>% 
      group_by(chr) %>% 
      mutate(ratio=(coverage/median(coverage+0.00001)))
    wg_pr_seg <- getAllSeg(wg_pr_rd)
    
    wg_pr_seg <- wg_pr_seg %>% 
      mutate(seg.mean=ifelse(seg.mean < -2.5,-2,seg.mean))
    
    temp <- wg_pr_seg %>% 
      group_by(chr) %>% 
      summarise(max_end = max(loc.end)) %>% 
      mutate(across("chr", str_replace, "chr", "")) %>% 
      arrange(as.numeric(chr)) %>% 
      mutate(loc_add = lag(cumsum(as.numeric(max_end)), default = 0)) %>% 
      mutate(chr = paste0("chr", chr))
    
    wg_pr_seg <- wg_pr_seg %>% 
      inner_join(temp, by = "chr") %>% 
      mutate(end_cum = loc_add + loc.end) 
    
    wg_pr_seg <- wg_pr_seg %>% 
      mutate(start_cum = end_cum- num.mark*1000)
    
    axis_set <- wg_pr_seg %>% 
      group_by(chr) %>% 
      summarize(center = mean(end_cum)) %>% 
      arrange(as.numeric(chr))
    
    label_seg <- wg_pr_seg %>% 
      filter(num.mark > 100) %>% 
      filter(seg.mean > 0.4 | dplyr::between(seg.mean,-1.5, -0.3))
    wg1 <- wg_pr_seg %>% 
      ggplot(aes(x = end_cum, y = seg.mean, color = chr))+
      geom_segment(aes(x = start_cum, y = seg.mean, xend = end_cum, yend = seg.mean+0.001))+
      # geom_point(shape = ".")+
      geom_point(data = label_seg, shape= 8, color = "red")+
      scale_x_continuous(label = axis_set$chr, breaks = axis_set$center)+
      theme_minimal() +
      theme( 
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_text(angle = 60, size = 8, vjust = 0.5)
      )+
      scale_rd+
      scale_size_continuous(range = c(0.5,3))+
      labs(x = NULL)+
      coord_cartesian(expand = F)
    w$hide()
    mod_plot_output_Server("wg_pr_rd", wg1, ranges, dnCNV_table)
    output$wg_rd_table <- renderTable({
      label_seg 
    })
  })
  
  
  ## Chr plot section
  ## Plots section
  
  
  ## RD plots
  observeEvent(input$btn_plot,{
    req(nrow(plots$pr_rd) != 0)
    
    showNotification("Plotting Read Depth Plot", duration = 12, type = "message")
    include_seg <- input$include_seg
    df <- rbindlist(list(plots$pr_seg,plots$m_seg,plots$f_seg))%>%
      filter(ID%in%include_seg)%>%
      mutate(ID=as.factor(ID))%>%
      mutate(seg.mean=ifelse(seg.mean < -2.5,-2,seg.mean))
    plots$xlabel=unique(df$chrom)[1]
    rd <- ggplot(plots$pr_rd, aes(x=V2, y=log2(ratio+0.00001))) +
      geom_point(shape=".")+
      #scattermore::geom_scattermore(shape=".",pixels=c(1024,1024))+
      geom_point(data = subset(plots$pr_rd, ratio < 0.7),aes(V2,log2(ratio+0.00001)),shape=".",color="green")+
      geom_point(data = subset(plots$pr_rd, ratio > 1.3),aes(V2,log2(ratio+0.00001)),shape=".",color="red")+
      geom_segment(data = df,aes(x=loc.start,y=seg.mean,xend=loc.end,yend=seg.mean,color=factor(ID)),size=1)+
      scale_color_manual(name="Segment",values = c("Proband"="blue","Mother"="#E69F00","Father"="#39918C"),)+
      ylim(-4,4)+
      xlab(plots$xlabel)+
      scale_rd+
      style_rd+
      scale_x_continuous(labels = scales::label_number())
    
 
    btnValrds <- mod_checkbox_Server("RD-static")
    ranges <- mod_plot_switch_Server("RD-static", btnValrds$box_state, rd, ranges, dnCNV_table, zoom= F)
    btnValrdd <- mod_checkbox_Server("RD-dynamic")
    ranges <- mod_plot_switch_Server("RD-dynamic", btnValrdd$box_state, rd, ranges, dnCNV_table)

  })
  
  #Baf-B plot
  observeEvent(input$btn_plot,{
    
    req(nrow(plots$snp_chr) != 0)
    
    noti_id <- showNotification("Plotting B-allele frequency plots", type = "message", duration = NULL)
    df <- plots$snp_chr%>%filter(likelyDN%in%c(input$include_dnSNV,"FALSE"))
    cols <- plots$SNPcols
    xlabel=unique(df$chrom)[1]
    
    
    snp_a <- ggplot(df, aes(x=start,y=pr_ALT_Freq,col=A_InhFrom))+
      geom_point(shape=".")+
      #scattermore::geom_scattermore(shape=".",pixels=c(1024,1024))+
      geom_point(data = subset(df, likelyDN %in%c("TRUE")),size = 2,shape=8,color="red")+
      scale_fill_manual("LikelyDN",limits=c("dnSNV"),values = "red")+
      xlab(xlabel)+
      scale_snp+
      style_snp+
      scale_colour_manual(values = cols)+
      guides(color = guide_legend(override.aes = list(size = 4)))+
      scale_x_continuous(labels = scales::label_number())
    
    btnVala <- mod_checkbox_Server("Baf-A_allele")
    ranges <- mod_plot_switch_Server("Baf-A_allele", btnVala$box_state, snp_a, ranges)
    
    snp_b <- ggplot(df, aes(x=start,y=pr_ALT_Freq,col=B_InhFrom))+
      geom_point(shape=".")+
      #scattermore::geom_scattermore(shape=".",pixels=c(1024,1024))+
      geom_point(data = subset(df, likelyDN %in%c("TRUE")),size = 2,shape=8,color="red")+
      scale_fill_manual("LikelyDN",limits=c("dnSNV"),values = "red")+
      xlab(xlabel)+
      scale_snp+
      style_snp+
      scale_colour_manual(values = cols)+
      guides(color = guide_legend(override.aes = list(size = 4)))+
      scale_x_continuous(labels = scales::label_number())
    
    removeNotification(noti_id)
    btnValb <- mod_checkbox_Server("Baf-B_allele")
    ranges <- mod_plot_switch_Server("Baf-B_allele", btnValb$box_state, snp_b, ranges)
    
  })
  
  ##Anno tracks
  observeEvent(input$btn_anno,{
    if (input$ref == "GRCh37"){
      p1_file = "NCBI_RefSeq_hg19_clean.bed.parquet"
      p2_file = "Claudia_hg19_MergedInvDirRpts_sorted.bed"
      p3_file = "SegDup_hg19_UCSC.bed"
      p4_file = "OMIM_gene2_hg19_UCSC_all.bed"
      p5_file = "gnomAD_allSV_hg19_UCSC.bed"
      p6_file = "rmsk_hg19_UCSC.parquet"
    } else {
      p1_file = "MANE.GRCh38.v1.0.refseq.gz.parquet"
      p2_file = NULL
      p3_file = NULL
      p4_file = NULL
      p5_file = NULL
      p6_file = NULL
    }
    
    if (nrow(values$pr_sv) != 0){
    pr_sv <- values$pr_sv %>% 
      filter(CHROM == chrn) %>% 
      filter(AVGLEN > 10000 & AVGLEN < 100000000)
    pr_sv <- pr_sv %>% 
      mutate(color = case_when(SVTYPE == "DEL" ~ "darkblue",
                               SVTYPE == "DUP" ~ "#8b0000",
                               SVTYPE == "INS" ~ "darkgreen", 
                               SVTYPE == "INV" ~ "darkorange", 
                               TRUE ~ "magenta3")) %>% 
      group_by(SVTYPE) %>% 
      mutate(idx = sample.int(n())/1000)
    pr_sv <- pr_sv %>% 
      mutate(start = POS, 
             end = as.integer(END)) %>% 
      relocate(CHROM, start, end) %>% 
      dplyr::select(-c(POS, END, REF, ALT, AVGLEN, MAPQ, RE, CIEND, CIPOS))
    pr_sv_plot <- ggplot(pr_sv, aes(x = POS, y = idx)) +
      annotate("rect", xmin = pr_sv$start, xmax = pr_sv$end, ymin = pr_sv$idx, ymax = pr_sv$idx+0.0001, color = pr_sv$color)+
      style_anno+
      scale_anno+
      ylab("pr_SV")
    btnVal_prsv <- mod_checkbox_Server("pr_sv")
    ranges <- mod_plot_switch_Server("pr_sv", btnVal_prsv$box_state, pr_sv_plot, ranges, dnCNV_table)
    anno_table_Server("pr_sv", pr_sv, ranges, chrn)
    }
    
    id <- showNotification("Pulling Data", type = "message", duration = NULL)
    chrn = input$chr
    path = "./data/"
    RefSeq <- read_parquet(paste0(path,p1_file))
    RefSeq <- RefSeq %>% 
      filter(seqname ==  chrn) %>%
      mutate(gene_num=round(as.numeric(as.factor(gene_id)),3)/100,
             strand=as.factor(strand))
      gene_x <- RefSeq%>%
        group_by(gene_num)%>%
        summarise(start=(min(start)+max(end))/2,
                  end=(min(start)+max(end))/2,
                  gene_id=unique(gene_id))
    p1 <- RefSeq %>% 
      ggplot(aes(xstart = start, xend = end, y = gene_num))+
      geom_intron(data = to_intron(RefSeq, "gene_num"), arrow.min.intron.length = 400)+
      geom_text(data=gene_x,aes(x=start,label=gene_id),vjust = -1.2,check_overlap = T,fontface="italic")+
      ylab("RefSeq")+
      style_anno+
      scale_genes+
      scale_fill_manual(values = c("+"="#E69F00","-"="#39918C"))+
      scale_x_continuous(labels = scales::label_number())
    
    IDR<-data.table::fread(paste0(path, p2_file)) %>% 
      dplyr::select(V1,V2,V3,V4,V5,V6) %>% 
      dplyr::rename("chrom" = V1, "start" = V2, "end" = V3, "name" = V4, "score" = V5, "strand" = V6) %>%
      dplyr::filter(chrom == chrn)
    IDR <- IDR %>%
      mutate(idx = sample(1:1, size = dim(IDR)[1], replace = T)/1000)
    IDR_pos <- IDR %>%
      filter(strand == "+")
    IDR_neg <- IDR %>%
      filter(strand == "-")
    p2 <- ggplot(IDR, aes(x = start, y = idx)) +
      geom_segment(data = IDR_pos, aes(x = start, y = idx, xend = end, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = "dodgerblue")+
      geom_segment(data = IDR_neg, aes(x = end, y = idx, xend = start, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = "dodgerblue3")+
      style_anno+
      ylab("IDR")
    SegDup <- data.table::fread(paste0(path, p3_file))
    SegDup <- SegDup %>%
      dplyr::rename("chrom" = "#chrom", "start" = "chromStart", "end" = "chromEnd") %>%
      filter(chrom ==chrn)
    SegDup <- SegDup %>%
      mutate(idx = sample(1:100, size = dim(SegDup)[1], replace = T)/1000) %>%
      mutate(color = case_when(level == "SD_low" ~ "gray40",
                               level == "SD_mid" ~ "yellow2",
                               level == "SD_high" ~ "darkorange"))
    SegDup_pos <- SegDup %>%
      filter(strand == "+")
    SegDup_neg <- SegDup %>%
      filter(strand == "-")
    p3 <- ggplot(SegDup, aes(x = start, y = idx)) +
      geom_segment(data = SegDup_pos, aes(x = start, y = idx, xend = end, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = SegDup_pos$color)+
      geom_segment(data = SegDup_neg, aes(x = end, y = idx, xend = start, yend = idx), arrow = arrow(length = unit(0.05, "inches")), color = SegDup_neg$color)+
      style_anno+
      scale_anno+
      ylab("SegDup")
    
    
    OMIM <- data.table::fread(paste0(path, p4_file))
    OMIM <- OMIM %>%
      filter(chrom == chrn) %>%
      mutate_at(vars(pheno_key), as.factor)
    OMIM <- OMIM %>%
      mutate(idx = sample(1:100, size = dim(OMIM)[1], replace = T)/1000) %>%
      mutate(color = case_when(pheno_key == "0" ~ "grey",
                               pheno_key == "1" ~ "lightgreen",
                               pheno_key == "2" ~ "green3",
                               pheno_key == "3" ~ "green4",
                               pheno_key == "4" ~ "purple", ))
    
    OMIM_label <- OMIM %>%
      filter(pheno_key  %in% c("3","4"))
    p4 <- ggplot(OMIM, aes(x = start, y = idx)) +
      annotate("rect", xmin = OMIM$start, xmax = OMIM$end, ymin = OMIM$idx, ymax = OMIM$idx+0.0001, color = OMIM$color)+
      annotate("text", x = OMIM_label$start, y = OMIM_label$idx, label = OMIM_label$gene_symbol, size = 3.5, hjust = 1.1)+
      style_anno+
      scale_anno+
      ylab("OMIM")
    
    gnomAD <- data.table::fread(paste0(path, p5_file))
    gnomAD <- gnomAD %>%
      dplyr::rename("chrom" = "#chrom", "start" = "chromStart", "end" = "chromEnd") %>%
      filter(chrom == chrn)
    gnomAD <- gnomAD %>%
      mutate(idx = sample(1:100, size = dim(gnomAD)[1], replace = T)/1000) %>%
      mutate(color = case_when(svtype == "BND" ~ "grey",
                               svtype == "DEL" ~ "red",
                               svtype == "DUP" ~ "blue",
                               svtype == "INS" ~ "darkorange",
                               TRUE ~ "cyan4"))
    p5 <- ggplot(gnomAD, aes(x = start, y = idx)) +
      annotate("rect", xmin = gnomAD$start, xmax = gnomAD$end, ymin = gnomAD$idx, ymax = gnomAD$idx+0.0001, color = gnomAD$color)+
      style_anno+
      scale_anno+
      ylab("gnomAD")
    
    
    rmsk <- read_parquet(paste0(path, p6_file))
    rmsk <- rmsk %>%
      dplyr::rename("chrom" = "#chrom") %>%
      filter(chrom ==  chrn)
    rmsk <- rmsk %>%
      mutate(idx = case_when(repClass == "SINE" ~ 0.014*7,
                             repClass == "LINE" ~ 0.014*6,
                             repClass == "LTR" ~ 0.014*5,
                             repClass == "DNA" ~ 0.014*4,
                             repClass == "Simple_repeat" ~ 0.014*3,
                             repClass == "Low_complexity" ~ 0.014*2,
                             TRUE ~ 0.014))
    p6 <- ggplot(rmsk, aes(x = start, y = idx)) +
      annotate("rect", xmin = rmsk$start, xmax = rmsk$end, ymin = rmsk$idx, ymax = rmsk$idx+0.0001, color = "black")+
      style_anno+
      scale_anno+
      ylab("RMSK")
    
    
    showNotification("Annotating", type = "message", duration = 8)
  
    btnVal1 <- mod_checkbox_Server("RefSeq")
    ranges <- mod_plot_switch_Server("RefSeq", btnVal1$box_state, p1, ranges, dnCNV_table)
    btnVal2 <- mod_checkbox_Server("IDR")
    ranges <- mod_plot_switch_Server("IDR", btnVal2$box_state, p2, ranges, dnCNV_table)
    btnVal3 <- mod_checkbox_Server("SegDup")
    ranges <- mod_plot_switch_Server("Segdup", btnVal3$box_state, p3, ranges, dnCNV_table)
    btnVal4 <- mod_checkbox_Server("OMIM")
    ranges <- mod_plot_switch_Server("OMIM", btnVal4$box_state, p4, ranges, dnCNV_table)
    btnVal5 <- mod_checkbox_Server("gnomAD")
    ranges <- mod_plot_switch_Server("gnomAD", btnVal5$box_state, p5, ranges, dnCNV_table)
    btnVal6 <- mod_checkbox_Server("RMSK")
    ranges <- mod_plot_switch_Server("RMSK", btnVal6$box_state, p6, ranges, dnCNV_table)
  
    removeNotification(id)
    anno_table_Server("RefSeq", RefSeq, ranges, chrn)
    anno_table_Server("IDR", IDR, ranges, chrn)
    anno_table_Server("SegDup", SegDup, ranges, chrn)
    anno_table_Server("OMIM", OMIM, ranges, chrn)
    anno_table_Server("gnomAD", gnomAD, ranges, chrn)
    anno_table_Server("rmsk", rmsk, ranges, chrn)
  })
  
  ## dnCNV table
  observeEvent(input$btn_dnCNV, {
    req(nrow(values$pr_rd)!=0)
    req(nrow(values$m_rd)!=0)
    req(nrow(values$f_rd)!=0)
    dnCNV_table$t <- mod_dnCNV_Server("dnCNV",plots$pr_seg, plots$m_seg, plots$f_seg)
    print("btn gen")
    print(dnCNV_table$t)
  })

  ## Show current ranges
  observe({
    output$cur_range <- renderText({
      req(!is.null(ranges$x))
      paste0("current range: ", round(ranges$x[1]), "-", round(ranges$x[2]), "    width: ", round(ranges$x[2])-round(ranges$x[1])+1)
      })
    })
  
  ## btn_goto
  observeEvent(input$btn_go,{
    req(!is.null(input$goto_reg))
    str <- stringr::str_trim(input$goto_reg)
    str <- strsplit(str,"-|_")
    if (length(str[[1]]) ==1) {
      showNotification("Looking up gene name", type = "message")
      path <- "./data/"
      p1_file <- "NCBI_RefSeq_hg19_clean.bed.parquet"
      RefSeq <- read_parquet(paste0(path,p1_file))
      search <- as.character(str[[1]])
      found <- RefSeq %>% 
        filter(seqname == input$chr) %>% 
        filter(grepl(search, gene_id, ignore.case = T))
      if (nrow(found) != 0){
        ranges$x <- c(as.numeric(min(found$start)-geneExtend),
                      as.numeric(max(found$end)+geneExtend))
      }
    } else if (length(str[[1]]) == 2){
      showNotification("Jumping to coordinates", type = "message")
      from <- as.numeric(str[[1]][1])
      to <- as.numeric(str[[1]][2])
      ranges$x <- c(from, to)
    }
  }, ignoreInit = T)


  
}


