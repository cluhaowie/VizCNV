

source("./mod/mod_plot_output.R")
source("./mod/mod_dnCNV.R")
source("./mod/mod_goto.R")

server <- function(input, output,session) {
  # Reavtive Values --------------------------
  values <- reactiveValues()
  
  values$data <- data.frame(stringsAsFactors = F)
  #values$work_data, store sv call, must include ALT 
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
  values$snp_gvcf_file_ref <- vector()
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
    list(values$data,values$pr_rd)
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
  output$ui_plot_anno <- shiny::renderUI({
    if(is.null(input$include_anno)){
      helpText("")
    } else {
      plotOutput(
        "plot_anno",
        height = 200,
        dblclick = "plot_anno_dblclick",
        brush = brushOpts(id = "plot_anno_brush",direction = "x",
                          resetOnNew = TRUE))
    }
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
  
  
  
 


  
  ## annotation panel
  # interactive plot regions-------
  output$brush_info <- renderPrint({
    brush <- input$plot1_brush
    chr <- input$chr
    if(!is.null(brush)){
      cat(paste0("select region: ",chr,":",
                 format(round(brush$xmin,0),big.mark=",",scientific = FALSE),
                 "-",
                 format(round(brush$xmax,0),big.mark=",",scientific = FALSE),
                 ", ",round(brush$xmax-brush$xmin,0),"bp"))
    }
  })
  
  
  ##work in progress
  ## old plot1
  # output$plot1 <- renderPlot({
  #   if(nrow(plots$pr_rd) == 0){
  #     return(NULL)
  #   }
  #   if(nrow(values$selected_record)==0){
  #     dup.df <- data.frame()
  #     del.df <- data.frame()
  #     trp.df <- data.frame()
  #     cpx.df <- data.frame()
  #   }else{
  #     dup.df <- subset(values$selected_record, ALT%in%c("<DUP>","DUP"))
  #     del.df <- subset(values$selected_record, ALT%in%c("<DEL>","DEL"))
  #     trp.df <- subset(values$selected_record, ALT%in%c("<TRP>","TRP"))
  #     cpx.df <- subset(values$selected_record, ALT%in%c("<CPX>","CPX"))
  #     dup.df$rowidx <- as.numeric(rownames(dup.df))
  #     del.df$rowidx <- as.numeric(rownames(del.df))
  #     trp.df$rowidx <- as.numeric(rownames(trp.df))
  #     cpx.df$rowidx <- as.numeric(rownames(cpx.df))
  #   }
  #   if(nrow(values$anno_rect)==0){
  #     anno_rect <- data.frame()
  #   }else{
  #     anno_rect <- values$anno_rect
  #     anno_rect$rowidx <- as.numeric(rownames(anno_rect))
  #   }
  #   plots$plot1 <- ext1()+
  #     coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)+
  #     ## pointer to the selected CNV call(suppressed)
  #     #annotate("point",x=ifelse(is.null(values$selected_record$POS),-100,values$selected_record$POS),y=rep(1.8,length(values$selected_record$POS)),shape=25,fill="orange",size=2)+ 
  #     geom_hline(yintercept=-2.5,col="black")+
  #     annotate("text",x=ifelse(is.null(ranges$x),-100000,sum(ranges$x)/2),y=-2.3,label="SVTYPE")+
  #     annotate("rect",xmin=dup.df$POS,xmax=dup.df$END,ymin=dup.df$rowidx%%4*0.1-3,ymax=(dup.df$rowidx%%4+1)*0.1-3,fill = CNVCOLOR6["DUP"])+
  #     annotate("rect",xmin=del.df$POS,xmax=del.df$END,ymin=del.df$rowidx%%4*0.1-3,ymax=(del.df$rowidx%%4+1)*0.1-3,fill = CNVCOLOR6["DEL"])+
  #     annotate("rect",xmin=trp.df$POS,xmax=trp.df$END,ymin=trp.df$rowidx%%4*0.1-3,ymax=(trp.df$rowidx%%4+1)*0.1-3,fill = CNVCOLOR6["TRP"])+
  #     annotate("rect",xmin=cpx.df$POS,xmax=cpx.df$END,ymin=cpx.df$rowidx%%4*0.1-3,ymax=(cpx.df$rowidx%%4+1)*0.1-3,fill = CNVCOLOR6["CPX"])+
  #     annotate("rect",xmin=anno_rect$xmin,xmax=anno_rect$xmax,fill=anno_rect$fill,ymin=anno_rect$rowidx*0-3,ymax=(anno_rect$rowidx-anno_rect$rowidx+2),alpha=0.3)+
  #     ggtitle(paste0(plots$xlabel,":",paste0(round(as.numeric(ranges$x)),collapse = "-")))
  #   plots$plot1
  # })
  # 
  # 
  
  # output$plot_anno <- renderPlot({
  #   if(length(plots$plot3) == 0){
  #     return(NULL)
  #   }
  #   if(nrow(values$anno_rect)==0){
  #     anno_rect <- data.frame()
  #   }else{
  #     anno_rect <- values$anno_rect
  #     anno_rect$rowidx <- as.numeric(rownames(anno_rect))
  #   }
  #   plots$plot3_dl <- plots$plot3+
  #     annotate("rect",xmin=anno_rect$xmin,xmax=anno_rect$xmax,fill=anno_rect$fill,ymin=anno_rect$rowidx*0,ymax=(anno_rect$rowidx-anno_rect$rowidx+2)+length(unique(plots$genelabel$gene_id)),alpha=0.3)+
  #     coord_cartesian(xlim = ranges$x, expand = FALSE)
  #   plots$plot3_dl
  # })

  
  #     if(nrow(values$data)!=0){
  #       values$work_data <- values$data%>%
  #         filter(CHROM==chr)%>%
  #         filter(POS>=brush$xmin,POS < brush$xmax)
  #     }
  #   } else {
  #     ranges$x <- NULL
  #     ranges$y <- NULL
  #     plots$genelabel <- data.frame(stringsAsFactors = F)
  #     plots$plot3 <- list()
  #     if(nrow(values$data)!=0){
  #       values$work_data <- values$data%>%filter(CHROM==chr)
  #     }
  #   }
  # })
# 
  
  ## goto btn
#   observeEvent(input$btl_goto,{
#     if(!is.null(input$goto_reg)){
#       gene <- as.character(input$goto_reg)
#       gene.exist <- genebase%>%filter(seqname==input$chr,
#                                       gene_id==gene,
#                                       type%in%c("exon"))%>%
#         dplyr::select(seqname,start,end,strand,transcript_id,gene_id,type)%>%
#         collect()
#       if(nrow(gene.exist)!=0){
#         ranges$x <- c(as.numeric(min(gene.exist$start)-geneExtend),
#                       as.numeric(max(gene.exist$end)+geneExtend))
#         gene.exist <- gene.exist%>%
#           mutate(gene_num=round(as.numeric(as.factor(gene_id)),3),
#                  strand=as.factor(strand))
#         gene_x <- gene.exist%>%
#           group_by(gene_num)%>%
#           summarise(start=(min(start)+max(end))/2,
#                     end=(min(start)+max(end))/2,
#                     gene_id=unique(gene_id))
#         plots$plot3 <- gene.exist%>%
#           ggplot(aes(xstart = start,xend = end,y = gene_num))+
#           ggtranscript::geom_range(aes(fill = strand)) +
#           ggtranscript::geom_intron(data = ggtranscript::to_intron(gene.exist, "gene_num"),aes(strand = strand))+
#           geom_text(data=gene_x,aes(x=start,label=gene_id),vjust = -1.2,check_overlap = T,fontface="italic")+
#           style_genes+scale_genes+
#           scale_fill_manual(values = c("+"="#E69F00","-"="#39918C"))+
#           scale_x_continuous(labels = scales::label_number())
#       }else{
#         goto_reg <- as.numeric(unlist(strsplit(input$goto_reg,"-|_")))
#         validate(need(length(goto_reg)==2,"Check at least two coordinates"))
#         ranges$x <- c(min(goto_reg)-geneExtend,max(goto_reg)+geneExtend)
#         plots$genelabel <- genebase%>%
#           filter(seqname==input$chr,
#                  start>(min(ranges$x)-geneExtend),
#                  end < (max(ranges$x)+geneExtend),
#                  type%in%c("exon"))%>%
#           dplyr::select(seqname,start,end,strand,transcript_id,gene_id,type)%>%
#           dplyr::collect()%>%
#           mutate(gene_num=round(as.numeric(as.factor(gene_id)),3),
#                  strand=as.factor(strand))
#         if(length(unique(plots$genelabel$gene_id))<maxtranscript){
#           gene_x <- plots$genelabel%>%
#             group_by(gene_num)%>%
#             summarise(start=(min(start)+max(end))/2,
#                       end=(min(start)+max(end))/2,
#                       gene_id=unique(gene_id))
#           plots$plot3 <- plots$genelabel%>%
#             ggplot(aes(xstart = start,xend = end,y = gene_num))+
#             ggtranscript::geom_range(aes(fill = strand)) +
#             ggtranscript::geom_intron(data = ggtranscript::to_intron(plots$genelabel, "gene_num"),aes(strand = strand))+
#             geom_text(data=gene_x,aes(x=start,label=gene_id),vjust = -1.2,check_overlap = T,fontface="italic")+
#             style_genes+scale_genes+
#             scale_fill_manual(values = c("+"="#E69F00","-"="#39918C"))+
#             scale_x_continuous(labels = scales::label_number())
#         }
#         
#       }
#       
#     }
#   })
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
    if(nrow(values$data) > 0){
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
  
  
  ## Plots section
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  observeEvent(input$btn_plot,{
    req(nrow(plots$pr_rd) != 0)
    
    showNotification("Plotting Read Depth Plot", duration = 8, type = "message")
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
    ranges <- mod_plot_switch_Server("RD-static", btnValrds$box_state, rd, ranges, zoom= F)
    btnValrdd <- mod_checkbox_Server("RD-dynamic")
    ranges <- mod_plot_switch_Server("RD-dynamic", btnValrdd$box_state, rd, ranges)
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
    btnVal2 <- mod_checkbox_Server("IDR")
    btnVal3 <- mod_checkbox_Server("SegDup")
    btnVal4 <- mod_checkbox_Server("OMIM")
    btnVal5 <- mod_checkbox_Server("gnomAD")
    btnVal6 <- mod_checkbox_Server("RMSK")
    ranges <- mod_plot_switch_Server("RefSeq", btnVal1$box_state, p1, ranges)
    ranges <- mod_plot_switch_Server("IDR", btnVal2$box_state, p2, ranges)
    ranges <- mod_plot_switch_Server("Segdup", btnVal3$box_state, p3, ranges)
    ranges <- mod_plot_switch_Server("OMIM", btnVal4$box_state, p4, ranges)
    ranges <- mod_plot_switch_Server("gnomAD", btnVal5$box_state, p5, ranges)
    ranges <- mod_plot_switch_Server("RMSK", btnVal6$box_state, p6, ranges)
  
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
    mod_dnCNV_Server("dnCNV",plots$pr_seg, plots$m_seg, plots$f_seg)
  })
  


  observe({
    output$cur_range <- renderText({
      req(!is.null(ranges$x))
      paste0("current range: ", round(ranges$x[1]), "-", round(ranges$x[2]), "    width: ", round(ranges$x[2])-round(ranges$x[1])+1)
      })
    })
  
  observeEvent(input$btn_go,{
    str <- stringr::str_trim(input$goto_reg)
    goto_coord <- as.numeric(unlist(strsplit(str,"-|_")))
    from <- goto_coord[1]
    to <- goto_coord[2]
    ranges$x <- c(from, to)
  }, ignoreInit = T)


  
}