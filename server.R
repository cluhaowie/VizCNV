source("./mod/mod_plot_output.R")
source("./mod/mod_dnCNV.R")
source("./mod/mod_upload.R")
source("./mod/mod_UCSC.R")
source("./helper/wg_plot.R")

server <- function(input, output,session) {
  # Reavtive Values --------------------------
  values <- reactiveValues()
  values$pr_sv <- data.frame(stringsAsFactors = F)
  values$m_sv <- data.frame(stringsAsFactors = F)
  values$f_sv <- data.frame(stringsAsFactors = F)
  values$pr_rd <- data.frame(stringsAsFactors = F)
  values$m_rd <- data.frame(stringsAsFactors = F)
  values$f_rd <- data.frame(stringsAsFactors = F)
  values$snp_gvcf_file_ref <- vector()
  values$snp_gvcf_file_path <- vector()
  values$snp_gvcf_file <- data.frame(stringsAsFactors = F)
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
  
  
  volumes <- c(Home="~/Downloads/","R installation" = R.home(),shinyFiles::getVolumes()())
  mod_rd_upload_Server("pr_rd",volumes=volumes,values) 
  mod_rd_upload_Server("m_rd",volumes=volumes,values) 
  mod_rd_upload_Server("f_rd",volumes=volumes,values) 
  mod_sv_upload_Server("pr_sv",volumes=volumes,values) 
  mod_sv_upload_Server("m_sv",volumes=volumes,values) 
  mod_sv_upload_Server("f_sv",volumes=volumes,values) 
  mod_snp_upload_Server("snp_file",volumes=volumes,values)

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
  output$ui_chkbox_RefSeq <- renderUI({
    if(!is.null(values$p1_file)){
      mod_checkbox_UI("RefSeq")
    }else{NULL}

  })
  output$ui_chkbox_IDR <- renderUI({
    if(!is.null(values$p2_file)){
      mod_checkbox_UI("IDR",value = F)
    }else{NULL}
  })
  output$ui_chkbox_SegDup <- renderUI({
    if(!is.null(values$p3_file)){
      mod_checkbox_UI("SegDup",value = F)
    }else{NULL}
  })
  output$ui_chkbox_OMIM <- renderUI({
    if(!is.null(values$p4_file)){
      mod_checkbox_UI("OMIM",value = F)
    }else{NULL}
  })
  output$ui_chkbox_gnomAD <- renderUI({
    if(!is.null(values$p5_file)){
      mod_checkbox_UI("gnomAD",value = F)
    }else{NULL}
  })

  output$ui_chkbox_RMSK <- renderUI({
    if(!is.null(values$p6_file)){
      mod_checkbox_UI("RMSK",value = F)
    }else{NULL}
  })


  
  # observe file uploaded and save in SQLdatabase---------

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
    if(input$ref=="hg38"){
      blacklist <- data.table::fread("GRCh38_unified_blacklist.bed.gz")%>%
        regioneR::toGRanges()
      values$ref_info <- data.table::fread("hg38.info.txt")
      values$p1_file <-  "hg38_MANE.v1.0.refseq.parquet"
      values$p2_file <-  NULL
      values$p3_file <-  "hg38_ucsc_sugdups.parquet"
      values$p4_file <-  NULL
      values$p5_file <-  NULL
      values$p6_file <-  "hg38_rmsk.parquet"
    }else if(input$ref=="hg19"){
      blacklist <- data.table::fread("ENCFF001TDO.bed.gz")%>%
        regioneR::toGRanges()
      values$ref_info <- data.table::fread("hg19.info.txt")
      values$p1_file <-  "NCBI_RefSeq_hg19_clean.bed.parquet"
      values$p2_file <-  "Claudia_hg19_MergedInvDirRpts_sorted.bed"
      values$p3_file <-  "hg19_ucsc_sugdups.parguet"
      values$p4_file <-  "OMIM_gene2_hg19_UCSC_all.bed"
      values$p5_file <-  "gnomAD_allSV_hg19_UCSC.bed"
      values$p6_file <-  "hg19_rmsk.parquet"
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
    norm_option <- input$norm_options
    chr <- input$chr
    if(nrow(values$pr_rd)==0){return(NULL)
    }else{
      w$show()
      plots$pr_rd <- normalization_method(values$pr_rd, chr, norm_option)
      plots$pr_seg <- SegNormRD(plots$pr_rd,id="Proband",seg.method = seg_option)
      w$hide()
    }
    if(nrow(values$m_rd)==0){return(NULL)
    }else{
      w$show()
      plots$m_rd <- normalization_method(values$m_rd, chr, norm_option)
      plots$m_seg <- SegNormRD(plots$m_rd,id="Mother",seg.method = seg_option)
      w$hide()
    }
    if(nrow(values$f_rd)==0){return(NULL)
    }else{
      w$show()
      plots$f_rd <- normalization_method(values$f_rd, chr, norm_option)
      plots$f_seg <- SegNormRD(plots$f_rd,id="Father",seg.method = seg_option)
      w$hide()
    }
  })
  observeEvent(input$btn_filter,{
    snp_gvcf_file_path=values$snp_gvcf_file_path
    ref_genome <- input$ref
    chr <- input$chr
    if(!chr%in%values$snp_gvcf_file_ref){
      if(chr%in%chrom_id){
        chr <- names(chrom_id)[which(chrom_id==chr)]
      }else{
        chr <- chrom_id[which(names(chrom_id)==chr)]
      }
    }
    if(length(snp_gvcf_file_path)==0){
      return(NULL)
    }else{
      w$show()
      loc.start <- 0
      loc.end <- values$ref_info%>%
        filter(chrom%in%c(chr,paste0("chr",chr)))%>%
        dplyr::select(seqlengths)%>%unlist
      range.gr <- GenomicRanges::GRanges(chr,ranges = IRanges(loc.start,loc.end))
      range.gr <- GenomicRanges::setdiff(range.gr, blacklist)
      plots$snp_chr <- ReadGVCF(snp_gvcf_file_path,ref_genome=input$ref,param = range.gr)%>%
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
  
  
  # Plots -----
  ### "Global" reactive values
  wg_ranges <- reactiveValues(x = NULL, pr = NULL, m = NULL, f = NULL)
  wg_dnCNV_table <- reactiveValues(t = data.frame(start_cum = c(0), end_cum = c(0), stringsAsFactors = F))
  
  ## reset plots upon changing chr
  observeEvent(input$btn_wg_rd, {
    wg_ranges$x <-  NULL
    wg_dnCNV_table$t <-  data.frame(start_cum = c(0), end_cum = c(0), stringsAsFactors = F)
  })
  
  ## WG Plot section
  observeEvent(input$btn_wg_rd, {
    req(nrow(values$pr_rd) != 0)
    w$show()
    rd <- values$pr_rd
    rd <- wg_norm(rd, input$wg_norm_options)
    wg_ranges$pr <- getAllSeg(rd)
    wg_pr <- wg_seg2plot(wg_ranges$pr)
    w$hide()
    wg_ranges <- mod_plot_wg_Server("wg_pr_rd", wg_pr, wg_ranges, wg_dnCNV_table)
  })
  observeEvent(input$btn_wg_rd, {
    req(nrow(values$m_rd) != 0)
    w$show()
    rd <- values$m_rd
    rd <- wg_norm(rd, input$wg_norm_options)
    wg_ranges$m <- getAllSeg(rd)
    wg_m <- wg_seg2plot(wg_ranges$m)
    w$hide()
    wg_ranges <- mod_plot_wg_Server("wg_m_rd", wg_m, wg_ranges, wg_dnCNV_table)
  })
  observeEvent(input$btn_wg_rd, {
    req(nrow(values$f_rd) != 0)
    w$show()
    rd <- values$f_rd
    rd <- wg_norm(rd, input$wg_norm_options)
    wg_ranges$f <- getAllSeg(rd)
    wg_f <- wg_seg2plot(wg_ranges$f)
    w$hide()
    wg_ranges <- mod_plot_wg_Server("wg_f_rd", wg_f, wg_ranges, wg_dnCNV_table)
  })

  ## dnCNV table
  observeEvent(input$btn_wg_dnCNV, {
    req(!is.null(wg_ranges$pr))
    req(!is.null(wg_ranges$m))
    req(!is.null(wg_ranges$f))
    seg_data <- mod_dnCNV_Server("wg_dnCNV",wg_ranges$pr, wg_ranges$m, wg_ranges$f)
    names(seg_data)[1] = "chr"
    print(seg_data)
    temp <- wg_ranges$pr %>% 
      group_by(chr) %>% 
      summarise(max_end = max(loc.end)) %>% 
      mutate(across("chr", str_replace, "chr", "")) %>% 
      arrange(as.numeric(chr)) %>% 
      mutate(loc_add = lag(cumsum(as.numeric(max_end)), default = 0)) %>% 
      mutate(chr = paste0("chr", chr))
    seg_data <- seg_data %>% 
      inner_join(temp, by = "chr") %>% 
      mutate(end_cum = loc_add + end) 
    seg_data <- seg_data %>% 
      mutate(start_cum = end_cum- width)
    wg_dnCNV_table$t <- seg_data
  })
  

  
  
  ## Plots section
  ranges <- reactiveValues(x = NULL, y = NULL,cur = NULL, click = NULL, pb = NULL)
  dnCNV_table <- reactiveValues(t = data.frame(start = c(0), end = c(0), stringsAsFactors = F),
                                hl = data.frame(start = c(0), end = c(0)), 
                                hl_col = c("white"))
  
  ## reset plots upon changing chr
  observeEvent(input$btn_plot, {
    ranges$x <-  NULL
    ranges$cur <-  NULL
    ranges$click <-  NULL
    ranges$pb <-  NULL
    dnCNV_table$t <-  data.frame(start = c(0), end = c(0), stringsAsFactors = F)
  })

  ## dynamic highlight
  mod_col_pick_Server("highlight", dnCNV_table, ranges)
  
  ## RD plots
  observeEvent(input$btn_plot,{
    req(nrow(plots$pr_rd) != 0)
    
    showNotification("Plotting Read Depth Plot", duration = 12, type = "message")
    include_seg <- input$include_seg
    df <- rbindlist(list(plots$pr_seg,plots$m_seg,plots$f_seg))%>%
      filter(ID%in%include_seg)%>%
      mutate(ID=as.factor(ID))%>%
      mutate(seg.mean=ifelse(seg.mean < -2.5,-2,seg.mean))
    plots$xlabel=input$chr
    rd <- ggplot(plots$pr_rd, aes(x=V2, y=log2(ratio+0.00001))) +
      geom_point(shape=".")+
      #scattermore::geom_scattermore(shape=".",pixels=c(1024,1024))+
      geom_point(data = subset(plots$pr_rd, ratio < 0.7),aes(V2,log2(ratio+0.00001)),shape=".",color="green")+
      geom_point(data = subset(plots$pr_rd, ratio > 1.3),aes(V2,log2(ratio+0.00001)),shape=".",color="red")+
      geom_segment(data = df,aes(x=loc.start,y=seg.mean,xend=loc.end,yend=seg.mean,color=factor(ID)),size=1)+
      scale_color_manual(name="Segment",values = c("Proband"="blue","Mother"="#E69F00","Father"="#39918C"),)+
      xlab(plots$xlabel)+
      scale_rd+
      style_rd
    
 
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
    ranges <- mod_plot_switch_Server("Baf-A_allele", btnVala$box_state, snp_a, ranges, dnCNV_table)
    
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
    ranges <- mod_plot_switch_Server("Baf-B_allele", btnValb$box_state, snp_b, ranges, dnCNV_table)
    
  })
  
  ##Anno tracks
  observeEvent(input$btn_anno,{
    id <- showNotification("Pulling Data", type = "message", duration = NULL)
    chrn = input$chr
    path = "./data/"
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
    
    if(!is.null(values$p1_file)){
      RefSeq_data <- read_parquet(paste0(path,values$p1_file),as_data_frame = F)
      RefSeq <- RefSeq_data %>% 
        filter(seqname ==  chrn) %>% 
        collect() %>%
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
      anno_table_Server("RefSeq", RefSeq, ranges, chrn)
    }
    if(!is.null(values$p2_file)){
      IDR<-data.table::fread(paste0(path, values$p2_file)) %>% 
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
      anno_table_Server("IDR", IDR, ranges, chrn)
    }
    if(!is.null(values$p3_file)){
      SegDup_data <- arrow::read_parquet(paste0(path, values$p3_file),as_data_frame = F)
      SegDup <- SegDup_data %>%
        dplyr::select(segdups.keep.col)%>%
        filter(chrom ==chrn)%>%
        dplyr::rename("start" = "chromStart", "end" = "chromEnd") %>%
        collect()
        
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
      anno_table_Server("SegDup", SegDup, ranges, chrn)
    }
    if(!is.null(values$p4_file)){
      OMIM <- data.table::fread(paste0(path, values$p4_file))
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
      anno_table_Server("OMIM", OMIM, ranges, chrn)
    }
    if(!is.null(values$p5_file)){
      gnomAD <- data.table::fread(paste0(path, values$p5_file))
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
      anno_table_Server("gnomAD", gnomAD, ranges, chrn)
      
    }
    if(!is.null(values$p6_file)){
      rmsk_data <- read_parquet(paste0(path, values$p6_file),as_data_frame = F)
      rmsk <- rmsk_data %>% 
        filter(chrom ==  chrn) %>% 
        collect()
      p6 <- rmsk %>% 
        mutate(idx = case_when(repClass == "SINE" ~ 0.014*7,
                               repClass == "LINE" ~ 0.014*6,
                               repClass == "LTR" ~ 0.014*5,
                               repClass == "DNA" ~ 0.014*4,
                               repClass == "Simple_repeat" ~ 0.014*3,
                               repClass == "Low_complexity" ~ 0.014*2,
                               TRUE ~ 0.014))%>%
        ggplot(., aes(x = start, y = idx)) +
        annotate("rect", xmin = rmsk$start, xmax = rmsk$end, ymin = rmsk$idx, ymax = rmsk$idx+0.0001, color = "black")+
        style_anno+
        scale_anno+
        ylab("RMSK")
      anno_table_Server("rmsk", rmsk, ranges, chrn)
    }
    showNotification("Annotating", type = "message", duration = 8)
    btnVal1 <- mod_checkbox_Server("RefSeq")
    btnVal2 <- mod_checkbox_Server("IDR")
    btnVal3 <- mod_checkbox_Server("SegDup")
    btnVal4 <- mod_checkbox_Server("OMIM")
    btnVal5 <- mod_checkbox_Server("gnomAD")
    btnVal6 <- mod_checkbox_Server("RMSK")
    ranges <- mod_plot_switch_Server("RefSeq", btnVal1$box_state, p1, ranges, dnCNV_table)
    ranges <- mod_plot_switch_Server("IDR", btnVal2$box_state, p2, ranges, dnCNV_table)
    ranges <- mod_plot_switch_Server("Segdup", btnVal3$box_state, p3, ranges, dnCNV_table)
    ranges <- mod_plot_switch_Server("OMIM", btnVal4$box_state, p4, ranges, dnCNV_table)
    ranges <- mod_plot_switch_Server("gnomAD", btnVal5$box_state, p5, ranges, dnCNV_table)
    ranges <- mod_plot_switch_Server("RMSK", btnVal6$box_state, p6, ranges, dnCNV_table)
  
    removeNotification(id)
  })
  
  
  ## dnCNV table
  observeEvent(input$btn_dnCNV, {
    req(nrow(values$pr_rd)!=0)
    req(nrow(values$m_rd)!=0)
    req(nrow(values$f_rd)!=0)
    dnCNV_table$t <- mod_dnCNV_Server("dnCNV",plots$pr_seg, plots$m_seg, plots$f_seg)
  })
  
  ## Show current ranges
  observe({
    output$cur_range <- renderText({
      req(!is.null(ranges$pb))
      paste0("range: ", round(ranges$pb[1]), "-", round(ranges$pb[2]), "    width: ", round(ranges$pb[2])-round(ranges$pb[1])+1)
    })
  })
  
  ## Show cur location
  observe({
    output$cur_loc <- renderText({
      req(!is.null(ranges$cur))
      paste0("location: ", round(ranges$cur))
    })
  })
  
  ## copy clicked location
  # observe({
  #   if(is.null(ranges$click)){ranges$click <- 0}
  #   clipr::write_clip(str_remove_all(as.character(round(ranges$click)), "[\r\n]"), allow_non_interactive = T)
  # })

  ## btn_goto
  observeEvent(input$btn_go,{
    req(!is.null(input$goto_reg))
    str <- stringr::str_trim(input$goto_reg)
    str <- strsplit(str,"-|_")
    if (length(str[[1]]) == 1) {
      if (!is.na(as.numeric(str[[1]]))){
        showNotification("Jumping to coordinates", type = "message")
        from <- as.numeric(str[[1]])-500000
        if (from < 0){from <- 0}
        to <- as.numeric(str[[1]])+500000
        ranges$x <- c(from, to)
      }else {
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
        }else{
          showNotification("Gene not found", type = "error")
        }

      }
    } else if (length(str[[1]]) == 2){
      showNotification("Jumping to coordinates", type = "message")
      from <- as.numeric(str[[1]][1])
      to <- as.numeric(str[[1]][2])
      ranges$x <- c(from, to)
    }
  }, ignoreInit = T)
  
  ##UCSC
  observe({
    mod_UCSC_Server("UCSC", input$ref, input$chr, ranges)
  })
  
  
  
  # ## buttons 
  # output$ui_dlbtn_tbl <- renderUI({
  #   if(nrow(values$pr_sv) > 0){
  #     tagList(shiny::actionButton("btl_select", "Select",icon("check")))
  #   }
  # })
  # output$ui_dlbtn_plt <- renderUI({
  #   if(length(plots$plot1) > 0){
  #     downloadButton("dl_plt", "Download")
  #   }
  # })
  # output$ui_clbtn_plt <- renderUI({
  #   if(length(plots$plot1) > 0){
  #     shiny::actionButton("cl_btn","Clear plot",icon("trash"))
  #   }
  # })
  # output$ui_dlbtn_dnsnv <- renderUI({
  #   if(length(plots$plot2) > 0){
  #     shiny::downloadButton("dl_btn_dnsnv","Download dnSNV")
  #   }
  # })
  # 
  
  # observeEvent(input$cl_btn,{
  #   plots$snp_chr <- data.frame(stringsAsFactors = F)
  #   plots$pr_rd <- data.frame(stringsAsFactors = F)
  #   input$filter_sv_table_rows_selected <- NULL
  # })
  # 
  # ## Download handler
  # output$dl_plt <- downloadHandler(
  #   filename = function(){
  #     paste0(input$chr,".pdf")
  #   },
  #   content = function(file){
  #     
  #     mylist <- list(plots$plot1,plots$plot3_dl,plots$plot2)
  #     mylist <- mylist[lengths(mylist)!= 0]
  #     n <- length(mylist)
  #     p <- cowplot::plot_grid(plotlist=mylist,ncol = 1,align = 'v',axis = 'lr')
  #     ggplot2::ggsave(filename =file, plot = p,device = "pdf",width =12 ,height = n*4,units = "in")
  #   }
  # )
  # output$dl_btn_dnsnv <- downloadHandler(
  #   filename = function(){paste("dnSNV_",input$chr,".csv")},
  #   content = function(file){
  #     df <- plots$snp_chr%>%filter(likelyDN%in%c(input$include_dnSNV))
  #     write.csv(df,file,row.names = F)
  #   }
  # )
  # 
  
}


