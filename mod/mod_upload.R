###WORK IN PROGRESS###
###UNLIKELY TO IMPLEMENT THIS##



mod_snp_upload_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"),label = paste0(id),accept=c("*.vcf","*.vcf.gz"),multiple = F,buttonLabel = "Browse...")
  )
}

mod_snp_upload_Server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observe({
        if(is.null(input$file)){return(NULL)}
        req(input$file)
        Rsamtools::indexTabix(input$file$datapath,format = "vcf")
        #read the header of p.VCF file to checking the chr id
        ref <- VariantAnnotation::scanVcfHeader(input$file$datapath)@reference
        showModal(modalDialog(
          title = "File upload",
          "The joint called SNP file has been uploaded and indexed"
        ))
      })
      return (ref)
    }
  )
}

mod_sv_upload_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"),label = paste0(id),accept=c("*.vcf","*.vcf.gz"),multiple = F,buttonLabel = "Browse...")
  )
}

mod_sv_upload_Server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observe({
        if(is.null(input$file)){return(NULL)}
        req(input$file)
        df <- bedr::read.vcf(input$file$datapath,split.info = T)$vcf
        showModal(modalDialog(
          title = "File upload",
          "The sv file has been uploaded"
        ))
      })
      return (df)
    }
  )
}

mod_rd_upload_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fileInput(ns("file"),label = paste0(id), accept=c("*.bed","*.bed.gz"),multiple = F, buttonLabel = "Browse...")
  )
}

mod_rd_upload_Server <- function(id) {
  moduleServer(
    id,
    function(input, output, session) {
      observeEvent(input$file,{ 
        df <- data.table::fread(input$file$datapath,header = F) %>% 
          mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))%>%
          regioneR::toGRanges()
        # df <- df[!df%over%blacklist,]
        df <-  df %>%
          as.data.table()%>%
          dplyr::select(-c("strand","width"))
        setnames(df,c("seqnames","start","end"),c("V1","V2","V3"))
        showModal(modalDialog(title = "File upload",
                              "The proband read depth file has been uploaded"))
        
      },ignoreInit = T)
      return (df) }
  )
}
