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
    fluidRow(
      column(width = 6,
             fileInput(ns("file"),label = NULL, accept=c("*.vcf","*.vcf.gz"),multiple = F, buttonLabel = "Browse...")
      ),
      column(width=1,h5("or")),
      column(width = 4,
             shinyFilesButton(id = ns("local_sv_file"), label = "Browse...", title = "Please select a file", multiple = F, viewtype = "detail")
      )
    )
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

#' shiny file input wrapper UI
#'
#' 
#'
#' 
#'
#' @return tagList of Shiny UI elements depends on type
#'
#' @examples
#' mod_rd_upload_UI("plot", type = "required")
#'
#' @export

mod_rd_upload_UI <- function(id) {
  ns <- NS(id)
    tagList(
      fluidRow(
        column(width = 6,
               fileInput(ns("file"),label = NULL, accept=c("*.bed","*.bed.gz"),multiple = F, buttonLabel = "Browse...")
        ),
        column(width=1,h5("or")),
        column(width = 4,
               shinyFilesButton(id = ns("local_rd_file"), label = "Browse...", title = "Please select a file", multiple = F, viewtype = "detail")
        )
      )
    )
  
  
}

#' shiny file input wrapper server
#'
#' 
#'
#' @param volumes a names vector with default local path to explore
#' @param values, a reactivevalues object to store read depth information
#' 
#'
#' @return no return, add a dataframe to values[[id]]
#'
#' @examples
#' mod_rd_upload_Server("proband", volumes, values)
#'
#' @export
mod_rd_upload_Server <- function(id,volumes,values) {
  moduleServer(
    id,
    function(input, output, session) {
      shinyFileChoose(input, "local_rd_file", roots = volumes, session=session)
      observeEvent(input$file, {
        values[[id]] <- data.table::fread(input$file$datapath,header = F) %>%
          mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))
        assign(id,values$rd_df)
        showModal(modalDialog(title = "File upload",
                              "The proband read depth file has been uploaded"))
      })
      observeEvent(input$local_rd_file,{
        if(is.integer(input$local_rd_file)){
          cat("no file\n")
        }else{
          local_rd_file <- parseFilePaths(volumes, input$local_rd_file)
          values[[id]]<- data.table::fread(local_rd_file$datapath,header = F)%>%
            mutate(V1=ifelse(!V1%like%"chr",paste0("chr",V1),V1))
          showModal(modalDialog(title = "File upload",
                                paste0("The ", id," read depth file has been uploaded")))
        }
      },ignoreInit = T)
    }
  )
}


