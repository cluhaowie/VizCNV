###WORK IN PROGRESS###
###UNLIKELY TO IMPLEMENT THIS##



mod_snp_upload_UI <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      column(width = 6,
             fileInput(ns("file"),label = NULL, accept=c("*.vcf","*.vcf.gz"),multiple = F, buttonLabel = "Browse...")
      ),
      column(width=1,h5("or")),
      column(width = 4,
             shinyFilesButton(id = ns("local_snp_file"), label = "Browse...", title = "Please select a file", multiple = F, viewtype = "detail")
      )
    )
  )
}

mod_snp_upload_Server <- function(id,volumes,values) {
  moduleServer(
    id,
    function(input, output, session) {
      shinyFileChoose(input, "local_snp_file", roots = volumes, session=session)
      observeEvent(input$file, {
        Rsamtools::indexTabix(input$file$datapath,format = "vcf")
        values[["snp_gvcf_file_ref"]] <- VariantAnnotation::scanVcfHeader(input$file$datapath)@reference
        values[["snp_gvcf_file_path"]] <- input$file$datapath
        showModal(modalDialog(title = "File upload",
                              "The joint called SNP file has been uploaded and indexed"))
      })
      observeEvent(input$local_snp_file,{
        if(is.integer(input$local_snp_file)){
          cat("no local snp file found \n")
        }else{
          local_snp_file <- parseFilePaths(volumes, input$local_snp_file)
          index.file <- paste0(local_snp_file$datapath,".tbi")
          values[["snp_gvcf_file_ref"]] <- VariantAnnotation::scanVcfHeader(local_snp_file$datapath)@reference
          values[["snp_gvcf_file_path"]] <-  local_snp_file$datapath
          if(!file.exists(index.file)){
            Rsamtools::indexTabix(local_snp_file$datapath,format = "vcf")
          }
          showModal(modalDialog(title = "File upload",
                                paste0("The joint called SNP file has been indexed")))
        }
      },ignoreInit = T)
    }
  )
}


#' shiny file input wrapper UI
#'
#' @return tagList of Shiny UI elements with either upload file or able 
#' to select for local path to file
#'
#' @examples
#' mod_sv_upload_UI("pr_sv")
#'
#' @export

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

mod_sv_upload_Server <- function(id,volumes,values) {
  moduleServer(
    id,
    function(input, output, session) {
      shinyFileChoose(input, "local_sv_file", roots = volumes, session=session)
      observeEvent(input$file, {
        values[[id]] <- bedr::read.vcf(input$file$datapath,split.info = T)$vcf
        showModal(modalDialog(title = "File upload",
                              "The SV file has been uploaded"))
      })
      observeEvent(input$local_sv_file,{
        if(is.integer(input$local_sv_file)){
          cat("no file\n")
        }else{
          local_sv_file <- parseFilePaths(volumes, input$local_sv_file)
          values[[id]]<- bedr::read.vcf(local_sv_file$datapath,split.info = T)$vcf
          showModal(modalDialog(title = "File upload",
                                paste0("The SV file has been uploaded")))
        }
      },ignoreInit = T)
    }
  )
}

#' shiny file input wrapper UI
#'
#' @return tagList of Shiny UI elements with either upload file or able 
#' to select for local path to file
#'
#' @examples
#' mod_rd_upload_UI("pr_rd")
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
        showModal(modalDialog(title = "File upload",
                              "The read depth file has been uploaded"))
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

#WIP
dyncheckboxUI <- function(id){
  ns <- NS(id)
  label <- switch(id,"pr_rd"="Proband","m_rd"="Mother","f_rd"="Father")
  shinyjs::useShinyjs()
  shinyjs::hidden(
    id = ns("chk"),
    tagList(
      checkboxInput(ns("id"),label = label)
    )
  )

}

dyncheckboxServer <- function(id,values){
  moduleServer(
    id,
    function(input, output, session){
      if(nrow(values[[id]])==0){
        shinyjs::hide(id = "chk")
      } else {
        shinyjs::show(id = "chk")
      }
    }
  )
}
