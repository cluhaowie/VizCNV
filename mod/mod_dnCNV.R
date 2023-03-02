
#' find denovo CNV function
#'
#' subtracting mother and father segments from proband segmends 
#'
#' 
#' @param pr proband read depth segments
#' @param mo mother read depth segments
#' @param fa father read depth segments
#'
#' @return GRanges object of remaining segments 
#'
#' @examples
#' get_dnCNV(pr, mo, fa)
#'
#' @export
get_dnCNV = function(pr, mo, fa, size_threshold = 150000){
  
  # Essentially bedtools subtract
  tmp = unlist(subtract(pr, mo))
  dnCNV = unlist(subtract(tmp,fa))
  mcols(dnCNV)$length = width(dnCNV)
  
  # Filtering for size above threshold (in kb)
  dnCNV_filtered = dnCNV %>% 
    filter(length > size_threshold)
  
  return (dnCNV_filtered)
  
}


#' find all denovo CNV 
#' 
#' 
#' preprocess the segments by separating them into gain and loss based on previously normalized log ratio
#' 
#' subtracting mother and father segments from proband segmends 
#'
#' labeling dnCNV type and combining them into one table
#' 
#' @param pr_seg proband read depth segments
#' @param mo_seg mother read depth segments
#' @param fa_seg father read depth segments
#'
#' @return data.frame of all dnCNV segments 
#'
#' @examples
#' get_dnCNV(pr_seg, mo_seg, fa_seg)
#'
#' @export

get_dnCNV_all = function(pr_seg, mo_seg, fa_seg){
  
  ## separating segments into two groups
  pr_dup <- pr_seg %>% filter(seg.mean > 0.3) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  mo_dup <- mo_seg %>% filter(seg.mean > 0.3) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  fa_dup <- fa_seg %>% filter(seg.mean > 0.3) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  pr_del <- pr_seg %>% filter(seg.mean < -0.2 & seg.mean > -5) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  mo_del <- mo_seg %>% filter(seg.mean < -0.2 & seg.mean > -5) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  fa_del <- fa_seg %>% filter(seg.mean < -0.2 & seg.mean > -5) %>% makeGRangesFromDataFrame(keep.extra.columns = T)

  
  dnCNV_gain = get_dnCNV(pr_dup, mo_dup, fa_dup)
  dnCNV_loss = get_dnCNV(pr_del, mo_del, fa_del)
  mcols(dnCNV_gain)$type = rep("gain", length(dnCNV_gain))
  mcols(dnCNV_loss)$type = rep("loss", length(dnCNV_loss))
  
  dnCNV_all = append(dnCNV_gain, dnCNV_loss) %>% 
    as.data.frame()
  
  return(dnCNV_all)
}


#' find denovo CNV UI
#'
#' rendering DT
#'
#'  
#'
#' @return reactive df table and 
#'
#' @examples
#' mod_dnCNV_UI("table")
#'
#' @export
mod_dnCNV_UI <- function(id) {
  ns <- NS(id)
  tagList(
    dataTableOutput(ns("table"))
  )
}


#' find dnCNV Server
#'
#' preprocess tables for rendering
#'
#' @param pr_seg proband read depth segments from segmentation algo
#' @param mo_seg mother read depth segments from segmentation algo
#' @param fa_seg father read depth segments from segmentation algo
#'
#' @return reactive df table 
#'
#' @examples
#' mod_dnCNV_Server("table", pr_seg, mo_seg, fa_seg)
#'
#' @export
mod_dnCNV_Server <- function(id, pr_seg, mo_seg, fa_seg) {
  moduleServer(
    id,
    function(input, output, session) {
      table <- reactive({
        get_dnCNV_all(pr_seg, mo_seg, fa_seg) 
      })
      output$table <- renderDataTable(table())
    }
  )
}


