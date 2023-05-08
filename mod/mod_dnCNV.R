
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
get_dnCNV = function(pr, mo, fa, size_threshold = 10000){
  
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
  min_seg_dup = 0.4
  max_seg_del = -0.5
  pr_dup_nrow <- pr_seg %>% filter(seg.mean > min_seg_dup)%>%nrow()
  pr_del_nrow <- pr_seg %>% filter(seg.mean < max_seg_del)%>%nrow()
  if(pr_dup_nrow==0&pr_del_nrow==0){
    print("No CNV calls on proband")
    return(NULL)
  }
  if(pr_dup_nrow==0&pr_del_nrow!=0){
    print("No dup calls on proband")
    pr_del <- pr_seg %>% filter(seg.mean < max_seg_del) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    mo_seg_nrow <- mo_seg %>% filter(seg.mean < max_seg_del)%>%nrow()
    fa_seg_nrow <- fa_seg %>% filter(seg.mean < max_seg_del)%>%nrow()
    if(mo_seg_nrow!=0){
      mo_del <- mo_seg %>% filter(seg.mean < max_seg_del) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    }else{
      mo_del <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
    }
    if(fa_seg_nrow!=0){
      fa_del <- fa_seg %>% filter(seg.mean < max_seg_del) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    }else{
      fa_del <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
    }
    dnCNV_loss = get_dnCNV(pr_del, mo_del, fa_del)
    mcols(dnCNV_loss)$type = rep("loss", length(dnCNV_loss))
    return(as.data.frame(dnCNV_loss))
  }
  if(pr_dup_nrow!=0&pr_del_nrow==0){
    print("No del calls on proband")
    pr_dup <- pr_seg %>% filter(seg.mean > min_seg_dup) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    mo_seg_nrow <- mo_seg %>% filter(seg.mean > min_seg_dup)%>%nrow()
    fa_seg_nrow <- fa_seg %>% filter(seg.mean > min_seg_dup)%>%nrow()
    if(mo_seg_nrow!=0){
      mo_dup <- mo_seg %>% filter(seg.mean > min_seg_dup) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    }else{
      mo_dup <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
    }
    if(fa_seg_nrow!=0){
      fa_dup <- fa_seg %>% filter(seg.mean > min_seg_dup) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    }else{
      mo_dup <- GRanges(c(seqnames=NULL,ranges=NULL,strand=NULL))
    }
    dnCNV_gain = get_dnCNV(pr_dup, mo_dup, fa_dup)
    mcols(dnCNV_gain)$type = rep("gain", length(dnCNV_gain))
    return(as.data.frame(dnCNV_gain))
  }
  if(pr_dup_nrow!=0&pr_del_nrow!=0){
    pr_dup <- pr_seg %>% filter(seg.mean > min_seg_dup) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    mo_dup <- mo_seg %>% filter(seg.mean > min_seg_dup) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    fa_dup <- fa_seg %>% filter(seg.mean > min_seg_dup) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    pr_del <- pr_seg %>% filter(seg.mean < max_seg_del) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    mo_del <- mo_seg %>% filter(seg.mean < max_seg_del) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    fa_del <- fa_seg %>% filter(seg.mean < max_seg_del) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    dnCNV_gain = get_dnCNV(pr_dup, mo_dup, fa_dup)
    dnCNV_loss = get_dnCNV(pr_del, mo_del, fa_del)
    mcols(dnCNV_gain)$type = rep("gain", length(dnCNV_gain))
    mcols(dnCNV_loss)$type = rep("loss", length(dnCNV_loss))
    dnCNV_all = append(dnCNV_gain, dnCNV_loss) %>% 
      as.data.frame()
    return(dnCNV_all)
  }
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
      output$table <- renderDataTable(table(),
                                      extensions=c("Responsive","Buttons"),
                                      server = T,
                                      editable = F,
                                      filter = list(position = 'top', clear = T),
                                      options = list(dom = 'Bfrtip',
                                                     buttons = c('txt','csv', 'excel')))
      return(table() %>% as.data.frame())
    }
  )
}


