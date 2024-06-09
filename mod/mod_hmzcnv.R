#' find homozygous cnv
#' @description There are two methods to identify homologous CNVs.
#' The first one (default) is to intersect proband segments
#' and Baf segments; the second approach is to
#' intersect proband seg with both parents segments.
#'
#' @importFrom GenomicRanges pintersect findOverlaps
#' @importFrom S4Vectors queryHits subjectHits mcols
#' @importFrom IRanges width
#' @param pr Grange object of proband read depth segments
#' @param pr_baf Grange object proband baf segments
#' @param mo Grange object mother read depth segments
#' @param fa Grange object father read depth segments
#'
#' @return GRanges object of proband seg with overlapp percentage annotate
#'
#' @examples
#' get_hmzCNV(pr, pr_baf);
#' get_hmzCNV(pr, mo, fa)

.get_hmzcnv <-  function(pr, mo, fa){
    # intersect mo and fa segment
    temp <- GenomicRanges::intersect(mo, fa, ignore.strand=T )
    overlapping_indices <- findOverlaps(pr,temp,ignore.strand=T)
    overlap_size <- GenomicRanges::pintersect(pr[queryHits(overlapping_indices)], temp[subjectHits(overlapping_indices)],ignore.strand=T)
    overlap_percentage <- round(width(overlap_size) / width(pr[queryHits(overlapping_indices)]) * 100,1)
    
    hmz = pr[queryHits(overlapping_indices)]
    mcols(hmz)$overlap_length = width(overlap_size)
    mcols(hmz)$overlap_percentage <- overlap_percentage
    hmz <- hmz[mcols(hmz)$overlap_percentage>50&mcols(hmz)$num.mark>1]
  return(hmz)
}

.get_hmz_proband <- function(pr, pr_baf){
  # Find overlapping intervals
  overlapping_indices <- findOverlaps(pr, pr_baf,ignore.strand=T)
  # Calculate the percentage of overlap compare to pr
  overlap_size <- GenomicRanges::pintersect(pr[queryHits(overlapping_indices)], pr_baf[subjectHits(overlapping_indices)],ignore.strand=T)
  overlap_percentage <- round(width(overlap_size) / width(pr[queryHits(overlapping_indices)]) * 100,1)
  # Annotate the percentage of overlap 
  hmz = pr[queryHits(overlapping_indices)]
  mcols(hmz)$overlap_percentage <- overlap_percentage
  hmz <- hmz[mcols(hmz)$overlap_percentage>50&mcols(hmz)$num.mark>1]
  return(hmz)
}

#' find all hmz CNV 
#' 
#' 
#' @description preprocess the segments by separating them into gain and loss 
#' based on previously normalized log ratio labeling cnv type and combining 
#' them into one table
#' 
#' @param pr_seg proband read depth segments
#' @param mo_seg mother read depth segments
#' @param fa_seg father read depth segments
#'
#' @return GRanges object of proband seg with overlapp percentage annotate
#'
#' @examples
#' get_hmzCNV_all(pr_seg, mo_seg, fa_seg)
#'
#' @export

get_hmzCNV_all <- function(pr_seg, pr_baf=NULL,mo_seg=NULL,fa_seg=NULL){
  if(is.null(pr_baf)&is.null(mo_seg)&is.null(fa_seg)) {stop("require at least one method")}
  if(is.null(mo_seg)|is.null(fa_seg)){
    pr_dup <- pr_seg %>% filter(seg.mean > 0.4) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    pr_del <- pr_seg %>% filter(seg.mean < -0.5) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    hmzCNV_gain <-  .get_hmz_proband(pr_dup,pr_baf)
    hmzCNV_loss <-  .get_hmz_proband(pr_del,pr_baf)
    mcols(hmzCNV_gain)$type = rep("gain", length(hmzCNV_gain))
    mcols(hmzCNV_loss)$type = rep("loss", length(hmzCNV_loss))
    hmzCNV_all <-  append(hmzCNV_gain, hmzCNV_loss)%>%base::as.data.frame()
  }else if (!is.null(mo_seg)&!is.null(fa_seg)){
    ## separating segments into two groups
    pr_dup <- pr_seg %>% filter(seg.mean > 0.4) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    mo_dup <- mo_seg %>% filter(seg.mean > 0.4) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    fa_dup <- fa_seg %>% filter(seg.mean > 0.4) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    pr_del <- pr_seg %>% filter(seg.mean < -0.5) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    mo_del <- mo_seg %>% filter(seg.mean < -0.5) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    fa_del <- fa_seg %>% filter(seg.mean < -0.5) %>% makeGRangesFromDataFrame(keep.extra.columns = T)
    hmzCNV_gain <-  .get_hmzcnv(pr=pr_dup, mo=mo_dup, fa=fa_dup)
    hmzCNV_loss <-  .get_hmzcnv(pr_del, mo=mo_del, fa=fa_del)
    mcols(hmzCNV_gain)$type = rep("gain", length(hmzCNV_gain))
    mcols(hmzCNV_loss)$type = rep("loss", length(hmzCNV_loss))
    hmzCNV_all <-  append(hmzCNV_gain, hmzCNV_loss)%>%base::as.data.frame()
  }
  return(hmzCNV_all)
}


#' find hmz CNVs UI
#'
#' rendering hmz CNV DT
#' 
#' @return reactive df table 
#'
#' @examples
#' mod_hmzcnv_UI("table")
#'
#' @export
mod_hmzcnv_UI <- function(id) {
  ns <- NS(id)
  tagList(
    dataTableOutput(ns("table"))
  )
}

#' find hmz CNVs Server
#'
#' preprocess tables for rendering
#'
#' @param pr_seg proband read depth segments from segmentation algo
#' @param pr_baf proband Baf segments processed with file
#' @param mo_seg mother read depth segments from segmentation algo
#' @param fa_seg father read depth segments from segmentation algo
#'
#' @return reactive df table 
#'
#' @examples
#' mod_hmzcnv_Server ("table", pr_seg, pr_baf, mo_seg, fa_seg)
#'
#' @export
mod_hmzcnv_Server <- function(id, pr_seg, pr_baf, mo_seg, fa_seg) {
  moduleServer(
    id,
    function(input, output, session) {
      if(nrow(pr_baf)==0){
        table <- reactive({
          get_hmzCNV_all(pr_seg, mo_seg = mo_seg, fa_seg = fa_seg) 
        })
      } else{
        table <- reactive({
          get_hmzCNV_all(pr_seg, pr_baf)
          })
      }

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


