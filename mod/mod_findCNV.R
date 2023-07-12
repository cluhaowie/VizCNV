getSeg = function(df, idx){
  df = df %>% 
    filter(chr == idx)
  slm = SLMSeg::SLM(
    log2(df$ratio + 0.00001),
    omega = 0.3,
    FW = 0,
    eta = 0.00001
  )
  res <- rle(slm[1, ])
  idx <- sapply(seq_along(res$lengths),function(i){
    if(i==1){return(1)}
    start.idx=1+sum(res$lengths[1:(i-1)])
    return(start.idx)
  })
  chr=df$chr[idx]
  start=df$start[idx]
  end=c(df$start[c(idx[-1],end(df$start)[1])])
  res.dt <- data.table(chr=chr,loc.start=start,loc.end=end,num.mark=res$lengths,seg.mean=res$values)
  return (res.dt)  
}

getAllSeg <-function(df){
  out = list()
  for (c in paste0("chr",c(1:22, "X"))){
    tmp = getSeg(df,c)
    out = rbind(out, tmp)  
  }
  return(out)
}

normalization = function(rd){
  out <- rd %>% 
    group_by(chr) %>% 
    mutate(ratio=coverage/median(coverage+0.00001))
  return(out)
}

get_cnv_all <- function(pr_seg, mo_seg, fa_seg){
  pr_dup <- pr_seg %>% 
    filter(seg.mean > 0.4 & seg.mean < 3) %>% 
    filter(num.mark >=10) %>% 
    mutate(type = "gain")
  
  pr_del <- pr_seg %>% 
    filter(seg.mean < -0.8 & seg.mean > -5) %>% 
    filter(num.mark >=10) %>% 
    mutate(type = "loss")
  
  x <- rbind(pr_dup, pr_del)
  type <- x$type
  x <- x %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  
  pr_idx <- findOverlaps(x, pr_seg %>% makeGRangesFromDataFrame(keep.extra.columns = T), minoverlap = 10000, select = "last")
  pr_log <- pr_seg[pr_idx,]$seg.mean
  
  y <- mo_seg %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  hits <- findOverlaps(x, y)
  ol <- pintersect(x[queryHits(hits)], y[subjectHits(hits)])
  pol <- width(ol)/width(x[queryHits(hits)])
  t <- hits %>% 
    as.data.frame()
  t$pol <- pol
  
  slice <- dplyr::slice
  t <- t %>% 
    group_by(queryHits) %>% 
    slice(which.max(pol))
  idx <- t$subjectHits
  mo_log <- mo_seg$seg.mean[idx]
  
  y <- fa_seg %>% makeGRangesFromDataFrame(keep.extra.columns = T)
  hits <- findOverlaps(x, y)
  ol <- pintersect(x[queryHits(hits)], y[subjectHits(hits)])
  pol <- width(ol)/width(x[queryHits(hits)])
  t <- hits %>% 
    as.data.frame()
  t$pol <- pol
  t <- t %>% 
    group_by(queryHits) %>% 
    slice(which.max(pol))
  idx <- t$subjectHits
  fa_log <- fa_seg$seg.mean[idx]
  mcols(x) <- cbind(pr_log, mo_log, fa_log, type)
  x <- x %>% as.data.frame()
  return(x)
}

get_overlap <- function(df, SegDup_merge){
  
  # create empty lists to store the results
  # cnt_rs <- vector("integer", length = nrow(df))
  # str_rs <- vector("character", length = nrow(df))
  # target_rs <- vector("character", length = nrow(df))
  # sl_rs <- vector("character", length = nrow(df))
  # ll_rs <- vector("character", length = nrow(df))
  # cnt_omim <- vector("integer", length = nrow(df))
  # str_omim <- vector("character", length = nrow(df))
  cnt_sd <- vector("integer", length = nrow(df))
  
  # iterate over each row of the 'df' data frame
  for (i in seq_len(nrow(df))) {
    
    # target <- target_list %>% 
    #   filter(BH == df$FAMILY[i])
    # target <- unlist(strsplit(unlist(target$Candidates), split = ", "))
    # create a GRanges object for the genomic region of interest
    
    gr <- GRanges(df[i, 1], IRanges(df[i, 2], df[i, 3]))
    
    # subset reference genes that overlap with the genomic region
    # ol <- subsetByOverlaps(RefSeq_gr %>% filter(seqnames == df[i, 1]), gr)
    # genes <- unique(mcols(ol)$gene_id)
    # # store the results in the corresponding list elements
    # cnt_rs[i] <- length(genes)
    # str_rs[i] <- paste(genes, collapse = ", ")
    
    # # extract unique gene IDs
    # 
    # if (length(which(genes %in%target, arr.ind = T)) == 0){
    #   target_rs[i] = "NA"
    # } else {
    #   target_rs[i] <- paste0(target[which(target %in% genes)], collapse = ", ")
    # }
    # if (length(which(genes %in%short_list, arr.ind = T)) == 0){sl_rs[i] = "NA"} else {sl_rs[i] = paste0(short_list[which(short_list %in% genes)], collapse = ", ")}
    # if (length(which(genes %in%long_list, arr.ind = T)) == 0){ll_rs[i] = "NA"} else {ll_rs[i] = paste0(long_list[which(long_list %in% genes)], collapse = ", ")}
    
    # 
    # # subset reference genes that overlap with the genomic region
    # ol <- subsetByOverlaps(makeGRangesFromDataFrame(OMIM %>% filter(chrom == df[i, 1]), keep.extra.columns = T), gr)
    # genes <- unique(mcols(ol)$gene_symbol)
    # cnt_omim[i] <- length(genes)
    # str_omim[i] <- paste(genes, collapse = ", ")
    
    cnt_sd[i] <- countOverlaps(gr, makeGRangesFromDataFrame(SegDup_merge %>% filter(chrom == df[i, 1]), keep.extra.columns = T), type = "any", minoverlap = as.integer((df[i,3]-df[i,2])*0.98))
    
  }
  
  
  
  # combine the results with the 'df' data frame using mutate()
  out <- df %>%
    mutate(
      # refseqCount = cnt_rs, refseqID = str_rs,
      # OMIMCount = cnt_omim, OMIMID = str_omim,
      SDOverlap = cnt_sd)
  return(out)
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
#' mod_findCNV_UI("table")
#'
#' @export
mod_findCNV_UI <- function(id) {
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
#' mod_findCNV_Server("table", pr_seg, mo_seg, fa_seg)
#'
#' @export
mod_findCNV_Server <- function(id, pr_rd, mo_rd, fa_rd, SegDup_merge) {
  moduleServer(
    id,
    function(input, output, session) {
      
      names(pr_rd) <- c("chr", "start", "end", "coverage")
      names(mo_rd) <- c("chr", "start", "end", "coverage")
      names(fa_rd) <- c("chr", "start", "end", "coverage")
      
      pr_rd <- normalization(pr_rd)
      mo_rd <- normalization(mo_rd)
      fa_rd <- normalization(fa_rd)
      
      id <- showNotification("Finding CNVs", type = "message", duration = NULL)
      pr_seg <- getAllSeg(pr_rd %>% as.data.frame())
      mo_seg <- getAllSeg(mo_rd %>% as.data.frame())
      fa_seg <- getAllSeg(fa_rd %>% as.data.frame())
      cnv_all <- get_cnv_all(pr_seg, mo_seg, fa_seg) 
      df <- get_overlap(cnv_all, SegDup_merge)
      removeNotification(id = id)
      df <- df %>% 
        mutate(pr_log = as.numeric(pr_log),
               mo_log = as.numeric(mo_log),
               fa_log = as.numeric(fa_log),)
      df <- df %>% 
        mutate(pr_lvl = case_when(pr_log <= -1.525 ~ "HOM_DEL", 
                                  pr_log >= log2(1*0.9/2) & pr_log <= log2(1*1.1/2) ~ "HET_DEL",
                                  pr_log >= log2(2*0.9/2) & pr_log <= log2(2*1.1/2) ~ "NML",
                                  pr_log >= log2(3*0.9/2) & pr_log <= log2(3*1.1/2) ~ "DUP",
                                  pr_log >= log2(4*0.9/2) & pr_log <= log2(4*1.1/2) ~ "TRP",
                                  pr_log >= 1.175 ~ "MUL_GAIN",
                                  TRUE ~ "UND"),
               mo_lvl = case_when(mo_log <= -1.525 ~ "HOM_DEL", 
                                  mo_log >= log2(1*0.9/2) & mo_log <= log2(1*1.1/2) ~ "HET_DEL",
                                  mo_log >= log2(2*0.9/2) & mo_log <= log2(2*1.1/2) ~ "NML",
                                  mo_log >= log2(3*0.9/2) & mo_log <= log2(3*1.1/2) ~ "DUP",
                                  mo_log >= log2(4*0.9/2) & mo_log <= log2(4*1.1/2) ~ "TRP",
                                  mo_log >= 1.175 ~ "MUL_GAIN",
                                  TRUE ~ "UND"),
               fa_lvl = case_when(fa_log <= -1.525 ~ "HOM_DEL", 
                                  fa_log >= log2(1*0.9/2) & fa_log <= log2(1*1.1/2) ~ "HET_DEL",
                                  fa_log >= log2(2*0.9/2) & fa_log <= log2(2*1.1/2) ~ "NML",
                                  fa_log >= log2(3*0.9/2) & fa_log <= log2(3*1.1/2) ~ "DUP",
                                  fa_log >= log2(4*0.9/2) & fa_log <= log2(4*1.1/2) ~ "TRP",
                                  fa_log >= 1.175 ~ "MUL_GAIN",
                                  TRUE ~ "UND"))
      
      
      df <- df %>% 
        mutate(class = case_when(!pr_lvl %in% c("UND", "NML")  & mo_lvl == "NML" & fa_lvl == "NML" ~ "de_novo",
                                 pr_lvl != "UND" & pr_lvl ==  mo_lvl & fa_lvl == "NML" ~ "inh", #inh from M
                                 pr_lvl != "UND" & pr_lvl == fa_lvl & mo_lvl == "NML" ~ "inh", #ihn from F
                                 pr_lvl == "TRP" & mo_lvl %in% c("DUP", "TRP") & fa_lvl %in% c("DUP", "TRP")  ~ "inh",
                                 pr_lvl %in% c("HET_DEL", "HOM_DEL") & mo_lvl %in% c("HET_DEL", "HOM_DEL") & fa_lvl %in% c("HET_DEL", "HOM_DEL") ~ "inh",
                                 TRUE ~ "UND"))
      
      df <- df %>% 
        mutate(inh = case_when(class == "inh" & pr_lvl != "TRP" & pr_lvl == mo_lvl & mo_lvl != fa_lvl ~ "mo",
                               class == "inh" & pr_lvl != "TRP" & pr_lvl == fa_lvl & mo_lvl != fa_lvl~ "fa",
                               class == "inh" & pr_lvl == "TRP" & mo_lvl == "DUP" & mo_lvl == fa_lvl~ "both",
                               TRUE ~ "UND"))
      table <- reactive({
        df
      })
      output$table <- renderDataTable(table(),
                                      extensions=c("Responsive","Buttons"),
                                      server = F,
                                      editable = F,
                                      filter = list(position = 'top', clear = T),
                                      options = list(dom = 'Bfrtip',
                                                     buttons = c('copy','csv', 'excel',)))
      return(table() %>% as.data.frame())
    }
  )
}


