##for wg plots

#' Obtain segmented copy number calls for a single chromosome
#'
#' This function performs segmentation analysis on the normalized coverage data for a single chromosome and returns segmented copy number calls.
#' @param df A data frame containing the normalized coverage data to segment. It must have columns named "chr", "start", "end", and "ratio".
#' @param idx A character string specifying the chromosome to segment.
#' @return A data table with columns named "chr", "loc.start", "loc.end", "num.mark", and "seg.mean", containing the segmented copy number calls for the specified chromosome.

# getSeg <- function(df, idx){
#   df = df %>% 
#     filter(chr == idx)
#   slm = SLMSeg::SLM(
#     log2(df$ratio + 0.00001),
#     omega = 0.3,
#     FW = 0,
#     eta = 0.00001
#   )
#   res <- rle(slm[1, ])
#   idx <- sapply(seq_along(res$lengths),function(i){
#     if(i==1){return(1)}
#     start.idx=1+sum(res$lengths[1:(i-1)])
#     return(start.idx)
#   })
#   chr=df$chr[idx]
#   start=df$start[idx]
#   end=c(df$start[c(idx[-1],end(df$start)[1])])
#   res.dt <- data.table(chr=chr,loc.start=start,loc.end=end,num.mark=res$lengths,seg.mean=res$values)
#   return (res.dt)  
# }
## replace getSeg with SegNormRD function which handles chrY segmentation
SegNormRD <- function(df, id, seg.method = "cbs") {
  # make sure df only has one chr in V1
  if(mad(df$ratio)==0){
    res <- data.table(ID=id,chrom=unique(df$V1),
                      loc.start=min(df$V2),
                      loc.end=max(df$V3),
                      num.mark=nrow(df),
                      seg.mean=log2(median(df$ratio)+0.001)
                      )
    return(res)
  }
  if (seg.method == "cbs") {
    print("segment with CBS")
    CNA.obj <-
      DNAcopy::CNA(
        log2(df$ratio + 0.001),
        df$V1,
        df$V2,
        data.type = "logratio",
        sampleid = id
      )
    seg <- DNAcopy::segment(CNA.obj)
    seg.output <- seg$output
    seg.output$ID <- id
    print("done")
    return(seg.output)
  }
  ##EDIT: include SLM segmentation
  if (seg.method == "slm") {
    print("segment with SLM")
    slm <-
      SLMSeg::SLM(
        log2(df$ratio + 0.001),
        omega = 0.3,
        FW = 0,
        eta = 1e-5
      )
    res <- rle(slm[1, ])
    idx <- sapply(seq_along(res$lengths),function(i){
      if(i==1){return(1)}
      start.idx=1+sum(res$lengths[1:(i-1)])
      return(start.idx)
    })
    chr=df$V1[idx]
    start=df$V2[idx]
    end=c(df$V2[c(idx[-1],end(df$V2)[1])])
    res <- data.table(ID=id,chrom=chr,loc.start=start,loc.end=end,num.mark=res$lengths,seg.mean=res$values)
    print("done")
    return(res)
  }
  
}


#' Obtain segmented copy number calls for all chromosomes
#'
#' This function performs segmentation analysis on the normalized coverage data for all chromosomes and returns segmented copy number calls for each chromosome.
#' @param df A data frame containing the normalized coverage data to segment. It must have columns named "chr", "start", "end", and "ratio".
#' @return A data table with columns named "chr", "loc.start", "loc.end", "num.mark", and "seg.mean", containing the segmented copy number calls for each chromosome.
# getAllSeg <-function(df){
#   
#   out <- rbindlist(lapply(paste0("chr", c(1:22,"X","Y")),function(c){
#     df <- df %>% filter(V1 == c)
#     SegNormRD(df,id="NA", seg.method = "slm")
#   }))
#   return(out)
# }

#' Normalize coverage data using median ratio
#'
#' This function normalizes coverage data by dividing each coverage value by the median coverage value of either the whole genome or each chromosome.
#' @param df A data frame containing the coverage data to normalize. It must have columns named "chr", "start", "end", and "coverage".
#' @param norm_option A character string specifying the normalization option. "chr_med" (default) normalizes by chromosome median, and "wg_med" normalizes by whole genome median.
#' @return A data frame with an additional "ratio" column containing the normalized coverage values.
# wg_norm <- function(df, norm_option = "chr_med"){
#   names(df) <- c("chr", "start", "end", "coverage")
#   if (norm_option == "chr_med"){
#     tmp <- df%>%
#       group_by(chr)%>%
#       mutate(ratio=coverage/median(coverage+0.00001))
#   } else if (norm_option == "wg_med"){
#     tmp <- df%>%
#       mutate(ratio=coverage/median(coverage+0.00001))
#   }
#   return(tmp)
# }



wg_seg2plot <- function(seg_data){
  seg_data <- seg_data %>% 
    mutate(seg.mean=ifelse(seg.mean < -2.5,-2,seg.mean))
  
  temp <- seg_data %>% 
    group_by(chrom) %>% 
    summarise(max_end = max(loc.end)) %>% 
    mutate(across("chrom", str_replace, "chr", "")) %>% 
    arrange(as.numeric(chrom)) %>% 
    mutate(loc_add = lag(cumsum(as.numeric(max_end)), default = 0)) %>% 
    mutate(chrom = paste0("chr", chrom))
  
  seg_data <- seg_data %>% 
    inner_join(temp, by = "chrom") %>% 
    mutate(end_cum = loc_add + loc.end) 
  seg_data <- seg_data %>% 
    mutate(start_cum = end_cum- num.mark*1000)
  axis_set <- seg_data %>% 
    group_by(chrom) %>% 
    summarize(center = mean(end_cum)) %>% 
    arrange((chrom))
  label_seg_gain <- seg_data %>% 
    filter(num.mark > 100) %>% 
    filter(seg.mean >0.4)
  label_seg_loss <- seg_data %>% 
    filter(num.mark > 100) %>% 
    filter(dplyr::between(seg.mean,-1.5, -0.3))
  wg <- seg_data %>% 
    ggplot(aes(x = end_cum, y = seg.mean, color = chrom))+
    geom_segment(aes(x = start_cum, y = seg.mean, xend = end_cum, yend = seg.mean+0.001), linewidth = 1.25)+
    geom_point(data = label_seg_gain, shape= 8, color = "red")+
    geom_point(data = label_seg_loss, shape= 8, color = "green")+
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
    scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center)+
    coord_cartesian(expand = F)
  return(wg)
}