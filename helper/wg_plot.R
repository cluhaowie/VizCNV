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
  if(mad(df$ratio,na.rm = TRUE)==0){
    res <- data.table(ID=id,chrom=unique(df$V1),
                      loc.start=min(df$V2),
                      loc.end=max(df$V3),
                      num.mark=nrow(df),
                      seg.mean=log2(median(df$ratio)+0.00001)
                      )
    return(res)
  }
  if (seg.method == "cbs") {
    print("segment with CBS")
    CNA.obj <-
      DNAcopy::CNA(
        log2(df$ratio + 0.00001),
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
    chr=df$V1[idx]
    start=df$V2[idx]
    end=c(df$V2[c(idx[-1],end(df$V2)[1])])
    res <- data.table(ID=id,chrom=chr,loc.start=start,loc.end=end,num.mark=res$lengths,seg.mean=res$values)
    print("done")
    return(res)
  }
  
  if (seg.method == "slm_wes") {
    print("segment with SLM for WES data")
    logratio <- log2(df$ratio + 0.001)
    logratio[is.na(logratio)] <- 0
    slm <-
      SLMSeg::HSLM(
        logratio,
        pos_data = (df$V2+df$V3)/2,
        omega = 0.7,
        FW = 0,
        eta = 1e-5,
        stepeta=1000 
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
    mean_ztpm <- sapply(seq_along(res$lengths),
                        function(i){
                          if(i==1){
                            return(mean(df$V4[1:res$lengths[i]]))
                          }
                          z <- df$V4[(1+sum(res$lengths[1:(i-1)])):(sum(res$lengths[1:i]))]
                          return(mean(z))
                        })
    res <- data.table(ID=id,chrom=chr,loc.start=start,loc.end=end,num.mark=res$lengths,seg.mean=res$values,ztpm.mean=mean_ztpm)
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



wg_seg2plot <- function(seg_data,ref=NULL,min.num.mark=100){
  del.log2width <- 0.1
  dup.log2width <- 0.15
  del.range <- c(log2(1*(1-del.log2width )/2),log2(1*(1+del.log2width )/2))
  dup.range <- c(log2(3*(1-dup.log2width )/2),log2(3*(1+dup.log2width )/2))
  
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
    mutate(start_cum = loc.start + loc_add)
  axis_set <- seg_data %>% 
    group_by(chrom) %>% 
    summarize(center = mean(end_cum),end=max(end_cum)) %>% 
    arrange((chrom))
  label_seg_gain <- seg_data %>% 
    filter(num.mark > min.num.mark) %>% 
    filter(seg.mean > dup.range[1])
  label_seg_loss <- seg_data %>% 
    filter(num.mark > min.num.mark) %>% 
    filter(dplyr::between(seg.mean,del.range[1],del.range[2]))
  wg <- seg_data %>% 
    ggplot(.,aes(x = end_cum, y = seg.mean, color = chrom))+
    geom_segment(data= subset(seg_data, num.mark > min.num.mark),aes(x = start_cum, y = seg.mean, xend = end_cum, yend = seg.mean), linewidth = 1.25)+
    geom_segment(data = label_seg_gain, aes(x = start_cum, y = seg.mean, xend = end_cum, yend = seg.mean),linewidth = 2,color = "red")+
   # geom_segment(data = label_seg_gain, aes(x = (start_cum+end_cum)/2, y = 1.6, xend = (start_cum+end_cum)/2, yend = 1.6+0.1),linewidth = 1,color = "purple")+
    geom_segment(data = label_seg_loss, aes(x = start_cum, y = seg.mean, xend = end_cum, yend = seg.mean),linewidth = 2,color = "green")+
   # geom_segment(data = label_seg_loss, aes(x = (start_cum+end_cum)/2, y = -2, xend = (start_cum+end_cum)/2, yend = -2+0.1),linewidth = 1,color = "orange")+
    geom_vline(xintercept=axis_set$end,lty=2,lwd=0.5,color="lightgrey")+
    theme_minimal() +
    theme( 
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.x = element_text(angle = 60, size = 12,vjust = 0.5,color = "black"),
      axis.text.y = element_text(size = 10,color = "black")
    )+
    scale_rd+
    scale_size_continuous(range = c(0.5,3))+
    labs(x = ref)+
    scale_x_continuous(label = axis_set$chrom, breaks = axis_set$center)+
    coord_cartesian(expand = F)
  return(wg)
}
