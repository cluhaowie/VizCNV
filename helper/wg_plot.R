##for wg plots
getSeg <- function(df, idx){
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
  for (c in paste0("chr", c(1:22,"X"))){
    tmp = getSeg(df,c)
    out = rbind(out, tmp)  
  }
  return(out)
}


wg_norm <- function(df, norm_option = "chr_med"){
  names(df) <- c("chr", "start", "end", "coverage")
  if (norm_option == "chr_med"){
    tmp <- df%>%
      group_by(chr)%>%
      mutate(ratio=coverage/median(coverage+0.00001))
  } else if (norm_option == "wg_med"){
    tmp <- df%>%
      mutate(ratio=coverage/median(coverage+0.00001))
  }
  return(tmp)
}



wg_seg2plot <- function(seg_data){
  seg_data <- seg_data %>% 
    mutate(seg.mean=ifelse(seg.mean < -2.5,-2,seg.mean))
  
  temp <- seg_data %>% 
    group_by(chr) %>% 
    summarise(max_end = max(loc.end)) %>% 
    mutate(across("chr", str_replace, "chr", "")) %>% 
    arrange(as.numeric(chr)) %>% 
    mutate(loc_add = lag(cumsum(as.numeric(max_end)), default = 0)) %>% 
    mutate(chr = paste0("chr", chr))
  
  seg_data <- seg_data %>% 
    inner_join(temp, by = "chr") %>% 
    mutate(end_cum = loc_add + loc.end) 
  seg_data <- seg_data %>% 
    mutate(start_cum = end_cum- num.mark*1000)
  axis_set <- seg_data %>% 
    group_by(chr) %>% 
    summarize(center = mean(end_cum)) %>% 
    arrange((chr))
  label_seg_gain <- seg_data %>% 
    filter(num.mark > 100) %>% 
    filter(seg.mean >0.4)
  label_seg_loss <- seg_data %>% 
    filter(num.mark > 100) %>% 
    filter(dplyr::between(seg.mean,-1.5, -0.3))
  wg <- seg_data %>% 
    ggplot(aes(x = end_cum, y = seg.mean, color = chr))+
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
    scale_x_continuous(label = axis_set$chr, breaks = axis_set$center)+
    coord_cartesian(expand = F)
  return(wg)
}