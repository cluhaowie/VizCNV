# process_sv <- function(sv){
#   sv <- sv %>% 
#     filter(AVGLEN > 10000 & AVGLEN < 100000000,
#            SVTYPE != "BND")
#   sv <- sv %>% 
#     mutate(color = case_when(SVTYPE == "DEL" ~ "darkblue",
#                              SVTYPE == "DUP" ~ "#8b0000",
#                              SVTYPE == "INS" ~ "darkgreen", 
#                              SVTYPE == "INV" ~ "darkorange", 
#                              TRUE ~ "white")) %>% 
#     mutate(idx = sample.int(size = n(), n = 980, replace = T)/10000)
#   sv <- sv %>% 
#     mutate(start = POS, 
#            end = as.integer(END)) %>% 
#     relocate(CHROM, start, end) %>% 
#     dplyr::select(-c(POS, END, REF, ALT, AVGLEN, MAPQ, RE, CIEND, CIPOS))
#   return(sv)
# }

## update the function to handle different field of SV length type
process_sv <- function(sv) {
  
  # Dynamically filter using the detected SV length column
  sv <- sv %>% 
    mutate(SVLEN=as.numeric(END-POS))%>%
    filter(SVLEN > 10000 & SVLEN < 100000000,
           SVTYPE != "BND")
  
  sv <- sv %>% 
    mutate(color = case_when(SVTYPE == "DEL" ~ "darkblue",
                             SVTYPE == "DUP" ~ "#8b0000",
                             SVTYPE == "INS" ~ "darkgreen", 
                             SVTYPE == "INV" ~ "darkorange", 
                             TRUE ~ "white")) %>% 
    mutate(idx = sample.int(size = n(), n = 980, replace = TRUE)/10000)
  
  sv <- sv %>% 
    mutate(start = POS, 
           end = as.integer(END)) %>% 
    relocate(CHROM, start, end) %>% 
    dplyr::select(-c(POS, END, REF, ALT))
  
  return(sv)
}