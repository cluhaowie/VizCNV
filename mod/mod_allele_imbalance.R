#' Prune SNVs
#' 
#' @param gr snv genomic range object, return from readGVCF
#' @param target_spacing density options target_spacing=0, no filtering
#' 
#' @return pruned snv


prune_snvs <- function(gr, target_spacing=0) {
  # Calculate the distance to the next SNV
  distances <- c(diff(start(gr)), NA) # NA for the last SNV
  
  # Initialize a vector to keep track of which SNVs to keep
  keep <- rep(TRUE, length(gr))
  
  # Loop through the SNVs, removing those too close to a previous one
  for (i in seq_along(distances)) {
    if (!is.na(distances[i]) && distances[i] < target_spacing) {
      keep[i + 1] <- FALSE
    }
  }
  
  # Return the pruned set of SNVs
  return(gr[keep])
}

# ReadGVCF <- function(path_to_gVCF,ref_genome=ref_genome,param = param,target_spacing){
#   print("Grabbing regions")
#   expBAF=0.5 # expected BAF for diploid genome
#   vcf<- VariantAnnotation::readVcf(file = path_to_gVCF,genome = ref_genome,param = param)
#   vcf.gr <- vcf@rowRanges
#   GT <- VariantAnnotation::geno(vcf)$GT
#   AD <- VariantAnnotation::geno(vcf)$AD
#   DP <- VariantAnnotation::geno(vcf)$DP
#   PR_ID=colnames(GT)[1]
#   P1_ID=colnames(GT)[2]
#   P2_ID=colnames(GT)[3]
#   G1=c('0/0',"0|0")
#   G2=c('1/1',"1|1")
#   G3=c('0/1',"0|1")
#   GT <- as.data.table(GT)
#   setnames(GT,colnames(GT),c("index","P1","P2"))
#   GT.anno <- GT %>% 
#     mutate(B_InhFrom=case_when(index %in% G3 & P1 %in% G1 & P2  %in% c(G2,G3) ~ P2_ID, #cases 11,12
#                                index %in% G3 & P1 %in% c(G2,G3) & P2  %in% G1 ~ P1_ID, #case 13,16
#                                index %in% G2 & P1 %in% G1 & P2  %in% c(G2,G3) ~ P2_ID, #cases 20,21
#                                index %in% G2 & P1 %in% c(G2,G3) & P2  %in% G1 ~ P1_ID, #case 22,25
#                                TRUE ~ "Notphased")) %>% 
#     mutate(A_InhFrom=case_when(index %in% G1 & P1 %in% c(G1,G3) & P2  %in% G2 ~ P1_ID, #cases 3,6
#                                index %in% G1 & P1 %in% G2 & P2  %in% c(G1,G3) ~ P2_ID, #cases 7,8
#                                index %in% G3 & P1 %in% c(G1,G3) & P2  %in% G2 ~ P1_ID, #cases 12,15
#                                index %in% G3 & P1 %in% G2 & P2  %in% c(G1,G3) ~ P2_ID, #cases 16,17
#                                TRUE ~ "Notphased")) %>% 
#     mutate(B_col = case_when(B_InhFrom == P1_ID ~ "#E69F00",
#                              B_InhFrom == P2_ID ~ "#39918C",
#                              TRUE ~ "#999999")) %>% 
#     mutate(A_col = case_when(A_InhFrom == P1_ID ~ "#E69F00",
#                              A_InhFrom == P2_ID ~ "#39918C",
#                              TRUE ~ "#999999"))
#   
#   AD <- as.data.table(AD)
#   setnames(AD,colnames(AD),c("index","P1","P2"))
#   AD.anno <- AD%>%
#     mutate(index_ale_count=stringr::str_count(as.character(index),",|:"),
#            p1_ale_count=stringr::str_count(as.character(P1),",|:"),
#            p2_ale_count=stringr::str_count(as.character(P2),",|:"))%>%
#     mutate(index_REF_RD=sapply(index,"[[",1),
#            index_ALT_RD=sapply(index,"[[",2),
#            p1_REF_RD=sapply(P1,"[[",1),
#            p1_ALT_RD=sapply(P1,"[[",2),
#            p2_REF_RD=sapply(P2,"[[",1),
#            p2_ALT_RD=sapply(P2,"[[",2),
#            pr_count=index_ALT_RD+index_REF_RD,
#            p1_count=p1_REF_RD+p1_ALT_RD,
#            p2_count=p2_REF_RD+p2_ALT_RD,
#            pr_ALT_Freq=index_ALT_RD/(index_ALT_RD+index_REF_RD))%>%
#     mutate(pr_absBAF=abs(pr_ALT_Freq-expBAF))%>%
#     mutate(likelyDN=ifelse(p1_ALT_RD<2&p2_ALT_RD<2&index_ALT_RD>5&p1_REF_RD>10&p2_REF_RD>10&pr_count>10&pr_ALT_Freq>0.2,"TRUE","FALSE"))
#   
#   
#   AD.anno <- AD.anno[,c("pr_count","p1_count","p2_count","pr_ALT_Freq","pr_absBAF","likelyDN")]
#   merged <- cbind(GT.anno ,AD.anno)%>%
#     mutate(P1_phased_BAF=case_when(B_InhFrom==P1_ID ~ pr_ALT_Freq,
#                                    A_InhFrom==P2_ID ~pr_ALT_Freq),
#            P2_phased_BAF=case_when(B_InhFrom==P2_ID ~ pr_ALT_Freq,
#                                    A_InhFrom==P1_ID ~ pr_ALT_Freq))
#   
#   setnames(merged,c("index","P1","P2"),c(PR_ID,P1_ID,P2_ID))
#   mcols(vcf.gr) <- merged
#   vcf.gr <- prune_snvs(vcf.gr,target_spacing=target_spacing)
#   return(vcf.gr)
# }

ReadGVCF <- function(path_to_gVCF,ref_genome=ref_genome,param = param,target_spacing){
  # Safe read of VCF file with error handling
  tryCatch({
    vcf <- suppressWarnings({
      VariantAnnotation::readVcf(file = path_to_gVCF, genome = ref_genome, param = param)
    })
  }, error = function(e) {
    stop("Failed to read VCF file: ", e$message)
  })
  expBAF=0.5 # expected BAF for diploid genome
  
  vcf.gr <- vcf@rowRanges
  GT <- VariantAnnotation::geno(vcf)$GT
  AD <- VariantAnnotation::geno(vcf)$AD
  DP <- VariantAnnotation::geno(vcf)$DP
  PR_ID=colnames(GT)[1]
  P1_ID=colnames(GT)[2]
  P2_ID=colnames(GT)[3]
  G1=c('0/0',"0|0")
  G2=c('1/1',"1|1")
  G3=c('0/1',"0|1")
  GT <- as.data.table(GT)
  setnames(GT,colnames(GT),c("index","P1","P2"))
  GT.anno <- GT %>% 
    mutate(B_Index=case_when(index %in% G3 & P1 %in% G1 & P2 %in% G2 ~ 1, 
                               index %in% G3 & P1 %in% G1 & P2 %in% G3 ~ 2, 
                               index %in% G3 & P1 %in% G3 & P2 %in% G2 ~ 3,
                               index %in% G2 & P1 %in% G1 & P2 %in% G2 ~ 4,
                               index %in% G2 & P1 %in% G1 & P2 %in% G3 ~ 5,
                               index %in% G3 & P1 %in% G2 & P2 %in% G1 ~ 6,
                               index %in% G3 & P1 %in% G3 & P2 %in% G1 ~ 7,
                               index %in% G3 & P1 %in% G2 & P2 %in% G3 ~ 8,
                               index %in% G2 & P1 %in% G2 & P2 %in% G1 ~ 9,
                               index %in% G2 & P1 %in% G3 & P2 %in% G1 ~ 10,
                               index %in% G3 & P1 %in% G3 & P2 %in% G3 ~ 11,
                               index %in% G2 & P1 %in% G3 & P2 %in% G3 ~ 12,
                               index %in% G2 & P1 %in% G2 & P2 %in% G2 ~ 13,
                               index %in% G2 & P1 %in% G3 & P2 %in% G2 ~ 14,
                               index %in% G2 & P1 %in% G2 & P2 %in% G3 ~ 15,
                               index %in% G1 & P1 %in% G1 & P2 %in% G3 ~ 16,
                               index %in% G1 & P1 %in% G3 & P2 %in% G1 ~ 17,
                               index %in% G1 & P1 %in% G3 & P2 %in% G3 ~ 18,
                               index %in% G3 & P1 %in% G1 & P2 %in% G1 ~ 19,
                               index %in% G3 & P1 %in% G2 & P2 %in% G2 ~ 20,
                               TRUE ~ 21)) %>% 
    mutate(B_col = case_when(B_Index %in% c(6:10) ~ "#E69F00",
                             B_Index %in% c(1:5) ~ "#39918C",
                             B_Index %in% c(19:20) ~ "red",
                             B_Index %in% c(12:18) ~ "black",
                             TRUE ~ "#999999"))%>%
    mutate(B_InhFrom=case_when(B_Index %in% c(6:10) ~ P1_ID,
                               B_Index %in% c(1:5) ~ P2_ID,
                               B_Index %in% c(19:20) ~ "ME",
                               B_Index %in% c(12:18) ~ "AOH_signal",
                               TRUE ~ "Notphased"))
  
  AD <- as.data.table(AD)
  setnames(AD,colnames(AD),c("index","P1","P2"))
  AD.anno <- AD%>%
    mutate(index_ale_count=stringr::str_count(as.character(index),",|:"),
           p1_ale_count=stringr::str_count(as.character(P1),",|:"),
           p2_ale_count=stringr::str_count(as.character(P2),",|:"))%>%
    mutate(index_REF_RD=sapply(index,"[[",1),
           index_ALT_RD=sapply(index,"[[",2),
           p1_REF_RD=sapply(P1,"[[",1),
           p1_ALT_RD=sapply(P1,"[[",2),
           p2_REF_RD=sapply(P2,"[[",1),
           p2_ALT_RD=sapply(P2,"[[",2),
           pr_count=index_ALT_RD+index_REF_RD,
           p1_count=p1_REF_RD+p1_ALT_RD,
           p2_count=p2_REF_RD+p2_ALT_RD,
           pr_ALT_Freq=index_ALT_RD/(index_ALT_RD+index_REF_RD))%>%
    mutate(pr_absBAF=abs(pr_ALT_Freq-expBAF))%>%
    mutate(likelyDN=ifelse(p1_ALT_RD<2&p2_ALT_RD<2&index_ALT_RD>5&p1_REF_RD>10&p2_REF_RD>10&pr_count>10&pr_ALT_Freq>0.2,"TRUE","FALSE"))
  
  
  AD.anno <- AD.anno[,c("pr_count","p1_count","p2_count","pr_ALT_Freq","pr_absBAF","likelyDN")]
  merged <- cbind(GT.anno ,AD.anno)%>%
    mutate(P1_phased_BAF=case_when(B_Index %in% c(6:10) ~ pr_ALT_Freq),
           P2_phased_BAF=case_when(B_Index %in% c(1:5) ~ pr_ALT_Freq))
  
  setnames(merged,c("index","P1","P2"),c(PR_ID,P1_ID,P2_ID))
  mcols(vcf.gr) <- merged
  vcf.gr <- prune_snvs(vcf.gr,target_spacing=target_spacing)
  return(vcf.gr)
}
