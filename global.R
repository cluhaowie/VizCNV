#  ------------------------------------------------------------------------
#
# Title : App - VizCNV
#    By : Haowei Du
#  Date : April 2022
#    
#  ------------------------------------------------------------------------
options(timeout = 6000)
options(shiny.maxRequestSize=2*1024^3) ## max file size 1 Gb
options(shiny.autoreload=TRUE)
options(shiny.reactlog=TRUE) 
library(BiocManager)
options(repos = BiocManager::repositories())

# Packages ----------------------------------------------------------------
# Install missing packages from CRAN
list.of.packages <- c("dplyr", "data.table", "knitr", "testthat", "shiny", "shinydashboard",
                      "tippy","DT","ggplot2","RSQLite","shinyWidgets","shinyFiles","waiter")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Install missing packages from Bioconductor
biocLitePackages <- c("DNAcopy", "GenomicRanges", "VariantAnnotation","BSgenome.Hsapiens.UCSC.hg38","BSgenome.Hsapiens.UCSC.hg19","bedr") 
new.biocLitePackage <- biocLitePackages[!(biocLitePackages %in% installed.packages()[,"Package"])]
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(length(new.biocLitePackage)) BiocManager::install(new.biocLitePackage)
# Install SLMSeg package from local
local.packages <- c("SLMSeg")
add.packages <- local.packages[!local.packages%in% installed.packages()[,"Package"]]
if(length(add.packages)) install.packages("SLMSeg_1.0.tar.gz",type = " source ", repos = NULL)
#Cleaning ----
detach_all <- function() {
  basic.pkg <- c("package:stats", "package:graphics", "package:grDevices",
                 "package:utils", "package:datasets", "package:methods", "package:base")
  
  pkg.list <- search()[ifelse(unlist(gregexpr("package:", search())) == 1 ,TRUE, FALSE)]
  
  pkg.list <- setdiff(pkg.list, basic.pkg)
  
  lapply(pkg.list, detach, character.only = TRUE)
}
detach_all()
rm(list = ls())

# Loading ----

library(dplyr)
library(data.table)
library(knitr)
library(testthat)
library(shiny)
library(shinydashboard)
library(shinyFiles)
library(shinyWidgets)
library(shinyjs)
library(waiter)
library(fs)
library(tippy)
#library(shinycssloaders)
#library(shinytest)
#library(htmltools)
library(DT)
# for file processing
library(ggplot2)
library(bedr)
library(RSQLite) 
library(Rsamtools)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)

# set up local database -------

sqlitePath="data/database.sqlite"

saveData <- function(data,table) {
  # Connect to the database
  db <- DBI::dbConnect(RSQLite::SQLite(), sqlitePath)
  # Construct the update query by looping over the data fields
  # query <- sprintf(
  #   "INSERT INTO %s (%s) VALUES ('%s')",
  #   table, 
  #   paste(names(data), collapse = ", "),
  #   paste(data, collapse = "', '")
  # )
  dbWriteTable(db, name = table, value = data,overwrite = TRUE)
  # Submit the update query and disconnect
  #DBI::dbGetQuery(db, query)
  DBI::dbDisconnect(db)
}
loadData <- function(table) {
  # Connect to the database
  db <- DBI::dbConnect(SQLite(), sqlitePath)
  # Construct the fetching query
  #query <- sprintf("SELECT * FROM %s", table)
  # Submit the fetch query and disconnect
  #data <- DBI::dbGetQuery(db, query)
  dbReadTable(db, table)
  DBI::dbDisconnect(db)
  data
}

#' A general function to quickly import a region from tabix-indexed tab-separated files into a data frame
#' modified based on the https://gist.github.com/kauralasoo/464c63085f0516b378ef6eca7b2ce28b
#' 
#' @param tabix_file Path to tabix-indexed text file
#' @param param A instance of GRanges, RangedData, or RangesList 
#' provide the sequence names and regions to be parsed. Passed onto Rsamtools::scanTabix()
#' @param ... Additional parameters to be passed on to data.table::fread()
#'
#' @return List of data_frames, one for each entry in the param GRanges object.
#' @export

scanTabixDataFrame <- function(tabix_file, param, format, ...){
  tabix_list = Rsamtools::scanTabix(tabix_file, param = param, format=format)
  df_list = lapply(tabix_list, function(x,...){
    if(length(x) > 0){
      if(length(x) == 1){
        #Hack to make sure that it also works for data frames with only one row
        #Adds an empty row and then removes it
        result = paste(paste(x, collapse = "\n"),"\n",sep = "")
        result = data.table::fread(text=result, header = F, ...)
      }else{
        result = data.table::fread(text=x, header = F, ...)
      }
    } else{
      #Return NULL if the nothing is returned from tabix file
      result = NULL
    }
    return(result)
  }, ...)
  return(df_list)
}


SegNormRD <- function(df, id, seg.method = "cbs") {
  #SegNormRD.file <- paste0(id, ".seg.rds")
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
    return(seg.output)
  }
  ##EDIT: include SLM segmentation
  if (seg.method == "slm") {
    print("segment with SLM")
    df.ls <- base::split(df, f="V1")
    res <- lapply(df.ls, function(df) {
      slm <-
        SLMSeg::SLM(
          log2(df$ratio + 0.00001),
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
      res.dt <- data.table(ID=id,chrom=chr,loc.start=start,loc.end=end,num.mark=res$lengths,seg.mean=res$values)
    })
    res <- data.table::rbindlist(res)
  }
  
}
#ref_genome="GRCh38"

hg38.info <- seqinfo(BSgenome.Hsapiens.UCSC.hg38::Hsapiens)%>%as.data.frame()
hg38.info <- hg38.info %>% mutate(chrom=rownames(hg38.info))
hg19.info <- seqinfo(BSgenome.Hsapiens.UCSC.hg19::Hsapiens)%>%as.data.frame()
hg19.info <- hg19.info %>% mutate(chrom=rownames(hg19.info))
ReadGVCF <- function(path_to_gVCF,ref_genome=ref_genome,param = param){
  print("scaning the region")
  vcf<- VariantAnnotation::readVcf(file = path_to_gVCF,genome = ref_genome,param = param)
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
  GT.anno <- GT %>% mutate(InhFrom=ifelse(index%in%G3&P1%in%G1&P2%in%c(G2,G3),P2_ID,
                                          ifelse(index%in%G3&P1%in%c(G2,G3)&P2%in%G1,P1_ID,"Notphased")))
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
           pr_ALT_Freq=index_ALT_RD/(index_ALT_RD+index_REF_RD))%>%
    mutate(likelyDN=ifelse(p1_ALT_RD<2&p2_ALT_RD<2&index_ALT_RD>5&p1_REF_RD>10&p2_REF_RD>10&pr_count>10&pr_ALT_Freq>0.2,"TRUE","FALSE"))
  AD.anno <- AD.anno[,c("pr_count","pr_ALT_Freq","likelyDN")]
  merged <- cbind(GT.anno ,AD.anno)
  setnames(merged,c("index","P1","P2"),c(PR_ID,P1_ID,P2_ID))
  mcols(vcf.gr) <- merged
  return(vcf.gr)
}


## plot parameter
style_rd <- theme_classic()+
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "top",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    panel.grid.minor.y = element_blank(),
    #panel.grid.minor.x = element_line(colour = "grey50"),
    panel.grid.major.y = element_line(linetype = 5,colour = "grey50"),
    panel.grid.major.x = element_line(linetype = 5,colour = "grey50"),
    panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 16,face = "bold"),
    #axis.line.x = element_blank(),
    axis.ticks = element_line(color = "black"))
style_snp <- theme_classic()+
  theme(
    plot.title = element_text(face = "bold", size = 12),
    legend.position = "top",
    legend.title = element_text(colour="black", size=12, face="bold"),
    legend.text = element_text(size = 12),
    panel.border = element_blank(),
    panel.grid.minor.y = element_blank(),
    #panel.grid.minor.x = element_line(colour = "grey50"),
    panel.grid.major.y = element_line(linetype = 5,colour = "grey50"),
    panel.grid.major.x = element_line(linetype = 5,colour = "grey50"),
    panel.background = element_blank(),
    axis.text = element_text(color = "black", size = 12),
    axis.title = element_text(color = "black", size = 16,face = "bold"),
    #axis.line.x = element_blank(),
    axis.ticks = element_line(color = "black"))

scale_rd <- scale_y_continuous(name="Log2 Ratio",
                               limits=c(-3, 2),
                               breaks = c(round(log2(1/2),2),
                                          round(log2(2/2),2),
                                          round(log2(3/2),2),
                                          round(log2(4/2),2),
                                          round(log2(5/2),2),
                                          round(log2(6/2),2)
                               ))
scale_snp <- scale_y_continuous(name="B-allele frequency",
                                breaks = c(round(1/2,2),
                                           round(1/3,2),
                                           round(2/3,2),
                                           round(1/4,2),
                                           round(3/4,2),
                                           round(2/5,2),
                                           round(3/5,2)
                                ))
SNPCOLOR2 <- c("#E69F00","#0072B2")


