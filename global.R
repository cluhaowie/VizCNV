#  ------------------------------------------------------------------------
#
# Title : App - VizCNV
#    By : Haowei Du, Cliff Lun
#  Date : April 2022
#    
#  ------------------------------------------------------------------------
options(timeout = 6000)
options(scipen=999)
options(shiny.maxRequestSize=3*1024^3) ## max file size 3 Gb
options(shiny.autoreload=TRUE)
#options(shiny.reactlog=TRUE) 
#library(BiocManager)
#options(repos = BiocManager::repositories())


# Packages ----------------------------------------------------------------
# Install missing packages from CRAN, 'arrow' may be a problem
list.of.packages <- c("dplyr", "data.table", "shiny", "shinydashboard", "shinyFeedback",
                      "tippy","DT","ggplot2","shinyWidgets","shinyFiles","waiter",
                      "cowplot","devtools","BiocManager","arrow","colourpicker", "shinyjs","rclipboard",
                      "shinydashboardPlus","bs4Dash", "colourpicker","tidyr","config") 

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Install missing packages from Bioconductor
biocLitePackages <- c("DNAcopy", "GenomicRanges", "VariantAnnotation","bedr","plyranges") 
new.biocLitePackage <- biocLitePackages[!(biocLitePackages %in% installed.packages()[,"Package"])]
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
if(length(new.biocLitePackage)) BiocManager::install(new.biocLitePackage)
# Install SLMSeg package from github
local.packages <- c("SLMSeg","regioneR","ggtranscript")
add.packages <- local.packages[!local.packages%in% installed.packages()[,"Package"]]
if(length(add.packages)) devtools::install_github(c("cluhaowie/VizCNV/SLMSeg","bernatgel/regioneR","dzhang32/ggtranscript"))
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
library(plyranges)
library(stringr)
library(dplyr)
library(data.table)
library(shiny)
library(shinydashboard)
library(shinydashboardPlus)
library(bs4Dash) ## support bootstrap 4
library(shinyFiles)
library(shinyWidgets)
library(waiter)
library(fs)
library(tippy)
library(DT)
# for file processing
library(ggplot2)
library(bedr)
library(Rsamtools)
library(VariantAnnotation)
library(arrow) ## read parquet data
library(shinyjs)
library(ggtranscript)
library(rclipboard)
library(colourpicker)
library(tidyr)


# set up local database -------
config <- config::get()
maxSize_anno <- config$maxSize_anno 
maxtranscript <- config$maxtranscript 
geneExtend <- config$geneExtend 
minseg <- config$minseg 
minsegmean <- config$minsegmean


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


normalization_method <- function(df, chr, norm_option = "chr_med"){
  if (norm_option == "chr_med"){
    df_chr <- df%>%filter(V1 == chr) ## this is needed for chrY normalization
    tmp_median <- df_chr$V4[df_chr$V4!=0]
    tmp <- df%>%
      filter(V1 == chr)%>%
      mutate(ratio=V4/median(tmp_median+0.00001))
  } else if (norm_option == "wg_med"){
    tmp <- df%>%
      mutate(ratio=V4/median(V4+0.00001)) %>% 
      filter(V1 == chr)
  }
  return(tmp)
}

## plot parameter


style_rd <- theme_classic()+
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        legend.title = element_text(colour="black", size=12),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_line(linetype = 4,colour = "grey85"),
        panel.grid.major.y = element_line(linetype = 5,colour = "grey70"),
        panel.grid.major.x = element_line(linetype = 5,colour = "grey50"),
        panel.background = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = c("black", "white"), size = 10),
        axis.title = element_text(color = "black", size = 12),
        axis.ticks = element_line(color = "black"))

style_snp <- theme_classic()+
  theme(plot.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        legend.title = element_text(colour="black", size=12),
        legend.text = element_text(size = 12),
        panel.border = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_line(linetype = 4,colour = "grey85"),
        panel.grid.major.y = element_line(linetype = 5,colour = "grey70"),
        panel.grid.major.x = element_line(linetype = 5,colour = "grey50"),
        panel.background = element_blank(),
        axis.text.y = element_text(color = "black", size = 10),
        axis.text.x = element_text(color = c("black", "white"), size = 10),
        axis.title.y = element_text(color = "black", size = 12),
        axis.title.x = element_blank(),
        axis.ticks = element_line(color = "black"))

scale_rd <- scale_y_continuous(name="Log2 Ratio",
                               limits=c(-2.5, 2),
                               breaks = c(-2,
                                 round(log2(1/2),2),
                                 round(log2(2/2),2),
                                 round(log2(3/2),2),
                                 round(log2(4/2),2),
                                 round(log2(5/2),2),
                                 round(log2(6/2),2)
                               ))
scale_snp <- scale_y_continuous(name="B-allele frequency",
                                limits = c(-0.05, 1.1),
                                breaks = c(0,
                                           round(1/2,2),
                                           round(1/3,2),
                                           round(2/3,2),
                                           round(1/4,2),
                                           round(3/4,2),
                                           round(2/5,2),
                                           round(3/5,2),
                                           1))





style_anno <- theme_classic()+
  theme(
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), #remove x axis ticks
        axis.title.x = element_blank(),
        panel.grid.major.x = element_line(linetype = 5,colour = "grey50"),
        axis.text.y=element_text(color = "white"),  #remove y axis labels
        axis.ticks.y=element_blank()
  )

scale_anno <- scale_y_continuous(limits = c(-0.01,.11))


style_genes <- style_rd+
  theme(panel.grid.major.y = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(color = "white"))
scale_genes <- scale_y_continuous(labels = scales::label_number(accuracy = 0.01))

SNPCOLOR2 <- c("#39918C","#E69F00")
CNVCOLOR6 <- c("#00468b","#00468b","#8b0000","#8b0000","#008b46","#008b46")

SKYCOLOR <- c("chr1"="#fdfa01","chr2"="#a70106","chr3"="#beafc6","chr4"="#06fafb",
              "chr5"="#d75000","chr6"="#d60e51","chr7"="#fca4a3","chr8"="#fa6104",
              "chr9"="#A9A9A9","chr10"="#038201","chr11"="#05a5ff","chr12"="#fd01fe",
              "chr13"="#fd0100","chr14"="#fe5579","chr15"="#87d4a7","chr16"="#fa9308",
              "chr17"="#0154d1","chr18"="#e60052","chr19"="#52fc56","chr20"="#0325ff",
              "chr21"="#fefda7","chr22"="#fd8afc","chrX"="#004800","chrY"="#01b600")

names(CNVCOLOR6) <- c("<TRP>","TRP","<DUP>","DUP","<DEL>","DEL")
scale_SVType <- scale_fill_manual(CNVCOLOR6)
chrom_id <- c(1:22,"X","Y") ## incorporate chrY
names(chrom_id) <- paste0("chr",chrom_id)

## keep column from ucsc track
segdups.keep.col <- c("chrom","chromStart","chromEnd","strand","name","uid","fracMatch","fracMatchIndel","level")
rmsk.keep.col <- c("genoName","genoStart","genoEnd","strand","repName","repClass","repFamily")
