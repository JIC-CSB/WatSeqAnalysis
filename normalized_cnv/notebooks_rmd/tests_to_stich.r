library(devtools)
library(knitr)
library(GenomicRanges)
library(GenomicFeatures)
library(RSQLite)
library(ggbio)
library(rtracklayer)
require(plyr)
require(dplyr)
library("parallel")


setwd("~/Documents/public_code/WatSeqAnalysis/normalized_cnv/notebooks_rmd/")
load_all("../bio.cnv/R", TRUE)

params<-list()
params$data_path <- "~/Documents/WatSeq/Deletions/20200724/200bp/"
params$gff <- "~/Dropbox/JIC/WatSeq/Annotation/IWGSCv1.0_UTR_ALL.parts.tidy.gff3.gz"
#params$gff <- "~/Dropbox/JIC/WatSeq/Annotation/Triticum_aestivum.IWGSC.41.parts.tidy.gff3.gz"


gff  <- import(params$gff)
cnvs_file <- paste0(params$data_path, "/long_tables/cnv_out_M10/M10_chr6D_part2.cnv.for.db.tsv") 
cnvs_df_M10 <- read.csv(cnvs_file, sep = "\t" )
cnvs_M10 <- makeGRangesFromDataFrame(cnvs_df_M10, keep.extra.columns=T)
gff6D <- gff[seqnames(gff) == "chr6D_part2" & gff$type == "gene"]


to_ignore <- c("WATDE0009", "WATDE0039",  "WATDE0056", "WATDE0060" ,  "WATDE0090")

cnvs_M10 <- cnvs_M10[cnvs_M10$line %in% setdiff(unique(cnvs_M10$line), to_ignore) ]

filled_M10  <- fillAllCNVsAroundGenes(cnvs_M10, gff6D, mc.cores = 4)

test_line_filled<-filled_M10[filled_M10$line == "WATDE0456"]

test_filled_subset <- test_line_filled[start(test_line_filled) > 4400000 & end( test_line_filled) < 4500000]

autoplot(test_filled_subset, aes(fill=cnv_level))
test_filled_subset

cnvs_M10[start(cnvs_M10) > 4400000 & end( cnvs_M10) < 4500000 & cnvs_M10$line == "WATDE0456" ]

mini_test <- filled_M10[start(filled_M10) > 4400000 & end( filled_M10) < 4500000]

stich_cnv_level <- function(filled, line = "WATDE0456"){
  filled_sorted <- filled[filled$line == line]
  print(length((filled_sorted)))
  print(head(filled_sorted))
  print("SORTING")
  filled_sorted<- sort(filled_sorted)
  print(length((filled_sorted)))
  print(head(filled_sorted))
  #print("TO DF")
  #filled_sorted <- as.data.frame(filled_sorted)
  last <- filled_sorted[1,]
  
  ret <- filled_sorted[0,]
  #print(length((filled_sorted)))
  #print(head(filled_sorted))
  for(i in seq(1, length(filled_sorted))){
    current = filled_sorted[i,]
    if(as.character(seqnames(current)[1]) != as.character(seqnames(last)[1]) || current$cnv_level != last$cnv_level){
    #if(current$seqnames != last$seqnames || current$cnv_level != last$cnv_level ){
      ret[length(ret)+1,] <- last
      last <- current
    }else{
      #last$end = current$end
      end(last) <- end(current)
    }
  }
  ret[nrow(ret)+1,] <- last
  ret$norm_cov <- NULL
  ret
  #makeGRangesFromDataFrame(ret, keep.extra.columns = T)
}



test_line <- stich_cnv_level(filled_M10, line="WATDE0456")
View(data.frame(test_line))

length(filled_M10)


test_stiched <- stich_cnv_level(mini_test)
mini_test
test_stiched
autoplot(test_stiched, aes(fill=cnv_level))
autoplot(mini_test, aes(fill=cnv_level))