#!/usr/bin/env Rscript
library(devtools)
library(GenomicRanges)
library(RSQLite)
library("optparse")
load_all("./R", TRUE)

options(echo=TRUE) 


option_list = list(
    make_option(c("-M", "--max_gap"), 
                type="integer",
                default=3,
                help="Maximum gap length"
    ),
    make_option(c("-d", "--database"), type="character", default="dels.db", 
                help="SQLite database containing the deletion details (see documentation)", 
                metavar="FILE"),
    make_option(c("-o", "--out"), type="character", default="./cnv_out", 
                help="output folder [default= %default]. This folder is created if it doesn't exist", 
                metavar="DIR"),
    make_option(c("-m", "--minimum_length"), type="integer", default=200, 
                help="Minimum length [default= %default]")
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt

output_folder<-opt$out
#dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)
setwd(output_folder)


covs_db = dbConnect(drv=RSQLite::SQLite(), dbname=opt$database)
lines <- getLinesFromDB(covs_db)
print(lines)

regions <- getFullChromosomes(covs_db)

print(regions)


missing_chunks_fn <- "missing_chunks.txt"

all<- NULL
#file.create(missing_chunks_fn, showWarnings = TRUE)
for(line in lines$line){
    output_folder<-paste0(opt$out, "/", line)
    print(line)
    for(i in seq(from=1, to=length(regions), by=1)){
         region = regions[i]
         in_f=paste0(line, "/", seqnames(region)[1], ".csv")
         #print(in_f)
         #in_f=paste0(region, ".csv")
         in_pair = paste0(c(line,"," , seqnames(region)[1]))
         if(file.exists(in_f)){
            tmp <- read.csv(in_f)
            if(is.null(all)){
                all <- tmp
            }else{
                all <- rbind(all, tmp)
            }
         }else{
             cat(in_pair,file=missing_chunks_fn,sep="\n",append=TRUE) 
         }
#  #       dels <- getCNVInRegion(covs_db, region=region, line=opt$line, max_gap=opt$max_gap, het=FALSE)
  #      write.csv(dels, file=out_f, quote=FALSE, append = !first, row.names=FALSE, col.names = first )
    }
}

write.csv(all, file="all_cnv.csv", quote=FALSE, row.names=FALSE )

