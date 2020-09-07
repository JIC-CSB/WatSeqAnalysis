#!/usr/bin/env Rscript
library(devtools)
library(GenomicRanges)
library(RSQLite)
library("optparse")
load_all("./R", TRUE)

options(echo=TRUE) 


option_list = list(
    make_option(c("-l", "--line"),
                type="character", 
                default=NULL, 
                help="Coverage file", metavar="CHARACTER"),
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

if (is.null(opt$line)){
  print_help(opt_parser)
  stop("Please select a line to extract.", call.=FALSE)
}

covs_db = dbConnect(drv=RSQLite::SQLite(), dbname=opt$database)

regions <- getFullChromosomes(covs_db)
output_folder<-paste0(opt$out, "/", opt$line)
#output_folder<-opt$out
dir.create(output_folder, showWarnings = FALSE, recursive = TRUE)

setwd(output_folder)
print(regions)
first <- TRUE
for(i in seq(from=1, to=length(regions), by=1)){
    region = regions[i] 
    #end(region) <- 10000000
    out_f=paste0(seqnames(region)[1], ".csv")
    #out_f=paste0(opt$line, ".csv")
    #print(out_f)
    #print(region)
    dels <- getCNVInRegion(covs_db, region=region, line=opt$line, max_gap=opt$max_gap, het=FALSE)
    #print(length(dels) )
    #dels <- dels[length(dels) > opt$minimum_length]
    write.csv(dels, file=out_f, quote=FALSE, append = !first, row.names=FALSE, col.names = first )
    first <- FALSE
}


