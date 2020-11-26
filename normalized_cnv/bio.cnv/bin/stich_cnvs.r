#!/usr/bin/env Rscript
library(devtools)
library(GenomicRanges)
library(rtracklayer)
library("optparse")
load_all("./R", TRUE)

options(echo=TRUE) 


option_list = list(
    make_option(c("-i", "--input"), 
                type="character",
                default="all_cnvs.csv.gz",
                metavar="FILE",
                help="File with the output from merge_cnv_tables.r [default= %default]" 
    ),
    make_option(c("-g", "--gff"),
                type="character", 
                default="annotation.gff.gz", 
                help="GFF file contaiong the genes. The regions used are the lines of type 'gene' and the 'flanking' size is added [default= %default]", 
                metavar="FILE"),
    make_option(c("-o", "--out"), 
                type="character", 
                metavar="FILE", 
                default="./stichedcnv.txt", 
                help="Output folder with the stiched [default= %default]"
                ),
    make_option(c("-f", "--flanking"), 
                type="integer", 
                metavar="number",
                default=2000, 
                help="Length [default= %default]"),
    make_option(c("-p", "--cores"), 
                type="integer", 
                metavar="number",
                default=4,
                help="Number of cores to use [default= %default]"), 
    make_option(c("-n","--njob" ), 
                type="integer", 
                metavar="number", 
                default=1,
                help="Total number of jobs"),
    make_option(c("-j", "--job"), 
                type="integer", 
                metavar="number", 
                default=1,
                help="Current job. The output will be suffixed with this number."),
    make_option(c("-c", "--chromosome"), 
                type="character", 
                default=NULL, 
                help="Chromosome to process.")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt

gff  <- import(opt$gff)
gff <- gff[gff$type == "gene"]
cnvs <- read.csv(gzfile(opt$input))
cnvs <- makeGRangesFromDataFrame(cnvs, keep.extra.columns=T)
if(! is.null(opt$chromosome )){
    gff  <- gff[seqnames(gff) == opt$chromosome]
    cnvs <- cnvs[seqnames(cnvs) == opt$chromosome]
    opt$out <- paste0(opt$out,".",opt$chromosome)
}
gc()

#exit()
filled <- fillAllCNVsAroundGenes(cnvs, gff, opt$flanking, mc.cores=opt$cores)
gz1 = gzfile(paste0(opt$out,".filled.csv.gz"),"w");
write.csv(filled, file=gz1, row.names=F);
close(gz1) 

stiched <- stichAllCnvLevel(filled,mc.cores=opt$cores)

gz1 = gzfile(paste0(opt$out,".stiched.csv.gz"),"w");
write.csv(stiched, file=gz1, row.names=F);
close(gz1) 



