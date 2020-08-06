#!/usr/bin/env Rscript
library("optparse")
options(echo=TRUE) 
require('bio.tilling')
library(grid)
library(ggplot2)

countHigherThan <- function(x, min_cov=0){
    length( x[  x >= min_cov] )
}

filterRegionsWithLowCoverage <- function(mat, min_lines=0.0, min_cov=0){
  counts<- apply(mat, 1, countHigherThan, min_cov = min_cov )
    pass_filter <- counts > min_lines * ncol(mat)
    mat[pass_filter, ]
}

plotHistogram<-function(table,
                        column="size_cds",
                        probs = c( 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                        binwidth=0.1, trim_range=T, title_prefix ="", force_range=c(), filter_zero = TRUE){
    table <- data.frame(table)
    if(filter_zero){
       table<-table[table[,column]>0,] 
    }
    
    quantiles <- data.frame(quantile(table[,column], prob=probs,na.rm=TRUE, include.lowest=TRUE), stringsAsFactors=FALSE)
    quantiles$quant<-rownames(quantiles)
    colnames(quantiles)<-c("value", "quant")
    values<-quantiles$values
    local_mean<-mean(table[,column],na.rm=TRUE)
    local_sd<-sd(table[,column],na.rm=TRUE)
    local_max <-  max(table[,column],na.rm = T)
    p <- ggplot(table, aes_string(column))
    if(nrow(table) > 100){
        table <- within(table,
           quantile <- as.integer(
               cut(table[,column],
                   unique(quantile(table[,column],
                    prob=probs,
                    na.rm=TRUE,
                    include.lowest=TRUE))
                   )
               ))
        table$quantile<-ifelse(is.na(table$quantile),0,table$quantile)
        table$quantile<-as.factor(table$quantile)
        
        if(length(force_range) == 2){
           xmin <- force_range[1]
           xmax <- force_range[2]
        }else{
            
            iq <- quantiles$value[4] - quantiles$value[2]
            xmax <- quantiles$value[3] + (iq * 2)
            xmin <- quantiles$value[3] - (iq * 2)
            if(xmin < 0) xmin <- 0
            if(xmax > local_max)  xmax <- local_max + 1
        }
        
        p <- ggplot(table, aes_string(column, fill="quantile"))
        p <- p + geom_vline(data=quantiles,aes(xintercept=quantiles$value) )
        for(i in seq(1,nrow(quantiles))){
            x_pos<-quantiles$value[i]
            gtext <- textGrob(quantiles$quant[i], y=0.02,  gp = gpar(fontsize = 6,col = "red"))
            p <- p + annotation_custom(gtext, xmin=x_pos, xmax=x_pos)
        }
       
        p <- p +  scale_fill_brewer(palette="Dark2")
        if(length(force_range) == 2){
           xmin <- force_range[1]
           xmax <- force_range[2]
        }
        if(trim_range){
           p <- p  + xlim(xmin, xmax) 
        }    
    }
    p <- p + geom_histogram(binwidth=binwidth, position = "identity") + theme_bw()
    p <- p + theme(legend.position="none")
    p <- p + ggtitle(paste0(title_prefix, "Mean: ", round(local_mean,2),
        " SD:", round(local_sd,2),
        " CV:", round(local_sd/local_mean, 2),
        " Median:", round(median(table[,column],na.rm=TRUE),2),
        " Max:", round(local_max,2),
        " N:", nrow(table)))
    p <- p + theme(plot.title = element_text(size=6))
    p
}


normalizeCovs<-function(counts, exonsDF, remove.zero = F) {
  
  #Filtering the rows with AVG=0 
  rowsToKeep<-rowSums(counts)>0
  counts<-counts[rowsToKeep,]
  exonsDF <- exonsDF[rowsToKeep,]

  #counts <- as.matrix(counts)
  if(remove.zero){
       counts<- apply(counts, 2 , function(x){ifelse(x == 0, NA, x)}) 
  }
  
  exonLengths<-exonsDF$ExonL
  totalReadsPerSample<-apply(counts,2,sum, na.rm=T)
  multiplier<-outer(exonLengths, totalReadsPerSample )
  multiplier<-1/multiplier
  multiplier<-1000000000*multiplier
  #"RPKM-like covs"
  multiplier<-counts*multiplier
  
  #Normalize by multiplying by the mean of each row.  
  mean.mult <- apply(multiplier, 1, mean, na.rm = T)
  covs<-sweep(multiplier,MARGIN=1,mean.mult,'/')
  covs <- ifelse(is.na(covs), 0, covs)
  if(remove.zero){
    covs<- apply(covs, 2 , function(x){ifelse(is.na(x), 0, x)}) 
  }
  covs
}

option_list = list(
  make_option(c("-f", "--coverage_file"), 
  	type="character", 
  	default=NULL, 
  	help="Coverage file", 
  	metavar="character"),
  make_option(c("-s", "--sampleSD"), 
  	type="double", 
  	default=0.3, 
  	help="Maximum Standard Deviation allowed in sample [default= %default]"),
  make_option(c("-w", "--windowSD"), 
  	type="double", 
  	default=0.3,
  	help="Maximum Standard Deviation allowed in each window [default= %default]"),
  make_option(c("-z", "--gzip"), 
  	action="store_true", 
  	default=FALSE,
  	help="The coverage file is gzipped [%default]"),
  make_option(c("-o", "--out"), 
  	type="character", 
  	default="./deletions_out", 
  	help="output folder [default= %default]", 
  	metavar="character"),
  make_option(c("-c", "--min_cov"),
    type="integer", 
    default=5, 
    help="Minimum coverage for the region filter [default=%default]"),
  make_option(c("-l","--min_lines"),
    type="double", 
    default=0.05, 
    help="Fraction of the samples needed to have a coverage over 'min_cov' to be kept [default= %default]"
  ),
  make_option(c("-0", "--remove_zero"),
    action="store_true",
    default=FALSE, 
    help="Exclude regions without coverage from normalisation. "
  )  
); 

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);
opt
output_folder<-opt$out

if (is.null(opt$coverage_file)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (coverage_file).",
  	call.=FALSE)
}


dir.create(output_folder, recursive = TRUE)
setwd(output_folder)

filename<-unlist(opt$coverage_file)

is_gz <- opt$gzip

covs<-readCoverageTable(filename, is_gz=is_gz)
print_pdfs <- FALSE
if(print_pdfs){
  pdf("cov_distributions.pdf")
  for(line in colnames(covs)){
    p<-plotHistogram(covs, column=line , binwidth = 1 )
    print(p)
  }
  dev.off()  
}

if(opt$sampleSD > 0){
  covs<-filterSamplesPerSD(covs, maxSD=opt$sampleSD)
}
covs <- filterRegionsWithLowCoverage(covs, min_cov=opt$min_cov, min_lines=opt$min_lines)

if(print_pdfs){
  pdf("cov_filter_distributions.pdf")
  for(line in colnames(covs)){
    p<-plotHistogram(covs, column=line , binwidth = 1 )
    print(p)
  }
  dev.off()
}

df<-getExonsDF(covs)
mat<-normalizeCovs(covs, df, remove.zero=opt$remove_zero)
rm(covs)


#mat<-filterLowQualityExons(mat, maxSD=opt$windowSD)
gc()
write.csv(mat, file='mat.csv')
#saveRDS(mat, "mat.rds")

df<-getExonsDF(mat)
write.csv(df, file='df.csv')
#saveRDS(df, "df.rds")

libSD<-getLibSD(mat)
write.csv(libSD, file='libSD.csv')
#saveRDS(libSD, "libSD.rds")

pdf("norm_cov_distributions.pdf")
for(line in colnames(mat)){
  p<-plotHistogram(mat, column=line , binwidth = 0.01, trim_range=F )
  p <- p + coord_cartesian(xlim=c(0,2)) 
  print(p)
}
dev.off()

