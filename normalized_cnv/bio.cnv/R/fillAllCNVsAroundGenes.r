fillAllCNVsAroundGenes <- function(cnvs, gff, flank_size=2000,  lines=unique(cnvs$line), mc.cores=4) {
	flanked <-  exetendGR(gff, size=flank_size)
	cnvs <- subsetByOverlaps(cnvs, flanked, ignore.strand = TRUE )
	flanked <- BiocGenerics::union(flanked, cnvs)
	ret <- mclapply(lines, function(x){fillCNVsAroundGenes(cnvs, flanked, flank_size = 0, line=x)},  mc.cores = mc.cores) 
	ret <- as(ret, "GRangesList")
	unlist(ret)
}