fillCNVsAroundGenes <- function(cnvs, gff, flank_size=2000,  line="WATDE1004"){
	flanked <- gff
	if(flank_size > 0){
		flanked <-  exetendGR(gff, size=flank_size)
		flanked <- union(flanked, cnvs)
	}
	cnvs_in_line <- cnvs[cnvs$line  == line]
	without_coverage <- BiocGenerics::setdiff(flanked, cnvs_in_line, ignore.strand=TRUE)
	without_coverage <- without_coverage[width(without_coverage) > 1] #This is a patch to remove tiny overhangs. 
	without_coverage$norm_cov <- 1
	without_coverage$cnv_level <- 1
	without_coverage$line <- line
	sort(c(without_coverage,cnvs_in_line))		
}