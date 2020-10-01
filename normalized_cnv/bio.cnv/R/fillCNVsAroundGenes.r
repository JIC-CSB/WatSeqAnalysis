fillCNVsAroundGenes <- function(cnvs, gff, flank_size=2000, line="WATDE0039"){
	flanked <-  flank(gff, flank_size) 
	cnvs_in_line <- cnvs[cnvs$line  == line]
	without_coverage <- BiocGenerics::setdiff(flanked, cnvs_in_line, ignore.strand=TRUE)
	without_coverage$norm_cov <- 1
	without_coverage$cnv_level <- 1
	without_coverage$line <- line
	c(without_coverage,cnvs_in_line)
}