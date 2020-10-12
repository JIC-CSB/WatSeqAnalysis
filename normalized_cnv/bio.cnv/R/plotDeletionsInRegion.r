plotDeletionsInRegion <- function(filled_cnvs, region, gff=NULL, lines=unique(filled_cnvs$line)){

	valuesToPlot <- subsetByOverlaps(filled_cnvs, region, type="any", ignore.strand=TRUE)
	valuesToPlot <- valuesToPlot[valuesToPlot$line %in% lines]
	ggplot(valuesToPlot, aes(fill=cnv_level), which=region) + 
	geom_alignment(facets = line ~ seqnames) + 
	theme(strip.text.y.right = element_text(angle = 0)) #+ xlim(region)
}