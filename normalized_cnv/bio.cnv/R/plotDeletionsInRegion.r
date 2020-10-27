plotDeletionsInRegion <- function(filled_cnvs, region, gff=NULL, lines=unique(filled_cnvs$line)){

	valuesToPlot <- subsetByOverlaps(filled_cnvs, region, type="any", ignore.strand=TRUE)
	valuesToPlot <- valuesToPlot[valuesToPlot$line %in% lines]

	valuesToPlot$copies <- ifelse( 
		valuesToPlot$cnv_level > 3,
		">3", 
		as.character(valuesToPlot$cnv_level))
	#valuesToPlot$copies <- as.character(valuesToPlot$cnv_level)

	ggplot(valuesToPlot, aes(fill=copies), which=region) + 
	geom_alignment(facets = line ~ seqnames) + 
	scale_fill_manual(values = c(
		"0"="#d73027", 
		"1"="#ffffbf",  
		"2"="#74add1", 
		"3"="#4575b4",
		">3"="#313695")) +
	theme(strip.text.y.right = element_text(angle = 0)) + xlim(region)
}