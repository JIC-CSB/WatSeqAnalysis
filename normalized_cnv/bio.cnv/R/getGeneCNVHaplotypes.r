getGeneCNVHaplotypes <- function(cnvs, gff, gene="TraesCS6D02G400000", flank_size=2000){
	gene_gff <- gff[gff$gene_id == gene]	
	flanked <-  flank(gene_gff, flank_size) 
	gene_gff 
}