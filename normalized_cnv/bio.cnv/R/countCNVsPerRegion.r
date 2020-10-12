countCNVsPerRegion <- function(cnvs_df){
	cnvs_df %>% count(seqnames, start, end, cnv_level)
}
