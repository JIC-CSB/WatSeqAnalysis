plotCNVDensity<-function(dels,tilewidth=100000, base=2, lens=NULL){
    #print(head(dels))
    #v2 <- dels[dels$type==type]
    v2 <- dels
    cov_v2 <- coverage(x = v2)
    bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(v2),
                                  tilewidth=tilewidth,
                                  cut.last.tile.in.chrom=TRUE)
    gr_bins <- binnedAverage(bins, cov_v2, "binned_cov")
    gr_bins <-fixRangesOrder(gr_bins, lens=lens, ordr=lens$seqname)
    gr_bins$log_cov <- ceiling(log(gr_bins$binned_cov+1,base=base))
    gr_bins$log_cov <- ifelse(gr_bins$log_cov < 1 ,0, gr_bins$log_cov )
    cov_ranges<-list()
    for(r in sort(unique(gr_bins$log_cov))){
        tmp_gr = gr_bins[gr_bins$log_cov == r]
        cov_ranges[as.character(r)] = paste0(round(min(tmp_gr$binned_cov),2),"-", round(max(tmp_gr$binned_cov),2))
    }
    gr_bins$cov_legend <- factor(as.character(cov_ranges[as.character(gr_bins$log_cov)]), levels=cov_ranges)
    title <- paste0(type, "  per bin size:", tilewidth, " Base: ", base)
    autoplot(gr_bins,  layout = "karyogram", aes(fill=cov_legend)) + ggtitle(title) + scale_fill_brewer(type="seq", palette = "Oranges") 
}