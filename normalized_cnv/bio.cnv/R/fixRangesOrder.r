fixRangesOrder<-function(gr,lens, ordr = c("chr1A","chr1B", "chr1D","chr2A","chr2B", "chr2D","chr3A",
                                           "chr3B", "chr3D","chr4A","chr4B", "chr4D","chr5A","chr5B", "chr5D","chr6A","chr6B", "chr6D","chr7A","chr7B", "chr7D","chrUn")){
    as_df <-data.frame(gr)
    as_df$seqnames <- factor(as.character(as_df$seqnames),levels = ordr )
    vect = lens$length
    names(vect) <-lens$seqname
    gr2 <- makeGRangesFromDataFrame(as_df,keep.extra.columns = T, seqinfo=vect)
    trim(gr2)
}

