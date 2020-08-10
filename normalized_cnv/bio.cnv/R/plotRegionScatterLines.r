plotRegionScatterLines <- function(covs_db,gff, region, lines=c(), window_size=250000, levels_count = 10 ){
    local_cov <- getWindowsInRange(covs_db,  ranges=region)
    gff_local <-  subsetByOverlaps(gff, region)
    gff_local <- gff_local[gff_local$type == "gene"]
    gff_local$pos_id <- seq(from = 0, to=length(gff_local)-1)
    gff_local$ymax <- 0.5 + (gff_local$pos_id %% levels_count)
    gff_local$ymin <- gff_local$pos_id %% levels_count

    p.norm.cov  <- ggplot(local_cov, aes(y=norm_cov, x=start, group= window_size * round(start / window_size) )) 
    p.norm.cov  <- p.norm.cov + theme_bw() + geom_boxplot(outlier.shape = NA) 
    p.norm.cov  <- p.norm.cov + theme(legend.position = "top")
    p.norm.cov  <- p.norm.cov + ylim(c(0,3))

    p.genes <- ggplot(gff_local) + geom_rect(aes(xmin=start, xmax=end , ymin = ymin, ymax = ymax ), stat = "identity", color = NA) 
    p.genes <- p.genes +  geom_text(aes(x=start, y =  ymax+0.2, label = gene_id), size=2, check_overlap=T ) 
    p.genes <- p.genes + theme_bw()

    heights = c()
    tracks_details <- list()

    for(line in lines){
        tmp_data <- data.frame(local_cov[local_cov$line == line])
        tracks_details[[line]] <-  p.norm.cov  + geom_point(data=(tmp_data), shape=".")
        heights <- append(heights,1)
    }

    tracks_details[["genes"]] <- p.genes
    heights <- append(heights,1)

    tks <- tracks(tracks_details, heights = heights )  + xlim(gr4D)  + scale_x_sequnit(unit = "Mb" ) 
    tks
    #gff_local
}