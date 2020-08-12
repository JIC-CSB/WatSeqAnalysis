plotRegionScatterLines <- function(covs_db,gff, region, show.ribons=F, genes=c(),  lines=c(), norm_cov_range = c(0,3), window_size=250000, levels_count = 10, show.coverage = FALSE  ){
    local_cov <- getWindowsInRange(covs_db,  ranges=region)
    gff_local <-  subsetByOverlaps(gff, region)
    gff_local <- gff_local[gff_local$type == "gene"]
    gff_local$pos_id <- seq(from = 0, to=length(gff_local)-1)
    gff_local$ymax <- 0.5 + (gff_local$pos_id %% levels_count)
    gff_local$ymin <- gff_local$pos_id %% levels_count
    gff_local$highligthed <- gff_local$gene_id %in% genes

    gff_h <- gff_local[gff_local$highligthed]


    p.norm.cov  <- ggplot() + geom_vline(xintercept=start(gff_h), col="red")
    p.norm.cov  <- p.norm.cov + theme_bw()

    #regions     <- data.frame( getRegionFromDB(covs_db, region=region, rename=FALSE) )
    if(show.ribons){
        regions     <- getRegionFromDB(covs_db, region=region, rename=TRUE) 
        regions$y_min = 1 - (2* regions$sd )
        regions$y_max = 1 + (2* regions$sd )

        #y_min = 0.5
        #y_max = 1.5
        p.norm.cov  <- p.norm.cov  + geom_ribbon(data = regions, mapping=aes(x = start ,ymin = y_min, ymax = y_max), fill = "grey70")
        #return(head(regions));
    }else{
         p.norm.cov <- p.norm.cov + geom_boxplot(data=data.frame(local_cov), mapping=aes(y=norm_cov, x=start, group= window_size * round(start / window_size) ), outlier.shape = NA) 
    }
    p.norm.cov  <- p.norm.cov + theme(legend.position = "top")
    p.norm.cov  <- p.norm.cov + ylim(norm_cov_range)

    p.genes <- ggplot(gff_local) + geom_rect(aes(xmin=start, xmax=end , ymin = ymin, ymax = ymax ), stat = "identity", color = NA) 
    p.genes <- p.genes +  geom_text(aes(x=start, y =  ymax+0.2, label = gene_id, color=highligthed), size=2, check_overlap=T ) 
    p.genes <- p.genes + theme_bw() + theme(legend.position = "none") + scale_colour_manual(values=c("TRUE"="red", "FALSE"="black"))

    heights = c()
    tracks_details <- list()

    for(line in lines){
        tmp_data <- data.frame(local_cov[local_cov$line == line])
        
        if(show.coverage){
            title_cov <- paste0(line, " cov")
             tracks_details[[title_cov]] <- ggplot()  + geom_point(data=(tmp_data), mapping = aes(x=start, y=cov),shape=".") + theme_bw()
             heights <- append(heights,1)
        }
       
        
        tracks_details[[line]] <-  p.norm.cov  + geom_point(data=(tmp_data), mapping = aes(x=start, y=norm_cov),shape=".")
        heights <- append(heights,1)
    }

    tracks_details[["genes"]] <- p.genes
    heights <- append(heights,1)

    tks <- tracks(tracks_details, heights = heights )  + xlim(gr4D)  + scale_x_sequnit(unit = "Mb" ) 
    tks
    #gff_local
    #regions
}