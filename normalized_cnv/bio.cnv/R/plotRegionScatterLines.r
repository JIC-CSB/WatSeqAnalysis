plotRegionScatterLines <- function(covs_db,gff, region, show.ribons=F,
    genes=c(),  lines=c(), norm_cov_range = c(0,3), window_size=250000, txdb=NULL ,
    gene.details=FALSE, levels_count = 10, 
    show.coverage = FALSE, sd.border=2.5, shape=".",
    stiched.cnv = NULL  ){

    local_cov <- getWindowsInRange(covs_db,  region=region,as.gr=TRUE)
    gff_local <-  subsetByOverlaps(gff, region)
    gff_local <- gff_local[gff_local$type == "gene"]
 
    p.norm.cov  <- ggplot() 
    if(length(gff_local)>0){
        gff_local$pos_id <- 1
        gff_local$pos_id <- seq(from = 0, to=length(gff_local)-1)
        gff_local$ymax <- 0.5 + (gff_local$pos_id %% levels_count)
        gff_local$ymin <- gff_local$pos_id %% levels_count
        gff_local$highligthed <- gff_local$gene_id %in% genes

        gff_h <- gff_local[gff_local$highligthed]

         p.norm.cov  <- p.norm.cov + geom_vline(xintercept=start(gff_h), col="red")
    }
    
    #regions     <- data.frame( getRegionFromDB(covs_db, region=region, rename=FALSE) )
    if(show.ribons){
        regions     <- getRegionFromDB(covs_db, region=region, rename=TRUE) 
        regions$y_min = 1 - (sd.border * regions$sd )
        regions$y_max = 1 + (sd.border * regions$sd )
        regions$y_min <- ifelse(regions$y_min < 0,  0 , regions$y_min)
        p.norm.cov  <- p.norm.cov  + geom_ribbon(data = regions, mapping=aes(x = start ,ymin = y_min, ymax = y_max), fill = "grey70")
    }else{
         p.norm.cov <- p.norm.cov + geom_boxplot(data=data.frame(local_cov), mapping=aes(y=norm_cov, x=start, group= window_size * round(start / window_size) ), outlier.shape = NA) 
    }
    p.norm.cov  <- p.norm.cov + theme_bw()  

    p.genes <- ggplot()
    if(length(gff_local)>0){
        if(is.null(txdb)){
            p.genes <- ggplot(gff_local) + geom_rect(aes(xmin=start, xmax=end , ymin = ymin, ymax = ymax ), stat = "identity", color = NA) 
            p.genes <- p.genes +  geom_text(aes(x=start, y =  ymax+0.2, label = gene_id, color=highligthed), size=2, check_overlap=T ) 
            p.genes <- p.genes + theme_bw() + theme(legend.position = "none") + scale_colour_manual(values=c("TRUE"="red", "FALSE"="black"))
        }else{
            p.genes <- autoplot(txdb, which = region) + theme_bw() 
        }
    }

    heights = c()
    tracks_details <- list()

    for(line in lines){
        tmp_data <- data.frame(local_cov[local_cov$line == line])
        
        if(show.coverage){
            title_cov <- paste0(line, " cov")
             tracks_details[[title_cov]] <- ggplot()  + geom_point(data=(tmp_data), mapping = aes(x=start, y=cov),shape=shape) + theme_bw()
             heights <- append(heights,1)
        }
       
        tmp_p <- p.norm.cov  + geom_point(data=(tmp_data), mapping = aes(x=start, y=norm_cov),shape=shape) 
        if(!is.null(stiched.cnv)){

            stich_local <-  subsetByOverlaps(stiched.cnv, region)
            tmp_st <- data.frame(stich_local[stich_local$line == line])
            tmp_p <- tmp_p + geom_segment(data=tmp_st, aes(x=start, xend=end , y=cnv_level, yend=cnv_level) , color="red")
        }
        
        tracks_details[[line]] <-  tmp_p + ylim(norm_cov_range) # + coord_cartesian(ylim = norm_cov_range)  
        heights <- append(heights,1)
    }

    if(length(gff_local)>0){
        tracks_details[["genes"]] <- p.genes
        heights <- append(heights,1)
    }
    tks <- tracks(tracks_details, heights = heights ) + scale_x_sequnit(unit = "Mb" )  + xlim(region)  
    tks
    #gff_local
    #regions
}