getRegionFromDB <- function(covs_db, rename=TRUE, region = null ){
    query <- "SELECT 
chrom, 
chromStart as start, 
chromEnd as end, 
reg_id, 
chromEnd - chromStart as length, 
sd as sd_region 
FROM 
regions "

    if(!is.null(region)){
        chr = seqnames(region)[1]
        first = start(region)[1]
        last = end(region)[1] 
        query <- paste0(query, "WHERE chrom = '", chr, "' and start BETWEEN ", first, " AND ",last )       
    }
    res <- dbSendQuery(covs_db,query)
    df  <- dbFetch(res, n = -1)
    if(rename){
        colnames(df) <- c("Scaffold", "Start", "Ends", "Exon", "ExonL", "sdExon")
    }
    df 
}

