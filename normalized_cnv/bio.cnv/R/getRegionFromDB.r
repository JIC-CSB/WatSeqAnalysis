getRegionFromDB <- function(covs_db, rename=TRUE, region = null ){
    query <- "SELECT 
*
FROM 
regions "

    if(!is.null(region)){
        chr = seqnames(region)[1]
        first = start(region)[1]
        last = end(region)[1] 
        query <- paste0(query, " WHERE chrom = '", chr, "' and chromStart BETWEEN ", first, " AND ",last )       
    }
    res <- dbSendQuery(covs_db,query)
    df  <- dbFetch(res, n = -1)
    if(rename){
        colnames(df) <- c("chrom", "start", "end", "reg_id", "sd")
    }
    df 
}

