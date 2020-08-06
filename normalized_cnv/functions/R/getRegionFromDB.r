getRegionFromDB <- function(covs_db){
    query <- "SELECT 
chrom as Scaffold, 
chromStart as Start, 
chromEnd as Ends, 
reg_id as Exon, 
chromEnd - chromStart as ExonL, 
sd as sdExon 
FROM 
regions "
    res <- dbSendQuery(covs_db,query)
    df  <- dbFetch(res, n = -1)
    df 
}

