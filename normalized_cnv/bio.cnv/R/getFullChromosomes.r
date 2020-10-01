getFullChromosomes <- function(db){
    query <- "SELECT chrom, 1  as chromStart,  MAX(chromEnd) as chromEnd FROM Regions GROUP BY chrom "
    res <- dbSendQuery(db, query)
    df  <- dbFetch(res, n = -1)
    makeGRangesFromDataFrame(df,
                    start.field="chromStart",
                    end.field="chromEnd", 
                    keep.extra.columns=T ) 
    # df
}