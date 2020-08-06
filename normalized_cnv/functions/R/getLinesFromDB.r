getLinesFromDB<-function(convs_db){
    query <- "SELECT * from  lines"
    res   <- dbSendQuery(covs_db,query)
    dbFetch(res, n = -1)
}