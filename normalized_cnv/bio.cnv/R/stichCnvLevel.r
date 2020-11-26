stichCnvLevel <- function(filled, line = "WATDE0456"){
    filled_sorted <- filled[filled$line == line]
    filled_sorted<- sort(filled_sorted)
    last <- filled_sorted[1,]  
    ret <- filled_sorted[0,]
    for(i in seq(1, length(filled_sorted))){
        current = filled_sorted[i,]
        if(as.character(seqnames(current)[1]) != as.character(seqnames(last)[1]) || current$cnv_level != last$cnv_level){
            ret[length(ret)+1,] <- last
            last <- current
        }else{
            end(last) <- end(current)
        }
    }
    ret[nrow(ret)+1,] <- last
    ret$norm_cov <- NULL
    ret
}