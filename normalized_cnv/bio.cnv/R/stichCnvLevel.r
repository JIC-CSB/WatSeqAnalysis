stichCnvLevel <- function(filled, line = "WATDE0456"){
    filled_sorted <- filled[filled$line == line]
    filled_sorted<- sort(filled_sorted)
    last <- filled_sorted[1,]  
    ret <- filled_sorted[0,]
    current <- NULL 
    i <- 0

    missing_cnv_level <- filled_sorted[is.na( filled_sorted$cnv_level) ]
    filled_sorted <- filled_sorted[!is.na( filled_sorted$cnv_level) ]
  
    if(length(missing_cnv_level)){
        print("MISSING:")
        print(missing_cnv_level)
    }    
#tryCatch({
    for(i in seq(1, length(filled_sorted))){
        current = filled_sorted[i,]
        if(as.character(seqnames(current)[1]) != as.character(seqnames(last)[1]) || current$cnv_level != last$cnv_level){
            ret[length(ret)+1,] <- last
            last <- current
        }else{
            end(last) <- end(current)
        }
    }
#},error=function(e){
#    print(i)
#    print(current)
#    print(last)
#  })
  ret[nrow(ret)+1,] <- last
  ret$norm_cov <- NULL
  ret
}