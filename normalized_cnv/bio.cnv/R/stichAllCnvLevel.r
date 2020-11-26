stichAllCnvLevel <- function(filled, mc.cores=4, lines=unique(filled$line)){
    ret <- mclapply(lines, function(x){
        stichCnvLevel(filled, line=x)
        },  mc.cores = mc.cores) 
	ret <- as(ret, "GRangesList")
    ret
}