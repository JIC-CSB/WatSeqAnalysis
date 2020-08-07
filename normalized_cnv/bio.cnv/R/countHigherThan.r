countHigherThan <- function(x, min_cov=0){
    length( x[  x >= min_cov] )
}