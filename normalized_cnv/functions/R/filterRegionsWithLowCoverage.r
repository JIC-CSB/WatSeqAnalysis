filterRegionsWithLowCoverage <- function(mat, min_lines=0.0, min_cov=0){
  counts<- apply(mat, 1, countHigherThan, min_cov = min_cov )
    pass_filter <- counts > min_lines * ncol(mat)
    mat[pass_filter, ]
}