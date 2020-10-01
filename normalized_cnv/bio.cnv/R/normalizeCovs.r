normalizeCovs<-function(counts, exonsDF, remove.zero = F) {
  
  #Filtering the rows with AVG=0 
  rowsToKeep<-rowSums(counts)>0
  counts<-counts[rowsToKeep,]
  exonsDF <- exonsDF[rowsToKeep,]

  #counts <- as.matrix(counts)
  if(remove.zero){
       counts<- apply(counts, 2 , function(x){ifelse(x == 0, NA, x)}) 
  }
  
  exonLengths<-exonsDF$ExonL
  totalReadsPerSample<-apply(counts,2,sum, na.rm=T)
  multiplier<-outer(exonLengths, totalReadsPerSample )
  multiplier<-1/multiplier
  multiplier<-1000000000*multiplier
  #"RPKM-like covs"
  multiplier<-counts*multiplier
  
  #Normalize by multiplying by the mean of each row.  
  mean.mult <- apply(multiplier, 1, mean, na.rm = T)
  covs<-sweep(multiplier,MARGIN=1,mean.mult,'/')
  covs <- ifelse(is.na(covs), 0, covs)
  if(remove.zero){
    covs<- apply(covs, 2 , function(x){ifelse(is.na(x), 0, x)}) 
  }
  covs
}