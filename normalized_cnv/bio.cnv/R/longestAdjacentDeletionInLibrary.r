longestAdjacentDeletionInLibrary <- function(library, tempDf, libValues, minForDelLib, maxForDel=0){
	tempExons <- rownames(tempDf)
	longestStretch<-0
	currentStretch<-0
	tempDf$minForDel <- 1 - ( 3 * tempDf$StdDev)
	prevDeletion<-F
	for (i in tempExons) {
		val<-libValues[i]
		minForDel <- tempDf[i,"minForDel"]
		if( minForDel < 0 ){
			next
		}
		if(val < minForDel && val < minForDelLib && val < maxForDel ){
			currentStretch <- currentStretch + 1
		}else{
			currentStretch <- 0
		}
		if(currentStretch > longestStretch){
			longestStretch<-currentStretch 
		}
			
	}
	longestStretch
}