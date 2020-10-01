#' Prepares the dataframe with the information per exon. This parses the rownames to get the scaffold
#' Start and end. It also calculates the exon length and the starad deviation of each exon. 
#' @param mat: The matrix with the data
#' @return A data frame with the following rows: Exon, Scaffold, Start, End, ExonL (exon length)
getExonsDF<-function(mat){
	names<-rownames(mat)
	exons<-strsplit(as.character(names), ":")
    df_names= as.data.frame(t(as.data.frame(exons,stringsAsFactors = F)),stringsAsFactors = F)
    colnames(df_names) <- c("Scaffold", "Start", "Ends")
    df_names$Exon<-names
    df_names$Start<-as.integer(df_names$Start)
    df_names$Ends<-as.integer(df_names$Ends)
	df_names$ExonL<-df_names$Ends-df_names$Start
    df_names$sdExon<-apply(mat,1,sd)
    df_names
}