plotLibraryOnChromosomeInterval <- function(covs_df,  lines=c("WATDE0027") , title="Normalised coverage", col = "cov" ){
    chrom <- unique(covs_df$chrom)
    from <- min(covs_df$chromStart)
    to <- max(covs_df$chromEnd)
	plot_Data <- ddply(covs_df, .(chromStart), mutate, Q1=quantile(cov, 1/4), Q3=quantile(cov, 3/4), IQR=Q3-Q1, upper.limit=Q3+1.5*IQR, lower.limit=Q1-1.5*IQR)
    plot_Data$position <- factor(plot_Data$chromStart)
    lines_plot_data<-subset(plot_Data, line %in% lines )
	gg <- ggplot(plot_Data, aes_string(x="position", y=col))
	gg <- gg + geom_boxplot(outlier.shape = ".")
	gg <- gg + geom_line(data=lines_plot_data,  aes_string(x="position", y=col,group="line", colour="line"))
	 gg <- gg + ggtitle(paste(title,"\n", chrom,":",from, "-", to))
     gg <- gg + theme_bw()
     gg <- gg + theme(axis.text.x = element_text(angle = 90, hjust = 1))
	 gg
    #head(plot_Data)
}
