library(ggplot2)
library(fields)
plotHistogram<-function(table,
                        column="size_cds",
                        probs = c( 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                        binwidth=0.1, trim_range=T, title_prefix ="", force_range=c(), filter_zero = TRUE){
    if(filter_zero){
       table<-table[table[,column]>0,] 
    }
    
    quantiles <- data.frame(quantile(table[,column], prob=probs,na.rm=TRUE, include.lowest=TRUE), stringsAsFactors=FALSE)
    quantiles$quant<-rownames(quantiles)
    colnames(quantiles)<-c("value", "quant")
    values<-quantiles$values
    local_mean<-mean(table[,column],na.rm=TRUE)
    local_sd<-sd(table[,column],na.rm=TRUE)
    local_max <-  max(table[,column],na.rm = T)
    p <- ggplot(table, aes_string(column))
    if(nrow(table) > 100){
        table <- within(table,
           quantile <- as.integer(
               cut(table[,column],
                   unique(quantile(table[,column],
                    prob=probs,
                    na.rm=TRUE,
                    include.lowest=TRUE))
                   )
               ))
        table$quantile<-ifelse(is.na(table$quantile),0,table$quantile)
        table$quantile<-as.factor(table$quantile)
        
        if(length(force_range) == 2){
           xmin <- force_range[1]
           xmax <- force_range[2]
        }else{
            
            iq <- quantiles$value[4] - quantiles$value[2]
            xmax <- quantiles$value[3] + (iq * 2)
            xmin <- quantiles$value[3] - (iq * 2)
            if(xmin < 0) xmin <- 0
            if(xmax > local_max)  xmax <- local_max + 1
        }
        
        p <- ggplot(table, aes_string(column, fill="quantile"))
        p <- p + geom_vline(data=quantiles,aes(xintercept=quantiles$value) )
        for(i in seq(1,nrow(quantiles))){
            x_pos<-quantiles$value[i]
            gtext <- textGrob(quantiles$quant[i], y=0.02,  gp = gpar(fontsize = 6,col = "red"))
            p <- p + annotation_custom(gtext, xmin=x_pos, xmax=x_pos)
        }
       
        p <- p +  scale_fill_brewer(palette="Dark2")
        if(length(force_range) == 2){
           xmin <- force_range[1]
           xmax <- force_range[2]
        }
        if(trim_range){
           p <- p  + xlim(xmin, xmax) 
        }    
    }
    p <- p + geom_histogram(binwidth=binwidth, position = "identity") + theme_bw()
    p <- p + theme(legend.position="none")
    p <- p + ggtitle(paste0(title_prefix, "Mean: ", round(local_mean,2),
        " SD:", round(local_sd,2),
        " CV:", round(local_sd/local_mean, 2),
        " Median:", round(median(table[,column],na.rm=TRUE),2),
        " Max:", round(local_max,2),
        " N:", nrow(table)))
    p <- p + theme(plot.title = element_text(size=6))
    p
}