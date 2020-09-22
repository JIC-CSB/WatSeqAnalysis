plotCnvsInLine<-function(validated_cnv, 
                                 line="WATDE0821", cnv_level=0){
    
    dels_to_print <- validated_cnv[
        validated_cnv$line == line &
        validated_cnv$cnv_level == cnv_level
    ]
    p<-autoplot(dels_to_print, layout = "karyogram")  
    p <- p + ggtitle(paste0(line, " CNV: ", cnv_level )) 
    #p +scale_color_manual(values = colors) + scale_fill_manual(values=colors) +
    p
}