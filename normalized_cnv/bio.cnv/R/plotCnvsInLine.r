plotCnvsInLine<-function(validated_cnv, 
                         line="WATDE0821", cnv_level=NULL, 
                         min_length =10000, remove_base = TRUE){
  title <- line
  dels_to_print <- validated_cnv[
    validated_cnv$line == line 
  ]
  if(!is.null(cnv_level)){
    dels_to_print <- dels_to_print[dels_to_print$cnv_level == cnv_level]
    title <- paste0(title, " CNV : ", cnv_level)
  }
  if(!is.null(min_length)){
    dels_to_print <- dels_to_print[width(dels_to_print)  > min_length]
    title <- paste0(title, " Min length: ", format(min_length, big.mark=",") )
  }
  if(remove_base){
    dels_to_print <- dels_to_print[dels_to_print$cnv_level != 1 ]
  }
  
  dels_to_print$CNV <- ifelse( dels_to_print$cnv_level >5 , "6+" , as.character(dels_to_print$cnv_level)  )

  colors = c ("0" = "#9e0142", 
              "1" = "#ffffff",
              "2" = "#f46d43",
              "3" = "#abdda4",
              "4" = "#66c2a5",
              "5" = "#3288bd",
              "6+"= "5e4fa2")
              
  p <- autoplot(dels_to_print, layout = "karyogram", aes(fill=CNV))  
  p <- p + ggtitle(title)
  p <- p + scale_color_manual(values = colors) + scale_fill_manual(values=colors)
  p 
}