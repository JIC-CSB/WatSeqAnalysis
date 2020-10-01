

fixRangesOrder<-function(gr,lens, ordr = c("chr1A","chr1B", "chr1D","chr2A","chr2B", "chr2D","chr3A","chr3B", "chr3D","chr4A","chr4B", "chr4D","chr5A","chr5B", "chr5D","chr6A","chr6B", "chr6D","chr7A","chr7B", "chr7D","chrUn")){
    as_df <-data.frame(gr)
    as_df$seqnames <- factor(as.character(as_df$seqnames),levels = ordr )
    vect = lens$length
    names(vect) <-lens$seqname
    gr2 <- makeGRangesFromDataFrame(as_df,keep.extra.columns = T, seqinfo=vect)
    trim(gr2)
}
deleted_bases_per_line<-function(dels){
    validated<-dels$validated
    validated$width<-width(validated)
    total_len <- sum(seqlengths(validated))
    w_hom <- vector()
    w_het <- vector()
    for(l in unique(validated$line)){
        v_hom <- validated[validated$line == l & validated$type == "Hom"]
        v_het <- validated[validated$line == l & validated$type == "Het"] 
        w_hom[l] <- sum(v_hom$width)
        w_het[l] <- sum(v_het$width)
        
    }
    df<-data.frame(line=names(w_hom), bases_hom = w_hom, bases_het=w_het, bases_tot=w_het+w_hom)
    df$per_hom <-100.0 * df$bases_hom / total_len
    df$per_het <-100.0 * df$bases_het / total_len
    df$per_tot <-100.0 * df$bases_tot / total_len
    df
    
}

plotHistogram<-function(table,
                        column="size_cds",
                        probs = c( 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                        binwidth=0.1, trim_range=T){
    table<-table[table[,column]>0,]
    quantiles <- data.frame(quantile(table[,column], prob=probs,na.rm=TRUE, include.lowest=TRUE), stringsAsFactors=FALSE)
    quantiles$quant<-rownames(quantiles)
    colnames(quantiles)<-c("value", "quant")
    values<-quantiles$values
    local_mean<-mean(table[,column])
    local_sd<-sd(table[,column])
    local_max <-  max(table[,column])
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

       iq <- quantiles$value[4] - quantiles$value[2]

       xmax <- quantiles$value[3] + (iq * 2)
       xmin <- quantiles$value[3] - (iq * 2)
       if(xmin < 0) xmin <- 0
       if(xmax > local_max)  xmax <- local_max + 1

       p <- ggplot(table, aes_string(column, fill="quantile"))
       p <- p + geom_vline(data=quantiles,aes(xintercept=quantiles$value) )
       for(i in seq(1,nrow(quantiles))){
        x_pos<-quantiles$value[i]
        gtext <- textGrob(quantiles$quant[i], y=0.02,  gp = gpar(fontsize = 6,col = "red"))
        p <- p + annotation_custom(gtext, xmin=x_pos, xmax=x_pos)
       }
           p <- p +  scale_fill_brewer(palette="Dark2")
       if(trim_range){
           p <- p  + xlim(xmin, xmax) 
       }    
    }




    p <- p + geom_histogram(binwidth=binwidth, position = "identity") + theme_bw()
    p <- p + theme(legend.position="none")
    p <- p + ggtitle(paste0("Mean: ", round(local_mean,2),
        " SD:", round(local_sd,2),
        " CV:", round(local_sd/local_mean, 2),
        " Median:", round(median(table[,column]),2),
        " Max:", round(local_max,2),
        " N:", nrow(table)))
    p <- p + theme(plot.title = element_text(size=6))

    stats_list<-list(mean=local_mean, sd = local_sd, cv = local_sd/local_mean,
        median =  median(table[,column],2),
        max = local_max,
        n = nrow(table)
        )
    p
}

plot_deletions_in_line<-function(deletions, 
                                lens=NULL,
                                line="J1.26_k80d50Mc5", 
                                chr_order = c("chr1A","chr1B", "chr1D","chr2A","chr2B", "chr2D","chr3A","chr3B", "chr3D","chr4A","chr4B", "chr4D","chr5A","chr5B", "chr5D","chr6A","chr6B", "chr6D","chr7A","chr7B", "chr7D","chrUn"), 
                                colors=    c("Hom"="#8c510a",
                                             "Het"="#01665e")
                                             ){
    
    select_line <- deletions$validated$line == line
    dels_to_print <- deletions$validated[ select_line ]
    dels_to_print<-data.frame(dels_to_print)
    dels_to_print$seqnames <- factor(as.character(dels_to_print$seqnames),levels = chr_order )
    lens2<-lens$length
    names(lens2) <- lens$seqname
    dels_to_print<-makeGRangesFromDataFrame(dels_to_print,keep.extra.columns = T, seqinfo=lens2)
    seqlengths(dels_to_print) <- lens[levels(seqnames(dels_to_print)),"length"]
    p<-autoplot(dels_to_print, layout = "karyogram", aes(color=type, fill=type)) + scale_color_manual(values = colors) + scale_fill_manual(values=colors) + ggtitle(line)
    p
}


plotDeletionDensity<-function(dels, type="Hom",tilewidth=100000, base=2){
    v2 <- dels[dels$type==type]
    cov_v2 <- coverage(x = v2)
    bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(v2),
                                  tilewidth=tilewidth,
                                  cut.last.tile.in.chrom=TRUE)
    gr_bins <- binnedAverage(bins, cov_v2, "binned_cov")
    gr_bins <-fixRangesOrder(gr_bins, lens=lens)
    gr_bins$log_cov <- ceiling(log(gr_bins$binned_cov+1,base=base))
    gr_bins$log_cov <- ifelse(gr_bins$log_cov < 1 ,0, gr_bins$log_cov )
    cov_ranges<-list()
    for(r in sort(unique(gr_bins$log_cov))){
        tmp_gr = gr_bins[gr_bins$log_cov == r]
        cov_ranges[as.character(r)] = paste0(round(min(tmp_gr$binned_cov),2),"-", round(max(tmp_gr$binned_cov),2))
    }
    gr_bins$cov_legend <- factor(as.character(cov_ranges[as.character(gr_bins$log_cov)]), levels=cov_ranges)
    title <- paste0(type, "  per bin size:", format(tilewidth, scientific=FALSE, big.mark=","), " Base: ", base)
    autoplot(gr_bins,  layout = "karyogram", aes(fill=cov_legend)) + ggtitle(title) + scale_fill_brewer(type="seq", palette = "Oranges") 
}

all_deletions_per_gene<-function(gff, deletions,gff_type="gene", id_col="gene_id"){
    chr_range<-gff[which(gff$type==gff_type)  ]
    hits     <- findOverlaps(deletions,chr_range)
    chr_range_df <- data.frame(chr_range)
    genes<-chr_range_df[subjectHits(hits), id_col]
    dels_df <- data.frame(deletions)
    dels_df <- dels_df[queryHits(hits), c("line","seqnames", "type")]
    cbind(dels_df, genes)
}

deletions_per_line<-function(dels){
    deletions<-dels$validated 
    deletions<-data.frame(deletions)
    deletions$grouped_windows <- NULL
    deletions$grouping <-NULL
    del_per_line <- sqldf("SELECT line, count(*) as no_del FROM deletions where type = 'Hom' GROUP BY line")
    quantiles <- data.frame(quantile(del_per_line[,"no_del"],probs = c( 0.1, 0.25, 0.5, 0.75, 0.9, 0.95),na.rm=TRUE, include.lowest=TRUE), stringsAsFactors=FALSE)
    quantiles$quant<-rownames(quantiles)
    colnames(quantiles)<-c("value", "quant")
    iq <- quantiles$value[4] - quantiles$value[2]
    m  <- mean(del_per_line$no_del)
    s_d <- sd(del_per_line$no_del)
    xmax <- quantiles$value[3] + (iq * 2)
    xmin <- quantiles$value[3] - (iq * 2)
    del_per_line$outlier <- del_per_line$no_del > xmax | del_per_line$no_del < xmin
    del_per_line$sample<- gsub(pattern = "_k80d50Mc5", replacement = "", x = del_per_line$line)
    del_per_line$sample <- gsub(pattern = ".", replacement="-", x=del_per_line$sample, fixed=TRUE)
    del_per_line
}

merge_deletions_with_metadata <- function(del_per_line){
    summary_ei <-read.csv(file = "../full_Set//ENQ-2360_PIP-1962_Skim-Seq_datasum_All Plates_31.07.18.csv")
    summary_ei$well <- gsub(pattern = "PRO1896_", replacement = "", x = summary_ei$Sample_ID)
    well_mapping <- read.csv( "../full_Set/well_mapping.csv")
    summary<-sqldf("SELECT summary_ei.*, well_mapping.sample FROM summary_ei LEFT JOIN well_mapping ON summary_ei.well = well_mapping.well")
    missing_samples <- sqldf("SELECT * from well_mapping WHERE well_mapping.well NOT IN (SELECT well from summary)")
    missing_samples <- sqldf("SELECT * from summary WHERE summary.well NOT IN (SELECT well from well_mapping )")
    cov_sum <- read.csv("../full_Set/coverage_summary.tsv", sep="\t")
    cov_sum <- cov_sum[cov_sum$Chromosome == "total",]
    summary_with_cov <- sqldf("SELECT summary.*, cov_sum.MappedBases, cov_sum.MappedReads, cov_sum.UnmappedReads, cov_sum.Coverage, cov_sum.MappedPercentage FROM summary LEFT JOIN cov_sum ON cov_sum.Sample = summary.well")
    head(summary_with_cov)
    names <- colnames(summary_with_cov)
    names[12] <- "LibQC"
    colnames(summary_with_cov) <- names 
    summary_with_cov
    
    del_per_line <- merge(del_per_line, summary_with_cov, by="sample")
    del_per_line
}

deleted_regions_categories<-function(dels, lens=lens, type="Hom", tilewidth=100000){
    v_dels <- dels$validated
    v1 <- v_dels[v_dels$type==type]
    
    del_per_line<-deletions_per_line(dels)
    outlier_lines <- del_per_line[del_per_line$outlier, ]


    
    bins <- GenomicRanges::tileGenome(GenomeInfoDb::seqlengths(v1),
                                  tilewidth=tilewidth,
                                  cut.last.tile.in.chrom=TRUE)
    
    v2 <- v1[!(v1$line  %in% outlier_lines$line)]
    cov_v2 <- coverage(x = v2)
    gr_bins <- binnedAverage(bins, cov_v2, "binned_cov")
    gr_bins <- fixRangesOrder(gr_bins, lens=lens)
    gr_bins$deleted <- ifelse(gr_bins$binned_cov > 0,"Deleted", "Non-deleted")
    
    
    v3 <- v1[v1$line %in% outlier_lines$line ]
    cov_v3 <- coverage(x = v3)
    gr_bins3 <- binnedAverage(bins, cov_v3, "binned_cov")
    gr_bins3 <- fixRangesOrder(gr_bins3, lens=lens)
    gr_bins3$deleted <- ifelse(gr_bins3$binned_cov > 0,"Deleted", "Non-deleted")
    
    gr_filter <- gr_bins[gr_bins$deleted == "Deleted"]
    gr_filter2 <- gr_bins3[gr_bins3$deleted == "Deleted"]
    gr_filter2 <- append(gr_filter2,gr_filter2)
    gr_filter_mix <- append(gr_filter, gr_filter2)
        
    cov_v4 <- coverage(x = gr_filter_mix)

    gr_bins2 <- binnedAverage(bins, cov_v4, "binned_cov")
    gr_bins2 <- fixRangesOrder(gr_bins2, lens=lens)
    gr_bins2$deleted <- ifelse(gr_bins2$binned_cov == 0,"Non-Deleted", 
                               ifelse(gr_bins2$binned_cov == 1, "Normal deletion line",
                                      ifelse(gr_bins2$binned_cov == 2,"Only in outlier line",
                                            "Outlier and normal deletion line")))
    
    gr_bins2$deleted <- factor(gr_bins2$deleted,
                levels = c("Non-Deleted", "Only in outlier line", "Outlier and normal deletion line","Normal deletion line" ))
 
     gr_bins2   
}


plotDeletedVsNoDeleted<-function(dels, lens=lens, type="Hom",tilewidth=100000){
    
    gr_bins2<-deleted_regions_categories(dels, lens=lens, type=type, tilewidth=tilewidth)
    
    title <- paste0(type, "  per bin size:", format(tilewidth, scientific=FALSE, big.mark=","))
    cols <- c("Non-Deleted" = "#e66101", "Only in outlier line" = "#fdb863", "Outlier and normal deletion line" = "#b2abd2" ,"Normal deletion line" = "#5e3c99" )
    autoplot(gr_bins2,  layout = "karyogram", aes(fill=deleted)) + ggtitle(title) + scale_colour_manual(values = cols) + scale_fill_manual(values = cols) 
}

doStatsForDataset <- function(path = "./Tables/all_dels_test_d1000_200k_to_10k_gap_3.rds",
                             seq_len="../WGAv1.0//sequence_lengths.txt", 
                             chr_order=c("chr1A","chr1B", "chr1D","chr2A","chr2B", "chr2D","chr3A","chr3B", "chr3D","chr4A","chr4B", "chr4D","chr5A","chr5B", "chr5D","chr6A","chr6B", "chr6D","chr7A","chr7B", "chr7D","chrUn"),
                             gff_path="../WGAv1.0//Triticum_aestivum.IWGSC.43.chr.gff3.gz", 
                              output_folder="./summary/dels_d100_200k_to_10k_gap3"){
    dels<-readRDS(path)
    lens <- read.csv(seq_len, sep="\t", stringsAsFactors=F)
    lens$seqname<-factor(lens$seqname, levels = chr_order)
    rownames(lens) <- lens$seqname
    
    gff_local <- import(gff_path)
    seqlengths(gff_local) <- lens[levels(seqnames(gff_local)),"length"]
    
    dels$validated <- fixRangesOrder(dels$validated , lens)
    dels$probs <- fixRangesOrder(dels$probs , lens)
    dels$dels <- fixRangesOrder(dels$dels , lens)
    
    deleted_bases<-deleted_bases_per_line(dels)
    dels_per_genes <- all_deletions_per_gene(gff_local, dels$validated)
    
    deleted_genes_per_line <- aggregate(dels_per_genes$genes, by=list(dels_per_genes$line), FUN=length)
    colnames(deleted_genes_per_line) <- c("line", "deleted_genes")
    deletion_lines_per_gene <- aggregate(dels_per_genes$line, by=list(dels_per_genes$genes), FUN=length)
    colnames(deletion_lines_per_gene) <- c("gene","deleted_in_lines")
    
    deletions<-dels$validated 
    deletions<-data.frame(deletions)
    deletions$grouped_windows <- NULL
    deletions$grouping <-NULL

    dels_per_line <- deletions_per_line(dels)
    dels_per_line <- merge_deletions_with_metadata(dels_per_line)
    
    dir.create(output_folder,recursive=T)
    wd <- getwd()
    setwd(output_folder)
    
    v <- dels$validated
    v$grouping <- NULL
    
    write.csv(deleted_genes_per_line, row.names=FALSE, file="deleted_genes_per_line.csv")
    write.csv(deletion_lines_per_gene,row.names=FALSE, file="deletion_lines_per_gene.csv")
    write.csv(deleted_bases,         row.names=FALSE, file="deleted_bases.csv")
    write.csv(v,         row.names=FALSE, file="validated_deletions.csv")
    write.csv(dels_per_line, row.names=FALSE, file="deletions_per_line.csv")
    
    pdf("all_plots.pdf",onefile = TRUE)
    
    print(plotHistogram(deleted_genes_per_line, column="deleted_genes", binwidth = 25, trim_range=F))
    print(plotHistogram(deleted_bases, column="per_hom", binwidth = 0.05, trim_range=F))
    print(plotHistogram(deletion_lines_per_gene, column="deleted_in_lines", binwidth = 1,trim_range = T))
    print(plotHistogram(deleted_genes_per_line, column="deleted_genes", binwidth = 25, trim_range=T))
    print(plotHistogram(deleted_genes_per_line, column="deleted_genes", binwidth = 25, trim_range=T))
    
    print(plotHistogram(dels_per_line,column = "no_del",binwidth = 5,trim_range = F))
    print(plotHistogram(dels_per_line,column = "no_del",binwidth = 1,trim_range = T))
    
    
    cols2 <- c("FALSE" = "#5e3c99", "TRUE" = "#e66101" )
    print(ggplot(dels_per_line, aes(x=no_del, y=Coverage, color=outlier)) + geom_point() + scale_colour_manual(values=cols2))
    print(ggplot(dels_per_line, aes(x=MappedPercentage, y=no_del, color=outlier)) + geom_point()+ scale_colour_manual(values=cols2))
    print(ggplot(dels_per_line, aes(x=MappedPercentage, y=Coverage, color=outlier)) + geom_point()+ scale_colour_manual(values=cols2))

    print(plotDeletedVsNoDeleted(dels, lens=lens,tilewidth=100000))
    print(plotDeletedVsNoDeleted(dels, lens=lens,tilewidth=1000000))
    
    print(plotDeletionDensity(dels$validated,tilewidth=100000, base=10))
    print(plotDeletionDensity(dels$validated,tilewidth=1000000,base=10))
    print(plotDeletionDensity(dels$validated,tilewidth=100000, base=2))
    print(plotDeletionDensity(dels$validated,tilewidth=1000000,base=2))
        
    for(l in sort(as.character(unique(dels$validated$line)))){
       print( plot_deletions_in_line(dels, line = l))
    }
    dev.off()
    setwd(wd)
}

# Function to extract legend
g_legend <- function(a.gplot){ 
    tmp <- ggplot_gtable(ggplot_build(a.gplot)) 
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
    legend <- tmp$grobs[[leg]] 
    legend
} 
                        
                     
plot_deletions_in_region <- function(deletions, 
                                     chr="chr4B", 
                                     start=58000000,
                                     end= 88000000,
                                     line="J7.28_k80d50Mc5", 
                                     gff=NULL,
                                     genes=c(),
                                     plot_validated=TRUE,
                                    window_sizes=c("200k", "100k", "075k","050k", "025k"),
                                    colors=c("200k"="#1b9e77FF",
                                             "100k"="#d95f02FF",
                                             "075k"="#7570b3FF",
                                             "050k"="#e7298aFF",
                                             "025k"="#66a61eFF",
                                             "010k"="#e6ab02FF",
                                             "Hom"="#8c510a",
                                             "Het"="#01665e", 
                                             "TRUE"="#cc8131",
                                             "FALSE"="#cdE5e5"
                                            )){
    validated<-trim(deletions$validated)
    validated <- validated[validated$line==line & seqnames(validated) == chr  ]
    
    
    dels<-trim(deletions$dels)
    dels <- dels[dels$line==line & seqnames(dels) == chr & dels$window %in% window_sizes ]
    
    validated_df <- data.frame(validated)
    validated_df$id = paste0(validated_df$seqnames,":", validated_df$start, "-", validated_df$end)
    validated_df$grouping<-NULL
    
    probs <- trim(deletions$probs)
    probs <- probs[probs$line==line & seqnames(probs) == chr & probs$window %in% window_sizes]
    probs_df<-data.frame(probs)
    
    scale_combined <- GRanges(chr, IRanges(start = start, end = end))
    
    p3 <- ggplot(probs , aes(y=norm_coverage, x=start,  color=window)) + geom_point( ) 
    p3 <- p3 + scale_color_manual(values = colors) 
    #p3 <- p3 + geom_segment(data=validated_df,aes(x=start, y=mean_cov,xend=end, yend=mean_cov), colour="black" ) 
    p3 <- p3 +  theme(legend.position='none') + ylim(0,2) 

    p4 <- autoplot(dels, aes(fill=window), colour = "#00000000")+ scale_fill_manual(values = colors)
    p4 <- p4 + theme(legend.position="bottom",legend.title = element_blank())
    
    alpha <-0
    if(plot_validated){
        alpha <- 1
    }
    p5 <- autoplot(validated, aes(fill=type), alpha=alpha)+ scale_fill_manual(values = colors) 
    
    
    
    p6 <- NULL
    
    t1 <-tracks(cov=p3, 
                del=p4, 
                val=p5, 
                heights = c(0.5,0.3,0.2), 
                title=line,
                xlim=scale_combined)
    if(!is.null(gff)){
        
        tmp_gff <- gff[gff$type == "gene"]
        tmp_gff<- subsetByOverlaps(tmp_gff, scale_combined)
        tmp_gff$select <- tmp_gff$gene_id %in% genes
        
        p6 <- autoplot(tmp_gff, aes(fill=select, color = select)) + scale_fill_manual(values = colors)+ scale_color_manual(values = colors)
        t1 <-tracks( 
                cov=p3, 
                del=p4, 
                val=p5,
                genes=p6,
                heights = c(0.5, 0.2,0.1,0.2), 
                title=line,
                xlim=scale_combined)
    }
    
    
    t1 <- t1 +
    scale_x_sequnit() + 
    theme_bw()              
    t1 
 }

 readDelsRDS<-function(path="./Tables/all_dels_test_d1000_200k_to_10k_gap_3.rds", lens=lens){
    dels<-readRDS(path)
    lens <- read.csv(lens, sep="\t", stringsAsFactors=F)
    dels$validated <- fixRangesOrder(dels$validated , lens)
    dels$probs <- fixRangesOrder(dels$probs , lens)
    dels$dels <- fixRangesOrder(dels$dels , lens)
    dels
}

merge_deletions<-function(left, right, name_left, name_right){
    lines <- unique(append(left$line, right$line))
    glst <- GRanges()
    for(l in lines){
        l_l <- left[left$line == l]
        r_l <- right[right$line == l]
        i_l <- intersect(l_l, r_l)
        
        ol_l <- setdiff(l_l, i_l)
        or_l <- setdiff(r_l, i_l)
        
        if(length(i_l) > 0){
            i_l$line = l
            i_l$set_name <- paste0(name_left, " ^ ", name_right)
        }
        
        if(length(ol_l) > 0){
            ol_l$line = l
            ol_l$set_name <- name_left
        }
       
        if(length(or_l) > 0){
            or_l$line = l
            or_l$set_name <- name_right
        }
        
        glst<-append(glst, i_l)
        glst<-append(glst, ol_l)
        glst<-append(glst, or_l)
    }
    glst
}

plot_merged_deletions<-function(merged_dels, line="J4.72_k80d50Mc5", colors=c("d1000_g3"="#a6611a", "d0200_g3"="#018571", "d1000_g3 ^ d0200_g3" = "#d01c8b", "d1000_g3 U d0200_g3" = "#377eb8"), title=line){
    to_plot <- merged_dels[merged_dels$line == line]
    autoplot(to_plot,  layout = "karyogram", aes(fill=set_name)) + scale_fill_manual(values = colors) + scale_colour_manual(values = colors) + ggtitle(title)
}


extend_merged_deletions<-function(merged_deletions, left, right, name_left, name_right){
    lines <- sort(unique(merged_deletions$line))
    name_intersect<-paste0(name_left, " ^ ", name_right)
    name_union<-paste0(name_left, " U ", name_right)
    gr_intersect <- merged_deletions[merged_deletions$set_name == name_intersect]
    glst <- GRanges()
    for(l in lines){
           
        l_l <- left[left$line == l]
        r_l <- right[right$line == l]
        i_l <- gr_intersect[gr_intersect$line == l]
        hits_l <- findOverlaps(i_l, l_l, type="any")
        hits_r <- findOverlaps(i_l, r_l, type="any")
        
        to_join <- union(l_l[subjectHits(hits_l)], r_l[subjectHits(hits_r)])
        if(length(to_join) > 0){
            to_join$line = l
            to_join$set_name = name_union
        }
        
        
        l_diff <- setdiff(l_l, to_join)
        if(length(l_diff) > 0){
          l_diff$line = l
            l_diff$set_name = name_left  
        }
        
        r_diff <- setdiff(r_l, to_join)
        if(length(r_diff) > 0){
            r_diff$line = l
            r_diff$set_name = name_right
        }
        
        glst<-append(glst, to_join)
        glst<-append(glst, l_diff)
        glst<-append(glst, r_diff)
        

    }
    glst
}