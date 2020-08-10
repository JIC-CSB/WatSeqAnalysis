getWindowsInRange <- function(db, chr="chr3A_part2", start=192557156 , end=192583886, ranges = null){

    if(!is.null(ranges)){
        chr=seqnames(ranges)[1]
        start=start(ranges)[1]
        end=end(ranges)[1]

    }
    inner_query <- paste0("SELECT *
	FROM REGIONS 
	WHERE
	regions.chrom = '", chr,"' 
	and ( regions.chromStart BETWEEN ", start," and ", end, "  OR
	regions.chromEnd BETWEEN  ", start," and ", end, ")")
    
    query <- paste0("SELECT chrom, chromStart, chromEnd, line, cov, norm_cov FROM 
(
SELECT * FROM covs 
WHERE covs.reg_id 
in (
	SELECT reg_id FROM (", inner_query,")
	
 ) 
) as  n_cov
JOIN 
lines on n_cov.line_id = lines.line_id
JOIN (
    ", inner_query,"
)as  regs ON regs.reg_id = n_cov.reg_id
ORDER BY chrom, chromStart, chromEnd, line
;")
    res <- dbSendQuery(db,query)
    df  <- dbFetch(res, n = -1)
    if(!is.null(ranges)){
        df <- makeGRangesFromDataFrame(df,
                         start.field="chromStart",
                         end.field="chromEnd", 
                         keep.extra.columns=T )
    }
    df 
}