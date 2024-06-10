mylodint <- function (results, chr, qtl.index, drop = 1.5, lodcolumn = 1, expandtomarkers = FALSE) 
### made a minor change to lodint function. Should perhaps put this to the authors
{
    if (!("scanone" %in% class(results))) {
        if (!("qtl" %in% class(results))) 
            stop("Input must have class \"scanone\" or \"qtl\".")
        else {
            if (!("lodprofile" %in% names(attributes(results)))) 
                stop("qtl object needs to be produced by refineqtl with keeplodprofile=TRUE.")
            else {
                if (lodcolumn != 1) {
                  warning("lod column ignored if input is a qtl object.")
                  lodcolumn <- 1
                }
                results <- attr(results, "lodprofile")
                if (missing(qtl.index)) {
                  if (length(results) == 1) 
                    results <- results[[1]]
                  else stop("You must specify qtl.index.")
                }
                else {
                  if (qtl.index < 1 || qtl.index > length(results)) 
                    stop("qtl.index misspecified.")
                  results <- results[[qtl.index]]
                }
                chr <- results[1, 1]
            }
        }
    }
    else {
        if (lodcolumn < 1 || lodcolumn + 2 > ncol(results)) 
            stop("Argument lodcolumn misspecified.")
        if (missing(chr)) {
            if (length(unique(results[, 1])) > 1) 
                stop("Give a chromosome ID.")
        }
        else {
            if (is.na(match(chr, results[, 1]))) 
                stop("Chromosome misspecified.")
            results <- results[results[, 1] == chr, ]
        }
    }
    if (all(is.na(results[, lodcolumn + 2]))) 
        return(NULL)
    maxlod <- max(results[, lodcolumn + 2], na.rm = TRUE)
    w <- which(!is.na(results[, lodcolumn + 2]) & results[, lodcolumn + 
        2] == maxlod)[1]### if two are identical, this is a problem, so I added[1]
    o <- range(which(!is.na(results[, lodcolumn + 2]) & results[, 
        lodcolumn + 2] > maxlod - drop))
    if (length(o) == 0) 
        o <- c(1, nrow(results))
    else {
        if (o[1] > 1) 
            o[1] <- o[1] - 1
        if (o[2] < nrow(results)) 
            o[2] <- o[2] + 1
    }
    if (expandtomarkers) {
        markerpos <- (1:nrow(results))[-grep("^c.+\\.loc-*[0-9]+(\\.[0-9]+)*$", 
            rownames(results))]
        if (any(markerpos <= o[1])) 
            o[1] <- max(markerpos[markerpos <= o[1]])
        if (any(markerpos >= o[2])) 
            o[2] <- min(markerpos[markerpos >= o[2]])
    }
    rn <- rownames(results)[c(o[1], w, o[2])]
    if (any(table(rn) > 1)) {
        rn[2] <- paste(rn[2], "")
        if (rn[1] == rn[3]) 
            rn[3] <- paste(rn[3], " ")
    }
    results <- results[c(o[1], w, o[2]), ]
    rownames(results) <- rn
    class(results) <- c("scanone", "data.frame")
    results
}
