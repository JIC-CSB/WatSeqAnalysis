### R script which modifies some of the R/qtl functions
### This is mainly about chromosome labelling at the axis
### 14/12/2009 Luzie U. Wingen, JIC
### 1/03/2018 improved mylodint Luzie U. Wingen, JIC

mylodint <- function (results, chr, qtl.index, drop = 1.5, qtlmarker='',lodcolumn = 1, expandtomarkers = FALSE) 
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
    if((!is.na(qtlmarker))&(qtlmarker!=''))
        w <- which(row.names(results)==qtlmarker)
    else
        w <- which(!is.na(results[, lodcolumn + 2]) & results[, lodcolumn + 2] == maxlod)[1]### if two are identical, this is a problem, so I added[1]
    if (is.na(w))
        print('qtl marker not found on chromosome')
    oindex <- which(!is.na(results[,lodcolumn+2]) & results[,lodcolumn+2] > maxlod - drop) ### if qtl goes up again
    gap <- which(diff(oindex)!=1) ### is there an up again
    gap <- c(1,gap,length(oindex))
    for (g in 2:length(gap))
    {
        onew <- oindex[gap[g-1]:gap[g]]
        if (w%in%onew)
            break
    }
###    o <- range(which(!is.na(results[,lodcolumn+2]) & results[,lodcolumn+2] > maxlod - drop))
    o <- range(onew)
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

plot.qtl <- function (x, chr, horizontal = FALSE, shift = TRUE, show.marker.names = FALSE, alternate.chrid = FALSE,labels.rotate=TRUE,  ...) 
###All marker names are shown if put to TRUE. 
{
    if (!("qtl" %in% class(x))) 
        stop("input should be a qtl object")
    if (length(x) == 0) 
        stop("  There are no QTL to plot.")
    map <- attr(x, "map")
    if (is.null(map)) 
        stop("qtl object doesn't contain a genetic map.")
    if (missing(chr)) 
        chr <- names(map)
    else {
        chr <- matchchr(chr, names(map))
        map <- map[chr]
        class(map) <- "map"
    }
    if (horizontal) 
        plot.map(map, horizontal = horizontal, shift = shift, 
            show.marker.names = show.marker.names, alternate.chrid = alternate.chrid, 
            ylim = c(length(map) + 0.5, 0), ...)
    else plot.map(map, horizontal = horizontal, shift = shift, 
        show.marker.names = show.marker.names, alternate.chrid = alternate.chrid, 
        xlim = c(0.5, length(map) + 1), ...)
    whchr <- match(x$chr, names(map))
    thepos <- x$pos
    thepos[is.na(whchr)] <- NA
    if (any(!is.na(thepos))) {
        whchr <- whchr[!is.na(whchr)]
        if (shift) 
            thepos <- thepos - sapply(map[whchr], min)
        if (is.matrix(map[[1]])) 
            whchr <- whchr - 0.3
        if (length(grep("^.+@[0-9\\.]+$", x$name)) == length(x$name)) 
            x$name <- x$altname
        if (horizontal) {
            arrows(thepos, whchr - 0.35, thepos, whchr, lwd = 2, 
                col = "red", len = 0.05)
            text(thepos, whchr - 0.4, x$name, col = "red", adj = c(0.5, 
                0))
        }
        else {
            arrows(whchr + 0.35, thepos, whchr, thepos, lwd = 2, 
                col = "red", len = 0.05)
            text(whchr + 0.4, thepos, x$name, col = "red", adj = c(0, 
                0.5))
        }
    }
    attr(x, "formula") <- NULL ### changed qtl to x
    attr(x, "pLOD") <- NULL### changed qtl to x
    invisible()
}

plot.scanone <- function (x, x2, x3, chr, lodcolumn = 1, incl.markers = TRUE, 
    xlim, ylim, lty = 1, col = c("black", "blue", "red"), lwd = 2, 
    add = FALSE, gap = 25, mtick = c("line", "triangle","none"), show.marker.names = FALSE, 
    alternate.chrid = FALSE, labels.rotate=TRUE, bandcol = NULL, ...) 
###All marker names are shown if put to TRUE. 
{
    if (!any(class(x) == "scanone") || (!missing(x2) && !any(class(x2) == 
        "scanone")) || (!missing(x3) && !any(class(x3) == "scanone"))) 
        stop("Input should have class \"scanone\".")
    if (!is.factor(x$chr)) 
        x$chr <- factor(x$chr, levels = unique(x$chr))
    dots <- list(...)
    mtick <- match.arg(mtick)
    if (length(dim(x)) != 2) 
        stop("Argument x must be a matrix or data.frame.")
    if (!missing(x2) && length(dim(x2)) != 2) 
        stop("Argument x2 must be a matrix or data.frame.")
    if (!missing(x3) && length(dim(x3)) != 2) 
        stop("Argument x3 must be a matrix or data.frame.")
    if (length(lodcolumn) == 1) 
        lodcolumn <- rep(lodcolumn, 3)[1:3]
    else if (length(lodcolumn) == 2) {
        if (missing(x2)) 
            x2 <- x
        lodcolumn <- lodcolumn[c(1, 2, 3)]
    }
    else {
        if (missing(x2)) 
            x2 <- x
        if (missing(x3)) 
            x3 <- x
    }
    lodcolumn <- lodcolumn + 2
    second <- third <- TRUE
    if (missing(x2) && missing(x3)) 
        second <- third <- FALSE
    if (missing(x3)) 
        third <- FALSE
    if (missing(x2)) 
        second <- FALSE
    if (lodcolumn[1] > ncol(x) || (second && lodcolumn[2] > ncol(x2)) || 
        (third && lodcolumn[3] > ncol(x3))) 
        stop("Argument lodcolumn misspecified.")
    out <- x[, c(1:2, lodcolumn[1])]
    if (second) 
        out2 <- x2[, c(1:2, lodcolumn[2])]
    if (third) 
        out3 <- x3[, c(1:2, lodcolumn[3])]
    if (length(lty) == 1) 
        lty <- rep(lty, 3)
    if (length(lwd) == 1) 
        lwd <- rep(lwd, 3)
    if (length(col) == 1) 
        col <- rep(col, 3)
    if (missing(chr) || length(chr) == 0) 
        chr <- unique(as.character(out[, 1]))
    else chr <- matchchr(chr, unique(out[, 1]))
    out <- out[!is.na(match(out[, 1], chr)), ]
    if (second) 
        out2 <- out2[!is.na(match(out2[, 1], chr)), ]
    if (third) 
        out3 <- out3[!is.na(match(out3[, 1], chr)), ]
    onechr <- FALSE
    if (length(chr) == 1) {
        gap <- 0
        onechr <- TRUE
    }
    temp <- out
    begend <- matrix(unlist(tapply(temp[, 2], temp[, 1], range)), 
        ncol = 2, byrow = TRUE)
    rownames(begend) <- unique(out[, 1])
    begend <- begend[as.character(chr), , drop = FALSE]
    len <- begend[, 2] - begend[, 1]
    if (!onechr) 
        start <- c(0, cumsum(len + gap)) - c(begend[, 1], 0)
    else start <- 0
    maxx <- sum(len + gap) - gap
    if (all(is.na(out[, 3]))) 
        maxy <- 1
    else maxy <- max(out[, 3], na.rm = TRUE)
    if (second) 
        maxy <- max(c(maxy, out2[, 3]), na.rm = TRUE)
    if (third) 
        maxy <- max(c(maxy, out3[, 3]), na.rm = TRUE)
    old.xpd <- par("xpd")
    old.las <- par("las")
    par(xpd = FALSE, las = 1)
    on.exit(par(xpd = old.xpd, las = old.las))
    if (missing(ylim)) 
        ylim <- c(0, maxy)
    if (missing(xlim)) {
        if (onechr) 
            xlim <- c(0, max(out[, 2]))
        else xlim <- c(-gap/2, maxx + gap/2)
    }
    if (!add) {
        if (onechr) {
            if ("ylab" %in% names(dots)) {
                if ("xlab" %in% names(dots)) {
                  plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
                    ...)
                }
                else {
                  plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
                    xlab = "Map position (cM)", ...)
                }
            }
            else {
                if ("xlab" %in% names(dots)) {
                  plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
                    ylab = dimnames(out)[[2]][3], ...)
                }
                else {
                  plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
                    xlab = "Map position (cM)", ylab = dimnames(out)[[2]][3], 
                    ...)
                }
            }
        }
        else {
            if ("ylab" %in% names(dots)) {
                if ("xlab" %in% names(dots)) {
                  plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
                    xaxt = "n", xaxs = "i", ...)
                }
                else {
                  plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
                    xaxt = "n", xlab = "Chromosome", xaxs = "i", 
                    ...)
                }
            }
            else {
                if ("xlab" %in% names(dots)) {
                  plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
                    xaxt = "n", ylab = dimnames(out)[[2]][3], 
                    xaxs = "i", ...)
                }
                else {
                  plot(0, 0, ylim = ylim, xlim = xlim, type = "n", 
                    xaxt = "n", xlab = "Chromosome", ylab = dimnames(out)[[2]][3], 
                    xaxs = "i", ...)
                }
            }
        }
    }
    if (!add && !onechr && !is.null(bandcol)) {
        u <- par("usr")
        for (i in seq(2, by = 2, length(chr))) {
            rect(min(out[out[, 1] == chr[i], 2]) + start[i] - 
                gap/2, u[3], max(out[out[, 1] == chr[i], 2]) + 
                start[i] + gap/2, u[4], border = bandcol, col = bandcol)
        }
        abline(h = u[3:4])
    }
    xtick <- NULL
    xticklabel <- NULL
    for (i in 1:length(chr)) {
        x <- out[out[, 1] == chr[i], 2] + start[i]
        y <- out[out[, 1] == chr[i], 3]
        if (length(x) == 1) {
            g <- max(gap/10, 2)
            x <- c(x - g, x, x + g)
            y <- rep(y, 3)
        }
        lines(x, y, lwd = lwd[1], lty = lty[1], col = col[1])
        if (!add && !onechr) {
            tloc <- mean(c(min(x), max(x)))
            xtick <- c(xtick, tloc)
            xticklabel <- c(xticklabel, as.character(chr[i]))
        }
        if (second) {
            x <- out2[out2[, 1] == chr[i], 2] + start[i]
            y <- out2[out2[, 1] == chr[i], 3]
            if (length(x) == 1) {
                g <- max(gap/10, 2)
                x <- c(x - g, x, x + g)
                y <- rep(y, 3)
            }
            lines(x, y, lty = lty[2], col = col[2], lwd = lwd[2])
        }
        if (third) {
            x <- out3[out3[, 1] == chr[i], 2] + start[i]
            y <- out3[out3[, 1] == chr[i], 3]
            if (length(x) == 1) {
                g <- max(gap/10, 2)
                x <- c(x - g, x, x + g)
                y <- rep(y, 3)
            }
            lines(x, y, lty = lty[3], col = col[3], lwd = lwd[3])
        }
        if (!add) {
            nam <- dimnames(out)[[1]][out[, 1] == chr[i]]
            wh.genoprob <- grep("^c.+\\.loc-*[0-9]+", nam)
            if (length(wh.genoprob) == 0) 
                wh.genoprob <- seq(along = nam)
            else wh.genoprob <- (seq(along = nam))[-wh.genoprob]
            pos <- out[out[, 1] == chr[i], 2][wh.genoprob] + 
                start[i]
            if (incl.markers) 
              if (mtick!="none"){
                if (mtick == "line") 
                  rug(pos, 0.02, quiet = TRUE)
                else {
                  a <- par("usr")
                  points(pos, rep(a[3] + diff(a[3:4]) * 0.01, 
                    length(pos)), pch = 17, cex = 1.5)
                }
            }
            if (show.marker.names) {
                a <- par("usr")
                text(pos, rep(a[3] + diff(a[3:4]) * 0.03, length(pos)), 
                  nam[wh.genoprob], srt = 90, adj = c(0, 0.5))
            }
        }
    }
    if (!add && !onechr) {
        if (!alternate.chrid || length(xtick) < 2) {
          axis(side = 1, at = xtick, labels = FALSE)
          if (labels.rotate)
            text(xtick,-0.6,labels=xticklabel,srt=45,adj=1,xpd=TRUE,cex=1) ### 45 degree labels
          else
            for (i in seq(along = xtick)) axis(side = 1, at = xtick[i],
                            labels = xticklabel[i])
        }
        else {
            odd <- seq(1, length(xtick), by = 2)
            even <- seq(2, length(xtick), by = 2)
            for (i in odd) {
                axis(side = 1, at = xtick[i], labels = "")
                axis(side = 1, at = xtick[i], labels = xticklabel[i], 
                  line = -0.4, tick = FALSE)
            }
            for (i in even) {
                axis(side = 1, at = xtick[i], labels = "")
                axis(side = 1, at = xtick[i], labels = xticklabel[i], 
                  line = +0.4, tick = FALSE)
            }
        }
    }
    invisible()
}

plot.rf <- function (x, chr, what = c("both", "lod", "rf"), alternate.chrid = FALSE, 
    zmax = 12, mark.diagonal = FALSE, labels.rotate=TRUE,...) 
{
    if (!any(class(x) == "cross")) 
        stop("Input should have class \"cross\".")
    what <- match.arg(what)
    if ("onlylod" %in% names(attributes(x$rf)) && attr(x$rf, 
        "onlylod")) {
        onlylod <- TRUE
        what <- "lod"
    }
    else onlylod <- FALSE
    if (!missing(chr)) 
        x <- subset(x, chr = chr)
    if (!("rf" %in% names(x))) {
        warning("Running est.rf.")
        x <- est.rf(x)
    }
    g <- x$rf
    old.xpd <- par("xpd")
    old.las <- par("las")
    par(xpd = TRUE, las = 1)
    on.exit(par(xpd = old.xpd, las = old.las))
    if (!onlylod) {
        if (any(is.na(g))) 
            g[is.na(t(g))] <- NA
        g[row(g) > col(g) & g > 0.5] <- 0.5
        g[row(g) > col(g)] <- -4 * (log2(g[row(g) > col(g)]) + 
            1)/12 * zmax
    }
    diag(g) <- zmax
    g[!is.na(g) & g > zmax] <- zmax
    g[is.na(g)] <- -1
    if (what == "lod" && !onlylod) {
        g[row(g) > col(g)] <- t(g)[row(g) > col(g)]
    }
    else if (what == "rf") {
        g[row(g) < col(g)] <- t(g)[row(g) < col(g)]
    }
    br <- c(-1, seq(-1e-06, zmax, length = 257))
    image(1:ncol(g), 1:nrow(g), t(g), ylab = "Markers", xlab = "Markers", 
        breaks = br, col = c("lightgray", rev(rainbow(256, start = 0, 
            end = 2/3))))
    if (mark.diagonal) {
        for (i in 1:ncol(g)) segments(i + c(-0.5, -0.5, -0.5, 
            +0.5), i + c(-0.5, +0.5, -0.5, -0.5), i + c(-0.5, 
            +0.5, +0.5, +0.5), i + c(+0.5, +0.5, -0.5, +0.5))
    }
    n.mar <- nmar(x)
    n.chr <- nchr(x)
    a <- c(0.5, cumsum(n.mar) + 0.5)
    abline(v = a, xpd = FALSE, col = "white")
    abline(h = a, xpd = FALSE, col = "white")
    abline(h = 0.5 + c(0, nrow(g)), xpd = FALSE)
    abline(v = 0.5 + c(0, nrow(g)), xpd = FALSE)
    a <- par("usr")
    wh <- cumsum(c(0.5, n.mar))
    chrnam <- names(x$geno)
    chrpos <- (wh[-1] + wh[-length(wh)])/2
    if (!alternate.chrid || length(chrnam) < 2) {
      for (i in seq(along = chrpos)) 
        axis(side = 4, at = chrpos[i], labels = chrnam[i], tick = FALSE, line = -0.8)
      if (labels.rotate)
        {
          axis(side = 3, at = chrpos, labels = FALSE)
          text(chrpos,a[4]+5,labels=chrnam,srt=315,adj=1,xpd=TRUE,cex=1) ### 315 degree labels
        }
      else
        for (i in seq(along = chrpos)) 
          axis(side = 3, at = chrpos[i], labels = chrnam[i], tick = FALSE, line = -0.8)
    }
    else {
        odd <- seq(1, length(chrpos), by = 2)
        even <- seq(2, length(chrpos), by = 2)
        for (i in odd) {
            axis(side = 3, at = chrpos[i], labels = chrnam[i], 
                line = -0.8, tick = FALSE)
            axis(side = 4, at = chrpos[i], labels = chrnam[i], 
                line = -0.8, tick = FALSE)

          }
        for (i in even) {
            axis(side = 3, at = chrpos[i], labels = chrnam[i], 
                line = 0, tick = FALSE)
            axis(side = 4, at = chrpos[i], labels = chrnam[i], 
                line = 0, tick = FALSE)
        }
    }
    dots <- list(...)
    if ("main" %in% names(dots)) 
        title(main = dots$main)
    else {
        if (what == "lod") 
            title(main = "Pairwise LOD scores")
        else if (what == "rf") 
            title(main = "Recombination fractions")
        else title("Pairwise recombination fractions and LOD scores")
    }
    invisible()
}

geno.image <- function (x, chr, reorder = FALSE, main = "Genotype data",
                        alternate.chrid = FALSE, labels.rotate=TRUE,  ...) 
{
    cross <- x
    if (!any(class(cross) == "cross")) 
        stop("Input should have class \"cross\".")
    if (!missing(chr)) 
        cross <- subset(cross, chr = chr)
    type <- class(cross)[1]
    if (type == "bc" || type == "f2") {
        chrtype <- sapply(cross$geno, class)
        if (any(chrtype == "X")) {
            for (i in which(chrtype == "X")) cross$geno[[i]]$data <- reviseXdata(type, 
                "simple", getsex(cross), geno = cross$geno[[i]]$data, 
                cross.attr = attributes(cross))
        }
    }
    Geno <- pull.geno(cross)
    if (type != "4way") {
        thecolors <- c("white", "#E41A1C", "#377EB8", "#4DAF4A", 
            "#984EA3", "#FF7F00")
        thebreaks <- seq(-0.5, 5.5, by = 1)
    }
    else {
        if (max(Geno, na.rm = TRUE) <= 5) {
            thecolors <- c("white", "#E41A1C", "#377EB8", "#4DAF4A", 
                "#984EA3", "#FF7F00")
            thebreaks <- seq(-0.5, 5.5, by = 1)
        }
        else {
            thecolors <- c("white", "#8DD3C7", "#FFFFB3", "#BEBADA", 
                "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", 
                "#D9D9D9", "#BC80BD")
            thebreaks <- seq(-0.5, 10.5, by = 1)
        }
    }
    o <- 1:nrow(Geno)
    if (reorder) {
        if (is.numeric(reorder)) {
            if (reorder < 1 || reorder > nphe(cross)) 
                stop("reorder should be TRUE, FALSE, or an integer between 1 and", 
                  nphe(cross))
            o <- order(cross$pheno[, reorder])
        }
        else {
            wh <- sapply(cross$pheno, is.numeric)
            o <- order(apply(cross$pheno[, wh, drop = FALSE], 
                1, sum))
        }
    }
    g <- t(Geno[o, ])
    g[is.na(g)] <- 0
    old.xpd <- par("xpd")
    old.las <- par("las")
    par(xpd = TRUE, las = 1)
    on.exit(par(xpd = old.xpd, las = old.las))
    image(1:nrow(g), 1:ncol(g), g, ylab = "Individuals", xlab = "Markers", 
        col = thecolors, breaks = thebreaks)
    n.mar <- nmar(cross)
    n.chr <- nchr(cross)
    a <- c(0.5, cumsum(n.mar) + 0.5)
    b <- par("usr")
    segments(a, b[3], a, b[4] + diff(b[3:4]) * 0.02)
    abline(h = 0.5 + c(0, ncol(g)), xpd = FALSE)
    a <- par("usr")
    wh <- cumsum(c(0.5, n.mar))
    x <- 1:n.chr
    for (i in 1:n.chr) x[i] <- mean(wh[i + c(0, 1)])
    thechr <- names(cross$geno)
    if (!alternate.chrid || length(thechr) < 2) {
      if (labels.rotate)
        {
          axis(side = 3, at = x, labels = FALSE)
          text(x,a[4]+3,labels=thechr,srt=315,adj=1,xpd=TRUE,cex=1) ### 315 degree labels
        }
      else
        for (i in seq(along = x)) axis(side = 3, at = x[i], thechr[i], 
                        tick = FALSE, line = -0.5)
    }
    
    else {
        odd <- seq(1, length(x), by = 2)
        even <- seq(2, length(x), by = 2)
        for (i in odd) axis(side = 3, at = x[i], labels = thechr[i], 
            line = -0.75, tick = FALSE)
        for (i in even) axis(side = 3, at = x[i], labels = thechr[i], 
            line = +0, tick = FALSE)
    }
    title(main = main)
    invisible()
}

plot.map <- function (x, map2, chr, horizontal = FALSE, shift = TRUE, show.marker.names = FALSE, 
    alternate.chrid = FALSE,labels.rotate=TRUE, ...) 
{
    dots <- list(...)
    if ("main" %in% names(dots)) {
        themain <- dots$main
        usemaindefault <- FALSE
    }
    else usemaindefault <- TRUE
    if ("xlim" %in% names(dots)) {
        xlim <- dots$xlim
        usexlimdefault <- FALSE
    }
    else usexlimdefault <- TRUE
    if ("ylim" %in% names(dots)) {
        ylim <- dots$ylim
        useylimdefault <- FALSE
    }
    else useylimdefault <- TRUE
    if ("xlab" %in% names(dots)) 
        xlab <- dots$xlab
    else {
        if (horizontal) 
            xlab <- "Location (cM)"
        else xlab <- "Chromosome"
    }
    if ("ylab" %in% names(dots)) 
        ylab <- dots$ylab
    else {
        if (horizontal) 
            ylab <- "Chromosome"
        else ylab <- "Location (cM)"
    }
    map <- x
    if (any(class(map) == "cross")) 
        map <- pull.map(map)
    if (!missing(map2) && any(class(map2) == "cross")) 
        map2 <- pull.map(map2)
    if (!any(class(map) == "map") || (!missing(map2) && !any(class(map2) == 
        "map"))) 
        warning("Input should have class \"cross\" or \"map\".")
    if (!missing(map2) && is.matrix(map[[1]]) != is.matrix(map2[[1]])) 
        stop("Maps must be both sex-specific or neither sex-specific.")
    if (!missing(chr)) {
        map <- map[matchchr(chr, names(map))]
        if (!missing(map2)) 
            map2 <- map2[matchchr(chr, names(map2))]
    }
    sex.sp <- FALSE
    if (is.matrix(map[[1]])) {
        one.map <- FALSE
        sex.sp <- TRUE
        if (!missing(map2)) {
            if (is.logical(map2)) {
                horizontal <- map2
                map2 <- lapply(map, function(a) a[2, ])
                map <- lapply(map, function(a) a[1, ])
            }
            else {
                Map1 <- lapply(map, function(a) a[1, , drop = TRUE])
                Map2 <- lapply(map, function(a) a[2, , drop = TRUE])
                Map3 <- lapply(map2, function(a) a[1, , drop = TRUE])
                Map4 <- lapply(map2, function(a) a[2, , drop = TRUE])
                old.mfrow <- par("mfrow")
                on.exit(par(mfrow = old.mfrow))
                par(mfrow = c(2, 1))
                class(Map1) <- class(Map2) <- class(Map3) <- class(Map4) <- "map"
                plot.map(Map1, Map3, horizontal = horizontal, 
                  shift = shift, show.marker.names = show.marker.names, 
                  alternate.chrid = alternate.chrid)
                plot.map(Map2, Map4, horizontal = horizontal, 
                  shift = shift, show.marker.names = show.marker.names, 
                  alternate.chrid = alternate.chrid)
                return(invisible(NULL))
            }
        }
        else {
            map2 <- lapply(map, function(a) a[2, ])
            map <- lapply(map, function(a) a[1, ])
        }
    }
    else {
        if (!missing(map2)) 
            one.map <- FALSE
        else one.map <- TRUE
    }
    if (one.map) {
        n.chr <- length(map)
        if (!show.marker.names) {
            chrpos <- 1:n.chr
            thelim <- range(chrpos) + c(-0.5, 0.5)
        }
        else {
            chrpos <- seq(1, n.chr * 2, by = 2)
            thelim <- range(chrpos) + c(-0.35, 2.35)
        }
        if (shift) 
            map <- lapply(map, function(a) a - a[1])
        maxlen <- max(unlist(lapply(map, max)))
        if (horizontal) {
            old.xpd <- par("xpd")
            old.las <- par("las")
            par(xpd = TRUE, las = 1)
            on.exit(par(xpd = old.xpd, las = old.las))
            if (usexlimdefault) 
                xlim <- c(0, maxlen)
            if (useylimdefault) 
                ylim <- rev(thelim)
            plot(0, 0, type = "n", xlim = xlim, ylim = ylim, 
                yaxs = "i", xlab = xlab, ylab = ylab, yaxt = "n")
            a <- par("usr")
            for (i in 1:n.chr) {
                segments(min(map[[i]]), chrpos[i], max(map[[i]]), 
                  chrpos[i])
                segments(map[[i]], chrpos[i] - 0.25, map[[i]], 
                  chrpos[i] + 0.25)
                if (show.marker.names) 
                  text(map[[i]], chrpos[i] + 0.35, names(map[[i]]), 
                    srt = 90, adj = c(1, 0.5))
            }
            if (!alternate.chrid || length(chrpos) < 2) {
                for (i in seq(along = chrpos)) axis(side = 2, 
                  at = chrpos[i], labels = names(map)[i])
            }
            else {
                odd <- seq(1, length(chrpos), by = 2)
                even <- seq(2, length(chrpos), by = 2)
                for (i in odd) {
                  axis(side = 2, at = chrpos[i], labels = "")
                  axis(side = 2, at = chrpos[i], labels = names(map)[i], 
                    line = -0.4, tick = FALSE)
                }
                for (i in even) {
                  axis(side = 2, at = chrpos[i], labels = "")
                  axis(side = 2, at = chrpos[i], labels = names(map)[i], 
                    line = +0.4, tick = FALSE)
                }
            }
        }
        else {
            old.xpd <- par("xpd")
            old.las <- par("las")
            par(xpd = TRUE, las = 1)
            on.exit(par(xpd = old.xpd, las = old.las))
            if (usexlimdefault) 
                xlim <- thelim
            if (useylimdefault) 
                ylim <- c(maxlen, 0)
            plot(0, 0, type = "n", ylim = ylim, xlim = xlim, 
                xaxs = "i", xlab = xlab, ylab = ylab, xaxt = "n")
            a <- par("usr")
            for (i in 1:n.chr) {
                segments(chrpos[i], min(map[[i]]), chrpos[i], 
                  max(map[[i]]))
                segments(chrpos[i] - 0.25, map[[i]], chrpos[i] + 
                  0.25, map[[i]])
                if (show.marker.names) 
                  text(chrpos[i] + 0.35, map[[i]], names(map[[i]]), 
                    adj = c(0, 0.5))
            }
            if (!alternate.chrid || length(chrpos) < 2) {
             if (labels.rotate)
                {
                  axis(side = 1, at = chrpos, labels =FALSE)
                  text(chrpos,a[3]+(a[3]/50),labels=names(map),srt=45,adj=1,xpd=TRUE,cex=1) ### 45 degree labels
                }
              else
                for (i in seq(along = chrpos)) axis(side = 1, 
                  at = chrpos[i], labels = names(map)[i])

              }
            else {
                odd <- seq(1, length(chrpos), by = 2)
                even <- seq(2, length(chrpos), by = 2)
                for (i in odd) {
                  axis(side = 1, at = chrpos[i], labels = "")
                  axis(side = 1, at = chrpos[i], labels = names(map)[i], 
                    line = -0.4, tick = FALSE)
                }
                for (i in even) {
                  axis(side = 1, at = chrpos[i], labels = "")
                  axis(side = 1, at = chrpos[i], labels = names(map)[i], 
                    line = +0.4, tick = FALSE)
                }
            }
        }
        if (usemaindefault) 
            title(main = "Genetic map")
        else if (themain != "") 
            title(main = themain)
    }
    else {
        if (is.matrix(map2[[1]])) 
            stop("Second map appears to be a sex-specific map.")
        if (length(map) != length(map2)) 
            stop("Maps have different numbers of chromosomes.")
        if (any(sapply(map, length) != sapply(map2, length))) 
            stop("Maps have different numbers of markers.")
        map1 <- map
        if (shift) {
            map1 <- lapply(map1, function(a) a - a[1])
            map2 <- lapply(map2, function(a) a - a[1])
        }
        n.chr <- length(map1)
        maxloc <- max(c(unlist(lapply(map1, max)), unlist(lapply(map2, 
            max))))
        if (!show.marker.names) {
            chrpos <- 1:n.chr
            thelim <- range(chrpos) + c(-0.5, 0.5)
        }
        else {
            chrpos <- seq(1, n.chr * 2, by = 2)
            thelim <- range(chrpos) + c(-0.4, 2.4)
        }
        if (!horizontal) {
          old.xpd <- par("xpd")
            old.las <- par("las")
            par(xpd = TRUE, las = 1)
            on.exit(par(xpd = old.xpd, las = old.las))
            if (usexlimdefault) 
                xlim <- thelim
            if (useylimdefault) 
                ylim <- c(maxloc, 0)
            plot(0, 0, type = "n", ylim = ylim, xlim = xlim, 
                xaxs = "i", xlab = xlab, ylab = ylab, xaxt = "n")
            a <- par("usr")
            for (i in 1:n.chr) {
                if (max(map2[[i]]) < max(map1[[i]])) 
                  map2[[i]] <- map2[[i]] + (max(map1[[i]]) - 
                    max(map2[[i]]))/2
                else map1[[i]] <- map1[[i]] + (max(map2[[i]]) - 
                  max(map1[[i]]))/2
                segments(chrpos[i] - 0.3, min(map1[[i]]), chrpos[i] - 
                  0.3, max(map1[[i]]))
                segments(chrpos[i] + 0.3, min(map2[[i]]), chrpos[i] + 
                  0.3, max(map2[[i]]))
                segments(chrpos[i] - 0.3, map1[[i]], chrpos[i] + 
                  0.3, map2[[i]])
                if (show.marker.names) 
                  text(chrpos[i] + 0.35, map2[[i]], names(map2[[i]]), 
                    adj = c(0, 0.5))
            }
            if (!alternate.chrid || length(chrpos) < 2) {
               if (labels.rotate)
                {
                  axis(side = 1, at = chrpos, labels =FALSE)
                  text(chrpos,a[3]+(a[3]/50),labels=names(map1),srt=45,adj=1,xpd=TRUE,cex=1) ### 45 degree labels
                }
              else
                for (i in seq(along = chrpos)) axis(side = 1, 
                                    at = chrpos[i], labels = names(map1)[i])
            }
            else {
                odd <- seq(1, length(chrpos), by = 2)
                even <- seq(2, length(chrpos), by = 2)
                for (i in odd) {
                  axis(side = 1, at = chrpos[i], labels = FALSE)
                  axis(side = 1, at = chrpos[i], labels = names(map1)[i], 
                    line = -0.4, tick = FALSE)
                }
                for (i in even) {
                  axis(side = 1, at = chrpos[i], labels = FALSE)
                  axis(side = 1, at = chrpos[i], labels = names(map1)[i], 
                    line = +0.4, tick = FALSE)
                }
            }
        }
        else {
            old.xpd <- par("xpd")
            old.las <- par("las")
            par(xpd = TRUE, las = 1)
            on.exit(par(xpd = old.xpd, las = old.las))
            if (usexlimdefault) 
                xlim <- c(0, maxloc)
            if (useylimdefault) 
                ylim <- rev(thelim)
            plot(0, 0, type = "n", xlim = xlim, ylim = ylim, 
                xlab = xlab, ylab = ylab, yaxt = "n", yaxs = "i")
            a <- par("usr")
            for (i in 1:n.chr) {
                if (max(map2[[i]]) < max(map1[[i]])) 
                  map2[[i]] <- map2[[i]] + (max(map1[[i]]) - 
                    max(map2[[i]]))/2
                else map1[[i]] <- map1[[i]] + (max(map2[[i]]) - 
                  max(map1[[i]]))/2
                segments(min(map1[[i]]), chrpos[i] - 0.3, max(map1[[i]]), 
                  chrpos[[i]] - 0.3)
                segments(min(map2[[i]]), chrpos[i] + 0.3, max(map2[[i]]), 
                  chrpos[[i]] + 0.3)
                segments(map1[[i]], chrpos[i] - 0.3, map2[[i]], 
                  chrpos[i] + 0.3)
                if (show.marker.names) 
                  text(map2[[i]], chrpos[i] + 0.35, names(map2[[i]]), 
                    srt = 90, adj = c(1, 0.5))
            }
            if (!alternate.chrid || length(chrpos) < 2) {
              if (labels.rotate)
                {
                  axis(side = 1, at = chrpos, labels = FALSE)
                  text(chrpos,-0.6,labels=names(map1),srt=45,adj=1,xpd=TRUE,cex=1) ### 45 degree labels
                }
              else
                for (i in seq(along = chrpos)) axis(side = 2, 
                                    at = chrpos[i], labels = names(map1)[i])
#            for (i in seq(along = xtick)) axis(side = 1, at = xtick[i],
#                labels = xticklabel[i])
      
                }
            else {
                odd <- seq(1, length(chrpos), by = 2)
                even <- seq(2, length(chrpos), by = 2)
                for (i in odd) {
                  axis(side = 2, at = chrpos[i], labels = "")
                  axis(side = 2, at = chrpos[i], labels = names(map1)[i], 
                    line = -0.4, tick = FALSE)
                }
                for (i in even) {
                  axis(side = 2, at = chrpos[i], labels = "")
                  axis(side = 2, at = chrpos[i], labels = names(map1)[i], 
                    line = +0.4, tick = FALSE)
                }
            }
        }
        if (usemaindefault) {
            if (!sex.sp) 
                title(main = "Comparison of genetic maps")
            else title(main = "Genetic map")
        }
        else if (themain != "") 
            title(main = themain)
    }
    invisible()
}

matchchr <- function (selection, thechr) 
{
    if (is.factor(thechr)) 
        thechr <- as.character(thechr)
    if (length(thechr) > length(unique(thechr))) 
        stop("Duplicate chromosome names.")
    if (is.logical(selection)) {
        if (length(selection) != length(thechr)) 
            stop("Logical vector to select chromosomes is the wrong length")
        return(thechr[selection])
    }
    if (is.numeric(selection)) 
        selection <- as.character(selection)
    if (length(selection) > length(unique(selection))) {
        warning("Dropping duplicate chromosomes")
        selection <- unique(selection)
    }
    g <- grep("^-", selection)
    if (length(g) > 0 && length(g) < length(selection)) 
        stop("In selecting chromosomes, all must start with '-' or none should.")
    if (length(g) > 0) {
        selectomit <- TRUE
        selection <- substr(selection, 2, nchar(selection))
    }
    else selectomit <- FALSE
    wh <- match(selection, thechr)
    if (any(is.na(wh))) {
        warning("Chr ", paste(selection[is.na(wh)], collapse = ", "), 
            " not found")
        wh <- wh[!is.na(wh)]
        if (length(wh) == 0) 
            return(thechr)
    }
    if (selectomit) 
        return(thechr[-wh])
    thechr[sort(wh)]
}
