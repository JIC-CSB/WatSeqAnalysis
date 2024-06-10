### R script
### QTL analysis with the R/qtl packaga
### Functions for
### -data upload and formatting
### -QTL analysis
### -plotting results into pdf files
### 1/12/2017 Luzie U. Wingen, JIC
### update 04/05/2020 LUW

library(qtl)

mywd <- getwd()
###print(mywd)
###setwd("~/projects/qtl_analysis/rqtl_analysis/rqtl_jic_vs1.7")
print("Welcome to QTL analysis for wheat trait data",quote=FALSE)
print('Loading further source files.',quote=FALSE)
source("service_functions_files.R")
source("my_rqtl_functions.R")
###setwd(mywd)

checkMaps <- function(gtfnames,outdir,printpdf=TRUE)
### run through map diagnostic for a bundle of maps. That can
{
  print('welcome to R/QTL check Map',quote=FALSE)
  maplist <- c()
  for (gtfn in gtfnames)
    {
      print(gtfn)
      mapname <- sub('_rqtl_.*','',gtfn)
      mapname <- sub('^.*/','',mapname)
      pfn <- sub('genotype','phenotype',gtfn)
      print(gtfn)
      cor=T ### this has to be made a bit more brightly.
      if (cor==T)
        corL <- startCorCalcRQTL(pfn)
      MyCross <- read.cross(format="csvs",genfile=gtfn,phefile=pfn,,genotypes=c("a","b")) ###treat it as backcross
      print('read.cross object was established',quote=FALSE)
      namesMyCross <- names(MyCross$geno)
      print('Names of Linkage Groups:',quote=FALSE)
      namesMyCross <- renameLGR(namesMyCross)
      print(namesMyCross)
      names(MyCross$geno) <- namesMyCross
      MyCross <- mapAnalysis(MyCross,gtfn,outdir,printpdf=T,diagnostic=TRUE)
      maplist[mapname] <- MyCross
    }
  return(maplist)
}

nameLinkagegroups <- function(GTmatrix)
### linkage groups in the genotype matrix are renamed according to the
### wheat chromosome name - given that the name is attached to the marker name
{
  nameLGR <- c()
  for (lg in unique(GTmatrix[,2]))
    {
      lgindex <- GTmatrix[,2]==lg
      lgmarker <- GTmatrix[lgindex,1]
      lgname <- c()
      for (lgm in lgmarker)
        {
          lgsub <-  grep("^[0-9][A-D]",sub("^.*-","",lgm),value=TRUE)
          if (length(lgsub)==0)
            lgsub <-  grep("^[0-9][A-D]",sub("-.*$","",lgm),value=TRUE)
          if (length(lgsub)>0)
            {
              if (length(grep("/",lgsub))>0)
                for (i in 1:grep("/",lgsub))
                  {
                    lgname <- c(lgname,sub("/.*","",lgsub))
                    lgsub <- sub(".*/","",lgsub)
                  }
              lgname <- c(lgname,lgsub)
            }
        }
      if (length(lgname)>0)
        {
          lgsum <- summary(as.factor(sub("/.*","",lgname)))
          if (length(lgsum)>1)
            lgname <- names(lgsum[lgsum>1])
          lgname <- unique(lgname)
          if (length(lgname)>1)
            {
              lnamestart <- ""
              for (l in lgname)
                {
                  if (lnamestart=="")
                    lnamestart <- l
                  else
                    lnamestart <- paste(lnamestart,l,sep="/")
                }
              lgname <- lnamestart
            }
        }
       nameLGR <- rbind(nameLGR,c(lg, lgname))
    }
  return(nameLGR)
}

addToMatrix <- function(phenotypematrix,treatment,line,phenotypevec)
### adding mean values and source values for phenotypes to a matrix
{
  if (length(phenotypevec)<3)
    phenotypevec <- c(phenotypevec,rep(NA,(3-length(phenotypevec))))
  if (length(phenotypevec)>3)
    {
      print(paste("too many results:",length(phenotypevec),"in",treatment,line))
      print(phenotypevec)
      phenotypevec <- phenotypevec[c(1,2,4)]
    }
  phenotypematrix <- rbind(phenotypematrix,c(paste("X",line,sep=""),round(mean(phenotypevec,na.rm=TRUE),2),phenotypevec,treatment))
  return(phenotypematrix)
}


renameLGR <- function(namesChromosomes)
### Add 'lg' at the start of  numeric linkage group names 
### If names end '_<lgnumber>', this number will be converted to 1, 2 ..
{

    if(any(grepl('^[0-9]*$',namesChromosomes))) ### just numeric names
        namesChromosomes <- paste('lg',namesChromosomes,sep='')
    if(any(grepl('_',namesChromosomes))) ### any underscore
    {
        for (nm in namesChromosomes) ### check for underscores
        {
            if (!(is.na(nm)))
            {
                chr <- nm
                chrL <- strsplit(nm,'_')[[1]]
                for (j in (1:length(chrL)))
                {
                    if (length(grep('[1-7][A-D]*',chrL[j]))>0) ### wheat chromosoes
                        chr <- chrL[j]
                    break
                }
            }
            chrPat=paste('^',chr,'$',sep='')
            if (length(grep(chrPat,namesChromosomes))>0)
            {
                if (length(grep(chrPat,namesChromosomes))==1)
                {
                    namesChromosomes[grep(chrPat,namesChromosomes)] <- chr
                }
                else ### several linkage groups from one chromosome
                {
                    appendix=1
                    chrlgs <- namesChromosomes[grep(chrPat,namesChromosomes)]
                    chrnos <- grep(chrPat,namesChromosomes)
                    unprocessed <- grep('_',chrlgs)
                    if (length(unprocessed)>0)
                        for (u in unprocessed)
                        {
                            namesChromosomes[chrnos[u]] <- paste(chrPat,letters[appendix],sep="")
                            appendix=appendix+1
                        }
                }
            }
        }
    }
  return(namesChromosomes)
}

modelformula <- function(qtlnumber,m)
{
  model <- "y~Q"
  for (qi in 1:qtlnumber)
    model <- paste(model,qi,m,"Q",sep="")
  model <- sub("[^0-9]*$","",model)
  return(model)
}


findAdditiveEffect <- function(qtleffect,chr,qpos)
###finds additive effect for qtl table
{
####  print('in routine to finding additive effect')
  chreffect <- qtleffect[qtleffect[,1]%in%chr,]
  chreffect <- chreffect[!is.na(chreffect$a),] ### remove NaN which are a nuisance
  markerrows <- 1:length(rownames(chreffect))
  simmarrows <- grep(paste('c',chr,sep=""),rownames(chreffect))
  if (length(simmarrows)>0)
    markerrows <- markerrows[-simmarrows]
  markereffect <- chreffect[markerrows,]
    nearestmarker <-markereffect[abs(markereffect[,2]-qpos)==min(abs(markereffect[,2]-qpos)),]
    nearestmarker <- nearestmarker[sample(dim(nearestmarker)[1],1),]
  #  print('nearestmarker')
  #  print(nearestmarker)
    exactposition <- chreffect[chreffect[,2]%in%qpos,]
  #  print('exactposition')
  #  print(exactposition)
  #  print(dim(exactposition))
  if (dim(exactposition)[1]>0)
    {
      psample <- sample(dim(exactposition)[1],1)### chose one position
      exactposition <- exactposition[psample,]
      print(exactposition)
      tablerow <- c(sapply(as.numeric(exactposition[3:4]),round,3),rownames(exactposition),round(nearestmarker[,2],3),rownames(nearestmarker))
  #    print(tablerow)
    }
  else
  {
 #     print('dim exacpositions[1] >0')
      qtlrange <- c(qpos-10,qpos+10)
      findposition <- chreffect[(chreffect[,2]>qtlrange[1])&(chreffect[,2]<qtlrange[2]),]
#      print('findposition for 1 dim')
#      print(findposition)
      if (dim(findposition)[2]>0)
        {
            fmin <- min(abs(findposition[,2]-qpos[1]),na.rm=TRUE)
   #         print(paste('fmin',fmin))
            if (is.finite(fmin))
            {
                posindex <- abs(findposition[,2]-qpos[1])==fmin
              findposition <- findposition[posindex,]
              psample <- sample(dim(findposition)[1],1)### chose one position
              findposition <- findposition[psample,]
              tablerow <- c(sapply(as.numeric(findposition[3:4]),round,3),rownames(findposition),round(nearestmarker[,2],3),rownames(nearestmarker))
              if(is.na(findposition[3]))
              {
                  ### Moved to effect nearby ### needs improving
                  nonNArows <- (1:nrow(chreffect))[!is.na(chreffect[,3])]
                  posno <- (1:nrow(chreffect))[row.names(chreffect)%in%row.names(findposition)]
                  findpos <- chreffect[nonNArows,][min(abs(nonNArows-posno))==nonNArows-posno,]

                  print('findpos when there are NAs')
                  print(findpos)
                  tablerow <- c(sapply(as.numeric(findposition[3:4]),round,3),rownames(findposition),round(nearestmarker[,2],3),rownames(nearestmarker))
              }
            }
          else
            {
              print('Do not find a min(abs(findposition[,2]))',quote=FALSE)
              print(findposition)
              tablerow <- c(NA,NA,NA,round(nearestmarker[,2],3),rownames(nearestmarker))
            }
        }
      else
        {
          print('findposition is not a matrix',quote=FALSE)
          print(findposition)
          tablerow <- c(NA,NA,NA,round(nearestmarker[,2],3),rownames(nearestmarker))
        }
    }
  return(tablerow)
}

QTLplotFilename <- function(plotparlist,plotno=1,plotname='',traitnames="")
{
### directs the plots to the right filename
  filename=""
#  print("QTLplotFilename")
#  print("plotparlist")
#  print(plotparlist)
#  print("plotno")
#  print(plotno)
#  print("plotname")
#  print(plotname)
#  print("traitnames")
#  print(traitnames)
  ext=plotparlist$pmethod
    if (ext=="png")###single QTL:Q<trait>_<chr>_psr_<pop+env>.png
    {
      plotno=1
      ext="png"
      filename <- paste(sub('[.][a-z]*$','',plotparlist$filenamevector[plotno]),'.',ext,sep="")
      if (traitnames!='')
        filename <- sub("201",paste(traitnames,"_201",sep=""),filename)
      if (sum(plotparlist$pmethod%in%names(dev.list()))>0) ### shutting off devices from last plot
          dev.off(dev.list()[names(dev.list())%in%plotparlist$pmethod])
        names(dev.list())
    png(filename,width=2000, height=1200, units="px",pointsize="20")
###      png(filename,width=1000, height=600, units="px",pointsize="10",)
    }
    if (ext=="pdf")###many QTL plotted
    {
        if(!is.na(plotname))
            dev.set(plotparlist$devicelist[plotname])
        else
            dev.set(plotparlist$devicelist[plotno])
    }
  return(filename)
}

QTLSingleFilename <- function(plotparlist,plotno,traitnames,chro)
{
### if plotted into single plots, e.g. png or svg
#### single QTL plot in png should be: QHD_2B_psr_ParW700_CFLN17.png
    ext=plotparlist$pmethod
    print(chro)
    filename <- plotparlist$filenamevector[1]
    filename <- sub('(^.*/)([^_]*).*','\\1Q_trait_chro_psr_\\2_env',filename)
    filename <- paste(filename,ext,sep='.')
    if (traitnames!='')
        for(tr in traitnames)
        {
            filename <- sub("_trait_",sub('_.*','_',tr),filename)
            filename <- sub("_env",sub('.*_','_',tr),filename)
            chro <- sub('[0-9]$','',chro)
            print(paste('chromosome:',chro))
            filename <- sub("chro",chro,filename) ### removing
####    filename <- sub("chro",sub('[0-9]$','',chro),filename)) ### removing 
            if (sum(plotparlist$pmethod%in%names(dev.list()))>0) ### shut off device from last plot
                dev.off(dev.list()[names(dev.list())%in%plotparlist$pmethod])
            print(plotparlist$pmethod)
            if(plotparlist$pmethod=='svg')
                svg(filename,width=14, height=14,pointsize="12")
            if(plotparlist$pmethod=='png')
            {
###                png(filename,width=1000, height=600, units="px",pointsize="20")
                png(filename,width=2000, height=1200, units="px",pointsize="20")
                ##
                help(par)
                par('cex'=0.4)
                par('cex.axis'=1.5)
                par('cex.lab'=1.5)

            }



            nf <- layout(matrix(c(1,2),nrow=2,ncol=1),heights=c(2.8,1.2))
        }
    return(filename)
}


peakPlotPosition <- function(scan1out,chr,gap=25,DART=T,drop=0,ypos=0,oneplot=F)
### prints name of peak marker near the  peak and prints other markers at
### the bottom of plot with an even distribution.
{
  onechr=FALSE
  markercol <- c()
  if (length(chr)==1)
    onechr=TRUE
  begend <- matrix(unlist(tapply(scan1out[, 2], scan1out[, 1], range)), ncol = 2, byrow = TRUE)
  rownames(begend) <- unique(scan1out[, 1])
  begend <- begend[as.character(chr), , drop = FALSE]
  len <- begend[, 2] - begend[, 1]
  if (!onechr)
    start <- c(0, cumsum(len + gap)) - c(begend[, 1], 0)
  else
    start <- 0
  for (mp in 1:length(chr)) {
    peak <- scan1out[scan1out[, 1] == chr[mp], ]
    lodv <- peak[,3]
    tloc <- peak[lodv==max(lodv),2]###
    if (drop>0)
      peak <- peak[lodv>(max(lodv)-drop),]### mark all real marker above CI
    realmark <- peak[!grepl('[.]loc',row.names(peak)),]### exclude calc. genetic predictor
    if (!DART)
      if (sum(!grepl('wPt',row.names(realmark)),na.rm=TRUE)>2)
        realmark <- realmark[!grepl('wPt',row.names(realmark)),]### exclude calculated and DArT
    if (nrow(realmark)>0)
      {
        markerpeak <- realmark[abs(realmark[,2]-tloc)==min(abs(realmark[,2]-tloc)),] ### marker closest to peak
        tloc <- markerpeak[1,2]
        markerpos <-  seq(begend[1],begend[2],by=begend[2]/(nrow(realmark))) ### not sure if markerpos is working in this setting.
        markercol <- rep('darkgreen',nrow(realmark))
        markercol[rownames(realmark)%in%rownames(markerpeak)] <- 'red'
        for (j in 1:nrow(realmark))
          {
            dist=-0.5
            pcol=markercol[j]
            pstr=90
            if (rownames(realmark[j,])%in%rownames(markerpeak))
              {
                dist=0.5
                pstr=0
              }
            markercol <- c(markercol,pcol)
            yadd <- ypos+par('usr')[4]/25
            if (pcol=='red')
              text(realmark[j,2]+start[mp],realmark[j,3]+dist,row.names(realmark)[j],col=pcol,srt=pstr,adj=0.5)
          }
        ypos=par('usr')[3]
        if (oneplot)
          plotMarkerSegments(realmark,markerpos,markercol,ypos,yadd,start[mp])
      }
  }
  return(list('realmark'=realmark,'markerpos'=markerpos,'markercol'=markercol))
}

plotMarkerSegments <- function(realmark,markerpos,markercol,ypos,yadd,start=0)
{
  for (j in 1:length(markerpos))
    {
      plotxpos <- markerpos[j]
      text(plotxpos+start,yadd,row.names(realmark)[j],col=markercol[j],srt=90,adj=0,cex=0.8)
      lwd=0.05
      segments(realmark[j,2]+start,ypos,plotxpos+start,yadd,col='black')
    }
}

testparameter <- function()
{
    ### just to test
    colvec=c("blue","red","green","orange","black")
    printlwd=3
    pmain='test'
    drop=0
    QTLtraceL <- QTLtraceLsel
    try(plotRQTLtraces(QTLtraceLsel,paste(pmain,'-',sub('_$','',tr)),colvec=colvec,plotparlist,printlwd=printlwd))
    devicelist <- dev.list()
    names(devicelist) <- 'qtls'
    pmethod='out'
    plotparlist <- list(filenamevector='test.pdf',devicelist=devicelist,pmethod=pmethod)
    dev.off()
}


plotRQTLtraces <- function(qtltraceL,pmain,colvec=c("blue","red","green","orange","black"),plotparlist,drop=0,printlwd=8)### printlwd=3
{
    findCIs <- function(traceChrDfr,chro,myqtltable,colvec)
    {
      CIV <- c()
      bordermarker <- c()
      qtlplottable <- c()
      myqtltab <- myqtltable[grep(chro,myqtltable$chr),]
      qtlTrait <- myqtltab[,'trait']
      if (length(qtlTrait)>0)
          for (trC in qtlTrait)
          {
              trCol <- sub('^.*_','',trC)
              if (trC%in%myqtltab$trait)
                  qtlplottable <- rbind(qtlplottable,cbind(myqtltab[myqtltab$trait==trC,c('CIstart','CIend','pos')],col=colvec[trCol],stringsAsFactors=FALSE))
              CIV <- c(CIV,myqtltab[,'nearest marker'])
              bordermarker <- c(bordermarker,myqtltab[,"start marker"],myqtltab[,"end marker"])
          }
      bordermarker <- unique(bordermarker)
      qtlplottable <- data.frame(qtlplottable,stringsAsFactors=FALSE)
      return(list(qtlplottable,bordermarker,CIV))
    }
  plotQTLs <- function(traceChrDfr,pmain,traitnames,chr,colvec,printlwd,ymax,thresholdV)
  {
      traitvec <- sub(' Eff.*','',traitnames)
      for (trC in 3:ncol(traceChrDfr))
      {
          trcol <-traitvec[trC-2]
          trCol <- sub('^.*_','',trcol)
          if (trC==3)
          {
###              pmain <- sub('PxBaj','Paragon x Baj',pmain)
              pmain=sub("- ",paste("Q",chro,"-",sep=''),pmain)
              print(paste('pmain',print(pmain)))
              plot(chrQTs[[trC-2]],chr=chro,main=pmain,xlab="",ylab="",col=colvec[trCol],lwd=printlwd,ylim=c(0,max(c(ymax+0.5,thresholdV))),mtick="none",axes=F,frame.plot=T)
              axis(2)
              title(ylab='LOD',line=2)
          }
          else
              lines(traceChrDfr[,2],traceChrDfr[,trC],col=colvec[trCol],lwd=printlwd)
          maxlod <- max(traceChrDfr[,trC],na.rm=T)
          maxtrace <- trC
          ypos <- par('usr')[3]
          yadd <- par('usr')[4]/100
          peakypos <- ypos+(ncol(traceDfr)-1)*yadd
###          if (maxlod==ymax)
              realmarkerL <- peakPlotPosition(chrQTs[[trC-2]],chro,DART=T,drop=0,ypos=peakypos)
      }
      legend("topleft",traitnames,col=colvec[sub('^.*_','',traitvec)],lty=1,lwd=printlwd)
      for (j in 1:length(thresholdV))
          abline(h=thresholdV[j],lty=2,col=colvec[sub('^.*_','',traitvec)][j])
      return(list(realmarkerL,chrQTs))
  }
### main
###    par("cex"=1)
    if (plotparlist$pmethod=="pdf")
        par("cex"=0.7)
    pln <- QTLplotFilename(plotparlist,plotname="qtls")
    nf <- layout(matrix(c(1,2,3,4),nrow=4,ncol=1),heights=c(2.5,1,2.5,1),width=c(1,1,1,1))
    traitenvnames <- names(qtltraceL)
    selected=1:length(qtltraceL)#### no selection here anymore
    if (length(selected)>0)
    {
        traceDfr <- c()
        chrs <- c()
        thresholdV <- c()
        chrQTs <- list()
        myqtltable <- c()
        traitnameV <- c()
        for (no in selected)
        {
            QT <- qtltraceL[[no]]
            traitnameV <- c(traitnameV,names(qtltraceL)[no])
            chrQTs[traitenvnames[no]] <- list(QT$trace)
            chrs <- sort(unique(c(chrs,as.character(qtltraceL[[no]]$summary$chr))))
            traceDfr <- cbind(traceDfr,QT$trace$lod)
            thresholdV <-c(thresholdV, QT$threshold)
            if (!is.null(QT$qtltable)) ### if empty qtltable, no qtl
            {
                myqtltable <- rbind(myqtltable,QT$qtltable)
            }
        }
        colnames(traceDfr) <- traitnameV
        traceDfr <- cbind(QT$trace$chr,QT$trace$pos,traceDfr)
        rownames(traceDfr) <- rownames(QT$trace)
        if (!is.null(myqtltable))
        {
            row.names(myqtltable) <- make.unique(as.character(row.names(myqtltable)))
            myqtltable <- as.data.frame(myqtltable,stringsAsFactors=F)
        }
        myqtltable <- myqtltable[sort(row.names(myqtltable)),]
        chro <- myqtltable$chr[1]
####for (chro in "7Db") put this in if a fixed chromosome to be printed
        for (chro in unique(myqtltable$chr))
        {
            traceChrDfr <- traceDfr[QT$trace[,'chr']==chro,]
            ymax <- max(traceChrDfr[,3:ncol(traceChrDfr)],na.rm=T)
###       ymax <- 8.0 ### put this in, if a fixed ymax wanted
            if (plotparlist$pmethod!="pdf")
                nf <- layout(matrix(c(1,2),nrow=2,ncol=1),heights=c(2.8,1.2))
            effs <-   as.numeric(myqtltable[myqtltable[,'chr']==chro,'add eff'])
            traitE <- myqtltable[myqtltable[,'chr']==chro,'trait']
            addEff <- rep('Eff: A(-)',length(effs))
            if (sum(effs>0,na.rm=T)>0)
                addEff[effs>0] <- 'Eff: B(+)'
            traitnames <-traitnameV
            
            if (plotparlist$pmethod=="png"|plotparlist$pmethod=="svg")
            {
                print('plotting')
                pln <- QTLSingleFilename(plotparlist,plotno,traitnames,chro)
###    pln <- QTLSingleFilename(plotparlist,plotno,traitnames,chro,ext='png')
            }
            traitnames[traitnames%in%traitE] <-  paste(traitnames[traitnames%in%traitE],addEff)
          par(mar=c(0,4,2,2))###c(b,l,t,r)
          resultL <- plotQTLs(traceChrDfr,pmain,traitnames,chro,colvec,printlwd,ymax,thresholdV)
          realmarkerL <- resultL[[1]]
          chrQTs <- resultL[[2]]
          par(mar=c(4,4,0,2))###c(b,l,t,r)
          plot(chrQTs[[1]],chr=chro,col="white",lwd=printlwd,ylim=c(0,1),xlab="",ylab="",axes=F,frame.plot=T,mtick="none")
          axis(1,xlab="")
          CIL <- findCIs(traceChrDfr,chro,myqtltable,colvec) ### does this work
          ypos <- par('usr')[3]
          yadd <- (par('usr')[4]-ypos)/50
          CIV <- c()
          if(dim(CIL[[1]])[1]!=0)
          {
              qtlplottable <- CIL[[1]]
              bordermarker <- CIL[[2]]
              CIV <- CIL[[3]]
              for (j in 1:nrow(qtlplottable))
              {
                  ciypos <- ypos+yadd*j*2
                  lines(c(qtlplottable[j,'CIstart'],qtlplottable[j,'CIend']),c(ciypos,ciypos),col=qtlplottable[j,'col'],lwd=printlwd)
                  points(qtlplottable[j,'pos'],ciypos,col="black",pch=20,lwd=3)
              }
              ypos <- ypos+ par('usr')[4]/4 ## 25% for CIs
              abline(h=ypos)
              yadd <- ypos+ par('usr')[4]/10
              realmarkerL$markercol[row.names(realmarkerL$realmark)%in%bordermarker] <- 'black'
          }
          plotMarkerSegments(realmarkerL$realmark,realmarkerL$markerpos,realmarkerL$markercol,ypos,yadd,start=0)
          title(xlab=paste('Chr',chro,'map positions (cM)'),line=2)
#          title(sub=paste("QTL ",paste(CIV,collapse="|")," (",Sys.Date(),")",sep=""),line=3) ### method ?
          mtext('CI',side=2,at=0.1,line=1,cex=1.5)
          mtext('marker',side=2,at=0.6,line=1,cex=1.5)
          }
      }
    }


CImarkers <- function(ci,scanout)
### find markers in CI, excluding the border markers
{
    start <- which(row.names(ci)[1]==row.names(scanout))
    end <- which(row.names(ci)[3]==row.names(scanout))
    CImarkers <- row.names(scanout)[(start+1):(end-1)]
    CImarkers <- CImarkers[-grep(paste(ci[1,'chr'],'.loc',sep=""),CImarkers)]### removing the calculated
  return(CImarkers)
}

###MyCrossL <- qtlAnalysis(MyCross,MyCrossTotal,gtfn,outdir,pmethod=pmethod,epistasis=epistasis,traitDfr=traitDfr,alpha=alpha,alphaCIM=alphaCIM,testmode=F,LODthreshold=LODthreshold)
###qtlmethods=c("hk")
qtlAnalysis <- function(MyCross,gtfn,outdir='.',pmethod='pdf',epistasis=F,qtlmethods=c("hk"),traitDfr="",alpha=alpha,alphaCIM=alphaCIM,testmode=F,printlwd=2,LODthreshold=NA,saveLODFlag=TRUE)
### Function to perform different qtl analysis steps.
### Summary will be written into a "qtl_table" csv table.
### Methods: Haley-Knott and imputation: qtlmethods=c("hk","imp")
{
    if (saveLODFlag)
        if (is.na(file.info(paste(outdir,'csv/lodfiles/',sep='/'))$isdir)) ### is lodfile  directory present present?
            dir.create(paste(outdir,'csv/lodfiles/',sep='/'))### create it.
  print('Start of the QTL scan')
  finalqtlT <- c() ### Container to collect all relevant qtl information
  finalEpistasisT <- c() ### Container to collect all relevant epistasis
  if (!is.null(ncol(MyCross$pheno)))
    if (ncol(MyCross$pheno)>1)
      MyCross$pheno <- MyCross$pheno[,sort(colnames(MyCross$pheno))] ### alphabetical order of traits
  pmain <- sub("^.*/","",gtfn)
  pmain <- sub("_rqtl_genotypes.*","",pmain) ### plot title
  print(pmain)
  devicelist <- c() ### plotting into several device in parallel
  plotparlist <- c()
  if (pmethod=='pdf'|pmethod=='png') ### plotting options 
    {
      filetype=paste('.',pmethod,sep="")
      printlwd=2
      if (pmethod=='png')
        printlwd=6
      if (length(dev.list())>0)
        for (i in 1:length(dev.list()))
          dev.off()
      filenamevector <- makeFilename(paste(outdir,"plots",sep="/"),paste(pmain,"_qtls",sep=""),filetype)### set plot file names
      fnexpansion <- c("qtl_genome_overview","qtl_additive_effects","qtl_map_location","phenotype_histograms","epistasis_2Dplots")
      if(epistasis==F)
          fnexpansion <- fnexpansion[-grep("epistasis",fnexpansion)]
      subfn <- function(new,filenamevector){return(sub("qtls",new,filenamevector))}
      filenamevector <- c(sapply(fnexpansion,subfn,filenamevector),filenamevector)
      names(filenamevector)[length(filenamevector)] <- "qtls"
      if (pmethod=='pdf')
        for (fn in filenamevector)
          openPDF(fn,n=1,p='a4',h=10,w=10,pointsize=10)
      devicelist <- dev.list()
      names(devicelist) <- c(fnexpansion,'qtls')
      plotparlist <- list(filenamevector=filenamevector,devicelist=devicelist,pmethod=pmethod)
    }
  MyCross <- calc.genoprob(MyCross, step=1, error.prob=0.01) ### markers at 1cM steps
  MyCross <- sim.geno(MyCross, step=1,n.draws=16, error.prob=0.01) ### simulate genotype probabilities, 1cM steps
###    MyCrossbak <- MyCross
    MyCross$pheno <- MyCross$pheno[,apply(!is.na(MyCross$pheno),2,sum)>5,drop=FALSE] ### remove phenotypes with less than 5 NA ### drop=FALSE allows a single column to be selected
  traitNos <- (1:ncol(MyCross$pheno))[round(apply(MyCross$pheno=='',2,sum,na.rm=T)/nrow(MyCross$pheno))!=1] ### remove empty trait columns
###  traitNos <-grep('GPD',colnames(MyCross$pheno)) ### if only GPD traits wanted
  if ('id'%in%names(MyCross$pheno[traitNos])) ### should be removed earlier - just in case
      names(MyCross$pheno)
    traitNos <- traitNos[!(names(MyCross$pheno)[traitNos])%in%"id"]
  if (sum(grepl('^EM_',names(MyCross$pheno[traitNos])))>0) ### to remove emergence data
      traitNos <- traitNos[-grep('^EM_',names(MyCross$pheno[traitNos]))]
#  if (sum(grepl('^AWNS',names(MyCross$pheno[traitNos])))>0) ### to remove emergence data
#    traitNos <- traitNos[-grep('^AWNS',names(MyCross$pheno[traitNos]))]
  test <- try(out.all <- scanone(MyCross, method="hk",pheno.col=1:length(traitNos)),silent=TRUE)
  if (class(test)[1]=="try-error")
    print("no early qtl summary possible",quote=FALSE)
  else
    {
      testThreshold=2.5 ### arbitrary threshold for an early QTL summary 
      if (!is.na(LODthreshold))
          testThreshold <- LODthreshold
      print(paste("Early QTL summary using Haley-Knott algorithm and an arbitrary threshold of ",testThreshold,".",sep=""),quote=FALSE)
      print(summary(out.all,threshold=2.5,format="allpeaks",quote=FALSE))
    }
  qmcolvec <- topo.colors(length(qtlmethods)+1)
  names(qmcolvec) <- c(qtlmethods,'cim')
  QTLtraceL <- list() ### list to collect all QTL LOD traces.
    trNo=2
  for (trNo in traitNos)
    {
      print("--------------------------------------------------------------------",quote=FALSE)
      print(paste('QTL of trait',names(MyCross$pheno[trNo])),quote=FALSE)
      print("--------------------------------------------------------------------",quote=FALSE)
      trait <- names(MyCross$pheno)[trNo]
      environment=pmain ### better definition of environment
      if (grepl('_',trait)) ### if each trait has the env attached, use this
      {	
          environment <- sub('^[^_]*_','',trait)
          pmain <- sub('_.*','',pmain) ### pmain will also have the env attached. discard this.
      }
      traitmodel='normal'
      if (length(levels(as.factor(MyCross$pheno[[trNo]])))<6)
        traitmodel='np' ### non-parametric for non-normals
      if (length(levels(as.factor(MyCross$pheno[[trNo]])))<2)
      {
          print(MyCross$pheno[[trNo]])
          print(levels(as.factor(MyCross$pheno[[trNo]])))
          if (levels(as.factor(MyCross$pheno[[trNo]]))%in%c(0,1))
              traitmodel='binary' ### for binary count data
      }
      print(paste('Assumed distribution for trait:',traitmodel),quote=F)
      print("Single QTL scan first.",quote=FALSE)
      if(sum(is.na(MyCross$pheno[,trNo]))<(nrow(MyCross$pheno)-2))
        {
            out <- list()
            for (m in qtlmethods) ### QTL scan can be done with different algorithms
                out[[m]] <- scanone(MyCross, method=m,pheno.col=trNo,model=traitmodel)
            alphaInitial=alpha
            if (traitmodel=='np') ### no CIM possible. Strict alpha for scan.
              alphaInitial=alphaCIM
            operm.a <- try(scanone(MyCross, method=qtlmethods[1], n.perm=1000,model=traitmodel,pheno.col=trNo))### permutation for LOD thresholds
            if (class(operm.a)[1]=="try-error")
              {
                 print("scanone fails on this dataset")
                 break
               }
            threshold <- round((summary(operm.a, alpha=alphaInitial))[1],2) ### alpha set above, default 0.4, not stringent
            if (!is.na(LODthreshold))
                if (!is.na(LODthreshold))
                    {
                        print(paste("The pre-set threshold of ",LODthreshold," will be used. It is lower than the threshold calculated by the distribution (",threshold,") and any QTLs detected should be used with caution.",sep=""),quote=FALSE)
                        threshold <- LODthreshold
                    }
            summary.a <- summary(out[[qtlmethods[1]]], perms=operm.a, threshold=threshold,pvalues=TRUE)
            outMax <- round(max(out[[m]]$lod),2)
            print(paste("LOD threshold: ",threshold,"; max LOD values (",paste(qtlmethods,collapse='/'),"): ",outMax," (alpha=",alphaInitial,")",sep=""),quote=FALSE)
          if (is.infinite(threshold)|(is.na(threshold)))
            print("Threshold is not finite or NA, no QTLs detected.",quote=FALSE)
          else
          {
              out.use <- out[[1]]
              pln <- QTLplotFilename(plotparlist,plotname="qtl_genome_overview") ### set device to print
              ymax <- max(c(threshold,out.use[,3]),na.rm=TRUE)
              plot(out.use,main=paste("QTL Overview 1st scan:",names(MyCross$pheno)[trNo],"-",pmain),sub=Sys.Date(),lwd=printlwd,col=qmcolvec[1],ylim=c(0,max(ymax,10)))
              if(length(qtlmethods)>1)
                for (p in (2:length(qtlmethods)))
                  plot(out[[p]],col=qmcolvec[p],lwd=printlwd,add=TRUE)
              legend("topleft",c(qtlmethods),col=qmcolvec[1:length(qtlmethods)],lty=1)
              abline(h=threshold,lty=2)
              text(x=200,y=threshold,labels=alphaInitial)
              #### overview plot finished
              if (dim(summary.a)[1]>0)### qtl present ..
                {
                  chrs <- as.character(summary.a$chr)
                  chrsunique <- chrs[!duplicated(chrs)]
                  qtlpos <- summary.a[,2][!duplicated(chrs)]
                  trialqtls <- makeqtl(MyCross,chr=chrsunique,pos=qtlpos)### stores the qtls
                  if (traitmodel!='np')
                    {
                      print(paste("CIM taking QTL:",paste(chrsunique,collapse=', '),'into account.'),quote=FALSE)
                      markern <- length(chrsunique)
                      out.use <- try(cim(MyCross,pheno.col=trNo,n.marcovar=markern,window=20,method='hk'))
                      if(any("try-error"%in%class(out.use)))  ### catch an error
                      {
                          print("an errer occured in cim")
                          QTLtraceL[names(MyCross$pheno[trNo])] <- list(list("summary"=summary.a,"trace"=out[[1]],"qtltable"=c(),"threshold"=threshold,"effects"=c()))
                          break
                      }
                      thresholdCIM <- round((summary(operm.a, alpha=alphaCIM))[1],2) ### alpha stringent, 0.02
                      outMaxCIM <- round(max(out.use$lod),2)
                      print(paste("CIM LOD threshold: ",thresholdCIM,"; max LOD values: ",outMaxCIM," (alpha=",alphaCIM,")",sep=""),quote=FALSE)
                      summary.cim <- summary(out.use,threshold=thresholdCIM-0.5)
                      myalpha=alphaCIM
                      if (nrow(summary.cim)==0)
                          for (k in 1:3)
                          {
                              myalpha=2*myalpha
                              thresholdCIM <- round((summary(operm.a, alpha=myalpha))[1],2)
                              summary.cim <- summary(out.use,threshold=thresholdCIM)
                              if (nrow(summary.cim)>0)
                              {
                                  print(paste('QTL found for lower significance level alpha:',myalpha),quote=F)
                                  print(summary.cim,quote=F)
                                  break
                              }
                          }
                      pln <- QTLplotFilename(plotparlist,plotname="qtl_genome_overview") ### direct to file
                      plot(out.use,out[[m]],main=paste("CIM QTL overview:",names(MyCross$pheno)[trNo],"-",pmain),sub=Sys.Date(),lwd=printlwd,col=qmcolvec[c('cim',m)],ylim=c(0,max(ymax,10)))
                      add.cim.covar(out.use)### print add-covars on plot
                      legend("topleft",c('cim',m),col=qmcolvec[c('cim',m)],lty=1)
                      abline(h=thresholdCIM,lty=2)
                      myx=max(out.use[,2])
                      myx <- myx-myx/20
                      text(x=myx,y=thresholdCIM+0.3,labels=paste('CIM',myalpha))### this needs to be at max x
                      if (threshold!=thresholdCIM)
                      {
                          abline(h=threshold,lty=2)
                          text(x=myx,y=threshold-0.2,labels=alphaInitial)
                      }
                      if(nrow(summary.cim)>0)
                      {
                          chrs <- as.character(summary.cim$chr)
                          chrsunique <- chrs[!duplicated(chrs)]
                          qtlpos <- summary.cim[!duplicated(chrs),2]
                          trialqtls <- makeqtl(MyCross,chr=chrsunique,pos=qtlpos)### stores the qtls
                      }
                      ymax <- max(c(out.use[,3],thresholdCIM),na.rm=TRUE)
                    }
                  pln <- QTLplotFilename(plotparlist,plotname="qtl_additive_effects") ### direct to file
                  pylim <- range(MyCross$pheno[trNo],na.rm=TRUE)
                  effectS <- effectscan(MyCross,pheno.col=trNo,chr=chrsunique,get.se=TRUE,show.marker=TRUE,add.legend=TRUE,sub=paste("(",Sys.Date(),")",sep=""),main=paste("QTL-effects:",names(MyCross$pheno)[trNo],"-", pmain),ylim=pylim) ### a=additive, d=dominace
### parent A allele increasing if effect is negative, B allele increasing if positive
                  pln <- QTLplotFilename(plotparlist,plotname="qtl_map_location") ### direct to file
                  plot.qtl(trialqtls,main=paste("QTL location:",names(MyCross$pheno)[trNo],"-", pmain))
                  qtlmodel <- modelformula(length(chrsunique),"+") ### pure additive model
                  print(paste("qtlmodel: ",qtlmodel))
                  MyCross <- sim.geno(MyCross, step=4,n.draws=16, error.prob=0.01)### step 2, maybe slow
                  if(traitmodel=='np')
                      traitmodel='normal'
                  refinedQTL <- refineqtl(MyCross,pheno.col=trNo,trialqtls,formula=formula(qtlmodel),model=traitmodel)
                  
                  
                  lod.model <- fitqtl(MyCross,pheno.col=trNo,refinedQTL,formula=formula(qtlmodel),model=traitmodel)### multiple scan - p-values very low!
                  print(summary(lod.model),quote=F)### fitqtl summary
                  trialqtls <- refinedQTL
###added functionality to save the LOD scores to file
                  if(saveLODFlag)
                  {
                      saveLOD <-out.use
                      saveLOD <- cbind(saveLOD,out.use)
                      colnames(saveLOD)[c(3,4)] <- c('lod','lod refined')
                      newqtl <- attributes(trialqtls)[[4]]
                      for(i in 1:length(newqtl))
                          saveLOD[row.names(newqtl[[i]]),'lod refined'] <-  newqtl[[i]][,3]
                      lodfilename <- makeFilename(paste(outdir,'csv/lodfiles/',sep='/'),paste('QTL_LOD_scores_',pmain,'_',trait,sep=""),'.csv')
                      write.csv(saveLOD,lodfilename,quote=F)
                  }
                  print('Writing qtl summary table.',quote=F)
                  qtlcolno=19
                  if (trialqtls$n.qtl==1) ### one qtl present
                    {
                        lodtable <- round(lod.model$result.full["Model",c("LOD","%var","Pvalue(F)")],3)
                      qtltable <- matrix(c(lodtable,rep(NA,16)),ncol=qtlcolno)
                      row.names(qtltable) <- trialqtls$name[1]
                    }
                  else ### more qtls
                    {
                      lodtable <-lod.model$result.drop[,c("LOD","%var","Pvalue(F)")]
                      lodtable <-cbind(round(lodtable[,1:2],1),round(lodtable[,3],4))
                      lodtable <- lodtable[rownames(lodtable)%in%trialqtls$name,]
                      qtltable <- matrix(lodtable,ncol=3)
                      rowno <- nrow(qtltable)
                      qtltable <- cbind(qtltable,matrix(NA,nrow=rowno,ncol=(qtlcolno-ncol(qtltable))))
                      row.names(qtltable) <- row.names(lodtable)
                    }
                  for (l in 1:nrow(qtltable)) ### each QTL at a time
                  {
                      addeff <- findAdditiveEffect(effectS,chrsunique[l],qtlpos[l])
                      qtltable[l,4:8] <- addeff
                      ci <- mylodint(out.use,chr=chrsunique[l],qtlmarker=qtltable[l,6],drop=1.5,expandtomarkers=TRUE)### CI with 1.5 LOD scores
                      CIm <- CImarkers(ci,out.use)
                      qtltablerow <- c(sub("/.*","",chrsunique[l]),sapply(ci$pos,round,3),rownames(ci)[c(1,3)],trait,environment,pmain,round(colMeans(MyCross$pheno[trNo],na.rm=T),3))
                      if (length(qtltablerow) > 10 )### QTL worked
                          qtltablerow <- c(sub("/.*","",chrsunique[l]),sapply(ci$pos,round,3)[c(1:3)],rownames(ci)[c(1,3)],trait,environment,pmain,round(colMeans(MyCross$pheno[trNo],na.rm=T),3))
                      if (length(qtltablerow) < 10 ) ### no QTL found
                        qtltablerow <- c(qtltablerow,rep(NA,9-length(qtltablerow)),round(colMeans(MyCross$pheno[trNo],na.rm=T),3))
                      qtltable[l,9:18] <- qtltablerow ###
                      qtltable[l,19] <- paste(CIm,collapse='|')
                    }
                  cns <- c("LOD","%var","Pvalue(F)","add eff","se.a","sim marker","pos nearest marker","nearest marker","chr","CIstart","pos","CIend","start marker","end marker","trait","env","population","mean","CI markers")
                  colnames(qtltable) <- cns
                  QTLtraceL[names(MyCross$pheno[trNo])] <-list(list("summary"=summary.a,"trace"=out.use,"qtltable"=qtltable,"threshold"=threshold,"effects"=qtltable[,'add eff']))
                  newqtltable <- data.frame(matrix(qtltable[,c("chr","pos","LOD","%var","mean","add eff","sim marker","nearest marker","pos nearest marker","CIstart","CIend","start marker","end marker","Pvalue(F)","trait","env","population","CI markers")],ncol=18),stringsAsFactors=F)
                  colnames(newqtltable) <- c("chr","pos","LOD","%var","mean","add eff","sim marker","nearest marker","pos nearest marker","CIstart","CIend","start marker","end marker","Pvalue(F)","trait","env","population","CImarkers")
                  print("--------------------------------------------------------------------",quote=FALSE)
                  print(paste("QTL summary table:",names(MyCross$pheno[trNo]),'in population',sub("_.*$","",pmain)),quote=F)### need marker in CI
                  print("",quote=FALSE)
                  print(newqtltable[,1:16],ncol=16,quote=F)
                  print("",quote=FALSE)
                  print("--------------------------------------------------------------------",quote=FALSE)
                  rownames(newqtltable) <- paste(rownames(qtltable),qtltable[,15],sep="_")
                  totalrownames <- make.unique(c(row.names(finalqtlT),rownames(qtltable)))
                  finalqtlT <- rbind(finalqtlT,newqtltable)
                  if((epistasis)&(nrow(newqtltable)>1))
                  {
                      print("Test for epistasis is in testing phase. Epistasis scanning takes as long time.",quote=FALSE)
                      out2.hk <- scantwo(MyCross, method="hk",pheno.col=trNo)
                      pln <- QTLplotFilename(plotparlist,plotname="epistasis_2Dplots") ### direct to file           
                      plot(out2.hk,main=paste("QTL 2D interaction:",names(MyCross$pheno)[trNo],"-",pmain),sub=Sys.Date())

                      dfr2D <- cbind(trait,summary(out2.hk),stringsAsFactors=FALSE)
                      dfr2D$chr1 <- as.character(dfr2D$chr1)
                      dfr2D$chr2 <- as.character(dfr2D$chr2)
                      {
                          operm2 <- scantwo(MyCross, method="hk", n.perm=2)
                          print("summary(operm2, alpha=0.05)")
                          print(summary(operm2, alpha=0.05))
                          sumperm2 <- summary(operm2,alpha=0.05)
                          sumperm2 <- sumperm2[-length(sumperm2)]
                          names(sumperm2) <- paste('lod',names(sumperm2),sep='.')
                          dfr2D[nrow(dfr2D)+1,] <- ''
                          dfr2D[nrow(dfr2D),names(sumperm2)] <- sumperm2
                          dfr2D[nrow(dfr2D),c("trait","chr1","chr2")] <- c(trait,'permut','permut')
                      }
                      epifilename <- makeFilename(paste(outdir,'csv',sep='/'),paste('epistasis_scores_',pmain,trait,sep=""),'.csv')
                      finalEpistasisT <- rbind(finalEpistasisT,dfr2D)
                      write.csv(dfr2D,epifilename,row.names=T)
                  }
                }
              else
                  QTLtraceL[names(MyCross$pheno[trNo])] <- list(list("summary"=summary.a,"trace"=out[[1]],"qtltable"=c(),"threshold"=threshold,"effects"=c()))
              rm(out.use)
            }
          }
      else
        print('not enough values for QTL analysis present',quote=FALSE)
    }
  pln <- QTLplotFilename(plotparlist,plotname="phenotype_histograms") ### direct to file           
  for (trNo in traitNos)
    {
      if ((names(MyCross$pheno[trNo])!="id")& (sum(!(is.na(MyCross$pheno[trNo])))))
          plotPheno(MyCross,pheno.col=trNo,sub=paste(pmain,' (',Sys.Date(),')',sep=""))
    }
  if(length(finalqtlT)>0)### add trait definitions to qtl summary table
      if (nrow(traitDfr)>0)### Was the controlled vocabulary uploaded?
          if (length(traitDfr)>0) ### 
              {
                  finalqtlT <- as.data.frame(cbind(finalqtlT,rep('',nrow(finalqtlT)),rep('',nrow(finalqtlT))),stringsAsFactors=FALSE)
                  colnames(finalqtlT)[(ncol(finalqtlT)-1):ncol(finalqtlT)] <- c('longtraitnames','units')
                  if(length(grep('_',finalqtlT$trait))>0)
                  {
                      ext <- sub('.*_','',finalqtlT$trait)
                      if (identical(ext,finalqtlT$env))
                          finalqtlT$trait <- sub('_.*','',finalqtlT$trait)
                  }
                  for (trait in unique(finalqtlT[,'trait']))
                      if (trait%in%traitDfr[,'short_name'])### is trait abbreviation in the vocabulary file?
                          {
                              finalqtlT[finalqtlT[,'trait']==trait,'longtraitnames']  <-  traitDfr[traitDfr[,'short_name']==trait,'long_trait_name']
                              finalqtlT[finalqtlT[,'trait']==trait,'units'] <- traitDfr[traitDfr[,'short_name']==trait,'unit_abbreviation'] ### add the units of measurements
                          }
              }
  increaseAll <- as.character(as.numeric(as.numeric(finalqtlT[,"add eff"])>0)) ### additive effect
  increaseAll[increaseAll=='1'] <- 'B' ### positive effect: increasing effect comes from  B
  increaseAll[increaseAll=='0'] <- 'A'### negative effect: increasing effect comes from A
  finalqtlT <- cbind(finalqtlT,increaseAll) ### add to qtl summary table
  pmain <- sub("^.*/","",gtfn)
  pmain <- sub("_rqtl_genotypes.*","",pmain) ###  plot title
  qtlfilename <- makeFilename(paste(outdir,'csv',sep='/'),paste('qtl_table_',pmain,sep=""),'.csv')
  if (length(dev.list())>0) ### close the printing devices
      for (dv in 1:length(dev.list()))
          dev.off()
  if (nrow(finalqtlT)>0) ### results output if QTL present.
  {
      print(paste('The QTL summary table will be written to file:',qtlfilename),quote=F)
      print(head(finalqtlT[,c(1,2,3,4,6,8,10,14,15,16,18,19)]),quote=F)
      write.csv(finalqtlT,qtlfilename,row.names=T) ### qtl summary
      MyCross$qtltable <- finalqtlT
      pmain <- sub("^.*/","",gtfn)
      pmain <- sub("_rqtl_genotypes.*","",pmain) ### plot title 
      filenamevector <- makeFilename(paste(outdir,"plots",sep="/"),paste(pmain,"_qtls",sep=""),filetype)
      names(filenamevector)[length(filenamevector)] <- "qtls"
      if (pmethod=='pdf')
          for (fn in filenamevector)
              pdf(fn,paper='a4',pointsize=12,height=14)
      devicelist=dev.list()
      names(devicelist) <- 'qtls'
      plotparlist <- list(filenamevector=filenamevector,devicelist=devicelist,pmethod=pmethod)
      if (pmethod=='pdf')
      {
          if (length(QTLtraceL)>0) ### plot the single QTL plots
              try(allRQTLplots(QTLtraceL,pmain,plotparlist,printlwd=printlwd,reducetraitF=FALSE))
          }
  }
  if (epistasis)
      {
          epifilename <- makeFilename(paste(outdir,'csv',sep='/'),paste('epistasis_scores_',pmain,sep=""),'.csv')
          write.csv(finalEpistasisT,epifilename,row.names=T)
      }

  return(list(MyCross,QTLtraceL))
}

wrapLoadFromData <- function()
{
    dir()
    rm(list=ls())
    source('analyse_qtl.R')
    ##fn='RobClaFieldData2019.RData'
    ##plotFromRData(fn)
###plotFromRData(fn=fn,pmethod=pmethod)
###    fn='panel5CFis.RData'
###    fn='ParWat_F4_2020.RData
###    fn='ParWat912_andParSta_F5_2020.RData'
###    fn='ParCha_ParJCa.RData'
####    fn='ParBaj_JIC.RData'
###    fn='parbaj.RData'
###        fn='../rqtl_data/ParW238_Fungal_Phylum.RData'
###    fn='ParW238.RData'
    fn='ParGar.RData'
    pmethod='svg'
##    pmethod='png'
    print(pmethod)
    plotFromDataSingle(fn=fn,pmethod=pmethod)
}
###plotFromRData(fn="ParW238.RData",pmethod='pdf') ###
###plotFromRData(fn="../rqtl_data/ParW238_Bacteria_Phylum.RData",pmethod='pdf') ### 
####wrapLoadFromData()

plotFromRData <- function(fn=".RData",pmethod='pdf')
{
    print('Function: plotFromRData')
###    fn='parbaj.RData'
###    fn='../rqtl_data/ParW238_Bacteria_Phylum.RData'
###    fn='ParW238.RData'
###    fn='ParGar.RData'
    load(fn)
    mylist=ls()
    mylist=mylist[!mylist%in%c('allqtlL','outdir')]
    rm(list=mylist)
    ls()
    source('analyse_qtl.R')
    pmethod='pdf'
    if (pmethod=="pdf")
        filetype=".pdf"
    printlwd=3
    i=1
    for(i in seq(1,length(allqtlL),by=2))
    {
        pmain <- allqtlL[[i]][[3]]$pop[1]
        filenamevector <- makeFilename(paste(outdir,"plots",sep="/"),paste(pmain,"_qtls",sep=""),filetype)
        names(filenamevector)[length(filenamevector)] <- "qtls"
        for (fnv in filenamevector)
        {
            if(pmethod=='pdf')
                pdf(fnv,paper='a4',pointsize=12,height=14)
        }
        devicelist=dev.list()
        names(devicelist) <- 'qtls'
        plotparlist <- list(filenamevector=filenamevector,devicelist=devicelist,pmethod=pmethod)
        QTLtraceL <- allqtlL[[i+1]]
### for ParGar drought trial
        print(names(QTLtraceL))
        names(QTLtraceL) <- sub('_D*2','_DT2',names(QTLtraceL))
        names(QTLtraceL) <- sub('([IN][RI])_DT','_\\1-DT',names(QTLtraceL))
        print(names(QTLtraceL))
###        QTLtraceL <- QTLtraceL[names(QTLtraceL)[grep('JIC',names(QTLtraceL))]] ### for priya's paper review
###        names(QTLtraceL) <-  sub('^([^_]*)_([^_]*)$','\\2_\\1',names(QTLtraceL)) ###reverse bacteria vs phylum name
        if (length(QTLtraceL)>0) ### plot the single QTL plots
          allRQTLplots(QTLtraceL,pmain,plotparlist,printlwd=printlwd,reducetraitF=FALSE)
    }
}

###plotFromRData(fn="ParGar.RData",pmethod='pdf')
    
plotFromDataSingle <- function(fn=fn,pmethod='svg',printlwd=5)
{
### For png and svg plots, e.g. for the QTL for CerealsDB
### using function QTLSingleFilename in allRQTLplots
    print('Function: plotFromDataSingle')
    for(d in dev.list())
        dev.off(d)
    load(fn)
    mylist=ls()
    mylist=mylist[!mylist%in%c('allqtlL','outdir')]
    rm(list=mylist)
    source('analyse_qtl.R')
    if (pmethod=="jpg")
        filetype=".jpg"
    if (pmethod=="svg")
        filetype=".svg"
    if (pmethod=="png")
        filetype=".png"
    printlwd=6
    i=1
    for(i in seq(1,length(allqtlL),by=2))
    {
        pmain <- allqtlL[[i]][[3]]$pop[1]
        pmain <- sub('_rqtl_JIC_genotypes.csv','',pmain)
        filenamevector <- makeFilename(paste(outdir,"plots",sep="/"),paste(pmain,"_qtls",sep=""),filetype)
        names(filenamevector)[length(filenamevector)] <- "qtls"
        filenamevector
        for (fnv in filenamevector)
        {
            if(pmethod=='svg')
                svg(fnv,width=14,height=14,pointsize=12)
            if(pmethod=='png')
                png(fnv)
            if(pmethod=='jpg')
                jpg(fnv,width=14,height=14,pointsize=12)
        
        }
        devicelist=dev.list()
        names(devicelist) <- 'qtls'
        plotparlist <- list(filenamevector=filenamevector,devicelist=devicelist,pmethod=pmethod)
        QTLtraceL <- allqtlL[[i+1]]
        print(names(QTLtraceL))
        if (length(QTLtraceL)>0) ### plot the single QTL plots
            allRQTLplots(QTLtraceL,pmain,plotparlist,printlwd=printlwd)
    }
}



allRQTLplots <- function(QTLtraceL,pmain,plotparlist,printlwd=printlwd,reducetraitF=FALSE)
{
    traitenvnames <- names(QTLtraceL)
    traitenvnames
    traitnames <- sub('_.*','_',traitenvnames)
    envs <- sort(unique(sub('.*_','',traitenvnames)))
    envs
    if (sum(traitnames!=traitenvnames)>0)
    {	
        print('The following environments are present:',quote=F)
        print(envs,quote=F)
    }
    traitnames <- unique(traitnames)
    traitnames <- traitnames[order(nchar(traitnames))]
    colvec <- rainbow(length(envs)*2+1)[seq(1,length(envs)*2+1,by=2)]
    colvec <- colvec[1:length(envs)]
    names(colvec) <- envs
###    barplot(rep(5,length(colvec)),col=colvec,beside=true)
    if(length(envs)==1&1==2)
    { ### I don't think this works
        colvec <- rainbow(length(traitnames)*2+1)[seq(1,length(traitnames)*2+1,by=2)]
        colvec <- colvec[1:length(traitnames)]
        names(colvec) <- traitnames
    }
    print('colvec')
    print(colvec)
    if(reducetraitF)### combine similar traits in one plot
    {
        reducetraits <- c()
        for (tr in traitnames)
            if (length(grep(tr,traitnames))>1)
                reducetraits <- rbind(reducetraits,cbind(tr,traitnames[grep(tr,traitnames)]))
        reducetraits <- data.frame(matrix(reducetraits,ncol=2),stringsAsFactors=FALSE)
        reducetraits <- reducetraits[reducetraits[,1]!=reducetraits[,2],]
        print('trait names plotted together')
        print(reducetraits)
        traitnames <- traitnames[!traitnames%in%reducetraits[,1]]
    }
    print('Plotting QTL for the following traits:',quote=F)
    for (tr in traitnames)
    {
        print(tr,quote=F)
        trgrep <- paste('^',tr,sep='')
        qtltraceL <- QTLtraceL[grep(trgrep,names(QTLtraceL))]
        print('qtltraceL names')
	print(names(qtltraceL))
        try(plotRQTLtraces(qtltraceL,paste(pmain,'-',sub('_$','',tr)),colvec=colvec,plotparlist,printlwd=printlwd))
    }
    dev.off()
}

findQTLforQTLmodel <- function(qtl,intsummary,MyCross)
### find the name of a qtl in a R/QTL qtl object
{
  nameQTL <- function(qtlobj,chr,pos)
### finds the name of the QTL in the qtl object
    {
      qtlname <- NA
      qtlobjpos <- qtlobj$pos[qtlobj$chr%in%chr]
      for (p in qtlobjpos)
        if ((pos>=p-10)&(pos<=p+10))
          qtlname <- qtlobj$altname[(qtlobj$chr%in%chr)&(qtlobj$pos%in%p)]
      return(qtlname)
    }
  findInteractions <- function(qtl,intsummary)
  {
    interactingQTL <- list()
    for (im in 1:nrow(intsummary))
      {
        firstQTL <- nameQTL(qtl,as.character(intsummary[im,1]),intsummary[im,3])
        secondQTL <- nameQTL(qtl,as.character(intsummary[im,2]),intsummary[im,4])
        if (is.na(firstQTL[1]))
          {
            firstQTL <- paste("Q",as.numeric(sub("Q","",max(additiveQTL)))+1,sep="")
            additiveQTL <- c(additiveQTL,firstQTL)
            qtl <- makeqtl(MyCross,chr=c(qtl$chr,as.character(intsummary[im,1])),pos=c(qtl$pos,intsummary[im,3]))### stores the found qtls as qtls
          }
        if (is.na(secondQTL[1]))
          {
            secondQTL <- paste("Q",as.numeric(sub("Q","",max(additiveQTL)))+1,sep="")
            qtl <- makeqtl(MyCross,chr=c(qtl$chr,as.character(intsummary[im,2])),pos=c(qtl$pos,intsummary[im,4]))### stores the found qtls as qtls
            additiveQTL <- c(additiveQTL,secondQTL)
          }
        interactingQTL <- c(interactingQTL,list(c(firstQTL,secondQTL)))
      }
    return(list(qtl,interactingQTL))
  }
  findOverlappingInteractions <- function(interactingQTL)
    {
      interactions <- c()
      if(sum(duplicated(unlist(interactingQTL)),na.rm=TRUE)>0)
        {
          duplicates <- unique(unlist(interactingQTL)[duplicated(unlist(interactingQTL))])
          interactions <- list()
          for (d in duplicates)### find out which ones interact
            {
              partner <- c()
              for (iq in interactingQTL)
                if (length(grep(d,iq)>0))
                  partner <- c(partner,iq[-grep(d,iq)])
              interactions <-c(interactions,list(c(d,partner)))
            }
        }
      return(interactions)
    }
  rebuildInteractions <- function(interactions)
    {
      startagain=FALSE
      for (iA in 1:(length(interactions)-1))
        {
          for (j in (iA+1):length(interactions))
            if (length(intersect(interactions[[iA]],interactions[[j]]))>0)
              {
                interactions[[iA]] <- union(interactions[[iA]],interactions[[j]])
                interactions[[j]] <- c()
                print(interactions[[iA]],quote=F)
                startagain=TRUE
                break()
              }
          if (startagain)
            break()
        }
      return(list(startagain=startagain,interactions=interactions))
    }

  additiveQTL <- qtl$altname ### these are all QTL names
  qtlmodel <- mymodelformula(start="y~",qtlvector=additiveQTL,m="+")
  result <- findInteractions(qtl,intsummary)
  qtl <- result[[1]]
  interactingQTL <- result[[2]]
  additiveQTL <- qtl$altname ### these are all QTL names
  interactingQTL <- findOverlappingInteractions(interactingQTL)
  repeat ### builds combined interactionstable
    {
      startagain=FALSE
      if (length(interactingQTL)>1) ### find out if interactions share QTLS
        {
          result <- rebuildInteractions(interactingQTL)
          interactingQTL <- result[[2]]
          if (length(interactingQTL)>1) ### new round only necessary if there are more layers
            startagain <- result[[1]] ### new round only necessary if there was a change to interactions
        }
      if (!startagain)
        break()
    }
  interactingQTLmodel="+"
  for (iq in interactingQTL)
      interactingQTLmodel <- paste(interactingQTLmodel,mymodelformula(start="",qtlvector=as.vector(iq),m="*"),"+",sep="")
  additiveQTL <- additiveQTL[!additiveQTL%in%unlist(interactingQTL)]
  qtlmodel <- mymodelformula(start="y~",qtlvector=additiveQTL,m="+")
  qtlmodel <- paste(qtlmodel,interactingQTLmodel,sep="")
  qtlmodel <- sub("[+]$","",qtlmodel)
  print(qtlmodel,quote=F)
  return(list(qtlmodel,qtl))
}

mymodelformula <- function(start="y~",qtlvector=c("Q1","Q2"),m="*")
{
  model <- start
  for (qi in 1:length(qtlvector))
    model <- paste(model,qtlvector[qi],m,sep="")
  model <- sub(".$","",model)
  return(model)
}

mapAnalysis <- function(MyCross,gtfn="rqtl_map",outdir,printpdf=FALSE,diagnostic=FALSE,detailed=FALSE,errorscore=TRUE)
### Map analysis first before QTL analysis starts
{
  print("Map analysis",quote=F)
###  class(MyCross)[1] <- "dh"
  pmain <- sub("[.]csv$","",gtfn)
  pmain <- sub("_rqtl_genotype.*","_map",pmain)
  pmain <- sub("^.*/","",pmain)
  wpath <- paste(outdir,"plots",sep="/")
  filename <- makeFilename(wpath,pmain,".pdf")
  if (printpdf)
    openPDF(filename,n=1,p='a4',h=42,w=10,pointsize=10)
  MyCross <- jittermap(MyCross) ### adds little noise
  plot.map(MyCross,show.marker.names=TRUE,horizontal=TRUE,main=paste("Genetic map:",pmain))
  MyCross <- drop.nullmarkers(MyCross)
  MyCross <- est.rf(MyCross)### estimate pairwise recombination frequency
  plot.rf(MyCross,sub=paste("Pairwise recombination: ",pmain),sub=Sys.Date()) ###marker alignment
  print('Allele order check:',quote=FALSE)
  chkAl <- checkAlleles(MyCross)
  if (length(chkAl)>0)
    {
      print('Marker order could possibly be optimised',quote=F)
      checkchr <- unique(as.character(chkAl[,2]))
      for (chr in checkchr)
        if (length(MyCross$geno[[chr]]$map)>3)
          {
            print(paste('Best possible marker orders on',chr),quote=F)
            ripplewindow <- min(length(MyCross$geno[[chr]]$map)-1,5)
            print(summary(ripple(MyCross,chr,window=ripplewindow)))
          }
    }
  if (diagnostic)
    geno.image(MyCross,sub=paste("Genotype data:",sub=pmain))###Graphical Genotype
###  pl.missing(MyCross)### A plot of missing genotypes
  if (detailed)
    for (pchr in names(MyCross$geno))
      {
        print(pchr,quote=F)
        plot.geno(MyCross,pchr,main=paste("Assumed crossovers:", pmain),sub=pchr)
        mtext(pchr,1,line=2)
      }
  if (errorscore==TRUE) ### will have to introduce this. calc.errorlod takes long
    {
      print('Genotype error check:',quote=FALSE)
      MyCross <- calc.errorlod(MyCross, error.prob=0.01)### this can take a while
      topErr <- top.errorlod(MyCross)
      if (length(topErr)> 0)
        {
          newmap <- est.map(MyCross, error.prob=0.01)###calculates new map with genotyping error 1%
          plot.map(MyCross, newmap,main=paste("Map comparison:", pmain),sub=Sys.Date())
          newMyCross <- replace.map(MyCross, newmap)
          newMyCross <- calc.errorlod(MyCross, error.prob=0.01)
          top.errorlod(newMyCross)
        }
    }
  if (diagnostic)
    try(plot.info(MyCross))### missing informations on chromosomes
  if (printpdf)
    dev.off()
  return(MyCross)
}

findPhenotypeGenotypeUnion <- function(phenotype,GenotypeCross,basicfilename,ptlabel)
###finds the union between genotype and phenotype file
###phenotype is a named vector with names being the line names (=genotype names)
{
  help(is.numeric)
  if (is.numeric(phenotype))
    {
      pnames <- as.character(names(phenotype))
      gnames <- colnames(GenotypeCross)
      ptindex <- pnames%in%gnames
      gtindex <- gnames%in%pnames
      porder <- orderNames(pnames[ptindex])
      gorder <- orderNames(gnames[gtindex])
      dfrlength <- sum(ptindex)
      filestem <- paste(basicfilename,"cross",ptlabel,sep="_")
      filename <- paste(filestem,".csv",sep="")
      phenotypematrix <- cbind(pnames[ptindex][porder],phenotype[ptindex][porder])
      colnames(phenotypematrix)=c("id",ptlabel)
      write.csv(phenotypematrix,row.names=FALSE,quote=FALSE,file.path("../R_out/csv",filename))
      genotypematrix <- t(cbind(GenotypeCross[,1:3],GenotypeCross[,gtindex][,gorder]))
      filename <- sub("_rqtl_","_rqtl_genotype_",filename)
      genotypematrix <- cbind(c("id","","",gnames[gtindex][gorder]),genotypematrix)
      write.table(genotypematrix,col.names=FALSE,row.names=FALSE,quote=FALSE,file.path("../R_out/csv",filename),sep=",")
      print(paste("wrote genotype and phenotype file for QTL analysis:",ptlabel))
    }
  else
    print("Data are not numeric. I will not proceed with them")
}

standardiseNames <- function(namevector)
{
  for (nii in 1:length(namevector))
    if (length(grep("x[a-z]*$",namevector[ni],ignore.case=TRUE))>0)
      namevector[ni] <- toupper(substr(namevector[ni],1,4))
  return(namevector)
}

orderNames <- function(namevector)
### sorts line names starting with "X" followed by a number in numeric order
### parentlines, not just numbers, are sorted in alphabetical order after
### the progeny lines
{
  nv <- sub("X","",namevector)
  nv <- as.numeric(nv)
  ordernv <- order(nv)
  ordernonnumnv <- order(namevector[is.na(nv)])
  ordernv[is.na(nv)] <- ordernv[is.na(nv)][ordernonnumnv]
  return(ordernv)
}

reduceMap <- function(dhcross,markerdist=0.5)
### Large Maps will slow down QTL mapping but will not improve resolution.
### Marker numbers will be reduced before analysis.
{
  for (chrno in length(dhcross$geno):1)
    {
      if(length(dhcross$geno[[chrno]]$map)<3)
        {
          print(paste('Will drop LG',names(dhcross$geno[chrno]),'with markers:',paste(names(dhcross$geno[[chrno]]$map),collapse=',')))
          dhcross <- drop.markers(dhcross,names(dhcross$geno[[chrno]]$map))
        }
      else
        {
          throw.markers <- c()
          j.old=-1
          for(j in 1:length(dhcross$geno[[chrno]]$map))
            {
              m <- dhcross$geno[[chrno]]$map[j]
              if (j.old>0)
                {
                  if (m<(m.old+markerdist))
                    {
                      if (sum(is.na(dhcross$geno[[chrno]]$data[,j]))>sum(is.na(dhcross$geno[[chrno]]$data[,j.old])))
                        throw.markers <- c(throw.markers,j)
                      else
                        {
                          throw.markers <- c(throw.markers,j.old)
                          j.old=j
                        }
                    }
                  else
                    j.old=j
                }
              else
                j.old=j
              m.old=m
            }
          if (length(throw.markers)>0)
            print(paste('will drop markers from LG',names(dhcross$geno[chrno]),':',paste(names(dhcross$geno[[chrno]]$map)[throw.markers],collapse=',')))
          dhcross <- drop.markers(dhcross,names(dhcross$geno[[chrno]]$map[throw.markers]))
        }
    }
  print(paste('Marker number in map reduced to::',totmar(dhcross)))
  return(dhcross)
}
###DHcrossL <- startQTLanalysis(datadir,gtfn,ptfn,genotypeCodes=genotypeCodes,outdir,pmethod=pmethod,epistasis=epistasis,map=map,redMap=redMap,BC.gen=BC.gen,F.gen=F.gen,traitDfr=traitDfr,alpha=alpha,alphaCIM=alphaCIM,LODthreshold=LODthreshold)
startQTLanalysis <- function(datadir,gtfn,ptfn,genotypeCodes,outdir,pmethod="pdf",epistasis=F,map=T,redMap=F,diagnostic=F,crosstype=crosstype,BC.gen=NA,F.gen=NA,traitDfr,alpha=0.4,alphaCIM=0.02,LODthreshold=NA)
### This is the central function of the QTL analysis suite.
{
    AZgrepl <- function(x){return(grepl('[A-Za-z]',x))} ### looks for alphabetic characters
    used.gt <- unique(unlist(read.csv(gtfn,skip=2,row.names=1,as.is=T))) ### get gt codes
    mygts <- toupper(genotypeCodes)[toupper(genotypeCodes)%in%used.gt] ### gt codes ordered
    if(sum(sapply(mygts,AZgrepl),na.rm=T)==0)
        mygts <- tolower(genotypeCodes)[tolower(genotypeCodes)%in%used.gt]### lower
    if(sum(sapply(mygts,AZgrepl),na.rm=T)==0)
      print('Check the genotype codes - No characters detected',quote=F)
    {if (is.na(BC.gen)|is.na(F.gen))
         MyCrossTotal <- read.cross(format="csvs",genfile=gtfn,phefile=ptfn,,genotypes=mygts,crosstype=crosstype)###treat as backcross
    else
        if ((!is.na(BC.gen))&(!is.na(F.gen))) ### SSD at generation?
            {
                print(paste('Generation: F',F.gen),sep='',quote=F)
                print(paste('Generation: BC',BC.gen),sep='',quote=F)
                MyCrossTotal <- read.cross(format="csvs",genfile=gtfn,phefile=ptfn,,genotypes=mygts,BC.gen=BC.gen,F.gen=F.gen)###treat as backcross
            }
        else
            {
                print('Cross type is not well defined. Check files')
                break()
            }
    }
    {if (redMap)
        MyCross <- reduceMap(MyCrossTotal) ### reduces size of map
    else
        MyCross <- jittermap(MyCrossTotal)
    }
    phenoMyCross <- names(MyCross$pheno) ### names of phenotypes
    MyCross$pheno <-subset(MyCross$pheno,  select=phenoMyCross[!phenoMyCross%in%'id'])### remove 'id'
    genoMyCross <- names(MyCross$geno)### genotype table
    genoMyCross <- renameLGR(genoMyCross) ### renames linkage groupse
    names(MyCross$geno) 
    names(MyCross$geno) <- genoMyCross
    if (map) ### map analysis checks for map consistency
        MyCross <- mapAnalysis(MyCross,gtfn,outdir,printpdf=T,diagnostic=TRUE,errorscore=FALSE)
    print("Start of QTL analysis",quote=FALSE) ### QTL analysis will be started
    summary(MyCross)
    MyCrossL <- qtlAnalysis(MyCross,gtfn,outdir,pmethod=pmethod,epistasis=epistasis,traitDfr=traitDfr,alpha=alpha,alphaCIM=alphaCIM,testmode=F,LODthreshold=LODthreshold)
    return(MyCrossL)
}
