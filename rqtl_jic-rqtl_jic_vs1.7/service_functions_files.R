### R script
### Provides some functions to create filenames easier and to give
### some graphic layout settings for plotting (PDF, PS, JGP and PNG).
### 15/07/2012 Luzie Wingen, JIC

createDir <- function(dirname)
{
  finf <- file.info(dirname)
  if (is.na(finf$isdir)|(!finf$isdir))
    {
      if (!is.na(finf$size))      
        {
          print(paste('Directory',dirname,' is a file already.'))
          dirname <- paste(dirname,'2',sep="")
        }
      dir.create(paste(dirname))
    }
  return(dirname)
}

makeFilename <- function(dirname, namevector, fext,ts=TRUE)
###constructs filenames from dirname, file basename, specname and file extension
{
  da <- Sys.Date()
  if (dirname != "")
    dirname <- paste(dirname,"/",sep="")
  filename <- ""
  for (nv in namevector)
    if (nv!= "")
      filename <- paste(filename,nv,"_",sep="")
  if (ts)
	  filename <- paste(dirname,filename,da,fext,sep="")
  else
	  filename <- paste(dirname,filename,fext,sep="")
  return(filename)
}

openPostscript <- function(plotfname,n=2,hz=FALSE)
### opens Postscript file with two or several plots on one A4 sheet
### Different page layouts could be addressed in future
{
  postscript(file=plotfname,paper="a4",horizontal=hz)
  if (n==2)
    nf <- layout(matrix(c(1,2),2,1,byrow=TRUE))
  if (n==4)
    nf <- layout(matrix(c(1,2,3,4),2,2,byrow=TRUE))
  if (n==6)
    nf <- layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE))
  if (n==8)
    nf <- layout(matrix(c(1,2,3,4,5,6,7,8),4,2,byrow=TRUE))
  print(plotfname)
}

openPNG<- function(plotfname,plwd=1,pcex=1,pwidth=1000,pheight=1000)
### opens PNG file
{
  png(file=plotfname,width=pwidth, height=pheight)
  print(plotfname)
	par('lwd'=plwd)	
	par('cex'=pcex)	
}

openJPG<- function(plotfname,n=2,hz=FALSE,plwd=1,pcex=1,pwidth=1000,pheight=1000)
### opens Postscript file with two or several plots on one A4 sheet
### Different page layouts could be addressed in future
{
  jpeg(file=plotfname,width=pwidth, height=pheight, quality=94)
  print(plotfname)
	par('lwd'=plwd)	
	par('cex'=pcex)	
}
 
openPDF <- function(plotfname,n=2,p="a4",w=7,h=7,pointsize=12)
### opens Postscript file with two or several plots on one A4 sheet
### Different page layouts could be addressed in future
{
  pdf(file=plotfname,paper=p,width=w,height=h,pointsize=pointsize)
  if (n==2)
    nf <- layout(matrix(c(1,2),2,1,byrow=TRUE))
  if (n==4)
    nf <- layout(matrix(c(1,2,3,4),2,2,byrow=TRUE))
  if (n==6)
    nf <- layout(matrix(c(1,2,3,4,5,6),3,2,byrow=TRUE))
  if (n==8)
    nf <- layout(matrix(c(1,2,3,4,5,6,7,8),4,2,byrow=TRUE))
  print(plotfname)
}
