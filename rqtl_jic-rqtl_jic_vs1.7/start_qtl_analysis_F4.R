####################################################
### R script
### QTL analysis with the R/qtl package
### -data upload (expects correct format)
### -upload of trait short names
### -QTL analysis
### -producing qtl result table
### -printing results into PDF or JPG files
### 25/12/2015 Luzie U. Wingen, JIC
####################################################

### This folder, Rqtl_jic, containing script 'start_qtl_analysis.R'
### should be the location where the script is executed.
### Copy datafiles into the folder below this: 'rqtl_data'

#############################################
### Sourcing R functions from other files ###
#############################################
rm(list=ls()) ### clear the dataspace
qtlFunctionDir="."
sourcefiles=c("analyse_qtl.R") ### one source file
source(paste(qtlFunctionDir,sourcefiles,sep="/"))

##########################
### Directory settings ###
##########################
###basdir <- '.' ### General directory where analysis starts
basdir <- '..' ### Directory in Luzie's system
datadir <- paste(basdir,'rqtl_data',sep='/') ### Directory with data
print('datadir')
print(datadir)
print(dir(datadir,pattern='_rqtl_'))
outdir <-paste(basdir,"rqtl_out",sep="/") ### Directory for output
mywd <- getwd() ### local working directory
print(paste('Working directory is:',mywd))
if (is.na(file.info(outdir)$isdir)) ### Is outdir present?
    dir.create(outdir)### No - create it.
plotdir <- paste(outdir,"plots",sep="/")
if (is.na(file.info(plotdir)$isdir)) ### Is plot directory present?
    dir.create(plotdir)### No - create it.
csvdir <- paste(outdir,"csv",sep="/")
if (is.na(file.info(csvdir)$isdir)) ### Is csv directory present present?
    dir.create(csvdir)### No - create it.
print(paste('Output will be written to outdir:',outdir))
#### Is a file with controlled vocabulary present?
{if (!is.na(file.info(paste(datadir,'field_trial_phenotypes_or_traits.csv',sep='/'))$size))
  traitDfr <- read.table(paste(datadir,'field_trial_phenotypes_or_traits.csv',sep='/'),as.is=T,sep=';',header=T,nrow=100)
else ### create an empty data.frame for trait names
     {
         traitDfr=data.frame(matrix(NA,ncol=4,nrow=0))
         colnames(traitDfr) <- c( "trait_name","short_name","unit_abbreviation","description")
     }
 print('Will used controlled trait vocabulary if file: field_trial_phenotypes_or_traits.csv is present ',quote=F)
 if(nrow(traitDfr) > 0)
   print(traitDfr[7:13,1:2],quote=F)
}


###################################
### Analysis Parameter settings ###
###################################
fnames <- dir(path=datadir,pattern=paste("^.*","_rqtl_.*",".csv",sep=""),full.names=TRUE)
ptindex <- grep('phenotype',fnames)
gtindex <- grep('genotype',fnames)
gpindex <- grep('_pg',fnames) #### For files containing genotype and phenotype together
pmethod="pdf"  ###  pmethod="jpeg"
epistasis=F ### epistasis between QTL. Very computational expensive. FALSE or F default.
crosstype='riself' ### can be 'dh','f2','bc','riself'. Will be overwritten by BC.gen and F.gen
BC.gen=0 ### Complicated cross types. Set generation number or BC.gen=NA if you want crosstype 'dh', 'bc', 'f2' or 'riself'
F.gen=5 ### Complicated cross types. Set generation number or BC.gen=NA if you want crosstype 'dh', 'bc', 'f2' or 'riself'
map=F ### analyse the quality of the map - not needed but interesting
redMap=T ### reduce the marker number in the map - great for large maps, e.g. from iSelect or axiom data
genotypeCodes=c('A','H','B','-')### genotype codes in order AA, AB, BB, missing
alpha=0.4 ### analysis significant level 0.4 = 40%
alphaCIM=0.05 ### analysis significant level 0.05 = 5%
LODthreshold=NA ### To test a lower threshold then calculated from the data set to a value, e.g. 2.0. Set to NA if you go with the programme threshold.


###################################
### Start Analysis  ###
###################################
### Data files needed with the following names
### '_rqtl_genotypes.csv' and
### '_rqtl_phenotypes.csv' 
### as the ONLY DIFFERENCE in the file names.

### This could be made interactive in future.

FindFilesAndStart <- function(datadir,map=map,redMap=redMap,pmethod=pmethod,epistasis=epistasis,crosstype=crosstype,BC.gen=BC.gen,F.gen=Fgen,genotypeCodes=genotypeCodes,traitDfr=traitDfr,LODthreshold=LODthreshold)
### Wrapper to start the analysis
{
    print(paste('Program checks for data files in the \'',datadir,'\' subdirectory',sep=""),quote=F)
    fnames <- dir(path=datadir,pattern=paste("^.*","_rqtl_.*",".csv",sep=""),full.names=TRUE)
    ptindex <- grep('_phenotypes',fnames)
    gtindex <- grep('_genotypes',fnames)
    gtptindex <- grep('_rqtl_pg',fnames)
    allqtlL <- list() ### List - container with results
    if ((length(ptindex)!=0) & (length(gtindex)!=0))
        {
            for (pt in ptindex)### Loops through all phenotype files
                for (gt in gtindex)### Loops through all genotype files
                    if (sub('_rqtl_.*','',fnames[pt])==sub('_rqtl_.*','',fnames[gt])) ### only works on data file with similar names
                        {
                            gtfn <- fnames[gt] ### filename for genotype
                            ptfn <- fnames[pt] ### filename for phenotype
                            print('RQTL data files found:',quote=F)
                            print(paste(gtfn,ptfn),quote=F)
                            DHcrossL <- startQTLanalysis(datadir,gtfn,ptfn,genotypeCodes=genotypeCodes,outdir,pmethod=pmethod,epistasis=epistasis,map=map,redMap=redMap,BC.gen=BC.gen,F.gen=F.gen,traitDfr=traitDfr,alpha=alpha,alphaCIM=alphaCIM,LODthreshold=LODthreshold)
                            allqtlL <-c(allqtlL,DHcrossL) ### add to result to a final list
                        }
            if (exists("gtfn"))
                {
                    rdatafile <- paste(sub('_rqtl.*$','',gtfn),'.RData',sep='')
                    save.image(file=rdatafile) ### save the R binary
                }
        }
    else
        print('Datafiles not found. Please put files into the folder: rqtl_datadir')
    return(allqtlL)
}

### Executing the FindFilesAndStartFunction
allqtlL <- FindFilesAndStart(datadir=datadir,map=map,redMap=redMap,pmethod=pmethod,epistasis=epistasis,crosstype=crosstype,BC.gen=BC.gen,F.gen=F.gen,genotypeCodes=genotypeCodes,traitDfr=traitDfr,LODthreshold=LODthreshold)
