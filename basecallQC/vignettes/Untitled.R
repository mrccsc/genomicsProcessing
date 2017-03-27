## ---- cleaning, fig.show='hold',warning=FALSE,message=FALSE--------------
library(basecallQC)
fileLocations <- system.file("extdata",package="basecallQC")
runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=bcl2fastqparams)
head(cleanedSampleSheet)

## ---- basemasks, fig.show='hold',warning=FALSE,message=FALSE-------------
fileLocations <- system.file("extdata",package="basecallQC")
runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=bcl2fastqparams)
baseMasks <- createBasemasks(cleanedSampleSheet,param=bcl2fastqparams)
baseMasks

## ---- submitCommand, fig.show='hold',warning=FALSE,message=FALSE---------
toSubmit <- createBCLcommand(bcl2fastqparams,cleanedSampleSheet,baseMasks)
toSubmit

