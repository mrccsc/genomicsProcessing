library(basecallQC)
context("Test 5th Run example")

fileLocations <- system.file("extdata",package="basecallQC")
runXML <- dir(file.path(fileLocations,"Runs/170303_D00467_0231_ACAK02ANXX/"),pattern="runParameters.xml",full.names=TRUE)
config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
sampleSheet <- dir(file.path(fileLocations,"Runs/170303_D00467_0231_ACAK02ANXX/"),pattern="*\\.csv",full.names=TRUE)
bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)

#'


expect_that(length(unique(bclQC@cleanedSampleSheet$Sample_Project)) ==
              length(levels(read.delim(sampleSheet,header=T,sep=",")$Project)),
            is_true()
)
