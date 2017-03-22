library(basecallQC)
context("Test 3rd Run example")

fileLocations <- system.file("extdata",package="basecallQC")
runXML <- dir(file.path(fileLocations,"Runs/170222_D00467_0228_ACAL6PANXX/"),pattern="runParameters.xml",full.names=TRUE)
config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
sampleSheet <- dir(file.path(fileLocations,"Runs/170222_D00467_0228_ACAL6PANXX/"),pattern="*\\.csv",full.names=TRUE)
bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)

#'


expect_that(length(unique(bclQC@cleanedSampleSheet$Sample_Project)) ==
              length(levels(read.delim(sampleSheet,header=T,sep=",")$Project)),
            is_true()
)
