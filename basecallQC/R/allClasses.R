#' The Parameters for BCL2FastQparamters object.
#'
#' Parameter class and accessors for use with basecallQC
#'
#' @aliases BCL2FastQparams BCL2FastQparams-BCL2FastQparams
#'
#' @references See \url{http://mrccsc.github.io} for more details on soGGi workflows
#' @rdname BCL2FastQparams
#' @docType class
#' @return A BCL2FastQparams object.
#' @export
#'
setClass("BCL2FastQparams", representation(RunDir="character",OutDir="character",RunParameters = "list"))

#' The basecallQC object.
#'
#' Objects and methods to handle Illumina BCL inputs and output files.
#' Provides sample sheet cleanup, basecall command summary QC statistics for basecalling.
#'
#' @aliases basecallQC basecallQC-basecallQC
#'
#' @references See \url{http://mrccsc.github.io} for more details on soGGi workflows
#' @rdname basecallQC
#' @docType class
#' @return A basecallQC object.
#' @export

setClass("basecallQC", representation(BCL2FastQparams="BCL2FastQparams",RunMetadata = "data.frame",
                                      cleanedSampleSheet="data.frame",BCLCommand="character",baseMasks= "data.frame",
                                      baseCallMetrics="list",demultiplexMetrics="list",fqQCmetrics="list"))


#' Set Parameters for BCL2FastQparamters object.
#'
#' Parameter class and accessors
#'
#' @aliases setBCL2FastQparams setBCL2FastQparams-setBCL2FastQparams
#'
#' @references See \url{http://mrccsc.github.io} for more details on soGGi workflows
#' @rdname setBCL2FastQparams
#' @docType methods
#' @param runXML file path to runParameters.xml ,if not specified
#' looks in run directory.
#' @param config file path to config.ini ,if not specified
#' looks in run directory.
#' @param runDir file path to run directory.
#' @param outDir file path to run directory.
#' @param verbose TRUE or FALSE. Messages on or off. Warnings/errors persist
#' @return A BCL2FastQparams object.
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#' @export

setBCL2FastQparams <- function(runXML=NULL,config=NULL,runDir=NULL,outDir=NULL,verbose=TRUE){
  if(is.null(runDir)) runDir <- getwd(); if(verbose) message("No runDir specified, run directory set to working directory");
  if(is.null(runXML)){
    if(verbose) message("No location for runParameters.xml specified")
    runParameters <- file.path(runDir,"runParameters.xml")
    if(!file.exists(runParameters)) stop("No runParameters.xml found in run directory")
  }
  if(is.null(config)){
    if(verbose) message("No location for config.ini specified")
    config <- file.path(runDir,"config.ini")
    if(!file.exists(config)) stop("No config.ini found in run directory")
  }

  runParameters = runParams(runXML,config)

  if(is.null(outDir)){
    if(verbose) message("No location for outDir specified")
    outDir <- file.path(runDir,runParameters$runParams$Barcode)
    message("outDir set to",outDir)
  }
  new("BCL2FastQparams",
      RunDir=runDir,
      OutDir=outDir,
      RunParameters=runParameters
  )

}
#' The bclCall function is a constructor for basecallQC objects.
#'
#' @name basecallQC
#' @rdname basecallQC
#' @param bcl2fastqparams A BCL2FastQparams object as created by setBCL2FastQparams.
#' @param RunMetaData Any run metadata to attach (data.frame)
#' @param sampleSheet The samplesheet for basecalling
#' @param doFQMetric TRUE or FALSE. Perform ShortRead FastQ quality assessment
#' using ShortRead's qa and report function
#' @return basecallQC a basecallQC object
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' @export
basecallQC <- function(bcl2fastqparams,RunMetaData=NULL,sampleSheet=NULL,doFQMetric=FALSE){

  cleanedSampleSheet <- validateBCLSheet(sampleSheet,bcl2fastqparams)
  baseMasks <- createBasemasks(cleanedSampleSheet,bcl2fastqparams)
  toSubmit <- createBCLcommand(bcl2fastqparams,cleanedSampleSheet,baseMasks)
  basecallmetrics <- baseCallMetrics(bcl2fastqparams)
  demultiplexmetrics <- demultiplexMetrics(bcl2fastqparams)
  if(doFQMetric==TRUE){
    fastqs <- dir("~/Documents/test2/fastq2/",pattern="*.fastq.gz",full.names=TRUE)[1:2]
    fqQCmetrics <- qcShortRead(fastqs)
  }else{
    fqQCmetrics <- list(FQQC_Table = NULL,ShortReadQC=NULL)
  }
  basecallQC <- new("basecallQC",
                    BCL2FastQparams = bcl2fastqparams,
                    cleanedSampleSheet = cleanedSampleSheet,
                    baseMasks = baseMasks,
                    BCLCommand=toSubmit,
                    baseCallMetrics = basecallmetrics,
                    demultiplexMetrics = demultiplexmetrics,
                    fqQCmetrics=fqQCmetrics)
  return(basecallQC)
}

runParams <- function(runXML=NULL,config=NULL){
  runParams <- runParameters(runXML)
  configParams <- configParams(config)
  return(list(runParams=runParams,configParams=configParams))
}

baseCallMetrics <- function(bcl2fastqparams){
  convStatsXML <- file.path(bcl2fastqparams@OutDir,"Stats","ConversionStats.xml")
  if(!file.exists(convStatsXML)) return(list(convStatsProcessed=NULL,summarisedConvStats=NULL))
  convStatsProcessed <- processConvStats(convStatsXML)
  summarisedConvStats <- summariseConvStats(convStatsProcessed)
  return(list(convStatsProcessed=convStatsProcessed,
                       summarisedConvStats=summarisedConvStats))
}


demultiplexMetrics <- function(bcl2fastqparams){
  demuxStatsXML <- file.path(bcl2fastqparams@OutDir,"Stats","DemultiplexingStats.xml")
  if(!file.exists(demuxStatsXML)) return(list(demuxStatsProcessed=NULL,summarisedDemuxStats=NULL))
  demuxStatsProcessed <- processDemultiplex(demuxStatsXML)
  summarisedDemuxStats <- summariseDemuxStats(demuxStatsProcessed)
  return(list(demuxStatsProcessed=demuxStatsProcessed,
                       summarisedDemuxStats=summarisedDemuxStats))
}






