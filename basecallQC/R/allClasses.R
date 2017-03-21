#' The Parameters for BCL2FastQparamters object.
#'
#' Parameter class and accessors
#'
#' @aliases BCL2FastQparams BCL2FastQparams-BCL2FastQparams
#'
#' @references See \url{http://mrccsc.github.io} for more details on soGGi workflows
#' @rdname BCL2FastQparams
#' @docType class
#' @return A BCL2FastQparams object.
#' @examples
#'
#' warning("Put example here!")
#' @export
#'
setClass("BCL2FastQparams", representation(RunDir="character",RunParameters = "list"))

#' The basecallQC object.
#'
#' Objects and methods to handle Illumina BCL inputs and output files.
#' Provides summary QC statistics for basecalling
#'
#' @aliases basecallQC basecallQC-basecallQC
#'
#' @references See \url{http://mrccsc.github.io} for more details on soGGi workflows
#' @rdname basecallQC
#' @docType class
#' @param BCL2FastQparams Parameters for BCL2FastQC program
#' and current analysis Run.
#' @param RunMetadata Dataframe containing Sample_ID column and user-defined metadata.
#' @param baseCallMetrics directory of illumina runfolder.
#' @param demultiplexMetrics lr.
#' @param sampleSheet fr
#' @return A basecallQC object.
#' @examples
#'
#' warning("Put example here!")
#' @export

setClass("basecallQC", representation(BCL2FastQparams="list",RunMetadata = "data.frame",
                                      sampleSheet="list",baseCallMetrics="list",demultiplexMetrics="list"))


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
#' @param verbose TRUE or FALSE. Messages on or off. Warnings/errors persist
#' @return A BCL2FastQparams object.
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#' @export

setBCL2FastQparams <- function(runXML=NULL,config=NULL,runDir=NULL,verbose=TRUE){
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
  new("BCL2FastQparams",
      RunDir=runDir,
      RunParameters = runParams(runXML,config))

}
#' The bclCall function is a constructor for basecallQC objects.
#'
#' @name basecallQC
#' @rdname basecallQC
#' @param Run The pun to process
#' @param RunMetadata Any run metadata to attach (sata.frame)
#' @export
basecallQC <- function(bcl2fastqparams,RunMetaData=NULL,sampleSheet=NULL,
                       baseCallMetrics=NULL,demultiplexMetrics=NULL){

  cleanedSampleSheet <- validateBCLSheet(sampleSheet,bcl2fastqparams)
  baseMasks <- createBasemasks(cleanedSampleSheet,bcl2fastqparams)
  toSubmit <- createBCLcommand(bcl2fastqparams,cleanedSampleSheet,baseMasks)
  basecallmetrics <- baseCallMetrics(bcl2fastqparams)
  demultiplexmetrics <- demultiplexMetrics(bcl2fastqparams)

  basecallQC <- new("basecallQC",
                    Run = Run,
                    runParameters = bcl2fastqparams,
                    cleanedSampleSheet = cleanedSampleSheet,
                    baseMasks = baseMasks,
                    baseCallMetrics = basecallmetrics,
                    demultiplexMetrics = demultiplexmetrics)
  return(basecallQC)
}

runParams <- function(runXML=NULL,config=NULL){
  runParams <- runParameters(runXML)
  configParams <- configParams(config)
  return(list(runParams=runParams,configParams=configParams))
}

basecallMetrics <- function(bcl2fastqparams){
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
  demuxStatsProcessed <- processDemuxStats(demuxStatsXML)
  summarisedDemuxStats <- summariseDemuxStats(demuxStatsProcessed)
  return(list(demuxStatsProcessed=demuxStatsProcessed,
                       summarisedDemuxStats=summarisedDemuxStats))
}






