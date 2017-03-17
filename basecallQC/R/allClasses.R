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
#' @param Run directory of illumina runfolder.
#' @param RunMatadata Dataframe containing Sample_ID column and user-defined metadata.
#' @param baseCallMetrics directory of illumina runfolder.
#' @param params directory of illumina runfolder.
#' @param demultiplexMetrics lr.
#' @param sampleSheet fr
#' @return A basecallQC object.
#' @examples
#'
#' warning("Put example here!")
#' @export
setClass("basecallQC", representation(Run = "character", RunMetadata = "data.frame",
                                      params="list",sampleSheet="list",baseCallMetrics="list",demultiplexMetrics="list"))

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
setClass("BCL2FastQparams", representation(RunParameters = "character"))

#' Set Parameters for BCL2FastQparamters object.
#'
#' Parameter class and accessors
#'
#' @aliases setBCL2FastQparams setBCL2FastQparams-setBCL2FastQparams
#'
#' @references See \url{http://mrccsc.github.io} for more details on soGGi workflows
#' @rdname setBCL2FastQparams
#' @docType methods
#' @param runXML file path to runParameters.xml
#' @return A BCL2FastQparams object.
#' @examples
#'
#' warning("Put example here!")
#' @export

setBCL2FastQparams <- function(runXML=NULL){
  new("BCL2FastQparams",
      runParameters = runParams(runXML))
}
#' The bclCall function is a constructor for basecallQC objects.
#'
#' @name bclCall
#' @rdname bclCall
#' @param Run The pun to process
#' @param RunMetadata Any run metadata to attach (sata.frame)
#' @export
basecallQC <- function(Run=NULL,RunMetadata=NULL,runXML=NULL,config=NULL,sampleSheet=NULL,
                       baseCallMetrics=NULL,demultiplexMetrics=NULL){
  basecallQC <- new("basecallQC",
                    Run = Run,
                    RunMetadata = RunMetadata,
                    runParameters = setBCL2FastQparams(runXML,config),
                    cleanedSampleSheet = validateBCLSheet(sampleSheet,runParameters),
                    baseMasks = createBasemasks(cleanedSampleSheet,runParameters),
                    baseCallMetrics = baseCallMetrics(runParameters),
                    demultiplexMetrics = demultiplexMetrics(runParameters))
  return(basecallQC)
}

runParams <- function(runXML=NULL,config=NULL){
  runParams <- runParameters(runXML)
  configParams <- configParams(config)
  return(list(runParams=runParams,configParams=configParams))
}

basecallMetrics <- function(params){
  convStatsProcessed <- processConvStats(params)
  summarisedConvStats <- summariseConvStats(params)
  return(c(params,list(convStatsProcessed=convStatsProcessed,
                       summarisedConvStats=summarisedConvStats)))
}

demultiplexMetrics <- function(params){
  demuxStatsProcessed <- processDemuxStats(params)
  summarisedDemuxStats <- summariseDemuxStats(demuxStatsProcessed)
  return(c(params,list(demuxStatsProcessed=demuxStatsProcessed,
                       summarisedDemuxStats=summarisedDemuxStats)))
}






