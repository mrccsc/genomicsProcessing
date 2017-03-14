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
#' @param runFolder directory of illumina runfolder.
#' @param sampleSheet
#' @return A basecallQC object.
#' @examples
#'
#' warning("Put example here!")
#' @export
setClass("basecallQC", representation(Run = "character", RunMetadata = "data.frame"))

#' The bclCall function is a constructor for basecallQC objects.
#'
#' @name bclCall
#' @rdname bclCall
#' @param Run The pun to process
#' @param RunMetadata Any run metadata to attach (sata.frame)
#' @export
basecallQC <- function(Run,RunMetadata=NULL,params=NULL,baseCallMetrics=NULL,demultiplexMetrics=NULL){
  basecallQC <- new("basecallQC",
               Run = Run,
               RunMetadata = RunMetadata,
               parameters=runParams(params),
               baseCallMetrics=baseCallMetrics(params),
               demultiplexMetrics=demultiplexMetrics(params))
  return(basecallQC)
}

runParams <- function(params=NULL){
  runParams <- runParameters(params)
  configParams <- configParams(params)
  return(c(params,list(runParams=runParams,configParams=configParams)))
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
