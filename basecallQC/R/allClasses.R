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
bclCall <- function(Run,RunMetadata,params=NULL,baseCallMetrics,demultiplexMetrics){
  bclQC <- new("basecallQC",
               Run = Run,
               RunMetadata = RunMetadata,
               params=checkParams(params))
  return(bclQC)
}

checkParams <- function(params){
  runParams <- runParameters(params)
  configParams <- configParams(params)
}
