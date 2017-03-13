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
setClass("basecallQC", representation(Run = "character", RunMetadata = "matrix"))

