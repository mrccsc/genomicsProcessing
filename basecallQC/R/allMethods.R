#' readlengths
#' Readlengths as defined runParamaeters.xml
#'
#' @usage
#' \S4method{readlengths}{BCL2FastQparams}(object)
#'
#' @docType methods
#' @name readlengths
#' @rdname readlengths
#' @aliases readlengths readlengths,BCL2FastQparams-method
#'
#' @author Thomas Carroll
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=T)
#' readlength <- readlengths(bcl2fastqparams)
#' @export
#' @param object A BCL2FastQparams object
#' @return readlengths <- Readlengths as defined runParamaeters.xml
readlengths.bcl2fastqparams <-  function (object)
{
  dplyr:::select(bcl2fastqparams@RunParameters$runParams,Read1,Read2) %>% mutate_all(as.numeric)
}

setGeneric("readlengths", function(object="BCL2FastQparams") standardGeneric("readlengths"))

#' @rdname readlengths
#' @export
setMethod("readlengths", signature(object="BCL2FastQparams"), readlengths.bcl2fastqparams)

indexlengths.bcl2fastqparams <-  function (object)
{
  dplyr:::select(bcl2fastqparams@RunParameters$runParams,IndexRead1,IndexRead2)  %>% mutate_all(as.numeric)
}

setGeneric("indexlengths", function(object="BCL2FastQparams") standardGeneric("indexlengths"))

#' @rdname indexlengths
#' @export
setMethod("indexlengths", signature(object="BCL2FastQparams"), indexlengths.bcl2fastqparams)



