#' Read lengths
#' Readlengths as defined runParameters.xml
#'
#' @usage
#' \S4method{readlengths}{BCL2FastQparams}(object)
#'
#' @docType methods
#' @name readlengths
#' @rdname readlengths
#' @aliases readlengths readlengths,BCL2FastQparams-method, readlengths.bcl2fastqparams
#'
#' @author Thomas Carroll
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#' readlength <- readlengths(bcl2fastqparams)
#' @param object A BCL2FastQparams object
#' @return readlengths <- Readlengths as defined runParamaeters.xml
#' @export
readlengths.bcl2fastqparams <-  function (object)
{
  dplyr:::select(object@RunParameters$runParams,Read1,Read2) %>% mutate_all(as.numeric)
}

setGeneric("readlengths", function(object="BCL2FastQparams") standardGeneric("readlengths"))

#' @rdname readlengths
#' @export
setMethod("readlengths", signature(object="BCL2FastQparams"), readlengths.bcl2fastqparams)

#' Index lengths
#' Indexlengths as defined runParameters.xml
#'
#' @usage
#' \S4method{indexlengths}{BCL2FastQparams}(object)
#'
#' @docType methods
#' @name indexlengths
#' @rdname indexlengths
#' @aliases indexlengths indexlengths,BCL2FastQparams-method, indexlengths.bcl2fastqparams
#'
#' @author Thomas Carroll
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#' indexlength <- indexlengths(bcl2fastqparams)
#' @param object A BCL2FastQparams object
#' @return readlengths <- Readlengths as defined runParamaeters.xml
#' @export
indexlengths.bcl2fastqparams <-  function (object)
{
  dplyr:::select(object@RunParameters$runParams,IndexRead1,IndexRead2)  %>% mutate_all(as.numeric)
}

setGeneric("indexlengths", function(object="BCL2FastQparams") standardGeneric("indexlengths"))

#' @rdname indexlengths
#' @export
setMethod("indexlengths", signature(object="BCL2FastQparams"), indexlengths.bcl2fastqparams)



