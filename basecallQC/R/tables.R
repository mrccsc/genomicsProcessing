#' Generate per sample summary statistics
#'
#' Creates per sample summary statistics from demultiplex results
#'
#'
#' @docType methods
#' @name summariseConvStats
#' @rdname summariseConvStats
#'
#' @author Thomas Carroll
#'
#' @param ConvStats Results from a call to processConvStats.
#' @return A datatable of summarised per sample results.
#' @import stringr XML RColorBrewer methods raster
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' summaryDemuxTable(bclQC,output="static")
#' @export
summaryDemuxTable <- function(BCLQC,output="static"){

  toTable <- BCLQC@demultiplexMetrics$summarisedDemuxStats$Summary
  if(output=="static"){
    kable(toTable)
  }
  if(output=="html"){
    DT:::datatable(toTable)
  }
}
