#' Generate per sample summary statistics
#'
#' Creates per sample summary statistics from demultiplex results
#'
#'
#' @docType methods
#' @name summaryDemuxTable
#' @rdname summaryDemuxTable
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
#' @export
reportBCL <- function(BCLQC,reportOut=file.path("/Users/tcarroll/genomicsProcessing/","report.html"),output="static",reportRMDfile=NULL){
  BCLQCreport <- BCLQC
  fileLocations <- system.file("extdata",package="basecallQC")
  if(is.null(reportRMDfile)){
    reportRMD <- file.path(fileLocations,"reportRMDs","basecallqcreport.Rmd")
  }
  if(!file.exists(reportRMD)) stop()
  #dir.create(BCLQCreport@BCL2FastQparams@OutDir,showWarnings = F)
  render(reportRMD,
         output_file = reportOut)
}
