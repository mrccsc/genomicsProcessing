#' Generate basecallQC report
#'
#' Creates a summary report for basecalling and demultiplexing statistics.
#'
#'
#' @docType methods
#' @name reportBCL
#' @rdname reportBCL
#'
#' @author Thomas Carroll
#'
#' @param BCLQC A basecall QC object as returned from basecallQC function
#' @param reportOut Name of report file
#' @param reportOutDir Directory for the report file
#' @param output Whether the report contains frozen or sortable tables. Options are "static" and "html"
#' @param reportRMDfile RMD to be used for reporting. (Default uses standard report template)
#' @param FQQC TRUE or FALSE, whether to run ShortRead fastq QC.
#' @return RenderedReport A rendered report.
#' @import stringr XML RColorBrewer methods raster ShortRead prettydoc
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' reportBCL(bclQC,"TestReport.html",output="html")
#' @export
reportBCL <- function(BCLQC,reportOut="report.html",reportOutDir=getwd(),
                      output="static",reportRMDfile=NULL,FQQC=FALSE){
  BCLQCreport <- BCLQC
  doFQQC <- FQQC
  fileLocations <- system.file("extdata",package="basecallQC")
  if(is.null(reportRMDfile)){
    reportRMD <- file.path(fileLocations,"reportRMDs","basecallqcreport.Rmd")
  }
  if(!file.exists(reportRMD)) stop()
  #dir.create(BCLQCreport@BCL2FastQparams@OutDir,showWarnings = F)
  render(reportRMD,
         output_file = reportOut,output_dir=reportOutDir,params=list(title=paste0("basecallQC")))
}
