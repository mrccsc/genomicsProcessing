#' Generate a table of per sample summary demultiplexing statistics
#'
#' Creates a table of per sample summary statistics from demultiplex results
#'
#'
#' @docType methods
#' @name summaryDemuxTable
#' @rdname summaryDemuxTable
#'
#' @author Thomas Carroll
#'
#' @param BCLQC A basecall QC object as returned from basecallQC function
#' @param output Whether the report contains frozen or sortable tables. Options are "static" and "html"
#' @return Table A table for reporting demultiplexing results in an HTML.
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
    return(kable(toTable))
  }
  if(output=="html"){
    return(DT:::datatable(toTable))
  }
}

##' Generate a table of per sample summary basecalling statistics
#'
#' Creates a table of per sample summary statistics from basecalling results
#'
#'
#' @docType methods
#' @name summaryConvStatsTable
#' @rdname summaryConvStatsTable
#'
#' @author Thomas Carroll
#'
#' @param BCLQC A basecall QC object as returned from basecallQC function
#' @param output Whether the report contains frozen or sortable tables. Options are "static" and "html"
#' @return Table A table for reporting demultiplexing results in an HTML.
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
summaryConvStatsTable <- function(BCLQC,output="static"){

  toTable <- BCLQC@baseCallMetrics$summarisedConvStats$Sample_Stats
  if(output=="static"){
    return(kable(toTable))
  }
  if(output=="html"){
    return(DT:::datatable(toTable))
  }
}


#' Generate a table of per sample summary fastq QC statistics from ShortRead
#'
#' Creates a table of per sample summary statistics from fastq QC statistics from ShortRead
#'
#'
#' @docType methods
#' @name makeFQTable
#' @rdname makeFQTable
#'
#' @author Thomas Carroll
#'
#' @param BCLQC A basecall QC object as returned from basecallQC function
#' @param output Whether the report contains frozen or sortable tables. Options are "static" and "html"
#' @return Table A table for reporting demultiplexing results in an HTML.
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
#' #makeFQTable(bclQC,output="static")
#' @export

makeFQTable <- function(BCLQC,output="static"){
  fqQCTable <- BCLQC@fqQCmetrics$FQQC_Table
    if(!is.null(fqQCTable)){
      if(output=="static"){
        table <- kable(fqQCTable,escape = FALSE)
      }
      if(output=="html"){
        table <- DT:::datatable(fqQCTable,escape=FALSE)
      }
      return(table)
    }else{
      return(NULL)
  }
}
