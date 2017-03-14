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

#' The bclCall function is a constructor for basecallQC objects.
#'
#' @name bclCall
#' @rdname bclCall
#' @param Run The pun to process
#' @param RunMetadata Any run metadata to attach (sata.frame)
#' @export
basecallQC <- function(Run,RunMetadata=NULL,params=NULL,sampleSheet=NULL,baseCallMetrics=NULL,demultiplexMetrics=NULL){
  basecallQC <- new("basecallQC",
               Run = Run,
               RunMetadata = RunMetadata,
               parameters=runParams(params),
               sampleSheet=validateBCLSheet(params),
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


## I have some dummy parameters which exists in central robect..ill rethink
validateBCLSheet <- function(sampleSheet,param){
  ss <- fread(sampleSheet,sep=",",header=T,stringsAsFactors=F) %>%
    tbl_df %>%
    mutate(Sample_Project = if (exists('Sample_Name', where = .)) SampleName else NA,
           Lane = if (exists('Labe', where = .)) Lane else NA,
           Sample_ID = if (exists('Sample_ID', where = .)) Sample_ID else NA,
           Sample_Name = if (exists('Sample_Name', where = .)) SampleName else NA,
           index = if (exists('index', where = .)) index else NA,
           index2 = if (exists('index2', where = .)) index2 else NA
           ) %>%
    dplyr::select(Sample_Project,Lane,Sample_ID,Sample_Name,index,index2) %>%
    mutate(Sample_ID=gsub("^X","Sample_",validNames(Sample_ID))) %>%
    mutate(Sample_ID=gsub("\\?|\\(|\\)|\\[|\\]|\\\\|/|\\=|\\+|<|>|\\:|\\;|\"|\'|\\*|\\^|\\||\\&|\\.","_",Sample_ID)) %>%
    mutate(index=str_trim(index,"both"),
           index2=str_trim(index2,"both")) %>%
    mutate(Index=str_sub(index,1,index1Length(param)),
           index2=str_sub(index2,1,index1Length(param)))  # Will use runParamsIndexLength

  message("Read samplesheet ",basename(sampleSheetName)," discovered for run ",currentRun)
  index1Lengths <- str_length(ss$Index)
}




