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


#' Illumina Basecalling functions.
#'
#' Parses the Illumina samplesheet from run folders to create standardised
#'
#'
#' @docType methods
#' @name validateBCLSheet
#' @rdname validateBCLSheet
#'
#' @author Thomas Carroll and Marian Dore
#' @param sampleSheet
#' @param param

#' @return Cleaned samplesheet.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#'

#' fileLocations <- system.file("extdata",package="basecallQC")
#'
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)

#' cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=NULL)
#'
#' @export
validateBCLSheet <- function(sampleSheet,param=NULL){
  ss <- fread(sampleSheet,sep=",",header=T,stringsAsFactors=F) %>%
    {if(exists('Project', where = .) & !exists('Sample_Project', where = .)) rename(.,Sample_Project = Project) else .} %>%
    {if(exists('SampleID', where = .) & !exists('Sample_ID', where = .)) rename(.,Sample_ID = SampleID) else .} %>%
    {if(exists('ID', where = .) & !exists('Sample_ID', where = .)) rename(.,Sample_ID = ID) else .} %>%
    {if(exists('SampleName', where = .) & !exists('Sample_Name', where = .)) rename(.,Sample_Name = SampleName) else .} %>%
    {if(exists('Name', where = .) & !exists('Name', where = .)) rename(.,Sample_Name = Name) else .} %>%
    {if(exists('Index', where = .) & !exists('index', where = .)) rename(.,index = Index) else .} %>%
    {if(exists('Index2', where = .) & !exists('index2', where = .)) rename(.,index2 = Index2) else .} %>%
    mutate(Sample_Project = if (exists('Sample_Project', where = .)) Sample_Project else NA,
           Lane = if (exists('Lane', where = .)) Lane else NA,
           Sample_ID = if (exists('Sample_ID', where = .)) Sample_ID else NA,
           Sample_Name = if (exists('Sample_Name', where = .)) SampleName else NA,
           index = if (exists('index', where = .)) index else NA,
           index2 = if (exists('index2', where = .)) index2 else NA
           ) %>%
    select(Sample_Project,Lane,Sample_ID,Sample_Name,index,index2,everything()) %>%
    mutate(Sample_ID=gsub("^X","Sample_",validNames(Sample_ID))) %>%
    mutate(Sample_ID=gsub("\\?|\\(|\\)|\\[|\\]|\\\\|/|\\=|\\+|<|>|\\:|\\;|\"|\'|\\*|\\^|\\||\\&|\\.","_",Sample_ID)) %>%
    mutate(index=str_trim(index,"both"),
           index2=str_trim(index2,"both")) #   %>%
    #mutate(Index=str_sub(index,1,index1Length(param)),    # Will use runParamsIndexLength
    #       index2=str_sub(index2,1,index2Length(param)))  # Will use runParamsIndexLength

  message("Read samplesheet ",basename(sampleSheet)," discovered for run ",currentRun)
  #index1Lengths <- str_length(ss$Index)
  ss
}




