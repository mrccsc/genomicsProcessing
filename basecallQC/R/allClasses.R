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
basecallQC <- function(Run,RunMetadata=NULL,params=NULL,baseCallMetrics=NULL,demultiplexMetrics=NULL){
  basecallQC <- new("basecallQC",
               Run = Run,
               RunMetadata = RunMetadata,
               parameters=runParams(params),
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

validateBCLSheet <- function(sampleSheet){
  ss <- fread(sampleSheet,sep=",",header=T,stringsAsFactors=F)
  ss$SampleID <- gsub("^X","Sample_",validNames(ss$SampleID))
  ss$SampleID <-gsub("\\?|\\(|\\)|\\[|\\]|\\\\|/|\\=|\\+|<|>|\\:|\\;|\"|\'|\\*|\\^|\\||\\&|\\.","_",ss$SampleID)
  ss$Index <- gsub(" ","",ss$Index)
  toBeTrimmed <- which(lapply(ss$Index,nchar) > as.numeric(currentRunParameters$IndexRead1))
  for(i in toBeTrimmed){
    ss$Index[toBeTrimmed] <-  substr(ss$Index[toBeTrimmed],0,as.numeric(currentRunParameters$IndexRead1))
  }

  message("Read samplesheet ",basename(sampleSheetName)," discovered for run ",currentRun)
  index1Lengths <- unlist(lapply(ss$Index,function(x)nchar(x)))

  if(any(colnames(ss) %in% "Index2")){
    ss$Index2 <- gsub(" ","",ss$Index2)
    index2Lengths <- unlist(lapply(ss$Index2,function(x)nchar(x)))
    index2NAs <- unlist(lapply(ss$Index2,function(x)is.na(x)))
    index2Lengths[index2NAs] <- 0
    allIndexTypes <- paste0(index1Lengths,"-",index2Lengths)
    uniqueIndexTypes <- unique(allIndexTypes)
    #ss$SampleID <- gsub("[[:punct:]]", "_", ss$SampleID).)
    #ss$SampleID <- gsub("[^[:alnum:]]", "_", ss$SampleID).)
    ss$Index <- gsub("-NA|-$|^-$","",paste(ss$Index,ss$Index2,sep="-"))
    ss <- ss[,-grep("Index2",colnames(ss))]
  }else{
    allIndexTypes <- paste0(index1Lengths)
    uniqueIndexTypes <- unique(allIndexTypes)
    #ss$SampleID <- gsub("[[:punct:]]", "_", ss$SampleID).)
    #ss$SampleID <- gsub("[^[:alnum:]]", "_", ss$SampleID).)
  }
}
