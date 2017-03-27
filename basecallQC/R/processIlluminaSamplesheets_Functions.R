
#' Illumina samplesheet cleaning
#'
#' Parses the Illumina samplesheet from run folders to create standardised samplesheet for BCL2Fastq > Version 2.18
#'
#'
#' @docType methods
#' @name validateBCLSheet
#' @rdname validateBCLSheet
#'
#' @author Thomas Carroll and Marian Dore
#' @param sampleSheet Data.frame of samplesheet for Illumina basecalling (see vignette for more details)
#' @param param A BCL2FastQparams object
#' @return cleanedSampleSheet A cleaned samplesheet.
#' @import stringr XML RColorBrewer methods raster BiocStyle lazyeval
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#' cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=bcl2fastqparams)
#'
#' @export
validateBCLSheet <- function(sampleSheet,param=bcl2fastqparams){
  #runParam <- runParameters(param)
  fread(sampleSheet,sep=",",header=T,stringsAsFactors=F) %>%
    tbl_df %>%
  {if(exists('Project', where = .) & !exists('Sample_Project', where = .)) dplyr:::rename(.,Sample_Project = Project) else .} %>%
  {if(exists('SampleID', where = .) & !exists('Sample_ID', where = .)) dplyr:::rename(.,Sample_ID = SampleID) else .} %>%
  {if(exists('ID', where = .) & !exists('Sample_ID', where = .)) dplyr:::rename(.,Sample_ID = ID) else .} %>%
  {if(exists('SampleName', where = .) & !exists('Sample_Name', where = .)) dplyr:::rename(.,Sample_Name = SampleName) else .} %>%
  {if(exists('Name', where = .) & !exists('Sample_Name', where = .)) dplyr:::rename(.,Sample_Name = Name) else .} %>%
  {if(exists('index', where = .) & !exists('Index', where = .)) dplyr:::rename(.,Index = index) else .} %>%
  {if(exists('index2', where = .) & !exists('Index2', where = .)) dplyr:::rename(.,Index2 = index2) else .} %>%
    mutate(Sample_Project = if (exists('Sample_Project', where = .)) Sample_Project else NA,
           Lane = if (exists('Lane', where = .)) Lane else NA,
           Sample_ID = if (exists('Sample_ID', where = .)) Sample_ID else NA,
           Sample_Name = if (exists('Sample_Name', where = .)) Sample_Name else NA,
           Index = if (exists('Index', where = .)) Index else NA,
           Index2 = if (exists('Index2', where = .)) Index2 else NA) %>%
    tbl_df %>%
    dplyr:::select(Sample_Project,Lane,Sample_ID,Sample_Name,Index,Index2,everything()) %>%
    mutate(Sample_ID=gsub("^X\\d+.\\.","Sample_",validNames(Sample_ID))) %>%
    mutate(Sample_ID=gsub("\\?|\\(|\\)|\\[|\\]|\\\\|/|\\=|\\+|<|>|\\:|\\;|\"|\'|\\*|\\^|\\||\\&|\\.","_",Sample_ID)) %>%
    mutate(Index=str_trim(Index,"both"),
           Index2=str_trim(Index2,"both"))    %>%
  mutate(Index=str_sub(Index,1,as.numeric(indexlengths(param)$IndexRead1)),    # Will use runParamsIndexLength
         Index2=str_sub(Index2,1,as.numeric(indexlengths(param)$IndexRead2)))  # Will use runParamsIndexLength
}

#' Functions to create basemasks for basecalling from Illumina samplesheet.
#'
#' Parses the Illumina samplesheet for versions > 2.18 and creates basemasks.
#'
#'
#' @docType methods
#' @name createBasemasks
#' @rdname createBasemasks
#'
#' @author Thomas Carroll and Marian Dore
#' @param cleanedSampleSheet Data.frame of cleaned samplesheet for Illumina basecalling (see vignette for more details)
#' @param param A BCL2FastQparams object
#' @return basemasks A data.frame containing basecall masks for reads and indexes.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#'

#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#'
#' cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=bcl2fastqparams)
#' basemasks <- createBasemasks(cleanedSampleSheet,param=bcl2fastqparams)
#'
#' @export
createBasemasks <- function(cleanedSampleSheet,param){
  indexCombinations <- cleanedSampleSheet %>%
    mutate(Index2=ifelse(is.na(Index2), "", Index2),Index=ifelse(is.na(Index), "", Index)) %>%
    mutate(indexLength=str_length(Index),indexLength2=str_length(Index2)) %>%
    dplyr:::group_by(Lane) %>% dplyr:::count(indexLength,indexLength2)

  if(nrow(indexCombinations) == length(unique(indexCombinations$Lane))){
    baseMasks <- indexCombinations %>%
      mutate(index1Mask = str_c(str_dup("I",indexLength),
                                str_dup("N",indexlengths(param)$IndexRead1-indexLength)),
             index2Mask = str_c(str_dup("I",indexLength2),
                                str_dup("N",indexlengths(param)$IndexRead2-indexLength2))) %>%
      mutate(read1Mask = str_c(str_dup("Y",as.numeric(readlengths(param)$Read1))),
             read2Mask = str_c(str_dup("Y",as.numeric(readlengths(param)$Read1)))) %>%
      mutate(read1Mask = str_replace(read1Mask,"Y$","N"),
             read2Mask = str_replace(read2Mask,"Y$","N")) %>%
      mutate(basemask = str_c(read1Mask,index1Mask,index2Mask,read2Mask,sep=",")) %>%
      mutate(basemask = str_c(Lane,":",basemask)) %>%
      mutate(basemask = str_replace(basemask,",,",",")) %>%
      tbl_df %>%
      dplyr:::select(Lane,basemask,read1Mask,index1Mask,index2Mask,read2Mask)
      }
}

#' Functions to create command for Illumina demultiplexing using fastq2BCL versions > 2.18 .
#'
#' Creates the command to be used for demultiplexing with fastq2BCL versions > 2.18
#'
#'
#' @docType methods
#' @name createBCLcommand
#' @rdname createBCLcommand
#'
#' @author Thomas Carroll and Marian Dore
#' @param bcl2fastqparams A BCL2FastQparams object
#' @param cleanedSampleSheet Data.frame of cleaned samplesheet for Illumina basecalling (see vignette for more details)
#' @param baseMasks A data.frame of basemasks as created by createBasemasks function
#' @return bclCommand A character vector containing the command for Illumina basecalling
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#'

#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
#'
#' cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=bcl2fastqparams)
#' baseMasks <- createBasemasks(cleanedSampleSheet,param=bcl2fastqparams)
#' toSubmit <- createBCLcommand(bcl2fastqparams,cleanedSampleSheet,baseMasks)
#' @export
createBCLcommand <- function(bcl2fastqparams,cleanedSampleSheet,baseMasks){
  sampleSheetLocation <- paste0(file.path(bcl2fastqparams@RunDir,bcl2fastqparams@RunParameters$runParams$Barcode),".csv")
  bclPath <- bcl2fastqparams@RunParameters$configParams[bcl2fastqparams@RunParameters$configParams$name == "configureBclToFastq","value"]
  write.table("[DATA]",file=sampleSheetLocation,sep="",quote=F,row.names=F)
  write.table(cleanedSampleSheet,file=sampleSheetLocation,sep=",",quote=F,row.names=F,append = T)
  baseMasksToUse <- str_c("--use-bases-mask ",dplyr:::select(tbl_df(baseMasks),basemask)$basemask,collapse = " ")
  bclcommand <- str_c(as.vector(bclPath$value),"--output-dir ",bcl2fastqparams@OutDir,"--sample-sheet",sampleSheetLocation,baseMasksToUse,sep=" ")
  return(bclcommand)
}
