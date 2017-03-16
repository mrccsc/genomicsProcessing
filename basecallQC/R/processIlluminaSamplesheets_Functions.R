
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
#' @param sampleSheet mm
#' @param param mn
#' @return Cleaned samplesheet.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runParameters <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=runParameters)
#'
#' @export
validateBCLSheet <- function(sampleSheet,param=NULL){
  runParam <- runParameters(param)
  fread(sampleSheet,sep=",",header=T,stringsAsFactors=F) %>%
    tbl_df %>%
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
           index2 = if (exists('index2', where = .)) index2 else NA) %>%
    tbl_df %>%
    dplyr:::select(Sample_Project,Lane,Sample_ID,Sample_Name,index,index2,everything()) %>%
    mutate(Sample_ID=gsub("^X\\d+.\\.","Sample_",validNames(Sample_ID))) %>%
    mutate(Sample_ID=gsub("\\?|\\(|\\)|\\[|\\]|\\\\|/|\\=|\\+|<|>|\\:|\\;|\"|\'|\\*|\\^|\\||\\&|\\.","_",Sample_ID)) %>%
    mutate(index=str_trim(index,"both"),
           index2=str_trim(index2,"both"))    %>%
  mutate(Index=str_sub(index,1,as.numeric(runParam$IndexRead1)),    # Will use runParamsIndexLength
         index2=str_sub(index2,1,as.numeric(runParam$IndexRead2)))  # Will use runParamsIndexLength
}

#' Functions to create basemasks for basecalling from Illumina samplesheet.
#'
#' Parses the Illumina samplesheet from run folders to create standardised
#'
#'
#' @docType methods
#' @name createBasemasks
#' @rdname createBasemasks
#'
#' @author Thomas Carroll and Marian Dore
#' @param sampleSheet mm
#' @param param mn
#' @return basemasks Basemasks.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#'

#' fileLocations <- system.file("extdata",package="basecallQC")
#'
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' runParameters <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' cleanedSampleSheet <- validateBCLSheet(sampleSheet,param=runParameters)
#' cleanedSampleSheet <- createBasemasks(cleanedSampleSheet,param=runParameters)
#'
#' @export
createBasemasks <- function(cleanedSampleSheet,param=NULL){
  runParam <- runParameters(param)
  indexCombinations <- cleanedSampleSheet %>%
    mutate(indexLength=str_length(index),indexLength2=str_length(index2)) %>%
    group_by(Lane) %>% count(indexLength,indexLength2)

  if(nrow(indexCombinations) == length(unique(indexCombinations$Lane))){
    baseMasks <- indexCombinations %>%
      mutate(index1Mask = str_c(str_dup("Y",indexLength),
                                str_dup("N",as.numeric(runParam$IndexRead1)-indexLength)),
             index2Mask = str_c(str_dup("Y",indexLength2),
                                str_dup("N",as.numeric(runParam$IndexRead2)-indexLength2))) %>%
      mutate(read1Mask = str_c(str_dup("Y",as.numeric(runParam$Read1))),
             read2Mask = str_c(str_dup("Y",as.numeric(runParam$Read2)))) %>%
      mutate(read1Mask = str_replace(read1Mask,"Y$","N"),
             read2Mask = str_replace(read2Mask,"Y$","N"))
      }
}


#' Illumina Basecalling functions.
#'
#' Parses the Illumina samplesheet from run folders to create standardised
#'
#'
#' @docType methods
#' @name readyBasecalling
#' @rdname readyBasecalling
#'
#' @author Thomas Carroll and Marian Dore
#' @param Run_folders_WithRTA Run folders with RTA.complete file.
#' @param subFoldersFull All run folders.
#' @param config Config file (ini format) - see manual for full details.
#' @param bclVersion Character vector of "New" or "Old". Refers to versions prior to 1.2.1 (?) as old.

#' @return Shell commands for basecalling.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#'
#' library(raster)
#' library(XML)

#' subFoldersFull <- "/ifs/data/Hiseq/Runs/161105_D00467_0205_AC9L0AANXX"
#' Run_folders_WithRTA <- "/ifs/data/Hiseq/Runs/161105_D00467_0205_AC9L0AANXX"
#' fileLocations <- system.file("extdata",package="basecallQC")
#'
#' configFile <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' config = data.frame(readIniFile(configFile))

#' #bclcommands <- readyBasecalling(Run_folders_WithRTA,subFoldersFull,config,bclVersion="New")
#'
#' @export
readyBasecalling <- function(Run_folders_WithRTA,subFoldersFull,config,bclVersion="New"){
  runParams <- vector("list",length=length(Run_folders_WithRTA))
  shellBCLs <- vector()

  # samplesheet functions should be into own method

  for(i in 1:length(Run_folders_WithRTA)){
    currentRun <- subFoldersFull[grepl(Run_folders_WithRTA[i],subFoldersFull)]
    xmlFromPresentRunFolder <- xmlParse(file.path(currentRun,"runParameters.xml"))
    currentRunParameters <- xmlToDataFrame(xmlFromPresentRunFolder)
    currentRunParameters <- currentRunParameters[!is.na(currentRunParameters$ExperimentName),,drop=F]
    sampleSheetName <- file.path(currentRun,paste0(currentRunParameters$Barcode,".csv"))
    if(file.exists(sampleSheetName)){
      ss <- read.delim(sampleSheetName,sep=",",quote=NULL,header=T,stringsAsFactors=F)
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
      indexLoop <- 1
      for(l in uniqueIndexTypes){
        tempss <- ss[allIndexTypes %in% l,]
        tempss[is.na(tempss)] <- ""
        if(!ncol(tempss) > 9){
          diff <- 10-ncol(tempss)
          addColumns <- matrix(nrow=nrow(tempss),ncol=diff)
          colnames(addColumns) <- paste0("Dummy",seq(1,diff))
          tempss <- cbind(tempss,addColumns)
        }
        if(!file.exists(gsub("\\.csv",paste0("_",l,"\\.csv"),sampleSheetName))){
          if(bclVersion == "old"){
            write.table(tempss,file=gsub("\\.csv",paste0("_",l,"\\.csv"),sampleSheetName),quote=F,sep=",",row.names=F,col.names=T)
          }else{
            Sample_Project <- tempss$Project
            Lane <- tempss$Lane
            Sample_ID <- tempss$SampleID
            Sample_Name <- tempss$SampleID
            index <- tempss$Index
            if(any(grepl("-",index))){
              indexMat <- matrix(unlist(strsplit(tempss$Index,"-")),ncol=2,byrow=T)
              index <- indexMat[,1]
              index2 <- indexMat[,2]
              tempss_New <- cbind(Sample_Project,Lane,Sample_ID,Sample_Name,index,index2)
            }else{
              tempss_New <- cbind(Sample_Project,Lane,Sample_ID,Sample_Name,index)
            }
            dataSectionString <- "[Data]"
            write.table(dataSectionString,file=gsub("\\.csv",paste0("_",l,"\\.csv"),sampleSheetName),quote=F,sep=",",row.names=F,col.names=F)
            write.table(tempss_New,file=gsub("\\.csv",paste0("_",l,"\\.csv"),sampleSheetName),quote=F,sep=",",row.names=F,col.names=T,append=T)
          }
        }
        #if(!dir.exists(gsub("\\.csv",paste0("_",l),sampleSheetName))){
        #  dir.create(gsub("\\.csv",paste0("_",l),sampleSheetName),showWarnings = F)
        #}
        if(!file.exists(gsub("\\.csv",paste0("_",l,"\\.sh"),sampleSheetName))){
          indexLengths <- unlist(strsplit(l,"-"))
          if(as.numeric(as.vector(currentRunParameters$IndexRead1)) != 0){
            numberOfNs <- as.numeric(as.vector(currentRunParameters$IndexRead1)) - as.numeric(indexLengths[1])
            if(indexLengths[1] == 0){
              indexMask1 <- ",n*"
            }else{
              indexMask1 <-  paste0(",I",indexLengths[1],paste0(rep("n",numberOfNs),collapse=""))
            }
          }else{
            indexMask1 <- ""
          }
          if(as.numeric(as.vector(currentRunParameters$IndexRead2)) != 0){
            numberOfNs <- as.numeric(as.vector(currentRunParameters$IndexRead2)) - as.numeric(indexLengths[2])
            if(indexLengths[2] == 0){
              indexMask2 <- ",n*"
            }else{
              indexMask2 <-  paste0(",I",indexLengths[2],paste0(rep("n",numberOfNs),collapse=""))
            }
          }else{
            indexMask2 <- ""
          }
          if(currentRunParameters$Read2 > 0){
            usebasemask <- paste0("y*n",indexMask1,indexMask2,",y*n")
          }else{
            usebasemask <- paste0("y*n",indexMask1,indexMask2)
          }

          if(bclVersion == "Old"){
            programToRun <- as.vector(config[config$section == "programs" & config$name == "configureBclToFastq","value"])
            runBCLcommand <- paste0(
              programToRun,
              " --input-dir ",
              file.path(currentRun,as.vector(config[config$section == "paths" & config$name == "inputdir","value"])),
              " --sample-sheet ",
              gsub("\\.csv",paste0("_",l,"\\.csv"),sampleSheetName),
              " --fastq-cluster-count=0 --use-bases-mask ",usebasemask
              ,
              " --output-dir ",
              gsub("\\.csv",paste0("_",l),sampleSheetName),
              "\n"
            )
          }else{
            programToRun <- as.vector(config[config$section == "programs" & config$name == "BclToFastq2","value"])
            runBCLcommand <- paste0(
              programToRun,
              " --input-dir ",
              file.path(currentRun,as.vector(config[config$section == "paths" & config$name == "inputdir","value"])),
              " --sample-sheet ",
              gsub("\\.csv",paste0("_",l,"\\.csv"),sampleSheetName),
              " --use-bases-mask ",usebasemask
              ,
              " -r 4 -p 4 -w 4 -d 4 --output-dir ",
              gsub("\\.csv",paste0("_",l),sampleSheetName),
              "\n"
            )
          }

          write.table("#!/bin/bash",file=gsub("\\.csv",paste0("_",l,"\\.sh"),sampleSheetName),col.names=F,row.names=F,append=F,quote=F)
          write.table(runBCLcommand,file=gsub("\\.csv",paste0("_",l,"\\.sh"),sampleSheetName),col.names=F,row.names=F,append=T,quote=F)
          write.table("#!/bin/bash",file=gsub("\\.csv",paste0("_",l,"ForQSUB\\.sh"),sampleSheetName),col.names=F,row.names=F,append=F,quote=F)
          write.table(paste0("",gsub("\\.csv",paste0("_",l,"\\.sh"),sampleSheetName)),
                      file=gsub("\\.csv",paste0("_",l,"ForQSUB\\.sh"),sampleSheetName),
                      col.names=F,row.names=F,append=T,quote=F)

          shellBCLs[p] <- gsub("\\.csv",paste0("_",l,"ForQSUB\\.sh"),sampleSheetName)
          p <- p+1
        }
      }
    }else{
      stop("No samplesheet ",basename(sampleSheetName)," discovered for run ",basename(currentRun))
    }
  }
  return(bclVersion)
}



makeQsubs <- function(bclcommands,launch=TRUE){
  qsubCommands <- list()
  for(i in 1:length(bclcommands)){
    dir.create(gsub("ForQSUB\\.sh","",basename(bclcommands[i])),showWarnings = F)
    qsubCommands[[i]] <- paste0("qsub -v PATH -cwd ",
                                "-o ",gsub("\\.sh","\\.out",bclcommands[i])," ",
                                "-e ",gsub("\\.sh","\\.error",bclcommands[i])," ",
                                "-N ",gsub("\\.sh","",basename(bclcommands[i]))," ",
                                bclcommands[i]
    )
    if(launch){
      system(qsubCommands[[i]],intern=T)
    }
  }
  return(qsubCommands)
}

#' Illumina Basecalling functions.
#'
#' Launchs Illumina BCL commands using appropriate BioCparallel (eventually).
#'
#'
#' @docType methods
#' @name launchBCLcommands
#' @rdname launchBCLcommands
#'
#' @author Thomas Carroll and Marian Dore
#' @param Run_folders_WithRTA Run folders with RTA.complete file.
#' @param bclcommands BCL comands from readyBasecalling function.
#' @param config Config file (ini format) - see manual for full details.
#' @param Launch TRUE or FALSE. TRUE makes QSub call to launch BCL command. TRUE and FALSE return BCL command as character string.

#' @return BCL command as character string.
#' @import stringr XML RColorBrewer methods raster
#' @examples
#'
#' library(raster)
#' library(XML)

#' subFoldersFull <- "/ifs/data/Hiseq/Runs/161105_D00467_0205_AC9L0AANXX"
#' Run_folders_WithRTA <- "/ifs/data/Hiseq/Runs/161105_D00467_0205_AC9L0AANXX"
#' fileLocations <- system.file("extdata",package="basecallQC")
#'
#' configFile <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' config = data.frame(readIniFile(configFile))

#' #bclcommands <- readyBasecalling(Run_folders_WithRTA,subFoldersFull,config,bclVersion="New")
#' #qsubcommands <- makeQsubs(bclcommands[1],launch=T)

#' @export
launchBCLcommands <- function(bclcommands,launch=TRUE){
  qsubCommands <- list()
  for(i in 1:length(bclcommands)){
    dir.create(gsub("ForQSUB\\.sh","",basename(bclcommands[i])),showWarnings = F)
    qsubCommands[[i]] <- paste0("qsub -v PATH -cwd ",
                                "-o ",gsub("\\.sh","\\.out",bclcommands[i])," ",
                                "-e ",gsub("\\.sh","\\.error",bclcommands[i])," ",
                                "-N ",gsub("\\.sh","",basename(bclcommands[i]))," ",
                                bclcommands[i]
    )
    if(launch){
      system(qsubCommands[[i]],intern=T)
    }
  }
  return(qsubCommands)
}

split_path <- function(path) {
  if (dirname(path) %in% c(".", path)) return(basename(path))
  return(c(basename(path), split_path(dirname(path))))
}

indev_QCshortRead <- function(fastqs) {
  # fl <- dir(system.file(package="ShortRead", "extdata", "E-MTAB-1147"),
  #                   pattern="fastq",full.names = T)
  # names(fl) <- c("Myc1","Myc2")
  # refqqc <- ShortRead::qa(fl)
  # refqqc <- ShortRead::qa(fastqs)
  dirPath <- system.file(package="ShortRead", "extdata", "E-MTAB-1147")
  fls <- dir(dirPath, "fastq.gz", full=TRUE)

  refqqc <- ShortRead::qa(fastqs)

  # coll <- ShortRead:::QACollate(ShortRead:::QAFastqSource(fls),
  #                               ShortRead:::QAReadQuality(),
  #                               ShortRead:::QAAdapterContamination(),
  #                               ShortRead:::QANucleotideUse(),
  #                               ShortRead:::QAQualityUse(),
  #                               ShortRead:::QASequenceUse(),
  #                   ShortRead:::QAFrequentSequence(n=10),
  #                   ShortRead:::QANucleotideByCycle(),
  #                   ShortRead:::QAQualityByCycle())
  # refqqc <- ShortRead:::qa2(coll, verbose=TRUE)
  #
  return(refqqc)
}




