library(raster)
library(XML)

subFoldersFull <- "/ifs/data/Hiseq/Runs/161105_D00467_0205_AC9L0AANXX"
Run_folders_WithRTA <- "/ifs/data/Hiseq/Runs/161105_D00467_0205_AC9L0AANXX"

########################Basecalling########################################
readyBasecalling <- function(Run_folders_WithRTA,subFoldersFull,config,bclVersion="New"){
  runParams <- vector("list",length=length(Run_folders_WithRTA))
  shellBCLs <- vector()
  for(i in 1:length(Run_folders_WithRTA)){
    currentRun <- subFoldersFull[grepl(Run_folders_WithRTA[i],subFoldersFull)]
    xmlFromPresentRunFolder <- xmlParse(file.path(currentRun,"runParameters.xml"))
    currentRunParameters <- xmlToDataFrame(xmlFromPresentRunFolder)
    currentRunParameters <- currentRunParameters[!is.na(currentRunParameters$ExperimentName),,drop=F]
    sampleSheetName <- file.path(currentRun,paste0(currentRunParameters$Barcode,".csv"))
    if(file.exists(sampleSheetName)){  # check whether the sample sheet exist
        ss <- read.delim(sampleSheetName,sep=",",quote=NULL,header=T,stringsAsFactors=F)
        # (1) check SampleID integrity
            ss$SampleID <- gsub("^X","Sample_",validNames(ss$SampleID))
            ss$SampleID <-gsub("\\?|\\(|\\)|\\[|\\]|\\\\|/|\\=|\\+|<|>|\\:|\\;|\"|\'|\\*|\\^|\\||\\&|\\.","_",ss$SampleID)
        # (2) check Project integrity
            # "UP": Unknown project
            ss[ss$Project==""|ss$Project==" "|is.na(ss$Project),]$Project<-"UP"
            
        # (3) a. check the index sequence in $Index column, i.e. remove the space character
            ss$Index <- gsub(" ","",ss$Index)
            # b. check the index sequence length in $Index column
            toBeTrimmed <- which(lapply(ss$Index,nchar) > as.numeric(currentRunParameters$IndexRead1))
            for(i4trim in toBeTrimmed){
              #ss$Index[i4trim] <-  substr(ss$Index[i4trim],0,as.numeric(currentRunParameters$IndexRead1))
              ss$Index[i4trim] <-  substr(ss$Index[i4trim],1,as.numeric(currentRunParameters$IndexRead1)-1)
            }
        # (4) after checking the SampleID and index, print out message
        message("Read samplesheet ",basename(sampleSheetName)," discovered for run ",currentRun)
        index1Lengths <- unlist(lapply(ss$Index,function(x)nchar(x)))
        # (5) check whether there is "Index2" column in the sample sheet
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
            # after merging Index2 information to Index, delete the Index2 column from ss
            ss <- ss[,-grep("Index2",colnames(ss))]
        }else{
            allIndexTypes <- paste0(index1Lengths)
            uniqueIndexTypes <- unique(allIndexTypes)
            #ss$SampleID <- gsub("[[:punct:]]", "_", ss$SampleID).)
            #ss$SampleID <- gsub("[^[:alnum:]]", "_", ss$SampleID).)
        }
        indexLoop <- 1  # unused parameter; for future development
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
  return(shellBCLs)
}


###########################makeqsubs################################
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



config = data.frame(readIniFile("~/genomicsProcessing/config.ini"))
p <- 1


bclcommands <- readyBasecalling(Run_folders_WithRTA,subFoldersFull,config,bclVersion="New")
qsubscommands <- makeQsubs(bclcommands,launch=TRUE)

qsubscommands
shellBCLs
bclcommands


      
            
            
            
            