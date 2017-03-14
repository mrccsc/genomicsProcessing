#' Demultiplex parser.
#'
#' Parses the DemultiplexingStats.xml from Illumina Basecalling
#'
#'
#' @docType methods
#' @name processDemultiplex
#' @rdname processDemultiplex
#'
#' @author Thomas Carroll
#'
#' @param Projects Demultiplex locations.
#' @return A datatable of XML results.
#' @import stringr XML RColorBrewer methods raster dplyr magrittr tidyr reshape2 ggplot2 lubridate
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#'
#' demuxStats <- dir(fileLocations,pattern="DemultiplexingStats.xml",full.names=TRUE)
#' demuxProcessed <-processDemultiplex(demuxStats)
#'
#' @export
processDemultiplex <- function(demuxStatsXML){
  if(is.null(demuxStatsXML)){
    demuxStatsXML <- dir(fileLocations,
        pattern="DemultiplexingStats.xml",
        full.names=TRUE,
        recursive = T)
  }
  demuxStatsXMLparse <- xmlTreeParse(demuxStatsXML)
  demuxStatsXML_root <- xmlRoot(demuxStatsXMLparse)
  Projects <- demuxStatsXML_root[[1]]
  flowcellID <- xmlAttrs(Projects)
  library(dplyr)
  Projects_Info <- list()
  for(p in 1:length(Projects)){
    Project <- Projects[[p]]
    Project_Name <- xmlAttrs(Project)
    Project_Sample_Info <- list()
    for(s in 1:length(Project)){
      Sample <- Project[[s]]
      Sample_Name <- xmlAttrs(Sample)
      Sample_BarcodeExpected <- Sample[[1]]
      Sample_BarcodeExpected_Name <- xmlAttrs(Sample[[1]])
      Project_Sample_BarcodeExpected_Lane_Info <- list()
      for(b in 1:length(Sample_BarcodeExpected)){
        Lane <- Sample_BarcodeExpected[[b]]
        Lane_Name <- xmlAttrs(Lane)
        Lane_BarcodeCount <- as.integer(xmlValue(Lane[["BarcodeCount"]]))
        Lane_PerfectBarcodeCount <- as.integer(xmlValue(Lane[["PerfectBarcodeCount"]]))
        Project_Sample_BarcodeExpected_Lane_Info[[b]] <- data.frame(BarcodeCount = Lane_BarcodeCount,
                                                                    PerfectBarcodeCount = Lane_PerfectBarcodeCount)
        names(Project_Sample_BarcodeExpected_Lane_Info)[b] <- paste0("Lane",xmlAttrs(Lane))

      }
      psbeliMat <- sapply(Project_Sample_BarcodeExpected_Lane_Info,function(x)x)
      psbeliMatDF <- data.frame(BarcodeStat=rownames(psbeliMat),psbeliMat)
      psbeliMatDF <- tidyr::gather(psbeliMatDF,Lane,Count,-BarcodeStat)

      Project_Sample_Info[[s]] <- data.frame(Barcode = rep(Sample_BarcodeExpected_Name,nrow(psbeliMatDF)),
                                             psbeliMatDF)
      names(Project_Sample_Info)[s] <- Sample_Name

    }
    psiMat <- do.call(rbind,Project_Sample_Info)
    psiMatDF <- data.frame(Sample=gsub("\\..*","",rownames(psiMat)),psiMat)
    Projects_Info[[p]] <- data.frame(Project = rep(Project_Name,nrow(psiMat)),
                                     psiMatDF)
    names(Projects_Info)[p] <- Project_Name

  }
  Projects_DF2 <- do.call(rbind,Projects_Info)
  rownames(Projects_DF2) <- NULL
  return(Projects_DF2)

}



#' Converstion stats parser.
#'
#' Parses the ConversionStats.xml from Illumina Basecalling
#'
#'
#' @docType methods
#' @name processConvStats
#' @rdname processConvStats
#'
#' @author Thomas Carroll
#'
#' @param Projects ConversionStats locations.
#' @return A datatable of XML results.
#' @import stringr XML RColorBrewer methods raster
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#'
#' convStats <- dir(fileLocations,pattern="ConversionStats.xml",full.names=TRUE)
#' convProcessed <-processConvStats(convStats)
#'
#' @export

processConvStats <- function(ConvStatsXML){

  if(is.null(ConvStatsXML)){
    ConvStatsXML <- dir(fileLocations,
                         pattern="ConversionStats.xml",
                         full.names=TRUE,
                         recursive = T)
  }
  convStatsXMLparse <- xmlTreeParse(ConvStatsXML)
  convStatsXML_root <- xmlRoot(convStatsXMLparse)
  Projects <- convStatsXML_root[[1]]
  flowcellID <- xmlAttrs(Projects)

  Projects <- Projects[names(Projects) == "Project"]
  library(dplyr)
  Projects_Info <- list()
  for(p in 1:length(Projects)){
    Project <- Projects[[p]]
    Project_Name <- xmlAttrs(Project)
    Project_Sample_Info <- list()
    for(s in 1:length(Project)){
      Sample <- Project[[s]]
      Sample_Name <- xmlAttrs(Sample)
      Sample_BarcodeExpected <- Sample[[1]]
      Sample_BarcodeExpected_Name <- xmlAttrs(Sample[[1]])
      Project_Sample_BarcodeExpected_Lane_Info <- list()
      Lane_Info <- list()
      for(b in 1:length(Sample_BarcodeExpected)){
        Lane <- Sample_BarcodeExpected[[b]]
        Lane_Name <- xmlAttrs(Lane)
        Tile_Info <- list()
        for(t in 1:length(Lane)){
          Tile <- Lane[[t]]
          Tile_Name <- xmlAttrs(Tile)
          FilterState_Info <- list()
          for(f in 1:length(Tile)){
            FilterState <- Tile[[f]]
            FilterState_Name <- xmlName(FilterState)
            FilterState_ClusterCount <- as.integer(xmlValue(FilterState[["ClusterCount"]]))
            ReadNumber_Info <- list()
            for(r in 2:length(FilterState)){

              Read <-  FilterState[[r]]
              ReadNumber <-  xmlAttrs(Read)
              ReadNumber_Yield <- as.integer(xmlValue(Read[["Yield"]]))
              ReadNumber_YieldQ30 <- as.integer(xmlValue(Read[["YieldQ30"]]))
              ReadNumber_QualityScoreSum <- as.numeric(xmlValue(Read[["QualityScoreSum"]]))
              ReadNumber_Info[[r-1]] <- data.frame(Yield = ReadNumber_Yield,
                                                   Yield30 = ReadNumber_YieldQ30,
                                                   QualityScoreSum = ReadNumber_QualityScoreSum
              )

              names(ReadNumber_Info)[r-1] <- paste0("Read_",ReadNumber)

            }
            rniMat <- do.call(rbind,ReadNumber_Info)
            rniMatDF <- data.frame(ReadNumber=rownames(rniMat),rniMat)
            FilterState_Info[[f]] <- data.frame(Filter = rep(FilterState_Name,nrow(rniMatDF)),
                                                rniMatDF)
            names(FilterState_Info)[f] <- FilterState_Name

          }
          fsiMat <- do.call(rbind,FilterState_Info)
          fsiMatDF <- data.frame(fsiMat)
          Tile_Info[[t]] <- data.frame(Tile = rep(Tile_Name,nrow(fsiMatDF)),
                                       fsiMatDF)
          names(Tile_Info)[t] <- paste0("Tile_",Tile_Name)

        }
        tiMat <- do.call(rbind,Tile_Info)
        tiMatDF <- data.frame(tiMat)
        Lane_Info[[b]] <- data.frame(Lane = rep(Lane_Name,nrow(tiMatDF)),
                                     tiMatDF)
        names(Lane_Info)[b] <- paste0("Lane_",Lane_Name)
      }
      liMat <- do.call(rbind,Lane_Info)
      liMatDF <- data.frame(liMat)
      Project_Sample_Info[[s]] <- data.frame(Sample = rep(Sample_Name,nrow(liMatDF)),
                                             liMatDF)
      names(Project_Sample_Info)[s] <- paste0(Sample_Name)

    }
    psi2Mat <- do.call(rbind,Project_Sample_Info)
    psi2MatDF <- data.frame(psi2Mat)
    Projects_Info[[p]] <- data.frame(Project = rep(Project_Name,nrow(psi2MatDF)),
                                     psi2MatDF)
    names(Projects_Info)[p] <- paste0(Project_Name)

  }
  Projects_DF <- do.call(rbind,Projects_Info)
  rownames(Projects_DF) <- NULL
  return(Projects_DF)

}

#' Generate per sample summary statistics
#'
#' Creates per sample summary statistics from demultiplex results
#'
#'
#' @docType methods
#' @name summariseDemuxStats
#' @rdname summariseDemuxStats
#'
#' @author Thomas Carroll
#'
#' @param DemuxStats Results from a call to processDemultiplex.
#' @return A datatable of summarised per sample results.
#' @import stringr XML RColorBrewer methods raster
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#'
#' demuxStats <- dir(fileLocations,pattern="DemultiplexingStats.xml",full.names=TRUE)
#' demuxProcessed <- processDemultiplex(demuxStats)
#' samplesDemuxStats <- summariseDemuxStats(demuxProcessed)
#'
#' @export
summariseDemuxStats <- function(demuxProcessed, plot=T){
  #Lane_Stats4 <- Test %>% filter(Sample != "all" & BarcodeStat == "BarcodeCount") %>% group_by(Project,Sample) %>% summarise(Count=sum(as.numeric(Count)))
  #Lane_Stats <- Projects_DF %>% filter(Sample == "all") %>% group_by(Lane,Filter) %>% summarise(sum(Yield))


  temp <- demuxProcessed %>% tbl_df %>% mutate(Count = as.numeric(Count)) %>%
            filter(Sample != "all" & BarcodeStat == "BarcodeCount") %>%
            filter(Project != "default") # Computing percent label text and position for pie chart
  # temp <- temp %>% group_by(Lane) %>% mutate(labelperc=round(Count/sum(Count),2)*100) %>% group_by(Lane) %>% mutate(pos = cumsum(labelperc)- labelperc/2)
  # p1 <- temp %>% filter(Project != "default") %>% ggplot(aes(x=Project,y=Count,fill=Project))+geom_violin(alpha=0.3,scale="width")+geom_jitter(alpha=0.6)+theme(legend.position="bottom")
  # p2 <- ggplot(temp,aes(x=Lane,y=Count,fill=Sample))+geom_bar(stat = "identity")+theme_bw()+theme(legend.position="bottom")
  # p3 <- ggplot(temp,aes(x=Sample,y=Count,fill=Lane))+geom_bar(stat = "identity")+theme_bw()+coord_flip()+theme(legend.position="bottom")
  # if(plot){
  #   print(p1)
  #   print(p2)
  #   print(p3)
  # }
  #return(list(Summary=temp,Boxplot=p1,StackedBar=p2,Bar=p3))
  Lane_perTileStats <- demuxProcessed %>% group_by(Lane,Tile,Filter) %>% filter(Sample != "all") %>% summarise(Yield=sum(as.numeric(Yield)))
  LaneSample_perTileStats <- demuxProcessed %>% group_by(Lane,Sample,Tile,Filter) %>% filter(Sample != "all") %>% summarise(Yield=sum(as.numeric(Yield)))
  Sample_Stats <- demuxProcessed %>% filter(Sample != "all") %>% group_by(Sample,Filter) %>% summarise(Yield=sum(as.numeric(Yield)))
  Lane_Stats <- demuxProcessed %>% filter(Sample != "all") %>% group_by(Lane,Filter) %>% summarise(Yield=sum(as.numeric(Yield)))

  return(list(Summary=temp))
}

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
#'
#' convStats <- dir(fileLocations,pattern="ConversionStats.xml",full.names=TRUE)
#' convStatsProcessed <- processConvStats(convStats)
#' summarisedConvStats <- summariseConvStats(convStatsProcessed)
#'
#' @export
summariseConvStats <- function(demuxProcessed, plot=T){
  #Lane_Stats4 <- Test %>% filter(Sample != "all" & BarcodeStat == "BarcodeCount") %>% group_by(Project,Sample) %>% summarise(Count=sum(as.numeric(Count)))
  #Lane_Stats <- Projects_DF %>% filter(Sample == "all") %>% group_by(Lane,Filter) %>% summarise(sum(Yield))


  Lane_perTileStats <- demuxProcessed %>% group_by(Lane,Tile,Filter) %>% filter(Sample != "all") %>% summarise(Yield=sum(as.numeric(Yield)))
  LaneSample_perTileStats <- demuxProcessed %>% group_by(Lane,Sample,Tile,Filter) %>% filter(Sample != "all") %>% summarise(Yield=sum(as.numeric(Yield)))
  Sample_Stats <- demuxProcessed %>% filter(Sample != "all") %>% group_by(Sample,Filter) %>% summarise(Yield=sum(as.numeric(Yield)))
  Lane_Stats <- demuxProcessed %>% filter(Sample != "all") %>% group_by(Lane,Filter) %>% summarise(Yield=sum(as.numeric(Yield)))
  #Lane_Stats <- Projects_DF %>% filter(Sample == "all") %>% group_by(Lane,Filter) %>% summarise(sum(Yield))
  #p3 <- ggplot(data=Lane_perTileStats,aes(x=Lane,y=Yield))+geom_violin()
  #ggplot(data=Sample_Stats,aes(x=Sample,y=Yield,fill=Filter))+geom_boxplot()
  # if(plot){
  #   print(p3)
  # }
  # return(list(Lane_perTileStats=Lane_perTileStats,
  #             LaneSample_perTileStats=LaneSample_perTileStats,
  #             Sample_Stats=Sample_Stats,
  #             Lane_Stats=Lane_Stats,
  #             ViolinPlot=p3))
}


#' Parses parameters for illumina basecalling.
#'
#' Creates data.frame of run parameters.
#'
#'
#' @docType methods
#' @name runParameters
#' @rdname runParameters
#'
#' @author Thomas Carroll
#'
#' @param paramXML file path to runParameters.xml .
#' @return A datatable of run parameter.
#' @import stringr XML RColorBrewer methods raster dplyr data.table stringr
#' @examples
#'
#' fileLocations <- system.file("extdata",package="basecallQC")
#'
#' runParameters <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' runParametersProcessed <- runParameters(runParameters)
#'
#' @export
runParameters <- function(runParameters = NULL){
  if(is.null(runParameters)){
    runParameters = dir(pattern="runParameters.xml",
                        recursive=T,
                        full.names=T)
  }
  xmlFromRunParameters <- xmlParse(runParameters)
  currentRunParameters <- xmlToDataFrame(xmlFromRunParameters)
  currentRunParameters <- currentRunParameters[!is.na(currentRunParameters$ExperimentName),,drop=F]
  return(currentRunParameters)
}


