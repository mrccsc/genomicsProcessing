library(XML)
demuxStatsXML <- xmlTreeParse("~/Downloads/Stats/DemultiplexingStats.xml")
demuxStatsXML_root <- xmlRoot(demuxStatsXML)

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

Lane_Stats4 <- Projects_DF2 %>% filter(Sample != "all" & BarcodeStat == "PerfectBarcodeCount") %>% group_by(Project,Sample) %>% summarise(Count=sum(as.numeric(Count)))
#Lane_Stats <- Projects_DF %>% filter(Sample == "all") %>% group_by(Lane,Filter) %>% summarise(sum(Yield))

Lane_Stats4 %>% filter(Project != "default") %>% ggplot(aes(x=Project,y=Count))+geom_violin(alpha=0.3,scale="width")+geom_jitter(alpha=0.6)

temp <- Projects_DF2 %>% tbl_df %>% mutate(Count = as.numeric(Count)) %>%
  filter(Sample != "all" & BarcodeStat == "PerfectBarcodeCount") %>% 
  filter(Project != "default")
ggplot(temp,aes(x=Lane,y=Count,fill=Sample))+geom_bar(stat = "identity")


convStatsXML <- xmlTreeParse("~/Downloads/Stats/ConversionStats.xml")
convStatsXML_root <- xmlRoot(convStatsXML)

Projects <- convStatsXML_root[[1]]
Projects <- Projects[names(Projects) == "Project"]
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

Lane_perTileStats <- Projects_DF %>% group_by(Lane,Tile,Filter) %>% filter(Sample != "all") %>% summarise(Yield=sum(as.numeric(Yield)))
Sample_Stats <- Projects_DF %>% filter(Sample != "all") %>% group_by(Sample,Filter) %>% summarise(Yield=sum(as.numeric(Yield)))
Lane_Stats4 <- Projects_DF %>% filter(Sample != "all") %>% group_by(Lane,Filter) %>% summarise(Yield=sum(as.numeric(Yield)))
#Lane_Stats <- Projects_DF %>% filter(Sample == "all") %>% group_by(Lane,Filter) %>% summarise(sum(Yield))

ggplot(data=Lane_perTileStats,aes(x=Lane,y=Yield))+geom_violin()

