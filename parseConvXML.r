

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
    Sample_BarcodeExpected_Name <- xmlAttrs(Sample[[i]])
    Project_Sample_BarcodeExpected_Lane_Info <- list()
    for(b in 1:length(Sample_BarcodeExpected)){
      Lane <- Sample_BarcodeExpected[[b]]
      Lane_Name <- xmlAttrs(Lane)
      for(t in 1:length(Lane)){
        Tile <- Lane[[t]]
        Tile_Name <- xmlAttrs(Tile)
        for(f in 1:length(Tile)){
          FilterState <- Tile[[f]]
          FilterState_Name <- xmlAttrs(FilterState)
          FilterState_ClusterCount <- as.integer(xmlValue(FilterState[["ClusterCount"]]))
          for(r in 2:length(FilterState)){
            
            Read <-  FilterState[[r]]
            ReadNumber <-  xmlAttrs(Read)
            ReadNumber_Yield <- as.integer(xmlValue(FilterState[["Yield"]]))
            ReadNumber_YieldQ30 <- as.integer(xmlValue(FilterState[["YieldQ30"]]))
            ReadNumber_QualityScoreSum <- as.integer(xmlValue(FilterState[["QualityScoreSum"]]))
          }
        }
      }
    }
  }
}
