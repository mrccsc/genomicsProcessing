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
    Sample_BarcodeExpected_Name <- xmlAttrs(Sample[[i]])
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
Projects_DF <- do.call(rbind,Projects_Info)
rownames(Projects_DF) <- NULL

