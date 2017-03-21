
#' Functions to create basemasks for basecalling from Illumina samplesheet.
#'
#' Parses the Illumina samplesheet from run folders to create standardised
#'
#'
#' @docType methods
#' @name passFilterBar
#' @rdname passFilterBar
#'
#' @author Thomas Carroll and Marian Dore
#' @param BCLQC baseCallQC A  basecall QC object
#' @param groupBy Character vector of lane and/or Sample
#' @return ggPlotObject A ggplot2 object.
#' @import stringr XML RColorBrewer methods raster BiocStyle
#' @examples
#' fileLocations <- system.file("extdata",package="basecallQC")
#' runXML <- dir(fileLocations,pattern="runParameters.xml",full.names=TRUE)
#' config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
#' sampleSheet <- dir(fileLocations,pattern="*\\.csv",full.names=TRUE)
#' outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
#' bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
#' bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#' plot <- passFilterBar(bclQC)

#' @export
passFilterBar <- function(BCLQC,groupBy=c("Lane")){
  groupByS <- unique(c(groupBy,"Filter"))
  groupByG <- unique(c(groupBy))
  toPlot <- bclQC@baseCallMetrics$convStatsProcessed %>% group_by_(.dots=as.list(groupByS)) %>% filter(Sample != "all") %>% summarise(Yield=sum(as.numeric(Yield)))
  toPlot <- toPlot %>% spread(Filter,Yield) %>% mutate(Ff=Raw-Pf) %>% dplyr:::select(-Raw) %>% tbl_df %>% gather(key=PassFilter,value=Yield,Ff,Pf)
  p <- ggplot(data=toPlot,aes_string(x=groupByG,y="Yield",fill="PassFilter"))+geom_bar(stat = "identity")+ coord_flip()
  return(p)
}
