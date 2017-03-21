
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
  toPlot <- BCLQC@baseCallMetrics$convStatsProcessed %>% group_by_(.dots=as.list(groupByS)) %>% filter(Sample != "all") %>% summarise(Yield=sum(as.numeric(Yield)))
  toPlot <- toPlot %>% spread(Filter,Yield) %>% mutate(Ff=Raw-Pf) %>% dplyr:::select(-Raw) %>% tbl_df %>% gather(key=PassFilter,value=Yield,Ff,Pf)
  p <- ggplot(data=toPlot,aes_string(x=groupByG,y="Yield",fill="PassFilter"))+geom_bar(stat = "identity")+ coord_flip()
  return(p)
}

#' Functions to create basemasks for basecalling from Illumina samplesheet.
#'
#' Parses the Illumina samplesheet from run folders to create standardised
#'
#'
#' @docType methods
#' @name passFilterBoxplot
#' @rdname passFilterBoxplot
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
#' plot <- passFilterBoxplot(bclQC,groupBy = "Sample")

#' @export
passFilterBoxplot <- function(BCLQC,groupBy=c("Lane")){
  groupByS <- unique(c("Lane","Sample","Tile","Filter"))
  groupByG <- unique(c(groupBy))
  toPlot <- BCLQC@baseCallMetrics$convStatsProcessed %>% group_by_(.dots=as.list(groupByS)) %>% filter(Sample != "all") %>% summarise(Yield=sum(as.numeric(Yield)))
  toPlot <- toPlot %>% spread(Filter,Yield) %>% mutate(Ff=Raw-Pf) %>% dplyr:::select(-Raw) %>% tbl_df %>% gather(key=PassFilter,value=Yield,Ff,Pf)
  p <- ggplot(data=toPlot,aes_string(x=groupByG,y="Yield",fill="PassFilter"))+geom_violin()+ coord_flip()+facet_grid(PassFilter~.)
  return(p)
}

# fileLocations <- system.file("extdata",package="basecallQC")
# runXML <- dir(file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/"),pattern="runParameters.xml",full.names=TRUE)
# config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
# sampleSheet <- dir(file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/"),pattern="*\\.csv",full.names=TRUE)
# outDir <- file.path(fileLocations,"Runs/161105_D00467_0205_AC9L0AANXX/C9L0AANXX/")
# bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),outDir,verbose=FALSE)
# bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
# plot <- passFilterBar(bclQC,groupBy = "Lane")
# plot <- passFilterBoxplot(bclQC,groupBy = "Sample")
#
# fileLocations <- system.file("extdata",package="basecallQC")
# runXML <- dir(file.path(fileLocations,"Runs/170217_D00467_0227_ACADN3ANXX/"),pattern="runParameters.xml",full.names=TRUE)
# config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
# sampleSheet <- dir(file.path(fileLocations,"Runs/170217_D00467_0227_ACADN3ANXX/"),pattern="*\\.csv",full.names=TRUE)
# bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
# bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#
# fileLocations <- system.file("extdata",package="basecallQC")
# runXML <- dir(file.path(fileLocations,"Runs/170222_D00467_0228_ACAL6PANXX/"),pattern="runParameters.xml",full.names=TRUE)
# config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
# sampleSheet <- dir(file.path(fileLocations,"Runs/170222_D00467_0228_ACAL6PANXX/"),pattern="*\\.csv",full.names=TRUE)
# bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
# bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#
#
# fileLocations <- system.file("extdata",package="basecallQC")
# runXML <- dir(file.path(fileLocations,"Runs/170303_D00467_0230_BCA7RDANXX/"),pattern="runParameters.xml",full.names=TRUE)
# config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
# sampleSheet <- dir(file.path(fileLocations,"Runs/170303_D00467_0230_BCA7RDANXX/"),pattern="*\\.csv",full.names=TRUE)
# bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
# bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
#
# fileLocations <- system.file("extdata",package="basecallQC")
# runXML <- dir(file.path(fileLocations,"Runs/170303_D00467_0231_ACAK02ANXX/"),pattern="runParameters.xml",full.names=TRUE)
# config <- dir(fileLocations,pattern="config.ini",full.names=TRUE)
# sampleSheet <- dir(file.path(fileLocations,"Runs/170303_D00467_0231_ACAK02ANXX/"),pattern="*\\.csv",full.names=TRUE)
# bcl2fastqparams <- setBCL2FastQparams(runXML,config,runDir=getwd(),verbose=FALSE)
# bclQC <- basecallQC(bcl2fastqparams,RunMetaData=NULL,sampleSheet)
