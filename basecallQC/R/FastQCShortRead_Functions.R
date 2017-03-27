qcShortRead <- function(fastqs,reportOutDir=getwd()){
  #qaSamples <- lapply(fastqs,ShortRead::qa)
  load("~/Documents/test2/tesr.RData")
  names(qaSamples) <- gsub("\\.fastq.*|fq.*","",basename(fastqs))
  reportOutFiles <- lapply(names(qaSamples),function(x){
    report(qaSamples[[x]],dest = file.path(reportOutDir,x))
  })
  links <- file.path(reportOutDir,names(qaSamples))
  links <- paste0("<a href='",
                  (file.path(links,"index.html")),"'>",
                  names(qaSamples),"</a>")
  qaSampleFrame <- do.call(ShortRead::rbind,qaSamples)
  fqQCTable <- data.frame(SampleNames=links,
                          as.data.frame(qaSampleFrame@.srlist$readCounts,
                          ))
  return(list(FQQC_Table = fqQCTable,ShortReadQC=qaSampleFrame))
}



