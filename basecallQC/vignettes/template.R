## ----rmarkdown, eval=FALSE-----------------------------------------------
#  install.packages("rmarkdown")

## ----code, echo = FALSE--------------------------------------------------
fileLocations <- system.file("extdata",package="basecallQC")

demuxStats <- dir(fileLocations,pattern="DemultiplexingStats.xml",full.names=TRUE)
processDemultiplex(demuxStats)

## ----macro, echo=FALSE---------------------------------------------------
macro <- function(name, pkg, description)
    sprintf('`` `r %s("%s")` `` %s %s', name, pkg,
            description, do.call(name, list(pkg)))

## ----figure, fig.width=4.0, fig.height=4.4-------------------------------
v = seq(0, 60i, length = 1000)
plot(abs(v)*exp(v), type = "l", col = "Royalblue")

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

