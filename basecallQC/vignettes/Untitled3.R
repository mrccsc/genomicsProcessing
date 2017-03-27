## ----rmarkdown, eval=FALSE-----------------------------------------------
#  install.packages("rmarkdown")

## ----code, echo = FALSE--------------------------------------------------
## The following redefinitions are only for printing the verbatim header below
doc_date = function() "`r doc_date()`"
pkg_ver = function(name) sprintf("`r pkg_ver('%s')`", name)

## ----macro, echo=FALSE---------------------------------------------------
macro <- function(name, pkg, description)
    sprintf('`` `r %s("%s")` `` %s %s', name, pkg,
            description, do.call(name, list(pkg)))

## ----figure, fig.width=4.0, fig.height=4.4-------------------------------
v = seq(0, 60i, length = 1000)
plot(abs(v)*exp(v), type = "l", col = "Royalblue")

## ----sessionInfo, echo=FALSE---------------------------------------------
sessionInfo()

