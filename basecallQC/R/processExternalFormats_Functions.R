#' Ini config file reader.
#'
#' @name configParams
#' @rdname configParams
#' @param Params Any parameters
#' @export
configParams <- function(params){
  data.frame(readIniFile(configFile))
}

