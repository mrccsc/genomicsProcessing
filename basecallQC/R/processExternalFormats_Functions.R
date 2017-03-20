#' Ini config file reader.
#'
#' @name configParams
#' @rdname configParams
#' @param Config Any parameters
#' @export
configParams <- function(Config){
  data.frame(readIniFile(config)) %>% tbl_df
}

