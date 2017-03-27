configParams <- function(config){
  data.frame(readIniFile(config)) %>% tbl_df
}

