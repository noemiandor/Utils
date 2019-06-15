source_https <- function(url, ...) {
  library(RCurl)
  
  script <- getURL(url, ssl.verifypeer = FALSE)
  
  eval(parse(text = script))
}