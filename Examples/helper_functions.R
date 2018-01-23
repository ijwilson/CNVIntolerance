



## Stolen from the install.load package
install.load <- function (package1, ...)  {   
  
  # convert arguments to vector
  packages <- c(package1, ...)
  
  # start loop to determine if each package is installed
  for(package in packages){
    
    # if package is installed locally, load
    if(package %in% rownames(installed.packages()))
      do.call('library', list(package))
    
    # if package is not installed locally, download, then load
    else {
      install.packages(package,repos="https://cran.r-project.org/")
      do.call("library", list(package))
    }
  } 
}


install.load.bioc <- function (package1, ...)  {   
  source("https://www.bioconductor.org/biocLite.R")
  # convert arguments to vector
  packages <- c(package1, ...)
  
  # start loop to determine if each package is installed
  for(package in packages){
    
    # if package is installed locally, load
    if(package %in% rownames(installed.packages()))
      do.call('library', list(package))
    
    # if package is not installed locally, download, then load
    else {
      biocLite(package)
      do.call("library", list(package))
    }
  } 
}