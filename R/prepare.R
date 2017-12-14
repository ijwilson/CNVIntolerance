## Prepare dropbox path

if (TRUE) {
  library(jsonlite)
if (.Platform$OS.type=="unix") {
  file <- file.path(Sys.getenv("HOME"), ".dropbox/info.json")
} else {
  file <- file.path(Sys.getenv("LOCALAPPDATA"), "Dropbox/info.json")
}
document <- fromJSON(readLines(file, warn = FALSE))
droppath <- document$personal$path

dropbox <- function(filename, dir="CNVIntolerance") {
  file.path(droppath, dir, filename)
}

rm(file, document)
}
HOME <- Sys.getenv("HOME")

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
## get filename for the posted directory with a timestamp
post <- function(filename) {
  timestring <- format(Sys.Date())
 # timestring <- gsub(" ", "_", timestring)
#  timestring <- gsub(":", "-", timestring)
  stem <- strsplit(filename, "\\.")[[1]][1]
  ext <- strsplit(filename, "\\.")[[1]][2]
  newfilename <- paste(stem, "_", timestring, ".",ext, sep="")
  here("posted", newfilename)
}

