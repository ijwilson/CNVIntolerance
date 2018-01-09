
## PRISM is available at http://compbio.cs.toronto.edu/prism

library(here)
source(here("R","prepare.R"))
install.load("data.table")
install.load.bioc("GenomicRanges", "GenomicAlignments")

del <- fread(dropbox("VUR/DEL_results.txt"))
dup <- fread(dropbox("VUR/DUP_results.txt"))

## There seem to be far too many of these

del_filter <- del[left_match>60 & right_match > 60]

conv <- function(xx) {
  g <- GRanges(seqnames=xx$X.chr, 
               IRanges(start=xx$gap_start, end=xx$gap_end)
               , ID=xx$ID
               , type=xx$type
               , support = xx$support_read_num
               , left_match = xx$left_match
               ,right_match = xx$right_match)
}


gdel <- conv(del_filter)
cov_gdel <- coverage(gdel)
## split by ID
sp <- split(gdel, gdel$ID)
## just get the overlaps.  
sp_reduce <- lapply(sp, reduce)
id_names <- names(sp_reduce)
spb <- mapply(function(x,y) {x$ID <- y; return(x)}, sp_reduce, id_names)


r <- Reduce(c, spb)
length(r)

cov2 <- coverage(r)
peaks <- slice(cov2$chr10, lower=24)
length(peaks)
peaksb <- peaks[width(peaks)>20]
gdup <- conv(dup)



table(dup$ID)

hist(table(del$ID))



## One approach here mght be to go along the chromosome and count the total number of individuals that 
## have deletions overlapping a particular point

## If you do this then we can think about just those that have an appreciable chance of being the 
## particular VUR cnv.