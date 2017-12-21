## VUR

## PRISM is available at http://compbio.cs.toronto.edu/prism

library(here)
source(here("R","prepare.R"))
install.load("data.table")

del <- fread(dropbox("VUR/DEL_results.txt"))
dup <- fread(dropbox("VUR/DUP_results.txt"))


table(dup$ID)

hist(table(del$ID))

plot(del$gap_start)


## One approach here mght be to go along the chromosome and count the total number of individuals that 
## have deletions overlapping a particular point

## If you do this then we can think about just those that have an appreciable chance of being the 
## particular VUR cnv.