## liftover
## scripts and code to convert cnv files from one to another coordincates

## Height
library(here)
source(here("R","prepare.R"))
install.load("data.table")
install.load.bioc("GenomicRanges", "rtracklayer")

height <- fread(dropbox("height_CNV_association_41467_2017_556_MOESM2_ESM.csv"))   ## standard CNV
height$ID <-  paste("ID", height$CHR, height$BP, sep="_")

cnvdat <- height

write.table(cbind(paste("chr",cnvdat$CHR,sep=""), cnvdat$BP, cnvdat$BP+1, cnvdat$ID,
            col.names = FALSE, row.names=FALSE, quote=FALSE, file=dropbox("height_cnv_hg18.bed"))


## liftOver cnv_hg18.bed hg18ToHg19.over.chain.gz cnv_hg19.bed unmapped

## Read the non-gene centric results.
scz.dup <- read.table("C:\\Users\\nijw\\Dropbox/CNVIntolerance/pgc_cnv/PGC_41K_QC_dup.cnv.results", header=TRUE)
scz.del <- read.table("C:\\Users\\nijw\\Dropbox/CNVIntolerance/pgc_cnv/PGC_41K_QC_del.cnv.results", header=TRUE)
scz.dup$ID <-  paste("ID", scz.dup$CHR, scz.dup$BP, sep="_"))
scz.dup$ID <-  paste("ID", scz.dup$CHR, scz.dup$BP, sep="_"))

cnvdat <- scz.dup

write.table(cbind(paste("chr",cnvdat$CHR,sep=""), cnvdat$BP, cnvdat$BP+1,cnvdat$ID,
            col.names = FALSE, row.names=FALSE, quote=FALSE, file=dropbox("scz.dup_cnv_hg18.bed"))

cnvdat <- scz.del

write.table(cbind(paste("chr",cnvdat$CHR,sep=""), cnvdat$BP, cnvdat$BP+1, cnvdat$ID,
            col.names = FALSE, row.names=FALSE, quote=FALSE, file=dropbox("scz.del_cnv_hg18.bed"))

if (FALSE) {
command <- "
  liftOver height_cnv_hg18.bed hg18ToHg19.over.chain.gz height_cnv_hg19.bed unmapped_height
  liftOver scz.dup_hg18.bed hg18ToHg19.over.chain.gz scz.dup_hg19.bed unmapped_scz.dup
  liftOver scz.del_hg18.bed hg18ToHg19.over.chain.gz scz.del_hg19.bed unmapped_scz.del
"
}

## now need to remap to the original data

height.remap <- fread(dropbox("height_cnv_hg19.bed"))
colnames(height.remap) <- c("chr19", "start19", "end19", "ID")
setkey(height.remap,  "ID")
setkey(height, "ID")

h <- height[height.remap]
h$start18 <- h$BP
h$BP <- h$start19
h$chr19 <- NULL
h$end19 <- NULL
h$ID <- NULL

scz.dup.remap <- fread(dropbox("scz.dup_hg19.bed"))
colnames(scz.dup.remap) <- c("chr19", "start19", "end19", "ID")
setkey(scz.dup.remap,  "ID")
setkey(scz.dup, "ID")

h <- scz.dup[scz.dup.remap]
h$start18 <- h$BP
h$BP <- h$start19
h$chr19 <- NULL
h$end19 <- NULL
h$ID <- NULL



