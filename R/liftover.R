## liftover
## scripts and code to convert cnv files from one to another coordincates

## Height
library(here)
source(here("R","prepare.R"))
install.load("data.table")
install.load.bioc("GenomicRanges") #, "rtracklayer")

height <- fread(dropbox("height_CNV_association_41467_2017_556_MOESM2_ESM.csv"))   ## standard CNV
height$ID <-  paste("ID", height$CHR, height$BP, sep="_")

## Read the non-gene centric results.
scz.dup <- fread(dropbox("pgc_cnv/PGC_41K_QC_dup.cnv.results"))
scz.del <- fread(dropbox("pgc_cnv/PGC_41K_QC_del.cnv.results"))
scz.dup$ID <-  paste("ID", scz.dup$CHR, scz.dup$BP, sep="_")
scz.del$ID <-  paste("ID", scz.del$CHR, scz.del$BP, sep="_")

###############################################################
## Gene centric scz analysis
scz.del.gene <- fread(dropbox("pgc_cnv/PGC_41K_QC_del_minimum8cnv.gene.results"))
scz.del.gene <- scz.del.gene[!is.na(scz.del.gene$BP_start_hg18),]
scz.all.gene <- fread(dropbox("pgc_cnv/PGC_41K_QC_all_minimum12cnv.gene.results"))
scz.all.gene <- scz.all.gene[!is.na(scz.all.gene$BP_start_hg18),]

if (FALSE) {
command <- "
  liftOver height_cnv_hg18.bed hg18ToHg19.over.chain.gz height_cnv_hg19.bed unmapped_height
  liftOver scz.dup_cnv_hg18.bed hg18ToHg19.over.chain.gz scz.dup_cnv_hg19.bed unmapped_scz.dup
  liftOver scz.del_cnv_hg18.bed hg18ToHg19.over.chain.gz scz.del_cnv_hg19.bed unmapped_scz.del

  liftOver scz.all_cnv_hg18.gene.bed hg18ToHg19.over.chain.gz scz.all_cnv_hg19.gene.bed unmapped_scz_all.dup
  liftOver scz.del_cnv_hg18.gene.bed hg18ToHg19.over.chain.gz scz.del_cnv_hg19.gene.bed unmapped_scz_del.dup
"
}

## now need to remap to the original data

height.remap <- fread(dropbox("height_cnv_hg19.bed"))
colnames(height.remap) <- c("chr19", "start19", "end19", "ID")
setkey(height.remap,  "ID")
setkey(height, "ID")

height_cnv19 <- height[height.remap]
height_cnv19$start18 <- height_cnv19$BP
height_cnv19$BP <- height_cnv19$start19
height_cnv19$chr19 <- NULL
height_cnv19$end19 <- NULL
height_cnv19$ID <- NULL

#=======================
scz.dup.remap <- fread(dropbox("scz.dup_cnv_hg19.bed"))
colnames(scz.dup.remap) <- c("chr19", "start19", "end19", "ID")
setkey(scz.dup.remap,  "ID")
setkey(scz.dup, "ID")

scz.dup19 <- scz.dup[scz.dup.remap]
scz.dup19$start18 <- scz.dup19$BP
scz.dup19$BP <- scz.dup19$start19
scz.dup19$chr19 <- NULL
scz.dup19$end19 <- NULL
scz.dup19$ID <- NULL
#=======================
scz.del.remap <- fread(dropbox("scz.del_cnv_hg19.bed"))
colnames(scz.del.remap) <- c("chr19", "start19", "end19", "ID")
setkey(scz.del.remap,  "ID")
setkey(scz.del, "ID")
scz.del19 <- scz.del[scz.del.remap]
scz.del19$start18 <- scz.del19$BP
scz.del19$BP <- scz.del19$start19
scz.del19$chr19 <- NULL
scz.del19$end19 <- NULL
scz.del19$ID <- NULL
#========================================
scz.del.gene.remap <- fread(dropbox("scz.del_cnv_hg19.gene.bed"))
colnames(scz.del.gene.remap) <- c("chr19", "start19", "end19", "Gene_symbol")
setkey(scz.del.gene.remap,  "Gene_symbol")
setkey(scz.del.gene, "Gene_symbol")
scz.del.gene19 <- scz.del.gene[scz.del.gene.remap]
#========================================
scz.all.gene.remap <- fread(dropbox("scz.all_cnv_hg19.gene.bed"))
colnames(scz.all.gene.remap) <- c("chr19", "start19", "end19", "Gene_symbol")
setkey(scz.all.gene.remap,  "Gene_symbol")
setkey(scz.all.gene, "Gene_symbol")
scz.all.gene19 <- scz.all.gene[scz.all.gene.remap]
#========================================
save(scz.del19, scz.dup19, height_cnv19, scz.all.gene19, scz.del.gene19, file=here("output", "remapped_cnv2.rda"))
