## liftover
## scripts and code to convert cnv files from one to another coordincates

## Height
library(here)
source(here("R","prepare.R"))
install.load("data.table")
install.load.bioc("GenomicRanges") #, "rtracklayer")

height <- fread(dropbox("height_CNV_association_41467_2017_556_MOESM2_ESM.csv"))   ## standard CNV
height$ID <-  paste("ID", height$CHR, height$BP, sep="_")
nrow(height)

write.table(cbind(paste("chr", height$CHR, sep=""), height$BP, height$BP+1, height$ID),
            col.names = FALSE, row.names=FALSE, quote=FALSE, file==here("output", "height_cnv_hg18.bed") )

## liftOver cnv_hg18.bed hg18ToHg19.over.chain.gz cnv_hg19.bed unmapped

## Read the non-gene centric results.

scz.dup <- fread(dropbox("pgc_cnv/PGC_41K_QC_dup.cnv.results"))
scz.del <- fread(dropbox("pgc_cnv/PGC_41K_QC_del.cnv.results"))
scz.dup$ID <-  paste("ID", scz.dup$CHR, scz.dup$BP, sep="_")
scz.del$ID <-  paste("ID", scz.del$CHR, scz.del$BP, sep="_")


write.table(cbind(paste("chr",scz.dup$CHR,sep=""), scz.dup$BP, scz.dup$BP+1, scz.dup$ID),
            col.names = FALSE, row.names=FALSE, quote=FALSE, file=here("output","scz.dup_cnv_hg18.bed") )

write.table(cbind(paste("chr",scz.del$CHR,sep=""), scz.del$BP, scz.del$BP+1, scz.del$ID),
            col.names = FALSE, row.names=FALSE, quote=FALSE, file=here("output","scz.del_cnv_hg18.bed"))

###############################################################
## Gene centric scz analysis


scz.del.gene <- fread(dropbox("pgc_cnv/PGC_41K_QC_del_minimum8cnv.gene.results"))
scz.del.gene <- scz.del.gene[!is.na(scz.del.gene$BP_start_hg18),]
scz.all.gene <- fread(dropbox("pgc_cnv/PGC_41K_QC_all_minimum12cnv.gene.results"))
scz.all.gene <- scz.all.gene[!is.na(scz.all.gene$BP_start_hg18),]


write.table(cbind(paste("chr",scz.all.gene$CHR, sep=""), scz.all.gene$BP_start_hg18, scz.all.gene$BP_start_hg18+1, scz.all.gene$Gene_symbol),
            col.names = FALSE, row.names=FALSE, quote=FALSE, file=here("output","scz.all_cnv_hg18.gene.bed") )

write.table(cbind(paste("chr",scz.del.gene$CHR, sep=""), scz.del.gene$BP, scz.del.gene$BP+1, scz.del.gene$Gene_symbol),
            col.names = FALSE, row.names=FALSE, quote=FALSE, file=here("output","scz.del_cnv_hg18.gene.bed"))






if (FALSE) {
command <- "
  liftOver height_cnv_hg18.bed hg18ToHg19.over.chain.gz height_cnv_hg19.bed unmapped_height
  liftOver scz.dup_cnv_hg18.bed hg18ToHg19.over.chain.gz scz.dup_cnv_hg19.bed unmapped_scz.dup
  liftOver scz.del_cnv_hg18.bed hg18ToHg19.over.chain.gz scz.del_cnv_hg19.bed unmapped_scz.del

  liftOver scz.all_cnv_hg18.gene.bed hg18ToHg19.over.chain.gz scz.all_cnv_hg19.gene.bed unmapped_scz_all.dup
  liftOver scz.del_cnv_hg18.gene.bed hg18ToHg19.over.chain.gz scz.del_cnv_hg19.gene.bed unmapped_scz_del.dup
"
}

