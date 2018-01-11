## Height
library(here)
source(here("R","prepare.R"))
install.load("data.table")
install.load.bioc("GenomicRanges")

scz <- fread(dropbox("pgc_cnv/PGC_41K_QC_all.cnv.results"))

gscz <- GRanges(seqnames=scz$CHR, IRanges(start=scz$BP, end=scz$BP+1)
                   ,  pscz=scz$z_pval )


if (!dir.exists(here("output")))
  dir.create(here("output"))
save(gscz, file=here("output", "scz_cnv.rda"))



install.load("qqman")

######################################################
opar <- par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(1.2,0.5,0))


scz$P <- scz$z_pval
manhattan(scz, main="CNV SCZ risk - Z test")
scz$P <- scz$exact_pval
manhattan(scz, main="CNV SCZ risk - exact test")



par(opar)

qq(scz$z_pval)

######################################################

scz <- scz[scz$NCNV>=20]
opar <- par(mfrow=c(2,1), mar=c(3,3,1,1), mgp=c(1.2,0.5,0))


scz$P <- scz$z_pval
manhattan(scz, main="CNV SCZ risk - Z test")
scz$P <- scz$exact_pval
manhattan(scz, main="CNV SCZ risk - exact test")



par(opar)
qq(scz$exact_pval)
