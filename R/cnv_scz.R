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


if (FALSE) {
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
}

### Intolerance

## look up the gene cnv tests
if (FALSE) {
  
  scz_gene <- fread(dropbox("pgc_cnv/PGC_41K_QC_all_minimum12cnv.gene.results"))
colnames(scz_gene)[1] <- "gene_symbol"

cnv_intolerance <- read.table(here("data", "Intolerance", "exac-final-cnv.gene.scores071316"), header=TRUE)

table(scz_gene$Gene_symbol %in% cnv_intolerance$gene_symbol)
nrow(cnv_intolerance)
w <- which(!(scz_gene$Gene_symbol %in% cnv_intolerance$gene_symbol))
scz_gene[w,]

a <- merge(scz_gene, cnv_intolerance, by="gene_symbol")

cor(cbind(a$glm_beta, a$del.score, a$dup.score, a$cnv.score, a$glm_pval))


#############################################################################################
## SCZ deletions

scz_del <- fread(dropbox("pgc_cnv/PGC_41K_QC_del_minimum8cnv.gene.results"))
colnames(scz_del)[1] <- "gene_symbol"

a <- merge(scz_del, cnv_intolerance, by="gene_symbol")

cor(cbind(a$AFF, a$UNAFF, a$glm_beta, a$glm_pval, a$dev_estimate, a$del.score, a$dup.score, a$cnv.score))


}