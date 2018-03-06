## prepare figures for presentation

library(here)

if (FALSE) {
################################################################################
## Manhattan CNV height plot 
source(here("R","prepare.R"))
install.load(c("data.table", "qqman"))
height <- fread(dropbox("height_CNV_association_41467_2017_556_MOESM2_ESM.csv"))
height[,SNP := paste("cnv",1:nrow(height),sep="")]
height$P <- height$`Pvalue Height`
png(file=here("figs","height_cnv_manhattan.png"),height=800, width=1200)
manhattan(height, main="Height")
dev.off()
################################################################################
## Manhattan CNV BMI plot 
height$P <- height$`Pvalue BMI`
png(file=here("figs","BMI_cnv_manhattan.png"),height=800, width=1200)
manhattan(height, main="BMI")
dev.off()

################################################################################
### Both

png(file=here("Documents","height_bmi_cnv_manhattan.png"),height=800, width=1400)
opar <- par(mfrow=c(2,1))
height[,SNP := paste("cnv",1:nrow(height),sep="")]
height$P <- height$`Pvalue Height`
manhattan(height, main="Height")
height$P <- height$`Pvalue BMI`
manhattan(height, main="BMI")
par(opar)
dev.off()

}
## Manhattan CNV BMI plot 
## Read the non-gene centric results.
scz.dup <- fread(dropbox("pgc_cnv/PGC_41K_QC_dup.cnv.results"))
scz.del <- fread(dropbox("pgc_cnv/PGC_41K_QC_del.cnv.results"))

png(file=here("Documents","scz_cnv_manhattan.png"),height=800, width=1400)

opar <- par(mfrow=c(2,1), mar=c(3,3,2,1))

scz.dup$P <- scz.dup$z_pval
manhattan(scz.dup, main="SCZ Duplications")

scz.del$P <- scz.del$z_pval
manhattan(scz.del, main="SCZ Deletions")

par(opar)
dev.off()



