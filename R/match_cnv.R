## get the joint CNV data set in hg19

library(here)
source(here("R","prepare.R"))

install.load("data.table")
install.load.bioc("GenomicRanges")
install.load.bioc("VariantAnnotation", "AnnotationHub", "TxDb.Hsapiens.UCSC.hg19.knownGene" ,"clusterProfiler", "")
install.load.bioc("ReactomePA")
load(here("output","remapped_cnv2.rda"))

gscz.dup19 <- GRanges(seqnames=scz.dup19$CHR, IRanges(start=scz.dup19$BP, width=1), 
                    z=scz.dup19$z, p=scz.dup19$z_pval, f=scz.dup19$NCNV/41321, 
                    fcontrol =scz.dup19$UNAFF/20227 , n=scz.dup19$NCNV, fcase = scz.dup19$AFF/21094)

gscz.del19 <- GRanges(seqnames=scz.del19$CHR, IRanges(start=scz.del19$BP, width=1), z=scz.del19$z, 
                    p=scz.del19$z_pval, f=scz.del19$NCNV/41321, 
                    fcontrol =scz.del19$UNAFF/20227 , n=scz.del19$NCNV, fcase = scz.del19$AFF/21094)

gheight_cnv19 <- GRanges(height_cnv19$CHR, IRanges(height_cnv19$BP, width=1), 
                       p_height = height_cnv19$`Pvalue Height`, beta_height= height_cnv19$`Beta Height`, 
                       p_bmi = height_cnv19$`Pvalue BMI`, beta_bmi = height_cnv19$`Beta BMI`, 
                       fdel=height_cnv19$F_DEL, fdup = height_cnv19$F_DUP)

nearest.scz <- nearest(gscz.del19, gscz.dup19)
all.scz <- cbind(data.frame(gscz.del19), data.frame(gscz.dup19[nearest.scz]) )[,-c(3:5,12,14:16)]

colnames(all.scz)[9:15] <- paste(colnames(all.scz)[2:8],"dup", sep="_")
colnames(all.scz)[2:8] <- paste(colnames(all.scz)[2:8],"del", sep="_")

## now join with height
all.scz2 <- all.scz[abs(all.scz$start_del-all.scz$start_dup)<=1,]
gall.scz19 <- GRanges(seqnames=all.scz2$seqnames, IRanges(start=all.scz2$start_del, width=1), 
                      z_del=all.scz2$z_del, p_del = all.scz2$p_del, f_del = all.scz2$f_del, fcontrol_del = all.scz2$fcontrol_del,
                      z_dup=all.scz2$z_dup, p_dup = all.scz2$p_dup, f_dup = all.scz2$f_dup, fcontrol_dup = all.scz2$fcontrol_dup)


nearest.h.scz <- nearest(gheight_cnv19, gall.scz19)
three <- cbind(data.frame(gall.scz19[nearest.h.scz]), data.frame(gheight_cnv19) )

three <- three[abs((three[,2] - three[,15]))<=10,]

## Should I add annotation to this?
hub <- AnnotationHub()


txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene
install.load.bioc("VariantAnnotation", "AnnotationHub", "TxDb.Hsapiens.UCSC.hg19.knownGene" ,"clusterProfiler")

FIQT <- function(z=z, min.p=10^-300){
  pvals<-2*pnorm(abs(z), low=F)
  pvals[pvals<min.p]<- min.p
  adj.pvals<-p.adjust(pvals,method="fdr")
  mu.z<-sign(z)*qnorm(adj.pvals/2,low=F)
  mu.z[abs(z)>qnorm(min.p/2,low=F)]<-z[abs(z)>qnorm(min.p/2,low=F)]
  mu.z
}

three$z_del2 <- FIQT(three$z_del, min.p = 1E-12)
three$z_dup2 <- FIQT(three$z_dup, min.p = 1E-12)
three$z_height <-   sign(three$beta_height)*qnorm(three$p_height/2,low=F)
three$z_height2 <- FIQT(three$z_height, min.p=1E-12)
three$z_BMI <-   sign(three$beta_bmi)*qnorm(three$p_bmi/2,low=F)
three$z_BMI2 <- FIQT(three$z_BMI, min.p=1E-12)

gthree <-  GRanges(Rle(paste0("chr", three$seqnames)), IRanges(start=three[,2], width=1), 
                   z_height = three$z_height2,
                   z_bmi = three$z_BMI2,
                   z_del = three$z_del2,
                   z_dup = three$z_dup2)
                   
loc.gthree <- locateVariants(gthree, txdb_hg19, AllVariants())
gthree.loc <- data.table(data.frame(gthree[loc.gthree$QUERYID], GENEID= loc.gthree$GENEID
                                        , LOCATION=loc.gthree$LOCATION, TXID = loc.gthree$TXID))

save(gthree.loc, file=here("output", "gthree.loc.rda"))



