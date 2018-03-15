## CNV annotation

## get the joint CNV data set in

library(here)
source(here("R","prepare.R"))
install.load("data.table")
install.load.bioc("GenomicRanges")
install.load.bioc("VariantAnnotation", "AnnotationHub", "TxDb.Hsapiens.UCSC.hg19.knownGene" ,"clusterProfiler", "")
install.load.bioc("ReactomePA")

hub <- AnnotationHub()


txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene


load(here("output", "remapped_cnv2.rda"))
####
FIQT <- function(z=z, min.p=10^-300){
  pvals<-2*pnorm(abs(z), low=F)
  pvals[pvals<min.p]<- min.p
  adj.pvals<-p.adjust(pvals,method="fdr")
  mu.z<-sign(z)*qnorm(adj.pvals/2,low=F)
  mu.z[abs(z)>qnorm(min.p/2,low=F)]<-z[abs(z)>qnorm(min.p/2,low=F)]
  mu.z
}

scz.del19$z2 <- FIQT(scz.del19$z, min.p = 1E-12)
scz.dup19$z2 <- FIQT(scz.dup19$z, min.p = 1E-12)
height_cnv19$z <-   sign(height_cnv19$`Beta Height`)*qnorm(height_cnv19$`Pvalue Height`/2,low=F)
height_cnv19$z2 <- FIQT(height_cnv19$z, min.p=1E-12)

if (FALSE) {
  library(ggplot2)
  ggplot(scz.del19, aes(x=z, y=z2)) + geom_point()
  ggplot(scz.dup19, aes(x=z, y=z2)) + geom_point()
  ggplot(height_cnv19, aes(x=z, y=z2)) + geom_point()
}

#---------------------------------------------
# Now make these into genomic ramges
#---------------------------------------------

gscz.del19 <-  GRanges(Rle(paste0("chr", scz.del19$CHR)), IRanges(start=scz.del19$start19, width=1), 
                       pval = scz.del19$z_pval, or = scz.del19$cmh_OR, z = scz.del19$z2)
loc_gscz.del19 <- locateVariants(gscz.del19, txdb_hg19, AllVariants())
gscz.del19.loc <- data.table(data.frame(gscz.del19[loc_gscz.del19$QUERYID], GENEID= loc_gscz.del19$GENEID, LOCATION=loc_gscz.del19$LOCATION, TXID = loc_gscz.del19$TXID))

gscz.dup19 <-  GRanges(Rle(paste0("chr", scz.dup19$CHR)), IRanges(start=scz.dup19$start19, width=1)
                       , pval = scz.dup19$z_pval, or = scz.dup19$cmh_OR,, z = scz.dup19$z2)
loc_gscz.dup19 <- locateVariants(gscz.dup19, txdb_hg19, AllVariants())
gscz.dup19.loc <- data.table(data.frame(gscz.dup19[loc_gscz.dup19$QUERYID], GENEID= loc_gscz.dup19$GENEID, LOCATION=loc_gscz.dup19$LOCATION, TXID = loc_gscz.dup19$TXID))

gheight19 <-  GRanges(Rle(paste0("chr", height_cnv19$CHR)), IRanges(start=height_cnv19$start19, width=1)
                      , p.height = height_cnv19$`Pvalue Height`, b = height_cnv19$`Beta Height`, z=height_cnv19$z2)
loc_gheight19 <- locateVariants(gheight19, txdb_hg19, AllVariants())
gheight19.loc <- data.table(data.frame(gheight19[loc_gheight19$QUERYID], GENEID= loc_gheight19$GENEID, LOCATION=loc_gheight19$LOCATION, TXID = loc_gheight19$TXID))

## Height genes with z > 3

height_path <- enrichPathway(unique(gheight19.loc$GENEID[abs(gheight19.loc$z)>2]),  pvalueCutoff=0.05, readable=T)
height_path
barplot(height_path)

height_path_down <- enrichPathway(unique(gheight19.loc$GENEID[gheight19.loc$z < -2 ]),  pvalueCutoff=0.05, readable=T)
height_path_down


## scz del genes

scz_path <- enrichPathway(unique(gscz.del19.loc$GENEID[abs(gscz.del19.loc$z)>2.5]),  pvalueCutoff=0.05, readable=T)
scz_path
barplot(height_path)

scz_del_path_up <- enrichPathway(unique(gscz.del19.loc$GENEID[gscz.del19.loc$z>3]),  pvalueCutoff=0.05, readable=T)
scz_del_path_up
scz_del_path_down <-  enrichPathway(unique(gscz.del19.loc$GENEID[gscz.del19.loc$z < -2]),  pvalueCutoff=0.05, readable=T)
scz_del_path_down

scz_dup_path_down <-  enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$z < -2]),  pvalueCutoff=0.05, readable=T)
scz_dup_path_down

scz_dup_path_up <-  enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$z > 3]),  pvalueCutoff=0.05, readable=T)
scz_dup_path_up



scz_del_path_down <- enrichPathway(unique(gscz.del19.loc$GENEID[gscz.del19.loc$pval<0.001 & gscz.del19.loc$or <1 ]),  pvalueCutoff=0.05, readable=T)
scz_del_path_down
plot(scz_del_path)

scz_dup_path <- enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001]),  pvalueCutoff=0.05, readable=T)
scz_dup_path
barplot(scz_dup_path )

scz_dup_path_up <- enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001 & gscz.dup19.loc$or > 1]),  pvalueCutoff=0.05, readable=T)
scz_dup_path_up

scz_dup_path_down <- enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001 & gscz.dup19.loc$or <1]),  pvalueCutoff=0.05, readable=T)
scz_dup_path_down
barplot(scz_dup_path_down )

scz_dup_path <- enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001]),  pvalueCutoff=0.05, readable=T)
scz_dup_path
barplot(scz_dup_path )

scz_both_path <- enrichPathway(unique(c(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001],gscz.del19.loc$GENEID[gscz.del19.loc$pval<0.01])) ,  pvalueCutoff=0.05, readable=T)
scz_both_path
plot(scz_both_path)

lh <- list(dup = unique(gscz.dup19.loc$GENEID[abs(gscz.dup19.loc$z) > 2]), 
           height = unique(gheight19.loc$GENEID[abs(gheight19.loc$z) > 2]),
           del = unique(gscz.del19.loc$GENEID[abs(gscz.del19.loc$z) > 2])
)

res_h <- compareCluster(lh, fun="enrichPathway")
plot(res_h)  


## Now I need to do something with matched CNVs to see if that is better
## Find the matched CNV dataset






