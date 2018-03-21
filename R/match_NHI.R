#' ## match Neuroticism, Height and Intelligence

library(here)
source(here("R","prepare.R"))
install.load("data.table")

## neuroticism summary results from my own computer
## Neuroticism paper is
## Association analysis in over 329,000 individuals identifies 116 independent variants influencing neuroticism
## Luciano et al. 2018.  Nature Geneticsvolume 50, pages6–11 (2018) doi:10.1038/s41588-017-0013-8
## UK biobank data

neuroticism <- fread("file:///D:/data/GWAS_Summary/Luciano_2017/SummaryStats.txt")
head(neuroticism)
colnames(neuroticism)[2] <- "snpid"
setkey(neuroticism, "snpid")

## This Height data is  
## Wood AR, Esko T, Yang J, Vedantam S, Pers TH, Gustafsson S et al. (2014). 
## Defining the role of common variation in the genomic and biological architecture 
## of adult human height. Nature Genetics 11:1173-86. 

GWAS_height_summary <- data.table(read.table(
  dropbox("GIANT_HEIGHT_Wood_et_al_2014_publicrelease_HapMapCeuFreq.txt.gz"), 
  header=TRUE, stringsAsFactors = FALSE))
head(GWAS_height_summary)
colnames(GWAS_height_summary)[1] <- "snpid"
setkey(GWAS_height_summary, "snpid")

#' Now join the first two before bringing another data set in

NH <- merge(neuroticism, GWAS_height_summary)

switch_snps <- (NH$a_1==NH$Allele1 & NH$a_0==NH$Allele2)
NH$N_res_beta[switch_snps] <- -NH$N_res_beta[switch_snps]
#` throw anything else away`
keep <- (NH$a_0==NH$Allele1 & NH$a_1==NH$Allele2) | (NH$a_0==NH$Allele2 & NH$a_1==NH$Allele1)
NH <- NH[keep]
#` and throw away any unneeded columns
NH[,a_0 := NULL]
NH[,a_1 := NULL]
NH[,N := NULL]
rm(GWAS_height_summary, neuroticism)
gc()




#' Genome-wide association meta-analysis of 78,308 individuals identifies new loci and genes influencing
#' human intelligence
#' Sniekers et al Nature Genetics volume 49, pages 1107–1112 (2017) doi:10.1038/ng.3869
#' https://www.nature.com/articles/ng.3869
#' 

GWAS_intelligence_summary <- fread("D:/data/GWAS_Summary/Intelligence/Sniekers.txt", stringsAsFactors = FALSE) 
head(GWAS_intelligence_summary)
colnames(GWAS_intelligence_summary)[2] <- "snpid"
setkey(GWAS_intelligence_summary, "snpid")


NHI <- NH[GWAS_intelligence_summary]





#' OK, so I've managed to merge this but If I look at the table of allele
#' names then it is clear that I need to swiotch some of the alleles around

table(snp_two$a_0, snp_two$a1)

#` This shoudl be fine by just changing the sign on the beta.  It should make no 
#` difference to the standard error or the p-value

switch_snps <- (snp_two$a_1==snp_two$a1 & snp_two$a_0==snp_two$a2)
snp_two$N_res_beta[switch_snps] <- -snp_two$N_res_beta[switch_snps]
#` throw anything else away`
keep <- (snp_two$a_1==snp_two$a1 & snp_two$a_0==snp_two$a2) | (snp_two$a_0==snp_two$a1 & snp_two$a_1==snp_two$a2)
snp_two <- snp_two[keep]
#` and throw away any unneeded columns
snp_two[,a_0 := NULL]
snp_two[,a_1 := NULL]
snp_two[,hg19chrc := NULL]
snp_two[,ngt := NULL]
rm(GWAS_scz_summary,neuroticism)
gc()

colnames(snp_two) <- c("snpid", "chrom", "pos","beta.N", "se.N", "p.N", "a1","a2", "bp","info", "or", "or.se", "p.scz")

setkey(snp_two, "snpid")
snp_three <- merge(snp_two, GWAS_height_summary)
rm(GWAS_height_summary)
gc()

snp_three
table(snp_three$a1==snp_three$Allele1)
#` these match much better so just throw mismatches away#
snp_three <- snp_three[snp_three$a1==snp_three$Allele1]
snp_three <- snp_three[snp_three$a2==snp_three$Allele2]
snp_three[,Allele1 := NULL]
snp_three[,Allele2 := NULL]
snp_three[,N := NULL]
snp_three[,bp := NULL]
colnames(snp_three) <- c("snpid", "chrom",  "pos", "beta.N","se.N", "p.N", "a1", "a2","info", "or", "or.se", "p.scz", "freq"
                         , "b.height", "se.height", "p.height")

rm(snp_two)

install.load.bioc("VariantAnnotation", "AnnotationHub", "TxDb.Hsapiens.UCSC.hg19.knownGene" )

snp_three$z.scz <- log(snp_three$or)/snp_three$or.se
snp_three$z.height <- snp_three$b.height/snp_three$se.height
snp_three$z.neu <- snp_three$beta.N/snp_three$se.N

snp_three$zb.scz <- FIQT(snp_three$z.scz)
snp_three$zb.height <- FIQT(snp_three$z.height)
snp_three$zb.neu <- FIQT(snp_three$z.neu)

snp_three$chrom <- paste0("chr", snp_three$chrom)

gthree <-  makeGRangesFromDataFrame(snp_three,
                                   keep.extra.columns=TRUE,
                                   ignore.strand=TRUE,
                                   seqnames.field=c("chrom"),
                                   start.field="pos",
                                   end.field=c("pos"))

seqnames(gthree)
rm(snp_three)
  


hub <- AnnotationHub()


txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

loc_three <- locateVariants(gthree, txdb_hg19, AllVariants())

newid <- paste(seqnames(loc_three), start(loc_three), loc_three$LOCATION, loc_three$GENEID, sep="_")
length(newid)
length(unique(newid))
loc_three_nodups <- loc_three[!duplicated(newid)]
table(loc_three$LOCATION)
table(loc_three_nodups$LOCATION)
#loc_both_hg19 <- loc_height_hg19
snp_shn <- data.table(data.frame(gthree[loc_three_nodups$QUERYID]
                                        , GENEID= loc_three_nodups$GENEID, 
                                        LOCATION=loc_three_nodups$LOCATION, 
                                        TXID = loc_three_nodups$TXID, 
                                        stringsAsFactors = FALSE))


save(snp_shn, file=here("output", "snps_SHN.rda"))


