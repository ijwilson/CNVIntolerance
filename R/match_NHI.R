#' ## match Neuroticism, Height and Intelligence

library(here)
source(here("R","prepare.R"))
install.load("data.table")

## neuroticism summary results from my own computer
## Neuroticism paper is
## Association analysis in over 329,000 individuals identifies 116 independent variants influencing neuroticism
## Luciano et al. 2018.  Nature Genetics volume 50, pages 6–11 (2018) doi:10.1038/s41588-017-0013-8
## UK biobank data

if (!file.exists(here("data","Luciano_2017","SummaryStats.txt"))) {
  if (!dir.exists(here("data")))
    dir.create(here("data"))
  download.file("http://www.psy.ed.ac.uk/ccace/downloads/Luciano_2017.zip", here("data", "luciano_2017.zip"))
  unzip(here("data", "luciano_2017.zip"), exdir=here("data"))
  file.remove(here("data", "luciano_2017.zip"))
  download.file("http://www.psy.ed.ac.uk/ccace/downloads/Hagenaars2017_UKB_MPB_summary_results.zip", here("data", "baldness.zip"))
  unzip(here("data", "baldness.zip"), exdir=here("data"))
  file.remove(here("data", "baldness.zip"))
}

neuroticism <- fread(here("data","Luciano_2017","SummaryStats.txt"))
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
#' 
 if (!file.exists(here("data", "sniekers_2017.txt.gz"))) {
   if (!dir.exists(here("data"))) dir.create(here("data"))
   download.file("http://ctg.cncr.nl/documents/p1651/sumstats.txt.gz", here("data", "sniekers_2017.txt.gz"))
 }
   
   
GWAS_intelligence_summary <- data.table(
  read.table(here("data", "sniekers_2017.txt.gz"), header=TRUE,  stringsAsFactors = FALSE))

#GWAS_intelligence_summary <- fread("D:/data/GWAS_Summary/Intelligence/Sniekers.txt", stringsAsFactors = FALSE) 
head(GWAS_intelligence_summary)
colnames(GWAS_intelligence_summary)[3] <- "snpid"

setkey(GWAS_intelligence_summary, "snpid")


NHI <- GWAS_intelligence_summary[NH]

#switch_snps <- (NHI$Allele1==NHI$alt & NHI$Allele2==NHI$ref)
#snp_two$N_res_beta[switch_snps] <- -snp_two$N_res_beta[switch_snps]
#` throw anything else away`
keep <- (NHI$Allele1==NHI$ref & NHI$Allele2==NHI$alt)
NHI <- NHI[keep]

#' OK, so I've managed to merge this but If I look at the table of allele
#' names then it is clear that I need to swiotch some of the alleles around


#` and throw away any unneeded columns
NHI[,Allele1 := NULL]
NHI[,Allele2 := NULL]
NHI[, BP := NULL]
NHI[,CHR:=NULL]

colnames(NHI) <- 
c("chrom","pos","snpid","ref","alt","MAF","beta.int","se.int","z.int","p.int","direction","beta.n","se.n","p.n"
  ,"f","beta.h","se.h","p.h")


rm(NH, GWAS_intelligence_summary)

install.load.bioc("VariantAnnotation", "AnnotationHub", "TxDb.Hsapiens.UCSC.hg19.knownGene" )

NHI[,z.h := beta.h/se.h]
NHI[,z.n := beta.n/se.n]

NHI$zb.h <- FIQT(NHI$z.h)
NHI$zb.n <- FIQT(NHI$z.n)
NHI$zb.i <- FIQT(NHI$z.int)

NHI$chrom <- paste0("chr", NHI$chrom)

gNHI <-  makeGRangesFromDataFrame(NHI,
                                   keep.extra.columns=TRUE,
                                   ignore.strand=TRUE,
                                   seqnames.field=c("chrom"),
                                   start.field="pos",
                                   end.field=c("pos"))

rm(NHI)
  


hub <- AnnotationHub()


txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

loc_NHI <- locateVariants(gNHI, txdb_hg19, AllVariants())

newid <- paste(seqnames(loc_NHI), start(loc_NHI), loc_NHI$LOCATION, loc_NHI$GENEID, sep="_")
length(newid)
length(unique(newid))
loc_NHI_nodups <- loc_NHI[!duplicated(newid)]
table(loc_NHI$LOCATION)
table(loc_NHI_nodups$LOCATION)
#loc_both_hg19 <- loc_height_hg19
snp_NHI <- data.table(data.frame(gNHI[loc_NHI_nodups$QUERYID]
                                        , GENEID= loc_NHI_nodups$GENEID, 
                                        LOCATION=loc_NHI_nodups$LOCATION, 
                                        TXID = loc_NHI_nodups$TXID, 
                                        stringsAsFactors = FALSE))


save(snp_NHI, file=here("output", "snp_NHI.rda"))



