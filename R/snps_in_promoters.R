#' ---
#' title: "SNPs in Promoters"
#' author: "Ian Wilson"
#' ---

#' Get a list of SNPs that are in gene promoters and then see what they look like

#'## Preparation
#'
do.plots=FALSE
sig_level <- 1E-6
#+ message=FALSE 
library(here)
source(here("R","prepare.R"))
install.load.bioc("GenomicRanges")
load(here("output", "genesGR.rda"))

#' get the promoter regions
promoters19 <- promoters(genesGR19)

#' Load Height GWAS results
#' 
load(here("output", "GWAS_height_summary.rda"))
#' 
#' Turn height into a Genomic Range
#' 
gheight <- GRanges(seqnames=Rle(GWAS_height_summary$seqnames), IRanges(start=GWAS_height_summary$pos, width=1, names=GWAS_height_summary$MarkerName)
                   , p = GWAS_height_summary$p, b=GWAS_height_summary$b, f=GWAS_height_summary$Freq.Allele1.HapMapCEU)
sig_height <- sum(GWAS_height_summary$p < sig_level)
n_height <- nrow(GWAS_height_summary)
#' Subset by those SNPs in promoter regions
promoterHeightSNPs <- subsetByOverlaps(gheight, promoters19)
geneHeightSNPs <-  subsetByOverlaps(gheight, genesGR19)
#' We have `r length(promoterHeightSNPs)` SNPs in promoter regions`
#' 

m_height <- matrix(c(sum(promoterHeightSNPs$p<sig_level), length(promoterHeightSNPs),
                     sum(geneHeightSNPs$p<sig_level), length(geneHeightSNPs),
                         sig_height, n_height)
                   , nrow=3, byrow=TRUE)

m_height[,1]/m_height[,2]


rm(GWAS_height_summary, gheight, genesGR18) # tidy
gc()
#' Load SCZ GWAS results
load(here("output", "GWAS_scz_summary.rda"))
gscz <- GRanges(seqnames=Rle(GWAS_scz_summary$chr), IRanges(start=GWAS_scz_summary$bp, width=1, names=GWAS_scz_summary$snpid)
                , p = GWAS_scz_summary$p, or=GWAS_scz_summary$or)
sig_scz <- sum(GWAS_scz_summary$p < sig_level)
n_scz <- nrow(GWAS_scz_summary)
#' Subset by those SNPs in promoter regions
geneSCZSNPs <-  subsetByOverlaps(gscz, genesGR19)
promoterSCZSNPs <- subsetByOverlaps(gscz, promoters19)
#' We have `r length(promoterSCZSNPs)` SNPs in promoter regions for SCZ`
#' 
m_scz <- matrix(c(sum(promoterSCZSNPs$p<sig_level), length(promoterSCZSNPs),
                     sum(geneSCZSNPs$p<sig_level), length(geneSCZSNPs),
                     sig_scz, n_scz)
                   , nrow=3, byrow=TRUE)

m_scz[,1]/m_scz[,2]



rm(GWAS_scz_summary, gscz)
gc()


if (do.plots==TRUE) {
  install.load("qqman")
  qq(promoterHeightSNPS$p)
  df <- as.data.frame(promoterHeightSNPs)
  df$idnames <- rownames(df)
  df$seqnames <- as.numeric(paste(df$seqnames))
  manhattan(df, chr="seqnames", bp="start", p="p", snp="idnames")

  qq(promoterSCZSNPs$p)
  df <- as.data.frame(promoterSCZSNPs)
  df$idnames <- rownames(df)
  df$seqnames <- as.numeric(paste(df$seqnames))
  manhattan(df, chr="seqnames", bp="start", p="p", snp="idnames")
  
  rm(df)
}

#' look at the CNV in genes

#' start with SCZ
o <- findOverlaps(promoterSCZSNPs, promoters19)

pp <- promoterSCZSNPs[queryHits(o)]
pp$Gene_symbol <- names(promoters19)[subjectHits(o)]


#' read SCZ cnv results

scz.del.gene <- read.table(dropbox("pgc_cnv/PGC_41K_QC_del_minimum8cnv.gene.results"), header=TRUE)
scz.all.gene <- read.table(dropbox("pgc_cnv/PGC_41K_QC_all_minimum12cnv.gene.results"), header=TRUE)

m <- match( scz.del.gene$Gene_symbol, pp$Gene_symbol)
#' get matching genes

scz.del.geneB <- scz.del.gene[!is.na(m),]
m <- m[!is.na(m)]
alld <- data.frame(scz.del.geneB, pp[m])
m <- match( scz.all.gene$Gene_symbol, pp$Gene_symbol)
#' get matching genes

scz.all.geneB <- scz.all.gene[!is.na(m),]
m <- m[!is.na(m)]
alld <- data.frame(scz.all.geneB, pp[m])

plot(alld$glm_beta, alld$or)
cor(alld$glm_beta, alld$or)


#'
#' Height
#' 
#' Label the promoter SNPS with the gene names 

o <- findOverlaps(promoterHeightSNPs, promoters19)   ## find overlaps between SNPs and promoters
pHeight <- promoterHeightSNPs[queryHits(o)]
pHeight$Gene_symbol <- names(promoters19)[subjectHits(o)]

#' read CNVs
height_cnv19 <- read.csv(here("data", "height_cnv_by_gene_hg19.csv"))
colnames(height_cnv19)[1] <- "Gene_symbol"

m <- match( height_cnv19$Gene_symbol, pHeight$Gene_symbol)
#' get matching genes

height_cnv19 <- height_cnv19[!is.na(m),]
m <- m[!is.na(m)]
alld <- data.frame(height_cnv19, pHeight[m])

plot(alld$mean_beta_height, alld$b)
plot(alld$mean_beta_height, alld$b)
