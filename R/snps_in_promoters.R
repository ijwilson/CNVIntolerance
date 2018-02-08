#' ---
#' title: "SNPs in Promoters"
#' author: "Ian Wilson"
#' ---

#' Get a list of SNPs that are in gene promoters and then see what they look like

#'## Data Preparation

library(here)
source(here("R","prepare.R"))
install.load.bioc("GenomicRanges")
load(here("output", "GWAS_height_summary.rda"))
load(here("output", "genesGR.rda"))

gheight <- GRanges(seqnames=Rle(GWAS_height_summary$seqnames), IRanges(start=GWAS_height_summary$pos, width=1, names=GWAS_height_summary$MarkerName)
                   , p = GWAS_height_summary$p, b=GWAS_height_summary$b, f=GWAS_height_summary$Freq.Allele1.HapMapCEU)

promoters19 <- promoters(genesGR19)
promoterHeightSNPS <- subsetByOverlaps(gheight, promoters19)
length(promoterHeightSNPS)
rm(GWAS_height_summary, gheight, genesGR18)
gc()

install.load("qqman")

qq(promoterHeightSNPS$p)

df <- as.data.frame(promoterHeightSNPS)
df$idnames <- rownames(df)
df$seqnames <- as.numeric(paste(df$seqnames))
manhattan(df, chr="seqnames", bp="start", p="p", snp="idnames")



load(here("output", "GWAS_scz_summary.rda"))
gscz <- GRanges(seqnames=Rle(GWAS_scz_summary$chr), IRanges(start=GWAS_scz_summary$bp, width=1, names=GWAS_scz_summary$snpid)
                   , p = GWAS_scz_summary$p, or=GWAS_scz_summary$or)
promoterSCZSNPS <- subsetByOverlaps(gscz, promoters19)
length(promoterSCZSNPS)

qq(promoterSCZSNPS$p)

df <- as.data.frame(promoterSCZSNPS)
df$idnames <- rownames(df)
df$seqnames <- as.numeric(paste(df$seqnames))
manhattan(df, chr="seqnames", bp="start", p="p", snp="idnames")


#' look at the CNVs


