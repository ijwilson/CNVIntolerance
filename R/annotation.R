

library(here)
source(here("R","prepare.R"))
install.load.bioc("VariantAnnotation", "AnnotationHub")


# Load height and 
load(here("output", "GWAS_height_summary.rda"))
load(here("output", "GWAS_scz_summary.rda"))

gheight <- GRanges(seqnames=Rle(paste0("chr", GWAS_height_summary$seqnames)), IRanges(start=GWAS_height_summary$pos, width=1, names=GWAS_height_summary$MarkerName)
                   , p = GWAS_height_summary$p, b=GWAS_height_summary$b, f=GWAS_height_summary$Freq.Allele1.HapMapCEU)
genome(gheight) <- "hg19"
gscz <- GRanges(seqnames=Rle(paste0("chr",GWAS_scz_summary$chr)), IRanges(start=GWAS_scz_summary$bp, width=1, names=GWAS_scz_summary$snpid)
                , p = GWAS_scz_summary$p, or=GWAS_scz_summary$or)
genome(gscz) <- "hg19"

rm(GWAS_height_summary, GWAS_scz_summary)
gc()
load(here("output", "GWAS_scz_summary.rda"))

hub <- AnnotationHub()


install.load.bioc("TxDb.Hsapiens.UCSC.hg19.knownGene")
txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

head(seqlevels(txdb_hg19))


loc_height_hg19 <- locateVariants(gheight, txdb_hg19, AllVariants())
loc_scz_hg19 <- locateVariants(gscz, txdb_hg19, AllVariants())

install.load.bioc("org.Hs.eg.db")
