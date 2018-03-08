
library(here)
source(here("R","prepare.R"))
install.load.bioc("VariantAnnotation", "AnnotationHub", "TxDb.Hsapiens.UCSC.hg19.knownGene")


if (!file.exists( here("output", "gwas_gr.rda"))) {
# Load height and scz
  load(here("output", "GWAS_height_summary.rda"))
  load(here("output", "GWAS_scz_summary.rda"))
  
  gheight <- GRanges(seqnames=Rle(paste0("chr", GWAS_height_summary$seqnames)), IRanges(start=GWAS_height_summary$pos, width=1, names=GWAS_height_summary$MarkerName)
                     , p = GWAS_height_summary$p, b=GWAS_height_summary$b, f=GWAS_height_summary$Freq.Allele1.HapMapCEU)
  genome(gheight) <- "hg19"
  gscz <- GRanges(seqnames=Rle(paste0("chr",GWAS_scz_summary$chr)), IRanges(start=GWAS_scz_summary$bp, width=1, names=GWAS_scz_summary$snpid)
                  , p = GWAS_scz_summary$p, or=GWAS_scz_summary$or)
  genome(gscz) <- "hg19"
  
  save(gheight, gscz, file = here("output", "gwas_gr.rda"))
  rm(GWAS_height_summary, GWAS_scz_summary )
  gc()
} else {
  load( here("output", "gwas_gr.rda"))
}


hub <- AnnotationHub()

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

head(seqlevels(txdb_hg19))

sig_height <- gheight[gheight$p<0.001]
sig_scz <- gscz[gscz$p<0.001]
  
loc_height_hg19 <- locateVariants(sig_height, txdb_hg19, AllVariants())
loc_scz_hg19    <- locateVariants(sig_scz, txdb_hg19, AllVariants())

## Now get up and down for height and risk scz by geneid

b <- sig_height[loc_height_hg19$QUERYID]$b

b <- data.frame(b, geneid = loc_height_hg19$GENEID)
or <- data.frame(or = sig_scz[loc_scz_hg19$QUERYID]$or, geneid = loc_scz_hg19$GENEID)

gupheight <- b$geneid[b$b>0.0]
gdownheight <- b$geneid[b$b<0.0]

gupscz <- or$geneid[or$or<1]
gdownscz <- or$geneid[or$or>1]

length(unique(gupheight))
length(unique(gdownheight))
length(unique(gupscz))
length(unique(gdownscz))

length(unique((gupheight[gupheight %in% gdownheight])))
length(unique((gupscz[gupscz %in% gdownscz])))
length(unique((gupheight[gupheight %in% gupscz])))
length(unique((gupheight[gupscz %in% gupheight])))
length(unique((gdownheight[gdownheight %in% gupscz])))
length(unique((gdownheight[gupscz %in% gdownheight])))

length(unique((gupheight[gupheight %in% gupscz])))
length(unique((gupheight[gupscz %in% gupheight])))
length(unique((gdownheight[gdownheight %in% gupscz])))
length(unique((gdownheight[gupscz %in% gdownheight])))
## Look at coding
b <- data.frame(b = sig_height[loc_height_hg19$QUERYID]$b, geneid = loc_height_hg19$GENEID)[loc_height_hg19$LOCATION=="coding",]
or <- data.frame(or = sig_scz[loc_scz_hg19$QUERYID]$or, geneid = loc_scz_hg19$GENEID)[loc_scz_hg19$LOCATION=="coding",]

gupheight <- b$geneid[b$b>0.0]
gdownheight <- b$geneid[b$b<0.0]

max(table(gdownscz))
gupscz <- or$geneid[or$or<1]
gdownscz <- or$geneid[or$or>1]



length(unique(gupheight))
length(unique(gdownheight))
length(unique(gupscz))
length(unique(gdownscz))

length(unique((gupheight[gupheight %in% gdownheight])))
length(unique((gupscz[gupscz %in% gdownscz])))
length(unique((gupheight[gupheight %in% gupscz])))
length(unique((gupheight[gupscz %in% gupheight])))
length(unique((gdownheight[gdownheight %in% gupscz])))
length(unique((gdownheight[gupscz %in% gdownheight])))

length(unique((gupheight[gupheight %in% gupscz])))
length(unique((gupheight[gupscz %in% gupheight])))
length(unique((gdownheight[gdownheight %in% gupscz])))
length(unique((gdownheight[gupscz %in% gdownheight])))

length(unique((gupscz[gupscz %in% gdownscz])))

install.load.bioc("org.Hs.eg.db")
