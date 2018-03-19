
library(here)
source(here("R","prepare.R"))

if (!dir.exists(here("output"))) {
  dir.create(here("output"))
}

library(data.table)
load(here("output", "snp_match.rda"))

### annotation

install.load.bioc("VariantAnnotation", "AnnotationHub", "TxDb.Hsapiens.UCSC.hg19.knownGene" ,"clusterProfiler")

snp_both$z.scz <- log(snp_both$or)/snp_both$or.se
snp_both$z.height <- snp_both$b/snp_both$b.se
snp_both$zb.scz <- FIQT(snp_both$z.scz)
snp_both$zb.height <- FIQT(snp_both$z.height)


gboth <- GRanges(seqnames=Rle(snp_both$chrom), IRanges(start=snp_both$pos, width=1, names=snp_both$snpid)
                 , p.height = snp_both$p.height, f=snp_both$freq, snpid = snp_both$snpid
                 , p.scz = snp_both$p.scz, or=snp_both$or, z.scz = snp_both$zb.scz, z.height = snp_both$zb.height )

rm(snp_both)
#OK then , want to set this up to look for coding variants, promoter variants and intronic variants 
#and then compre things.  Buit we could start with variants of estimated large effect

# I would probably also like to have a gene score - perhaps by separating out the different types
# of variation


hub <- AnnotationHub()


txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

loc_both_hg19 <- locateVariants(gboth, txdb_hg19, AllVariants())
## get rid of those lines that are in the same genes, but are duplicates (hit different transcripts) 
newid <- paste(seqnames(loc_both_hg19), start(loc_both_hg19), loc_both_hg19$LOCATION, loc_both_hg19$GENEID, sep="_")
length(newid)
length(unique(newid))
loc_both_hg19_nodups <- loc_both_hg19[!duplicated(newid)]
table(loc_both_hg19$LOCATION)
table(loc_both_hg19_nodups$LOCATION)
#loc_both_hg19 <- loc_height_hg19
both_annotated <- data.table(data.frame(gboth[loc_both_hg19_nodups$QUERYID]
                                        , GENEID= loc_both_hg19_nodups$GENEID, 
                                        LOCATION=loc_both_hg19_nodups$LOCATION, 
                                        TXID = loc_both_hg19_nodups$TXID, 
                                        stringsAsFactors = FALSE))

save(both_annotated, file=here("output", "both_annotated.rda"))

q("no")