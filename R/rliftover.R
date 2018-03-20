## rliftover
## biodonductor workflow
# https://www.bioconductor.org/help/workflows/liftOver/

library(here)
source(here("R","prepare.R"))
install.load.bioc("rtracklayer")
# , "rtracklayer")

chainfile <- here("data", "hg18ToHg19.over.chain")
if (!file.exists(chainfile)) {
  download.file("http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz", here("data", "hg18ToHg19.over.chain.gz"))
  lines <- readLines(here("data", "hg18ToHg19.over.chain.gz"))
  con <- file(chainfile)
  writeLines(lines, con=con)
  close(con)
}
chain <- import.chain(chainfile)


install.load("data.table")
install.load.bioc("GenomicRanges") #, "rtracklayer")

height <- fread(dropbox("height_CNV_association_41467_2017_556_MOESM2_ESM.csv"))   ## standard CNV
height$chrom <- paste0("chr", height$CHR)

gheight <- GRanges(seqnames=paste0("chr",height$CHR), IRanges(start=height$BP, width=1, names = height$ID)
                   , F_DEL=height$F_DEL)
gheight <- makeGRangesFromDataFrame(height,
                         keep.extra.columns=TRUE,
                         ignore.strand=TRUE,
                         seqinfo=NULL,
                         seqnames.field=c("chrom"),
                         start.field="BP",
                         end.field=c("BP"))

gheight19 <- unlist(liftOver(gheight, chain))

## Read the non-gene centric results.
scz.dup <- fread(dropbox("pgc_cnv/PGC_41K_QC_dup.cnv.results"))
scz.dup$chrom <- paste0("chr",scz.dup$CHR)
scz.dup$chrom[scz.dup$CHR=="23"] <- "chrX"
gscz.dup <- makeGRangesFromDataFrame(scz.dup,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqnames.field=c("chrom"),
                                     start.field="BP",
                                     end.field=c("BP"))
gscz.dup19 <- unlist(liftOver(gscz.dup, chain))
scz.del <- fread(dropbox("pgc_cnv/PGC_41K_QC_del.cnv.results"))
scz.del$chrom <- paste0("chr",scz.del$CHR)
scz.del$chrom[scz.del$CHR=="23"] <- "chrX"
gscz.del <- makeGRangesFromDataFrame(scz.del,
                                     keep.extra.columns=TRUE,
                                     ignore.strand=TRUE,
                                     seqnames.field=c("chrom"),
                                     start.field="BP",
                                     end.field=c("BP"))

gscz.del19 <- unlist(liftOver(gscz.del, chain))



###############################################################
## Gene centric scz analysis
## This is a biut more complicated to liftover but I'll try
scz.del.gene <- fread(dropbox("pgc_cnv/PGC_41K_QC_del_minimum8cnv.gene.results"))
scz.del.gene <- scz.del.gene[!is.na(scz.del.gene$BP_start_hg18),]
scz.all.gene <- fread(dropbox("pgc_cnv/PGC_41K_QC_all_minimum12cnv.gene.results"))
scz.all.gene <- scz.all.gene[!is.na(scz.all.gene$BP_start_hg18),]




