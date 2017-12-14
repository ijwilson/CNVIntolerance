library(here)
source(here("R","prepare.R"))
install.load.bioc("GenomicRanges", "rtracklayer")

#file <- download.file(
#"ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/cnv/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed"
#, "output/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed")

exac <- import.bed(dropbox("exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed.txt"))
exac[[1]]$thick <- NULL
exac[[1]]$itemRgb <- NULL
exac[[2]]$thick <- NULL
exac[[2]]$itemRgb <- NULL

exac.del <- exac[[1]]
exac.del$type <- "deletion"
  
exac.dup <- exac[[2]]
exac.dup$type="duplication"
exac <- c(exac.del, exac.dup)

  ## We load in the exac list of genes with their new notation
#genefile <- download.file(
#  "ftp://ftp.broadinstitute.org/pub/ExAC_release/release1/cnv/exac-final-cnv.gene.scores071316"
#  , "output/exac-final-cnv.gene.scores071316")

df.exac.genes <- read.table(dropbox("exac-final-cnv.gene.scores071316.txt"), header=TRUE)
df.exac.genes$dup.score <- as.numeric(df.exac.genes$dup.score)
head(df.exac.genes)
  #
exac.genes <- GRanges(seqnames=paste("chr",df.exac.genes$chr,sep=""),
                        IRanges(start=df.exac.genes$start, end=df.exac.genes$end, names=df.exac.genes$gene),
                        gene_symbol=df.exac.genes$gene_symbol,
                        complexity=df.exac.genes$complexity,
                        dip=df.exac.genes$dip, del=df.exac.genes$del, dup=df.exac.genes$dup,
                        del.score=df.exac.genes$del.score, dup.score=df.exac.genes$dup.score,
                        flag=df.exac.genes$flag, num.targ=df.exac.genes$num_targ)
  
  
  
  
  
  
save(exac, exac.genes, file=here("output","exacbed.rda"))
  
library("here")
source(here("R","prepare.R"))
load(here("output", "rare.rda"))
load(here("output", "exacbed.rda"))

install.load.bioc("GenomicRanges")
## exac scores


getexacscore <- function(cnv) {
  o <- findOverlaps(cnv, exac.genes)
  if (length(o) ==0) return(0);
  if (cnv$cn<2) {
    scores <- exac.genes[subjectHits(o)]$del.score
    return(sum(scores))
  } else {
    scores <- exac.genes[subjectHits(o)]$dup.score
    return(sum(scores))
  }
}

grare1.exac.scores <- sapply(grare1, getexacscore)  
fqs.exac.scores <- sapply(fqsGR, getexacscore)
save(grare1.exac.scores,fqs.exac.scores, file=here("output","exacscores.rda"))

