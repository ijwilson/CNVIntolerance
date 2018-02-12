#'## Plots
library(here)
source(here("R","prepare.R"))
install.load.bioc("GenomicRanges", "ggbio")

bands <- read.table(dropbox("cytoBand.txt"))
gbands <- GRanges(substring(bands$V1,4), 
                  IRanges(bands$V2, end=bands$V3, name=paste(substring(bands$V1, 4), bands$V4, sep="")))

b22q11.2 <- GRanges("chr22", IRanges(18660553,21455556))


gbands["22q11.2"]
#` A Good place to start is http://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf

p.ideo <- Ideogram(genome="hg19")
p.ideo + xlim(GRanges("chr22", IRanges()))
load(system.file("data", "hg19IdeogramCyto.rda", package="biovizBase", mustWork=TRUE))
pl <- plotIdeogram(hg19IdeogramCyto, "chr22") + xlim(b22q11.2)

pl
autoplot(Homo.sapiens, which = b22q11.2, stat = "reduce")
install.load.bioc("Homo.sapiens")
class(Homo.sapiens)
##


tx <- autoplot(Homo.sapiens, which=b22q11.2)
tx
autoplot(Homo.sapiens, which  = wh, label.color = "black", color = "brown",
         fill = "brown")

wh <- b22q11.2
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
tx2 <- autoplot(txdb, which = wh)

tracks(pl, tx)

data(genesymbol, package = "biovizBase")
wh <- genesymbol[c("BRCA1", "NBR1")]
wh <- range(wh, ignore.strand = TRUE)
p.txdb <- autoplot(Homo.sapiens, which  = wh)
p.txdb
autoplot(Homo.sapiens, which  = wh, label.color = "black", color = "brown",
         fill = "brown")

library(biovizBase)
gr.txdb <- crunch(txdb, which = wh)
## change column to 'model'
colnames(values(gr.txdb))[4] <- "model"
grl <- split(gr.txdb, gr.txdb$tx_id)

## fake some randome names
names(grl) <- sample(LETTERS, size = length(grl), replace = TRUE)
grl

autoplot(grl, aes(type = model))
ggplot() + geom_alignment(grl, type = "model")


################
install.load.bioc("Gviz")
data(geneModels)


