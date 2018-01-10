library(here)
source(here("R","prepare.R"))

install.load("DBI","RSQLite")
install.load.bioc("org.Hs.eg.db","TxDb.Hsapiens.UCSC.hg19.knownGene")

txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene
genes <-  transcriptsBy(txdb, "gene")
genes.gr <- reduce(genes)
genes.df <- as(genes.gr, "data.frame")

x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
xx <- as.list(x[mapped_genes])

## ------------------------------------------------------------------------
m <- match(genes.df$group_name, names(xx))
genes.df <- cbind(genes.df[!is.na(m),], name = unlist(xx[m]))
genesGR <- GRanges(seqnames=substring(genes.df$seqnames,4),ranges = IRanges(genes.df$start,genes.df$end, names=genes.df$name), strand=genes.df$strand)

seqlevels(genesGR) <- paste("chr", seqlevels(genesGR), sep="")

if (!dir.exists(here("output"))) 
  dir.create("output")
save(genesGR, file="output/genesGR.rda")
