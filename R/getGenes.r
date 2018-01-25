library(here)
source(here("R","prepare.R"))

install.load("DBI","RSQLite")
install.load.bioc("org.Hs.eg.db"
                  ,"TxDb.Hsapiens.UCSC.hg18.knownGene", "TxDb.Hsapiens.UCSC.hg19.knownGene")

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
genesGR19 <- GRanges(seqnames=substring(genes.df$seqnames,4),ranges = IRanges(genes.df$start,genes.df$end, names=genes.df$name), strand=genes.df$strand)


txdb <- TxDb.Hsapiens.UCSC.hg18.knownGene
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
genesGR18 <- GRanges(seqnames=substring(genes.df$seqnames,4),ranges = IRanges(genes.df$start,genes.df$end, names=genes.df$name), strand=genes.df$strand)


if (!dir.exists(here("output"))) dir.create("output")
save(genesGR19, genesGR18, file="output/genesGR.rda")
