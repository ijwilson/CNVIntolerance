## Gene centric comparison of ExAc and snz CNV results
## Light on my extra code

source("helper_functions.R")

## read SCZ cnv results
scz.del <- read.table("C:\\Users\\nijw\\Dropbox/CNVIntolerance/pgc_cnv/PGC_41K_QC_del_minimum8cnv.gene.results", header=TRUE)
head(scz)

## read in ExEc scores
exac.scores <- read.table("C:\\Users\\nijw\\Dropbox/CNVIntolerance/exac-final-cnv.gene.scores071316.txt", header=TRUE)
head(exac.scores)
colnames(scz)[1] <- "gene_symbol"  ## changing the column name to match exac.scores

## We can now merge the two data.tables

merged.table <- merge(scz, exac.scores, by="gene_symbol")
head(merged.table)

nrow(scz)
nrow(exac.scores)
nrow(merged.table)

plot(-log10(merged.table$dev_pval), merged.table$del.score)
## Does not seem to be much of a relationship.

sig <- merged.table$dev_pval<0.05
mean(merged.table$del.score[sig==TRUE])
mean(merged.table$del.score[sig==FALSE])
## No sign of any difference.  What happens if we just look at highly significant genes
sig <- merged.table$dev_pval<1E-5
mean(merged.table$del.score[sig==TRUE])
mean(merged.table$del.score[sig==FALSE])
t.test(merged.table$del.score[sig==TRUE],merged.table$del.score[sig==FALSE] )
## now look at the single del score

mean(merged.table$del.sing.score[sig==TRUE])
mean(merged.table$del.sing.score[sig==FALSE])
t.test(merged.table$del.sing.score[sig==TRUE],merged.table$del.sing.score[sig==FALSE] )
## This is highly significant
sig <- merged.table$dev_pval<0.05

mean(merged.table$del.sing.score[sig==TRUE])
mean(merged.table$del.sing.score[sig==FALSE])
t.test(merged.table$del.sing.score[sig==TRUE],merged.table$del.sing.score[sig==FALSE] )

plot(-log10(merged.table$dev_pval), merged.table$del.sing.score)

### Now we can do the same for height

## read SCZ cnv results
height_cnv <- read.csv("C:\\Users\\nijw\\Dropbox/CNVIntolerance/height_CNV_association_41467_2017_556_MOESM2_ESM.csv")
head(height_cnv)

## Need to get the genes associated with each CVN

install.load.bioc("GenomicRanges")


gheight <- GRanges(seqnames=height_cnv$CHR, IRanges(start=height_cnv$BP, end=height_cnv$BP+1)
                   , fdel = height_cnv$F_DEL, fdup=height_cnv$F_DUP, p=height_cnv$Pvalue.Height, 
                   beta = height_cnv$Beta.Height )

## Need teh locations of some genes


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

## find the gene for each of these CNVs
co <- countOverlaps(gheight, genesGR)
gheight2 <- gheight[co>0]
gg <- findOverlaps(gheight2[co>0], genesGR)
dups <- duplicated(queryHits(gg))
gg <- gg[!dups]

gheight2$gene_symbol <- names(genesGR)[subjectHits(gg)] 

height_score <- tapply(gheight2$beta, gheight2$gene_symbol, mean)
height_positive <- tapply(gheight2$beta, gheight2$gene_symbol, function(x) sum(x>0))
height_positive <- tapply(gheight2$beta, gheight2$gene_symbol, function(x) sum(x>0))


p_score <- tapply(gheight2$p, gheight2$gene_symbol, function(x) mean(x<0.05))
n <- tapply(gheight2$p, gheight2$gene_symbol, length)

del_change <- by(as.data.frame(gheight2)[,c(6,8,9)], gheight2$gene_symbol, function(x) sum(x[,1]*abs(x[,3])), simplify = TRUE )



d <- data.frame(height_score, p_score, n, height_positive, del.change =as.vector(del_change))
d$gene_symbol <- rownames(d)

merged.height <- merge(d, exac.scores, by="gene_symbol")

install.load("ggplot2")

ggplot(merged.height, aes(x=del.change, y=del.score)) + geom_point()
