#' ## 3 way SNP analysis

library(here)
source(here("R","prepare.R"))
install.load("data.table")
load(here("output", "snp_NHI.rda"))

#' ## Simplfy things to get some simple statistics and plots

z <- data.frame(intelligence = snp_shn$zb.i, height=snp_shn$zb.h, neuroticism = snp_shn$zb.n) 

gene_z <- split(z, snp_shn$GENEID)
gene_mad <- by(z, snp_shn$GENEID, function(x) colMeans(abs(x)))
gene_m <- matrix(unlist(gene_mad), ncol=3, byrow=TRUE)
rownames(gene_m) <- names(gene_mad)
colnames(gene_m) <- c("intelligence", "height", "neuroticism")
gene_n <- tapply(snp_shn$GENEID, snp_shn$GENEID, length)
gene_cor <- lapply(gene_z, cor)

n_agree <- t(sapply(gene_z, function(x) apply(x, 2, function(z) sum(abs(z)>2))))

o <- as.vector(t(outer(c("i","","I"), c("h", "", "H"), paste0)))
levels <- as.vector(t(outer(as.vector(t(o)), c("n","","N"), paste0)))


classify <- function(zz) {
  zc <- apply(zz, 2, cut, c(-Inf, -2, 2, Inf), labels=FALSE)
  1+zc[,1]-1 + 3*(zc[,2]-1) + 9*(zc[,3]-1)
}

z_cut <- classify(z)
z_cut <- factor(z_cut, labels=levels)
table( scz=cut(z$scz, c(-Inf, -2, 2, Inf)),height= cut(z$height, c(-Inf,-2,2,Inf)))
table( scz=cut(z$scz, c(-Inf, -2, 2, Inf)), neuroticism=cut(z$neur, c(-Inf,-2,2,Inf)))
table( height=cut(z$height, c(-Inf, -2, 2, Inf)), neuroticism=cut(z$neur, c(-Inf,-2,2,Inf)))

## This needs to be done by gene or there is little point.


tb <- table(z_cut)
tb2 <- data.frame(tb[1:13],rev(tb[15:27]))
 tb4 <- data.frame( changes = c(paste(tb2[,3],paste(tb2[,1]),sep="/"), "none"), count= c(rowSums(tb2[,c(2,4)]), tb[14]))



tbtbsum((x$scz > 2 & x$height > 2 & x$neur > 2) | (x$scz < -2 & x$height < -2 & x$neur < -2)
  ))

length(gene_n)
colMeans(z)

