#' ## 3 way SNP analysis

library(here)
source(here("R","prepare.R"))
install.load("data.table")
load(here("output", "snp_NHI.rda"))

#' ## Simplfy things to get some simple statistics and plots
#' 
#' 

#+ original_z, fig.width=8, fig.height=8
plot(density(snp_NHI$z.h), xlim=c(-5,5), ylim=c(0,0.39), axes=FALSE, main="Z score densities", lwd=3, col=4, ylab="", xlab="Z")
lines(density(snp_NHI$z.int), col=2, lwd=3)
lines(density(snp_NHI$z.n), col=3, lwd=3)
curve(dnorm(x), add=T, col=1)
axis(1)


#+ transformed_z, fig.width=8, fig.height=8
plot(density(snp_NHI$zb.h), xlim=c(-5, 5), ylim=c(0,2.9), axes=FALSE, main="Transformed Z score densities", lwd=3, col=4, ylab="", xlab="Z")
lines(density(snp_NHI$zb.i), col=2, lwd=3)
lines(density(snp_NHI$zb.n), col=3, lwd=3)
curve(dnorm(x), add=T, col=1)
axis(1)



z <- data.frame(intelligence = snp_NHI$z.int, height=snp_NHI$z.h, neuroticism = snp_NHI$z.n) 

gene_z <- split(z, snp_NHI$GENEID)
gene_mad <- t(sapply(gene_z, function(x) colMeans(abs(x))))
gene_means <- t(sapply(gene_z, colMeans))
gene_ss <- t(sapply(gene_z, function(x) colSums(x^2)))

gene_m <- matrix(unlist(gene_mad), ncol=3, byrow=TRUE)
rownames(gene_m) <- names(gene_mad)
colnames(gene_m) <- c("intelligence", "height", "neuroticism")
gene_n <- tapply(snp_NHI$GENEID, snp_NHI$GENEID, length)
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
table( intelligence=cut(z$intelligence, c(-Inf, -2, 2, Inf)),height= cut(z$height, c(-Inf,-2,2,Inf)))
table( intelligence=cut(z$intelligence, c(-Inf, -2, 2, Inf)), neuroticism=cut(z$neur, c(-Inf,-2,2,Inf)))
table( height=cut(z$height, c(-Inf, -2, 2, Inf)), neuroticism=cut(z$neur, c(-Inf,-2,2,Inf)))

## This needs to be done by gene or there is little point.


tb <- table(z_cut)
tb2 <- data.frame(tb[1:13],rev(tb[15:27]))
 tb4 <- data.frame( changes = c(paste(tb2[,3],paste(tb2[,1]),sep="/"), "none"), count= c(rowSums(tb2[,c(2,4)]), tb[14]))

install.load("pander")
pander(tb4)

#' ##Gene Centic Analysis

hist(log10(gene_n), main="SNPs per gene", axes=FALSE, xlab="Number of SNPs")
axis(1, labels=c(1,10,100,1000), at=c(0,1,2,3))

pairs(gene_means[gene_n>=10,])
cor(gene_means[gene_n>=20,])

pairs(gene_mad[gene_n>=20,])
cor(gene_mad)

summary(lm(intelligence~neuroticism+height, weights=1/sqrt(gene_n), data=data.frame(gene_means)))

install.load.bioc("aroma.light")
pca <- wpca(gene_means, w=sqrt(gene_n), swapDirections = TRUE)
pca$vt

#` ### K means clustering of genes`
k <- kmeans(gene_means, centers=4)

ncpca <- prcomp(gene_means)

gene_means <- data.frame(gene_means)
pairs(gene_means, col=k$cluster)
gene_z2 <- gene_z[gene_n>10]

gene_t <- t(
  sapply(gene_z2, function(y) sapply(y, function(x) t.test(x,mu=0)$st))
  )




#' ## four sets
#' Look for the contribution of the different ways that a set of three can be related.
#' or look at the simple 2 way statistics





