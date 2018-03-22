#' ## 3 way SNP analysis

library(here)
source(here("R","prepare.R"))
install.load("data.table")
install.load.bioc("ReactomePA","clusterProfiler")
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


make_labels <- function(names=c("iI", "hH", "nN")) {
  as.vector(t(outer(as.vector(t(o)), c(substr(names[3],1,1),"-",substr(names[3],2,2)), paste0)))

}

make_folded_labels <- function(names=c("iI", "hH", "nN")) {
  levels <- make_labels(names)
  c(paste(levels[1:13],rev(levels[15:27]),sep="/"), "---")
}

make_labels()
make_folded_labels()

gene_test <- gene_z[gene_n>=10]
classifydf <- function(zz, cutval=2 ) {
  zc <- sapply(zz,  cut, c(-Inf, -cutval, cutval, Inf), labels=FALSE)
  res <- 1+zc[,1]-1 + 3*(zc[,2]-1) + 9*(zc[,3]-1)
  res[res>14] <- 28-res[res>14]
  tabulate(res, nbin=14)
}

class_NHI <- t(sapply(gene_test, classifydf, cutval=5))
colnames(class_NHI) <- make_folded_labels()
IHN2 <- rownames(class_NHI)[class_NHI[,1]>=2]
IhN2 <- rownames(class_NHI)[class_NHI[,7]>=2]


res2.2 <- compareCluster(list(IHN=IHN2, IhN=IhN2), fun="enrichPathway")
plot(res2.2)
res2.1 <- compareCluster(list(IHN=rownames(class_NHI)[class_NHI[,1]>=1], IhN=rownames(class_NHI)[class_NHI[,7]>=1]), fun="enrichPathway")
plot(res2.1)


IHN2 <- rownames(class_NHI)[class_NHI[,1]>=2]
IhN2 <- rownames(class_NHI)[class_NHI[,7]>=2]


res2.2 <- compareCluster(list(IHN=IHN2, iHN=IhN2), fun="enrichPathway")
plot(res2.2)
res2.1 <- compareCluster(list(IHN=rownames(class_NHI)[class_NHI[,1]>=2], iHN=rownames(class_NHI)[class_NHI[,9]>=2]), fun="enrichPathway")
plot(res2.1)

res2.1 <- compareCluster(list(IHN=rownames(class_NHI)[class_NHI[,1]>=2], 
                              iHN=rownames(class_NHI)[class_NHI[,9]>=2],
                              IhN = rownames(class_NHI)[class_NHI[,7]>=2],
                              IHn = rownames(class_NHI)[class_NHI[,3]>=2]), fun="enrichPathway")
plot(res2.1)

try_pathway <- function(counts=class_NHI, comparisons = c("ihn", "iHn", "ihN", "Ihn"), nmin=2) {
  cols <- lapply(comparisons, grep, colnames(counts))
  print(cols)
  l <- lapply(cols, function(x) rownames(counts)[rowSums(counts[,x, drop=FALSE]) >= nmin])
  print(sapply(l, length))
  names(l) <- comparisons
  res <- compareCluster( l, fun="enrichPathway")
  plot(res)
}
  
try_pathway(comparisons=c("I[Hh]N", "I[Hh]n", "IHn", "IHN", "/IN", "/In"))

try_pathway(comparisons=c("[Ii]HN", "[Ii]hN", "iHN", "IHN", "/IN", "/In"))

try_pathway(comparisons=c("[Ii]HN", "[Ii]hN", "iHN", "IHN", "/IN", "/In"))

class_NHI5 <- t(sapply(gene_test, classifydf, cutval=5))
colnames(class_NHI5) <- make_folded_labels()

try_pathway(class_NHI5, comparisons=c("I", "H", "N", "IN", "IH","HN" ))
