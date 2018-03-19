
#+ prepare_libraries
library(here)
source(here("R","prepare.R"))
install.load.bioc("clusterProfiler")
install.load("data.table", "pander")
install.load.bioc("ReactomePA")

#+ prepare_data, cache=TRUE
load(here("output", "both_annotated.rda"))

#' ### What do these data mean
#'
#' I have tried various versions of the pathway analysis and a gene centric analysis and it seems
#' that there is something here - in the correlation and anti-correlation of the SCZ and height SNPs - and
#' I feel that there is information about biology that can be glenaed from these data but I need to know how to 
#' get at it.
#' 
#' * What it is 
#' * what can I learn about biology from it
#' * Is there extra power that comes from looking at two independent data sets.
#' 
#'  I suppose there must be more evidence that you are looking at actual biuological pathways rather than just noise.
#'   
#' 
# 
# allg <- both_annotated[LOCATION!="intergenic"]
# allc <- both_annotated[LOCATION=="coding"]
# alli <- both_annotated[LOCATION=="intron"]
# allp <- both_annotated[LOCATION=="promoter"]
# allr <- both_annotated[LOCATION%in% c("fiveUTR", "threeUTR","promoter")]
# # 
# val <- c(-2, 2)
# gene_stats <-  by(both_annotated, both_annotated$GENEID, function(x) c(
#   sum(x$z.height < val[1] & x$z.scz < val[1]), 
#   sum(x$z.height < val[1] & x$z.scz > val[2]),
#   sum(x$z.height > val[2] & x$z.scz < val[1]),
#   sum(x$z.height > val[2] & x$z.scz > val[2]))
# , simplify = TRUE )

## A quick investigation of 

#+ acf_plot, fig.height=8, fig.width=8
opar <- par(mfrow=c(2,2))
acf(both_annotated$z.height, main="z height", lag.max=20)
acf(both_annotated$z.scz, main="z height", lag.max=20)

acf(abs(both_annotated$z.height), main="|z height|", lag.max=20)
acf(abs(both_annotated$z.scz), main="|z SCZ|", lag.max=20)
par(opar)

#' But what about the allele frequencies
#' Can I look at just the rare alleles and common alleles

rare <- both_annotated[pmin(both_annotated$f, 1-both_annotated$f) < 0.02]
nrow(rare)

rest <-  both_annotated[pmin(both_annotated$f, 1-both_annotated$f) > 0.1 & pmin(both_annotated$f, 1-both_annotated$f) < 0.3]
nrow(rest)
common <- both_annotated[pmin(both_annotated$f, 1-both_annotated$f) >0.45]
nrow(common)

plot(density(rare$z.height), xlim=c(-2,2), ylim=c(0,1.6)
     , main="Marginal densities of z height and z scz", xlab="z", axes=FALSE)
lines(density(common$z.height), lty=3)
lines(density(rare$z.scz), col="blue")
lines(density(common$z.scz), lty=3, col="blue")
lines(density(rest$z.height), lty=2)
lines(density(rest$z.scz), lty=2, col="blue")


both_annotated$fc <- cut(pmin(both_annotated$f, 1-both_annotated$f),c(0, 0.02,0.05,0.1,0.2,0.3,0.4, 0.5))
table(both_annotated$fc)
tb <- table(both_annotated$fc, both_annotated$LOCATION)
p <- t(t(tb)/colSums(tb))
colSums(p)
barplot(p)


both_annotated$fc <- cut(pmin(both_annotated$f, 1-both_annotated$f),c(0,0.05,0.1,0.2,0.3,0.4, 0.5))




ms <- tapply(abs(both_annotated$z.scz), both_annotated$fc, function(x) mean(x>1.5))
mh <- tapply(abs(both_annotated$z.height), both_annotated$fc, function(x) mean(x>1.5))
mbc <- by(both_annotated, both_annotated$fc,
  function(x) mean((x$z.height>1.5 & x$z.scz>1.5) | (x$z.height < -1.5 & x$z.scz < -1.5)))
mbd <- by(both_annotated, both_annotated$fc,
         function(x) mean((x$z.height>1.5 & x$z.scz < -1.5) | (x$z.height < -1.5 & x$z.scz > 1.5)))



ms <- tapply(abs(both_annotated$z.scz), both_annotated$LOCATION, function(x) mean(x>2.5))
mh <- tapply(abs(both_annotated$z.height), both_annotated$LOCATION, function(x) mean(x>2.5))
barplot(t(cbind(ms, mh)))

mbc <- by(both_annotated, both_annotated$LOCATION,
          function(x) mean((x$z.height>2.5 & x$z.scz>2.5) | (x$z.height < -2.5 & x$z.scz < -2.5)))
mbd <- by(both_annotated, both_annotated$LOCATION,
          function(x) mean((x$z.height>2.5 & x$z.scz < -2.5) | (x$z.height < -2.5 & x$z.scz > 2.5)))

plot(mbc, mbd)
abline(0,1)
plot(ms, mh)

round(cbind(concordant=(2*mbc)/(mh*ms), discordant = 2*mbd/(ms*mh))[-1,],2)


#' Maybe I need these sats by gene
#' 
#' 

z.height <- tapply(abs(both_annotated$z.height), both_annotated$GENEID, mean)
z.scz <- tapply(abs(both_annotated$z.scz), both_annotated$GENEID, mean)
z.diff<- tapply(abs(both_annotated$z.scz-both_annotated$z.height), both_annotated$GENEID, mean)
z.sum<- tapply(abs(both_annotated$z.scz+both_annotated$z.height), both_annotated$GENEID, mean)


diff <- (z.diff/z.sum > 0.8) & z.sum > 6
table(diff)
cond <- (z.diff/z.sum < 0.4) & z.sum > 6
table(cond)
table(diff, cond)


res <- compareCluster(list(con = names(cond)[cond],
                           dis = names(diff)[diff]), fun="enrichPathway")

#+ fig.height=6, fig.width=10
plot(res)




plot(density(both_annotated$z.height), xlim=c(-4,4), ylim=c(0,1.6)
     , main="Marginal densities of z height and z scz", xlab="z", axes=FALSE)
axis(1);axis(2)
lines(density(both_annotated$z.scz), col="green")
curve(dnorm(x), add=T, col="blue")

probs <- c(0.005,0.01,0.025,0.05,0.1,0.5,0.9,0.95,0.975,0.99, 0.995)
tbh <- quantile(both_annotated$z.height, probs=probs)
tbz <- quantile(both_annotated$z.scz, probs=probs)

install.load("pander")
pander(cbind("z height"=tbh,"z scz"=tbz))


#' ## Let's have a look at regression
#+ regression_by_gene

reg <- by(both_annotated, both_annotated$GENEID, function(x) {
  lm(x$z.height~x$z.scz)
}
)



corr.gene <- by(both_annotated, both_annotated$GENEID, function(x) {
    if (nrow(x)>=5) 
      return( list( cor.test(x$z.height, x$z.scz), c(max(abs(c(x$z.scz))), max(abs(x$z.height))))  )
    else return(NULL)
  }
)

corr.geneb <- Filter(Negate(is.null), corr.gene)
corr.p <- sapply(corr.geneb, function(x) x[[1]]$p.value)
corr.x <- sapply(corr.geneb, function(x) x[[1]]$estimate)
corr.df <- sapply(corr.geneb,function(x) x[[1]]$parameter)
corr.maxscz <- sapply(corr.geneb, function(x) x[[2]][1])
corr.maxh <- sapply(corr.geneb, function(x) x[[2]][2])

concordant <- corr.x > 0.75 & corr.df>=5 &   (abs(corr.maxscz) > 3 | abs(corr.maxh) > 3) 
discordant <- corr.x < -0.75 & corr.df>=5 &  (abs(corr.maxscz) > 3 | abs(corr.maxh) > 3) 
table(concordant, discordant)

res <- compareCluster(list(con = names(corr.geneb[concordant]),
                           dis = names(corr.geneb[discordant])), fun="enrichPathway")

#plot(res)




reg.coef <- sapply(reg,  coefficients)
pval <-  sapply(reg, function(x) anova(x)$`Pr(>F)`[1])   
  
  

reg0 <- by(both_annotated, both_annotated$GENEID, function(x) {
  if (nrow(x)>=5) 
    return(lm(x$z.height~x$z.scz-1))
  else return(NULL)
}
)

reg0b <- Filter(Negate(is.null), reg0)

length(reg0)
length(reg0b)

reg0.coef <- sapply(reg0b, coefficients)
pval0 <-  sapply(reg0b, function(x) anova(x)$`Pr(>F)`[1])   


both_annotated[,r:=sqrt(z.height^2 + z.scz^2)]
both_annotated[, theta:=atan2(z.height, z.scz)]
plot(both_annotated$theta, both_annotated$r)




plot(pval0, pval)

#' ### Let's try everything at once to see what we can get

#+ calculate_statistics, cache=TRUE
which(pval0<0.001 & abs(reg0.coef) < 0.05)


gene_stats <- lapply(list(
  a=c(-2,2), b=c(-2.2,2.2), c=c(-2.5,2.5), d=c(-2.8,2.8), e=c(-3,3), f=c(-4,4))
                     , function(val) 
                       by(both_annotated, both_annotated$GENEID, function(x) c(
                         sum(x$z.height < val[1] & x$z.scz < val[1]), 
                         sum(x$z.height < val[1] & x$z.scz > val[2]),
                         sum(x$z.height > val[2] & x$z.scz < val[1]),
                         sum(x$z.height > val[2] & x$z.scz > val[2])))
)

dis <- sapply(gene_stats, function(xx) sapply(xx, function(yy) yy[2]+yy[3]))
con <- sapply(gene_stats, function(xx) sapply(xx, function(yy) yy[1]+yy[4]))

#' and comments to break the chunk up
#+ get_summaries

n <- tapply(both_annotated$GENEID, both_annotated$GENEID, length)
plot(n, dis[,1]+con[,1])

plot(dis[,1], con[,1])

nsig <- dis+con
p <- nsig/matrix(n, ncol=6, nrow=length(n))

#+ plot_means, fig.height=5, fig.width=5
plot(colMeans(p))
opar <- par(mfrow=c(2,2))
for (i in 1:4) plot(dis[,i], p[,i])
par(opar)
#` if we take a significant level of 5%, then the probability of at least one (if they are uncorrelated)
#` is 1 - .975^2 ~ approx 0.1

# sig4 <-  by(both_annotated, both_annotated$GENEID, function(x) c(
#   sum(x$z.height < val[1]), 
#   sum(x$z.height > val[2]),
#   sum(x$z.scz < val[1]),
#   sum(x$z.scz > val[2])
# )
# )

#+ cluster_comparison, cache=TRUE

disgenes10 <- apply(dis, 2, function(x) names(x[x>=10]))
congenes10 <- apply(con, 2,function(x) names(x[x>=10]))

disgenes5 <- apply(dis,2, function(x) names(x[x>=5]))
congenes5 <- apply(con, 2,function(x) names(x[x>=5]))

disgenes2 <- apply(dis,2, function(x) names(x[x>=2]))
congenes2 <- apply(con, 2, function(x) names(x[x>=2]))

names(disgenes2) <- paste("dis", names(disgenes2), sep="_")
names(congenes2) <- paste("con", names(congenes2), sep="_")
names(disgenes5) <- paste("dis", names(disgenes5), sep="_")
names(congenes5) <- paste("con", names(congenes5), sep="_")
names(disgenes10) <- paste("dis", names(disgenes10), sep="_")
names(congenes10) <- paste("con", names(congenes10), sep="_")


res2 <- compareCluster(c(congenes2, disgenes2), fun="enrichPathway")
res5 <- compareCluster(c(congenes5, disgenes5), fun="enrichPathway")
res10 <- compareCluster(c(congenes10, disgenes10), fun="enrichPathway")

#+ plot_res2, fig.width=10, fig.height=6
plot(res2)

#+ plot_res5, fig.width=10, fig.height=6
plot(res5)

#+ plot_res10, fig.width=10, fig.height=6
plot(res10)


#` Find another way to rearrange the results to check for the optimum.

#+ create_list
lista <- list(d2=disgenes2[[1]], d5=disgenes5[[1]], d10=disgenes10[[1]], c2=congenes2[[1]], c5=congenes5[[1]], c10=congenes10[[1]])
listb <- list(d2=disgenes2[[2]], d5=disgenes5[[2]], d10=disgenes10[[2]], c2=congenes2[[2]], c5=congenes5[[2]], c10=congenes10[[2]])
listc <- list(d2=disgenes2[[3]], d5=disgenes5[[3]], d10=disgenes10[[3]], c2=congenes2[[3]], c5=congenes5[[3]], c10=congenes10[[3]])
listd <- list(d2=disgenes2[[4]], d5=disgenes5[[4]], d10=disgenes10[[4]], c2=congenes2[[4]], c5=congenes5[[4]], c10=congenes10[[4]])
liste <- list(d2=disgenes2[[5]], d5=disgenes5[[5]], d10=disgenes10[[5]], c2=congenes2[[5]], c5=congenes5[[5]], c10=congenes10[[5]])
listf <- list(d2=disgenes2[[6]], d5=disgenes5[[6]], d10=disgenes10[[6]], c2=congenes2[[6]], c5=congenes5[[6]], c10=congenes10[[6]])


#` Try to compare the clusters@?'`
#+ get_lista-f_comparisions, cache=TRUE
resa <- compareCluster(lista, fun="enrichPathway")
resb <- compareCluster(listb, fun="enrichPathway")
resc <- compareCluster(listc, fun="enrichPathway")
resd <- compareCluster(listd, fun="enrichPathway")
rese <- compareCluster(liste, fun="enrichPathway")
resf <- compareCluster(listf, fun="enrichPathway")

#+ resa, fig.width=10, fig.height=6
plot(resa)
#+ resb, fig.width=10, fig.height=6
plot(resb)
#+ resc, fig.width=10, fig.height=6
plot(resc)
#+ fig.width=10, fig.height=6
plot(resd)
#+ fig.width=10, fig.height=6
plot(rese)
#+ fig.width=10, fig.height=6
plot(resf)
#+ tryfinal, cache=TRUE,  fig.width=10, fig.height=6
tryfinal <- list(concordantgenes = congenes2[[2]], discordant_genes = disgenes2[[2]])

plot(compareCluster(tryfinal,  fun="enrichPathway"))



