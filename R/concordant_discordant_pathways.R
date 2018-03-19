
#+ prepare_data, cache=TRUE
library(here)
source(here("R","prepare.R"))
install.load.bioc("clusterProfiler")
install.load("data.table", "pander")
install.load.bioc("ReactomePA")

load(here("output", "both_annotated.rda"))


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


#+ calculate_statistics, cache=TRUE

#' ### Let's try everything at once to see what we can get

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

n <- tapply(both_annotated$GENEID, both_annotated$GENEID, length)
plot(n, dis[,1]+con[,1])

nsig <- dis+con
p <- nsig/matrix(n, ncol=6, nrow=length(n))

#+ plot_means, fig.height=5, fig.width=5
plot(colMeans(p))


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

dis <- sapply(gene_stats, function(xx) sapply(xx, function(yy) yy[2]+yy[3]))
con <- sapply(gene_stats, function(xx) sapply(xx, function(yy) yy[1]+yy[4]))


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

#+ plot_res2, fig.width=8, fig.height=6

plot(res2)

#+ plot_res5, fig.width=8, fig.height=6
plot(res5)

#+ fig.width=8, fig.height=6
plot(res10)


#` Find another way to rearrange the results to check for the optimum.

lista <- list(d2=disgenes2[[1]], d5=disgenes5[[1]], d10=disgenes10[[1]], c2=congenes2[[1]], c5=congenes5[[1]], c10=congenes10[[1]])
listb <- list(d2=disgenes2[[2]], d5=disgenes5[[2]], d10=disgenes10[[2]], c2=congenes2[[2]], c5=congenes5[[2]], c10=congenes10[[2]])
listc <- list(d2=disgenes2[[3]], d5=disgenes5[[3]], d10=disgenes10[[3]], c2=congenes2[[3]], c5=congenes5[[3]], c10=congenes10[[3]])
listd <- list(d2=disgenes2[[4]], d5=disgenes5[[4]], d10=disgenes10[[4]], c2=congenes2[[4]], c5=congenes5[[4]], c10=congenes10[[4]])
liste <- list(d2=disgenes2[[5]], d5=disgenes5[[5]], d10=disgenes10[[5]], c2=congenes2[[5]], c5=congenes5[[5]], c10=congenes10[[5]])
listf <- list(d2=disgenes2[[6]], d5=disgenes5[[6]], d10=disgenes10[[6]], c2=congenes2[[6]], c5=congenes5[[6]], c10=congenes10[[6]])


#` Try to compare the clusters@?'`
#+ cache=TRUE` 
resa <- compareCluster(lista, fun="enrichPathway")
resb <- compareCluster(listb, fun="enrichPathway")
resc <- compareCluster(listc, fun="enrichPathway")
resd <- compareCluster(listd, fun="enrichPathway")
rese <- compareCluster(liste, fun="enrichPathway")
resf <- compareCluster(listf, fun="enrichPathway")

#+ resa, fig.width=8, fig.height=6
plot(resa)
#+ fig.width=8, fig.height=6
plot(resb)
#+ fig.width=8, fig.height=6
plot(resc)
#+ fig.width=8, fig.height=6
plot(resd)
#+ fig.width=8, fig.height=6
plot(rese)
#+ fig.width=8, fig.height=6
plot(resf)

tryfinal <- list(concordantgenes = congenes2[[2]], discordant_genes = disgenes2[[2]])

plot(compareCluster(tryfinal,  fun="enrichPathway"))



