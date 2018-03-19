
library(here)
source(here("R","prepare.R"))
install.load.bioc("clusterProfiler")
install.load("data.table", "pander")

load(here("output", "both_annotated.rda"))



allg <- both_annotated[LOCATION!="intergenic"]
allc <- both_annotated[LOCATION=="coding"]
alli <- both_annotated[LOCATION=="intron"]
allp <- both_annotated[LOCATION=="promoter"]
allr <- both_annotated[LOCATION%in% c("fiveUTR", "threeUTR","promoter")]

val <- c(-3, 3)
by(both_annotated, both_annotated$GENEID, function(x) c(
  sum(x$z.height< -val &  x$z.scz < -val, x$z.height > val & x$z.scz > val))

tapply(both_annotated)



plot(allg$f2[allg$p.height<0.001], allg$b[allg$p.height<0.001], col=factor(allg$LOCATION[allg$p.height<0.001]))


rm(gboth)
gc()

install.load.bioc("ReactomePA")
##########################################

paired_levels <- (gboth$z.scz< -3) + 2*(gboth$z.scz > 3) + 4*(gboth$z.height< -3) + 8*(gboth$z.height>3)
fpaired_levels <- factor(paired_levels, levels=c(0,1,2,4,5,6,8,9,10), 
                         labels=c("none", "scz_decrease", "scz_increase"
                                  , "height_decrease", 
                                  "scz_height_decrease", "scz_decrease_height_increase",
                                  "height_increase",
                                  "height_increase_scz_decrease", "height_scz_increase")
                         )

getlists <- function(xx, extra="") {
  l <- list(
    unique(
      c(xx$GENEID[xx$z.scz    < -3 & xx$z.height>3],
      xx$GENEID[xx$z.height < -3 & xx$z.scz > 3])
      ),
     unique(c(xx$GENEID[xx$z.scz     > 3 & xx$z.height>3],
    xx$GENEID[xx$z.height < -3 & xx$z.scz < -3]))
  )
  names(l) =paste(c("up_down", "up_up"),extra,sep="")
  l
}


res <- compareCluster(getlists(allg), fun="enrichPathway")
plot(res)  


res <- compareCluster(c(getlists(alli, "_intron"), getlists(allc,"_coding"), getlists(allp,"_regulatory")), fun="enrichPathway")
plot(res)  


resp <- compareCluster(getlists(allp), fun="enrichPathway")
plot(resp)  


resc <- compareCluster(getlists(allc), fun="enrichPathway")
plot(resc)  


res <- compareCluster(l, fun="enrichPathway")
plot(res)  


l <- list(
  scz_down = unique(allg$GENEID[allg$z.scz< -3]),
  scz_up = unique(allg$GENEID[allg$z.scz > 3]),
  height_down = unique(allg$GENEID[allg$z.height< -3]),
    height_up= unique(allg$GENEID[allg$z.height>3])
)
res <- compareCluster(l, fun="enrichPathway")
plot(res)  
ld <- list(
  h_up = unique(allg$GENEID[allg$z.scz< -3 & allg$z.height > 3]),
  s_up = unique(allg$GENEID[allg$z.scz > 3 & allg$z.height < -3])
)

res <- compareCluster(l, fun="enrichPathway")
plot(res)  

res_discordant <- compareCluster(ld, fun="enrichPathway")
plot(res_discordant) 


get_genes <- function(xx, minpheight, minpscz, minb=-1E10, maxb=1E10, minor=0, maxor=1E9,minf2=0, maxf2=1) {
  unique(xx$GENEID[xx$p.height<=minpheight & xx$p.scz <= minpscz &
                              xx$b >= minb & xx$b <= maxb & xx$or >= minor & xx$or <= maxor & xx$f2 >= minf2 & xx$f2 <= maxf2])
}
get_path <- function(xx, minpheight=0.01, minpscz=1) {
  genes <- unique(xx$GENEID[xx$p.height<=minpheight & xx$p.scz <= minpscz])
  print(paste("Using" , length(genes), "genes"))
  enrichPathway(gene=genes, pvalueCutoff=0.05, readable=T)
}
##########################################
all_height_genes <- get_genes(allg, 0.001, 1)
all_height_path <- enrichPathway(gene=all_height_genes, pvalueCutoff=0.05, readable=T)
all_height_path
dotplot(all_height_path, showCat=15)
all_height_increase <- get_genes(allg, 0.001, 1, 0.04, 1E10, maxf2 = 0.05)
all_height_increase_path <- enrichPathway(gene=all_height_increase, pvalueCutoff=0.05, readable=T)
dotplot(all_height_increase_path, showCat=15)
all_height_decrease <- get_genes(allg, 0.001, 1, -10, -0.04, maxf2 = 0.05)
all_height_decrease_path <- enrichPathway(gene=all_height_decrease, pvalueCutoff=0.05, readable=T)
dotplot(all_height_decrease_path, showCat=15)

res_h <- compareCluster(list(increase=all_height_increase, decrease=all_height_decrease), fun="enrichPathway")
plot(res_h)
##########################################
## Going to do something a bit different
## Look at genes with a higher count of significant results
tb <- table(allg$GENEID[allp$p.height<0.001], cut(allg$b[allp$p.height<0.001], c(-10,-0.02,0.02,10)))

plot(tb[,1], tb[,3])
abline(0,1)
height_genes <- names(tb[tb>20])


tb <- table(allg$GENEID[allp$p.scz<0.001], cut(allg$or[allp$p.scz<0.001], c(0,0.97,0.99,1.01,1.03,100)))

allg <- allg[!is.na(f)]
ch <- allg$f > 0.5
allg$f2 <- allg$f
allg$f2[ch] <- 1.0-allg$f[ch]
allg$b[ch] <- -allg$b[ch]


plot(allg$f2[allg$p.height<0.001], allg$b[allg$p.height<0.001])

##########################################
coding_height <- get_path(allc, 0.001, 1)
dotplot(coding_height)

promoters_height <-  get_path(allp, 0.001, 1)
dotplot(promoters_height)
genes_height_promoters <- unique(allp$GENEID[allp$p.height<1E-3])
  
  height_promoters <- enrichPathway(gene=genes_height_promoters, pvalueCutoff=0.05, readable=T)
  head(as.data.frame(height_promoters))
  barplot(height_promoters, showCategory=12)  
  
  genes_height_introns <- unique(alli$GENEID[alli$p.height<1E-3])
  
  height_introns <- enrichPathway(gene=genes_height_introns, pvalueCutoff=0.05, readable=T)
  head(as.data.frame(height_introns))
  barplot(height_introns, showCategory=12)  
  
  gene_scz_coding <- unique(allc$GENEID[allc$p.scz<1E-3])
  scz_coding <- enrichPathway(gene=gene_scz_coding, pvalueCutoff=0.05, readable=T)
  head(as.data.frame(scz_coding))
  barplot(scz_coding, showCategory=12)
  



  
height_coding <- get_path(allc, 1E-3, 1)
barplot(height_coding)
height_promoters <- get_path(allp, 1E-3, 1)
barplot(height_promoters)
height_introns <- get_path(alli, 1E-3, 1)
barplot(height_introns)


scz_coding <- get_path(allc, 1, 1E-3)
barplot(scz_coding)
scz_introns <- get_path(alli, 1, 1E-3)
barplot(scz_introns)
scz_promoters <- get_path(allp, 1, 1E-3)
barplot(scz_promoters)

enrichMap(scz_introns, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
enrichMap(scz_promoters, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)

  enrichMap(height_promoters, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
  enrichMap(height_introns, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
  enrichMap(ydown, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
  
  cnetplot(height_promoters, categorySize="pvalue", foldChange=geneList)

  cnetplot(xdown, categorySize="pvalue", foldChange=geneList)
  cnetplot(yup, categorySize="pvalue", foldChange=geneList, fixed=TRUE)
  
  
#install.load.bioc("clusterProfiler")
  
## compare scz

ls <- list(scz_c = get_genes(allc, 1, 0.001), scz_p=get_genes(allp,1 , 0.001), scz_i =get_genes(alli, 1, 0.001))
res_s <- compareCluster(ls, fun="enrichPathway")
plot(res_s)  

## compare height
lh <- list(height_c = get_genes(allc, 0.001, 1), height_p=get_genes(allp,0.001 , 1), height_i =get_genes(alli, 0.001, 1))
res_h <- compareCluster(lh, fun="enrichPathway")
plot(res_h)  

la <- c(lh,ls)

res <- compareCluster(la, fun="enrichPathway")
plot(res)  

## Now just look up and down
## compare scz

s_dir <- list(scz_down = get_genes(allg, 1, 0.001,-1E10, 1E10, 0, 1 ), 
              scz_up=get_genes(allg ,1 , 0.001,-1E10, 1E10,1,1E10))
res_s_dir <- compareCluster(s_dir, fun="enrichPathway")
plot(res_s_dir)  
##

h_dir <- list(h_down = get_genes(allg, 0.001, 1, -1E10,    0, 0, 1E10), 
                  h_up=get_genes(allg ,0.001, 1,     0, 1E10, 0, 1E10))
res_h_dir <- compareCluster(h_dir, fun="enrichPathway")
plot(res_h_dir)

dir_all <- c(s_dir, h_dir)
res_dir_all <- compareCluster(dir_all, fun="enrichPathway")
plot(res_dir_all)
## compare height
lh <- list(height_c = get_genes(allc, 0.001, 1), height_p=get_genes(allp,0.001 , 1), height_i =get_genes(alli, 0.001, 1))
res_h <- compareCluster(lh, fun="enrichPathway")
plot(res_h)  

## Now what happens when I add cnvs
load(here("output", "remapped_cnv2.rda"))
gscz.del19 <-  GRanges(Rle(paste0("chr", scz.del19$CHR)), IRanges(start=scz.del19$start19, width=1), pval = scz.del19$z_pval, or = scz.del19$cmh_OR)
loc_gscz.del19 <- locateVariants(gscz.del19, txdb_hg19, AllVariants())
gscz.del19.loc <- data.table(data.frame(gscz.del19[loc_gscz.del19$QUERYID], GENEID= loc_gscz.del19$GENEID, LOCATION=loc_gscz.del19$LOCATION, TXID = loc_gscz.del19$TXID))

gscz.dup19 <-  GRanges(Rle(paste0("chr", scz.dup19$CHR)), IRanges(start=scz.dup19$start19, width=1), pval = scz.dup19$z_pval, or = scz.dup19$cmh_OR)
loc_gscz.dup19 <- locateVariants(gscz.dup19, txdb_hg19, AllVariants())
gscz.dup19.loc <- data.table(data.frame(gscz.dup19[loc_gscz.dup19$QUERYID], GENEID= loc_gscz.dup19$GENEID, LOCATION=loc_gscz.dup19$LOCATION, TXID = loc_gscz.dup19$TXID))

gheight19 <-  GRanges(Rle(paste0("chr", height_cnv19$CHR)), IRanges(start=height_cnv19$start19, width=1)
                      , p.height = height_cnv19$`Pvalue Height`, b = height_cnv19$`Beta Height`)
loc_gheight19 <- locateVariants(gheight19, txdb_hg19, AllVariants())
gheight19.loc <- data.table(data.frame(gheight19[loc_gheight19$QUERYID], GENEID= loc_gheight19$GENEID, LOCATION=loc_gheight19$LOCATION, TXID = loc_gheight19$TXID))

## Height genes

height_path_up <- enrichPathway(unique(gheight19.loc$GENEID[gheight19.loc$p.height<0.01 & gheight19.loc$b>0]),  pvalueCutoff=0.05, readable=T)
height_path_up
barplot(height_path)

height_path_down <- enrichPathway(unique(gheight19.loc$GENEID[gheight19.loc$p.height<0.01 & gheight19.loc$b<0]),  pvalueCutoff=0.05, readable=T)
height_path_down


barplot(height_path)

scz_del_path <- enrichPathway(unique(gscz.del19.loc$GENEID[gscz.del19.loc$pval<0.001]),  pvalueCutoff=0.05, readable=T)
scz_del_path
plot(scz_del_path)


scz_del_path_up <- enrichPathway(unique(gscz.del19.loc$GENEID[gscz.del19.loc$pval<0.001 & gscz.del19.loc$or >=1 ]),  pvalueCutoff=0.05, readable=T)
scz_del_path_up
plot(scz_del_path)


scz_del_path_down <- enrichPathway(unique(gscz.del19.loc$GENEID[gscz.del19.loc$pval<0.001 & gscz.del19.loc$or <1 ]),  pvalueCutoff=0.05, readable=T)
scz_del_path_down
plot(scz_del_path)

scz_dup_path <- enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001]),  pvalueCutoff=0.05, readable=T)
scz_dup_path
barplot(scz_dup_path )

scz_dup_path_up <- enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001 & gscz.dup19.loc$or > 1]),  pvalueCutoff=0.05, readable=T)
scz_dup_path_up

scz_dup_path_down <- enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001 & gscz.dup19.loc$or <1]),  pvalueCutoff=0.05, readable=T)
scz_dup_path_down
barplot(scz_dup_path_down )

scz_dup_path <- enrichPathway(unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001]),  pvalueCutoff=0.05, readable=T)
scz_dup_path
barplot(scz_dup_path )

scz_both_path <- enrichPathway(unique(c(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001],gscz.del19.loc$GENEID[gscz.del19.loc$pval<0.01])) ,  pvalueCutoff=0.05, readable=T)
scz_both_path
plot(scz_both_path)

lh <- list(dup = unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001]), 
                        height = unique(gheight19.loc$GENEID[gheight19.loc$p.height<0.001]),
                        del = unique(gscz.del19.loc$GENEID[gscz.del19.loc$pval<0.001])
                        )
                             
res_h <- compareCluster(lh, fun="enrichPathway")
plot(res_h)  

lh <- list(dup = unique(gscz.dup19.loc$GENEID[gscz.dup19.loc$pval<0.001]), 
           scz_down = get_genes(allg, 1, 0.001,-1E10, 1E10, 0, 1 ), 
           scz_up=get_genes(allg ,1 , 0.001,-1E10, 1E10,1,1E10),
           del = unique(gscz.del19.loc$GENEID[gscz.del19.loc$pval<0.001]))

res_h <- compareCluster(lh, fun="enrichPathway")
plot(res_h)  

