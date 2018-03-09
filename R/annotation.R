
library(here)
source(here("R","prepare.R"))
install.load.bioc("VariantAnnotation", "AnnotationHub", "TxDb.Hsapiens.UCSC.hg19.knownGene")


if (!file.exists( here("output", "gwas_gr.rda"))) {
# Load height and scz
  load(here("output", "GWAS_height_summary.rda"))
  load(here("output", "GWAS_scz_summary.rda"))
  
  gheight <- GRanges(seqnames=Rle(paste0("chr", GWAS_height_summary$seqnames)), IRanges(start=GWAS_height_summary$pos, width=1, names=GWAS_height_summary$MarkerName)
                     , p = GWAS_height_summary$p, b=GWAS_height_summary$b, f=GWAS_height_summary$Freq.Allele1.HapMapCEU)
  genome(gheight) <- "hg19"
  gscz <- GRanges(seqnames=Rle(paste0("chr",GWAS_scz_summary$chr)), IRanges(start=GWAS_scz_summary$bp, width=1, names=GWAS_scz_summary$snpid)
                  , p = GWAS_scz_summary$p, or=GWAS_scz_summary$or)
  genome(gscz) <- "hg19"
  
  save(gheight, gscz, file = here("output", "gwas_gr.rda"))
  rm(GWAS_height_summary, GWAS_scz_summary )
  gc()
} else {
  load( here("output", "gwas_gr.rda"))
}


#OK then , want to set this up to look for coding variants, promoter variants and intronic variants 
#and then compre things.  Buit we could start with variants of estimated large effect

# I would probably also like to have a gene score - perhaps by separating out the different types
# of variation

table(b=abs(gheight$b>0.05), gheight$p < 0.001)
sig_height <- gheight[gheight$p<0.001 & abs(gheight$b) > ]

table(or=abs(gscz$or-1)>0.2, gscz$p < 0.001)

sig_scz <- gscz[gscz$p<0.001 & abs(gscz$or-1)>0.2]




#sig_height <- gheight[gheight$p<0.001 & abs(gheight$b>0.05)]
#sig_scz <- gscz[gscz$p<0.001 & abs(gscz$or-1)> 0.05]

loc_height_hg19 <- locateVariants(sig_height, txdb_hg19, AllVariants())
table(loc_height_hg19$LOCATION)
loc_scz_hg19    <- locateVariants(sig_scz, txdb_hg19, AllVariants())
table(loc_scz_hg19$LOCATION)


hub <- AnnotationHub()

txdb_hg19 <- TxDb.Hsapiens.UCSC.hg19.knownGene

head(seqlevels(txdb_hg19))


## Now get up and down for height and risk scz by geneid


b <- data.frame(b = sig_height[loc_height_hg19$QUERYID]$b, geneid = loc_height_hg19$GENEID)
or <- data.frame(or = sig_scz[loc_scz_hg19$QUERYID]$or, geneid = loc_scz_hg19$GENEID)

gupheight <- b$geneid[b$b>0.0]
gdownheight <- b$geneid[b$b<0.0]
gupscz <- or$geneid[or$or<1]
gdownscz <- or$geneid[or$or>1]




install.load.bioc("GO.db", "org.Hs.eg.db")
#x <- org.Hs.egGO
#mapped_genes <- mappedkeys(x)

#table(unique(gdownheight) %in% mapped_genes)
#table(unique(gupheight) %in% mapped_genes)

#ugupheight <- unique(gupheight)
#ugupheight <- ugupheight[!is.na(ugupheight)]
#ugupheight <- ugupheight[ugupheight %in% mapped_genes]
#xx <- as.list(x[ugupheight])
if(length(xx) > 0) {
  # Try the first one
  got <- xx[[1]]
  got[[1]][["GOID"]]
  got[[1]][["Ontology"]]
  got[[1]][["Evidence"]]
}


mapped_genes <- mappedkeys(unique(gupheight))




install.load.bioc("ReactomePA")

if (FALSE) {
 
  #' https://bioconductor.org/packages/release/bioc/vignettes/ReactomePA/inst/doc/ReactomePA.html
  data(geneList)
  de <- names(geneList)[abs(geneList) > 1.5]
  head(de)
  x <- enrichPathway(gene=de, pvalueCutoff=0.05, readable=T)
  head(as.data.frame(x))

  xup <- enrichPathway(gene=unique(gupheight), pvalueCutoff=0.05, readable=T)
  head(as.data.frame(xup))
  barplot(xup, showCategory=12)
  
  xdown <- enrichPathway(gene=unique(gdownheight), pvalueCutoff=0.05, readable=T)
  head(as.data.frame(xdown))
  barplot(xdown, showCategory=12)
  
  yup <- enrichPathway(gene=unique(gupscz), pvalueCutoff=0.05, readable=T)
  head(as.data.frame(yup))
  barplot(yup, showCategory=12)
  
  ydown <- enrichPathway(gene=unique(gdownscz), pvalueCutoff=0.05, readable=T)
  head(as.data.frame(ydown))
  barplot(ydown, showCategory=12)
 # dotplot(xup, showCategory=15)

  enrichMap(xup, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
  
  enrichMap(xdown, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
  enrichMap(yup, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
  enrichMap(ydown, layout=igraph::layout.kamada.kawai, vertex.label.cex = 1)
  
  cnetplot(xup, categorySize="pvalue", foldChange=geneList)

  cnetplot(xdown, categorySize="pvalue", foldChange=geneList)
  cnetplot(yup, categorySize="pvalue", foldChange=geneList, fixed=TRUE)
  
  
  install.load.bioc("clusterProfiler")
  data(gcSample)
  res <- compareCluster(gcSample, fun="enrichPathway")
  plot(res)
  
  require(clusterProfiler)
  data(gcSample)
  res <- compareCluster(gcSample, fun="enrichPathway")
  plot(res)

  l <- list(down_scz=unique(gdownscz), up_height = unique(gupheight), up_scz = unique(gupscz), down_height=unique(gdownheight))
  res <- compareCluster(l, fun="enrichPathway")
  plot(res)  
  res <- compareCluster(l, fun="enrichKEGG")
  plot(res)  
  l <- list(up_height = unique(gupheight), down_height=unique(gdownheight))
  res <- compareCluster(l, fun="enrichPathway")
  plot(res)  
}


cnetplot(x, categorySize="pvalue", foldChange=geneList)length(unique(gupheight))
length(unique(gdownheight))
length(unique(gupscz))
length(unique(gdownscz))

length(unique((gupheight[gupheight %in% gdownheight])))
length(unique((gupscz[gupscz %in% gdownscz])))
length(unique((gupheight[gupheight %in% gupscz])))
length(unique((gupheight[gupscz %in% gupheight])))
length(unique((gdownheight[gdownheight %in% gupscz])))
length(unique((gdownheight[gupscz %in% gdownheight])))

length(unique((gupheight[gupheight %in% gupscz])))
length(unique((gupheight[gupscz %in% gupheight])))
length(unique((gdownheight[gdownheight %in% gupscz])))
length(unique((gdownheight[gupscz %in% gdownheight])))
## Look at coding
b <- data.frame(b = sig_height[loc_height_hg19$QUERYID]$b, geneid = loc_height_hg19$GENEID)[loc_height_hg19$LOCATION=="coding",]
or <- data.frame(or = sig_scz[loc_scz_hg19$QUERYID]$or, geneid = loc_scz_hg19$GENEID)[loc_scz_hg19$LOCATION=="coding",]

gupheight <- b$geneid[b$b>0.0]
gdownheight <- b$geneid[b$b<0.0]

max(table(gdownscz))
gupscz <- or$geneid[or$or<1]
gdownscz <- or$geneid[or$or>1]



length(unique(gupheight))
length(unique(gdownheight))
length(unique(gupscz))
length(unique(gdownscz))

length(unique((gupheight[gupheight %in% gdownheight])))
length(unique((gupscz[gupscz %in% gdownscz])))
length(unique((gupheight[gupheight %in% gupscz])))
length(unique((gupheight[gupscz %in% gupheight])))
length(unique((gdownheight[gdownheight %in% gupscz])))
length(unique((gdownheight[gupscz %in% gdownheight])))

length(unique((gupheight[gupheight %in% gupscz])))
length(unique((gupheight[gupscz %in% gupheight])))
length(unique((gdownheight[gdownheight %in% gupscz])))
length(unique((gdownheight[gupscz %in% gdownheight])))

length(unique((gupscz[gupscz %in% gdownscz])))

install.load.bioc("org.Hs.eg.db")
