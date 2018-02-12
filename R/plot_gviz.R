#` # Plots using Gviz
#`  [Vignette](http://www.bioconductor.org/packages/release/bioc/vignettes/ggbio/inst/doc/ggbio.pdf)
#`
#` [Talk](https://jmonlong.github.io/MonBUG17_Gviz/#1)
#`

#` ## Data
#` I want the SCZ CNV data by Gene and the overall 
#+ message=FALSE
library(here)
source(here("R", "prepare.R"))
install.load.bioc("GenomicRanges", "Gviz") 

chr <- "chr22"

b22q11.2 <- GRanges("chr22", IRanges(18660553,21455556))
genome(b22q11.2) <- "hg19"

gen <- genome(b22q11.2)
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = gen, chromosome = chr)

## Get the gene models
load(here("output", "genesGR.rda"))
seqnames(genesGR19) <- paste("chr", seqnames(genesGR19), sep="")
geneModel <- subsetByOverlaps(genesGR19, GRanges("22", IRanges(18660553,21455556)))  
grtrack <- GeneRegionTrack(geneModels, genome = gen, chromosome = chr, name = "Gene Model")
## Get the CNVs for this region

load(here("output", "remapped_cnv.rda"))

gheight_cnv <- GRanges(paste("chr",height_cnv19$CHR,sep=""), IRanges(height_cnv19$BP, width=1), 
                       p_height = height_cnv19$`Pvalue Height`, beta_height= height_cnv19$`Beta Height`, 
                       p_bmi = height_cnv19$`Pvalue BMI`, beta_bmi = height_cnv19$`Beta BMI`, 
                       fdel=height_cnv19$F_DEL, fdup = height_cnv19$F_DUP)
genome(gheight_cnv) <- "hg19"

gheight_cnv <- subsetByOverlaps(gheight_cnv, b22q11.2)
dtrack_height <- DataTrack(data = gheight_cnv$beta_height, start = start(gheight_cnv),
                        end = end(gheight_cnv), chromosome = chr, genome = gen,
                        name = "beta_height")


gSCZ_del_cnv <- GRanges(paste("chr",scz.del19$CHR,sep=""), IRanges(scz.del19$BP, width=1), 
                       p = scz.del19$cmh_pval, or= scz.del19$cmh_OR, 
                       unaff=scz.del19$UNAFF, aff = scz.del19$AFF)


gSCZ_dup_cnv <- GRanges(paste("chr",scz.dup19$CHR,sep=""), IRanges(scz.dup19$BP, width=1), 
                        p = scz.dup19$cmh_pval, or= scz.dup19$cmh_OR, 
                        unaff=scz.dup19$UNAFF, aff = scz.dup19$AFF)

gSCZ_del_cnv <- subsetByOverlaps(gSCZ_del_cnv, b22q11.2)
gSCZ_dup_cnv <- subsetByOverlaps(gSCZ_dup_cnv, b22q11.2)

dtrack_scz_del <- DataTrack(data = gSCZ_del_cnv$or, start = start(gSCZ_del_cnv),
                           end = end(gSCZ_del_cnv), chromosome = chr, genome = gen,
                           name = "OR Del", col="red")

dtrack_scz_dup <- DataTrack(data = gSCZ_dup_cnv$or, start = start(gSCZ_dup_cnv),
                            end = end(gSCZ_dup_cnv), chromosome = chr, genome = gen,
                            name = "OR Dup", col="green")


plotTracks(list( itrack, gtrack, dtrack_height, dtrack_scz_dup, dtrack_scz_del))
