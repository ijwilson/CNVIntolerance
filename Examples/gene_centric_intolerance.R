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

