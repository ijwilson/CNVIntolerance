### Height
library(here)
source(here("R","prepare.R"))
install.load("data.table")
install.load.bioc("GenomicRanges")


load(here("output", "height_cnv19.rda"))

height_cnv19 <- height_cnv19[,c(1,2,4,5,8:15)]
colnames(height_cnv19)[3:4] <- c("fdel", "fdup")


gheight_cnv <- GRanges(height_cnv19$CHR, IRanges(height_cnv19$BP, width=1), 
                       p_height = height_cnv19$`Pvalue Height`, beta_height= height_cnv19$`Beta Height`, 
                       p_bmi = height_cnv19$`Pvalue BMI`, beta_bmi = height_cnv19$`Beta BMI`, 
                       fdel=height_cnv19$fdel, fdup = height_cnv19$fdup)
genome(gheight_cnv) <- "hg19"
rm(height_cnv19)  ##keep the workspace tidy


load(here("output", "genesGR.rda"))
rm(genesGR18)

genesGR19 <- genesGR19[-1675] ## remove 2nd NPHP3-ACAD11
genesGR19 <- subsetByOverlaps(genesGR19, gheight_cnv)
o <- findOverlaps(gheight_cnv, genesGR19)

df <- data.frame(gene = names(genesGR19)[subjectHits(o)],
                  as.data.frame(gheight_cnv[queryHits(o)])
                                                )
df$strand <- NULL
df$width <- NULL



prop_sig_height <- tapply(df$p_height, df$gene, function(x) mean(x<0.001)) 
mean_sig_height_log10 <- tapply(df$p_height, df$gene, function(x) mean(log10(x))) 
m_beta_height <- tapply(df$beta_height, df$gene, mean) 
abs_beta_height <- tapply(df$beta_height, df$gene, function(x) mean(abs(x))) 
ch_del_height <- as.vector(by(df, df$gene, function(xx) mean(-xx$fdel*xx$beta_height/100)))
ch_dup_height <- as.vector(by(df, df$gene, function(xx) mean(xx$fdup*xx$beta_height/100)))

m_fdel <- tapply(df$fdel/100, df$gene, mean) 
m_fdup <- tapply(df$fdup/100, df$gene, mean) 

s_fdel <- tapply(df$fdel, df$gene, sd) 
s_fdup <- tapply(df$fdup, df$gene, sd) 

n<- tapply(df$fdel, df$gene, length)
ch_height <- ch_del_height+ch_dup_height

m <- match(names(genesGR19), names(s_fdup))


res <- data.frame(as.data.frame(genesGR19), mean_fdel=m_fdel[m],
                   mean_fdup=m_fdup[m], mean_beta_height=m_beta_height[m], 
                   proportion_sig_height=prop_sig_height, mean_sig_height_log10 = mean_sig_height_log10[m],
                   n=n[m]  )

plot(res$mean_beta_height, res$mean_fdel)

plot(res$mean_fdel, res$mean_fdup)

height_cnv_by_gene_hg19 <- res

write.csv(height_cnv_by_gene_hg19, file = "height_cnv_by_gene_hg19.csv", quote=FALSE)
