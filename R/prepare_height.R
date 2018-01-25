
library(here)
source(here("R","prepare.R"))
install.load("data.table")
#
height <- fread(dropbox("height_CNV_association_41467_2017_556_MOESM2_ESM.csv"))
height_hg19 <- fread(dropbox("cnv_hg19.bed"))
height_hg19[,V3:=NULL]
colnames(height_hg19) <- c("CHR", "BP", "ID")
height$ID <- paste("ID", height$CHR, height$BP, sep="_")
setkey(height, ID)
setkey(height_hg19, ID)

h <- height_hg19[height]
h[, ID := NULL]
h[,CHR := i.CHR]
h[,i.CHR := NULL]
h[,i.BP := NULL]


height_raw <- fread(dropbox("cnv_catalogue_CNVs_2018-01-24_03_51_11.csv"))

#check <- cbind(height_raw$BP, height$BP, height_raw$BP -  height$BP )
#check[1:10,]

h3 <- merge(height_raw, h, by=c("CHR", "BP"))
height_cnv19 <- h3

save(height_cnv19, file=here("output", "height_cnv19.rda"))
