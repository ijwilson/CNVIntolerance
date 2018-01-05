## Height
library(here)
source(here("R","prepare.R"))
install.load("data.table")

height <- fread(dropbox("height_CNV_association_41467_2017_556_MOESM2_ESM.csv"))


colnames(height)[5] <- "pvalue_BMI" 

## Simple of Manhattan plot of BMI
##plot(-log10(height$`Pvalue BMI`), col=height$CHR%%2+1)
## Simple plot of Height
##plot(-log10(height$pvalue_BMI), col=height$CHR%%2+1, main="Height")

install.load("qqman")
##qq(height$pvalue_BMI)

##acf(height$F_DUP)
######################################################

opar <- par(mfrow=c(4,1), mar=c(3,3,2,1))

height$P <- height$pvalue_BMI
manhattan(height, main="BMI")



height$P <- height$`Pvalue Weight`
manhattan(height, main="Weight")


height$P <- height$`Pvalue Height`
manhattan(height, main="Height")

height$P <- height$`Pvalue Waist-Hip ratio`
manhattan(height, main="Waist-Hip Ratio")
par(opar)

opar <- par(mfrow=c(2,2), mar=c(3,3,2,1))
qq(height$pvalue_BMI, main="BMI")
qq(height$`Pvalue Weight`, main="Weight")
qq(height$`Pvalue Height`, main="Height")
qq(height$`Pvalue Waist-Hip ratio`, main="Waist Hip Ratio")

par(opar)
