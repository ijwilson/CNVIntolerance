



library(here)
source(here("R","prepare.R"))
### SNP intolerance scores
if (FALSE) {
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt",
                dropbox("Intolerance/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt"))
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt",
                dropbox("Intolerance/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"))  
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt",
                dropbox("Intolerance/fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt"))
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/README_fordist_cleaned_exac_r03_z_data_pLI_2016_01_13.txt",
                dropbox("Intolerance/README_fordist_cleaned_exac_r03_z_data_pLI_2016_01_13.txt")) 
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/README_fordist_cleaned_nonTCGA_z_data_pLI_2016_01_13.txt",
              dropbox("Intolerance/README_fordist_cleaned_nonTCGA_z_data_pLI_2016_01_13.txt"))
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/functional_gene_constraint/README_fordist_cleaned_nonpsych_z_data_pLI_2016_01_13.txt",
                dropbox("Intolerance/README_fordist_cleaned_nonpsych_z_data_pLI_2016_01_13.txt"))        
}
a1 <- read.table(dropbox("Intolerance/fordist_cleaned_exac_nonTCGA_z_pli_rec_null_data.txt"), header=TRUE)
a2 <- read.table(dropbox("Intolerance/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt"), header=TRUE)
a3 <- read.table(dropbox("Intolerance/fordist_cleaned_nonpsych_z_pli_rec_null_data.txt"), header=TRUE)

# Columns:
# transcript - Ensembl transcript ID
# gene - Ensembl gene symbol
# chr - chromosome
# n_exons - number of exons in the transcript
# cds_start - beginning of the transcript's coding sequence if + strand, end if - strand
# cds_end - end of the transcript's coding sequence if + strand, beginning if - strand
# bp - number of coding base pairs
# mu_syn - probability of a synonymous mutation across the transcript
# mu_mis - probability of a missense mutation across the transcript
# mu_lof - probability of a loss-of-function mutation across the transcript
# n_syn - number of rare (MAF < 0.1%) synonymous variants found in ExAC r0.3
# n_mis - number of rare (MAF < 0.1%) missense variants found in ExAC r0.3
# n_lof - number of rare (MAF < 0.1%) loss-of-function variants found in ExAC r0.3
# exp_syn - depth adjusted number of expected rare (MAF < 0.1%) synonymous variants
# exp_mis - depth adjusted number of expected rare (MAF < 0.1%) missense variants
# exp_lof - depth adjusted number of expected rare (MAF < 0.1%) loss-of-function variants
# syn_z - corrected synonymous Z score
# mis_z - corrected missense Z score
# lof_z - corrected loss-of-function Z score
# pLI - the probability of being loss-of-function intolerant (intolerant of both heterozygous and homozygous lof variants)
# pRec - the probability of being intolerant of homozygous, but not heterozygous lof variants
# pNull - the probability of being tolerant of both heterozygous and homozygous lof variants


## CNV Intolerance Scores
if (FALSE) {
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/cnv/exac-final-cnv.gene.scores071316",
                dropbox("Intolerance/exac-final-cnv.gene.scores071316"))
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/cnv/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed",
              dropbox("Intolerance/exac-final.autosome-1pct-sq60-qc-prot-coding.cnv.bed"))  
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/cnv/README.cnv_bed", 
              dropbox("Intolerance/README.cnv_bed"))
  download.file("ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/cnv/README.cnv_gene_scores",
              dropbox("Intolerance/README.cnv_gene_scores"))
}

