
#VPATH=output

#%.rda: %.R
#		R CMD BATCH $<

%.pdf: %.R
		R CMD BATCH $<

all: alldata  #docs

#docs figs

alldata:  output/genesGR.rda  output/height_cnv.rda output/scz_cnv.rda output/GWAS_height_summary.rda \
          output/GWAS_scz_summary.rda output/snp_match.rda output/both_annotated.rda output/snp_NHI.rda

  
## liftover will only work on linux machines        
liftover: output/remapped_cnv2.rda output/hg18ToHg19.over.chain.gz output/scz.all_cnv_hg18.gene.bed 

#output/snpinfo.rda output/exacscores.rda output/dbCNV.rda  \
          output/exacbed.rda output/ignore.regions.rda output/individuals.rda output/rare.rda output/relationships.rda

docs: alldata  concordant_discordant_pathways.html
        
figs: Figs/fig1.pdf alldata

## intermediate data files

output/genesGR.rda: R/getGenes.R
	R CMD BATCH --no-save R/getGenes.R

output/height_cnv.rda: R/height.R
	R CMD BATCH --no-save R/height.R

output/scz_cnv.rda: R/cnv_scz.R
	R CMD BATCH --no-save R/cnv_scz.R
	
output/snp_NHI.rda: R/match_NHI.R
	R CMD BATCH --no-save R/match_NHI.R
	
output/GWAS_height_summary.rda: R/GWAS_height_summary.R
	R CMD BATCH R/GWAS_height_summary.R
	
output/GWAS_scz_summary.rda: R/GWAS_scz_summary.R
	R CMD BATCH --no-save R/GWAS_scz_summary.R 
	
output/snp_match.rda: R/match_snps.R 
	R CMD BATCH --no-save R/match_snps.R 
	
output/both_annotated.rda: R/match_SCZ_height_snps_annotate.R output/snp_match.rda
	R CMD BATCH --no-save R/match_SCZ_height_snps_annotate.R

concordant_discordant_pathways.html: R/concordant_discordant_pathways.R
	Rscript -e 'library(rmarkdown); library(here); render(here("R","concordant_discordant_pathways.R"), output_format="html_document", output_dir=here("output"))'

output/gthree.loc.rda: R/match_cnv.R output/remapped_cnv2.rda
	R CMD BATCH --no-save R/match_cnv.R 
	

## Files and stuff for liftover
output/hg18ToHg19.over.chain.gz:
	wget http://hgdownload.soe.ucsc.edu/goldenPath/hg18/liftOver/hg18ToHg19.over.chain.gz
	mv -f hg18ToHg19.over.chain.gz output

output/scz.all_cnv_hg18.gene.bed: R/prepare_liftover.R
	R CMD BATCH --no-save R/prepare_liftover.R

output/scz.all_cnv_hg19.gene.bed: output/scz.all_cnv_hg18.gene.bed R/prepare_liftover.R output/hg18ToHg19.over.chain.gz
	liftOver output/height_cnv_hg18.bed output/hg18ToHg19.over.chain.gz output/height_cnv_hg19.bed output/unmapped_height
	liftOver output/scz.dup_cnv_hg18.bed output/hg18ToHg19.over.chain.gz output/scz.dup_cnv_hg19.bed output/unmapped_scz.dup
	liftOver output/scz.del_cnv_hg18.bed output/hg18ToHg19.over.chain.gz output/scz.del_cnv_hg19.bed output/unmapped_scz.del
	liftOver output/scz.all_cnv_hg18.gene.bed output/hg18ToHg19.over.chain.gz output/scz.all_cnv_hg19.gene.bed output/unmapped_scz_all.dup
	liftOver output/scz.del_cnv_hg18.gene.bed output/hg18ToHg19.over.chain.gz output/scz.del_cnv_hg19.gene.bed output/unmapped_scz_del.dup

output/remapped_cnv2.rda: R/liftover.R output/scz.all_cnv_hg19.gene.bed 
	R CMD BATCH --no-save R/liftover.R

	
output/exacbed.rda:  R/exacbed.R output/genesGR.rda 
	R CMD BATCH R/exacbed.R

output/ignore.regions.rda: R/excludeRegions.R
	R CMD BATCH R/excludeRegions.R

output/rare.rda:	R/rare.R output/dbCNV.rda output/filteredqs.rda
	R CMD BATCH R/rare.R

output/exacscores.rda: R/exacscores.R output/exacbed.rda output/rare.rda output/filteredqs.rda
	R CMD BATCH R/exacscores.R

output/snpinfo.rda: R/snpinfo.R
	R CMD BATCH R/snpinfo.R

output/relationships.rda: R/relationships.R output/individuals.rda
	R CMD BATCH R/relationships.R
	
output/A_1358_g20.rda: R/A_1358_g20.R R/genotypes.R
	R CMD BATCH R/A_1358_g20.R

## Figures

Figs/fig1.pdf: Figs/fig1.R output/rare.rda output/individuals.rda output/exacbed.rda
	R CMD BATCH Figs/fig1.R	

## Documents



Reports/check_individuals.html: Reports/check_individuals.Rmd R/start.R
	Rscript -e 'library(rmarkdown); library(here); render(here("Reports","check_individuals.Rmd"), output_format="html_document")'

Reports/2017_06_22_results_paragraph.html: Reports/2017_06_22_results_paragraph.Rmd output/rare.rda 
	Rscript -e 'library(rmarkdown); library(here); render(here("Reports","2017_06_22_results_paragraph.Rmd"), output_format="html_document")'
	Rscript -e 'library(rmarkdown); library(here); render(here("Reports","2017_06_22_results_paragraph.Rmd"), output_format="word_document")'

Reports/2017-08-01_results_paragraph2.html: Reports/2017-08-01_results_paragraph2.Rmd output/rare.rda 
	Rscript -e 'library(rmarkdown); library(here); render(here("Reports","2017-08-01_results_paragraph2.Rmd"), output_format="html_document")'
	Rscript -e 'library(rmarkdown); library(here); render(here("Reports","2017-08-01_results_paragraph2.Rmd"), output_format="word_document")'

Reports/associated_cnvs.html: Reports/associated_cnvs.Rmd
	Rscript -e 'library(rmarkdown); library(here); render(here("Reports","associated_cnvs.Rmd"), output_format="html_document")'

Reports/2017-06-29_Relationships.html: Reports/2017-06-29_Relationships.Rmd output/relationships.rda output/individuals.rda 
	Rscript -e 'library(rmarkdown); library(here); render(here("Reports","2017-06-29_Relationships.Rmd"), output_format="html_document")'





## Housekeeping

clean:
	rm *.Rout cytohg19.tar.gz humanChrLens.txt Reports/*.log R/*.Rout
