
#VPATH=output

#%.rda: %.R
#		R CMD BATCH $<

%.pdf: %.R
		R CMD BATCH $<

all: alldata 
#docs figs
alldata:  output/genesGR.rda 

#output/snpinfo.rda output/exacscores.rda output/dbCNV.rda output/genesGR.rda \
          output/exacbed.rda output/ignore.regions.rda output/individuals.rda output/rare.rda output/relationships.rda

docs: alldata Reports/associated_cnvs.html Reports/check_individuals.html \
        Reports/2017_06_22_results_paragraph.html Reports/2017-06-29_Relationships.html \
        Reports/2017-08-01_results_paragraph2.html 
        
figs: Figs/fig1.pdf alldata

## intermediate data files

output/genesGR.rda: R/getGenes.R
	R CMD BATCH R/getGenes.R

output/height_cnv.rda: R/height.R
	R CMD BATCH R/height.R



output/exacbed.rda:  R/exacbed.R output/genesGR.rda 
	R CMD BATCH R/exacbed.R

output/ignore.regions.rda: R/excludeRegions.R
	R CMD BATCH R/excludeRegions.R

output/individuals.rda:	R/alternative_individuals.R
	R CMD BATCH 	R/alternative_individuals.R

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
