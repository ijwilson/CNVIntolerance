---
output:
  word_document: default
  html_document: default
---
# Genome Wide Association Studies, mutation intolerance and Copy Number Variants


## Description

This project aims to investigate constraints to CNV variation by
considering at the frequency of CNVs within genes that are thought to be
under strong and weaker selective pressures. Resources such as
phenoscanner (<http://www.phenoscanner.medschl.cam.ac.uk/phenoscanner>)
are cataloging the range of genotype phenotype associations on the human
genome. Databases of exonic CNV variation in large numbers of
individuals are available, for example ExAC
<http://exac.broadinstitute.org/>), as are attempts to measure how
tolerant genes are to mutations (<http://genic-intolerance.org>). We can
hypothesise that genes that are intolerant to mutation are also
intolerant to gene losses and duplications. This bioinformatic project
aims to integrate gene specific information about phenotype, gene
intolerance and CNV population frequency.

------------

## Databases

Note that the browser may not be the most appropriate way for us to
access these data. Downloading the data and working on
statistics/spreadsheet software could be better.

GWAS Catalog: [https://www.ebi.ac.uk/gwas/](https://www.ebi.ac.uk/gwas/)

*A resource for associations between SNPs and disease/phenotype. Others
are available, for example Grasp: 
[https://grasp.nhlbi.nih.gov/Overview.aspx](https://grasp.nhlbi.nih.gov/Overview.aspx)*

ExAc: <http://exac.broadinstitute.org/>

* Has information on CNV counts and intolerance scores for CNVs*

1000 genomes: <http://www.internationalgenome.org/>

*    Human genetic variation, both structural and at single sites.*

WTCCC1

data

<https://www.ebi.ac.uk/ega/studies/EGAS00000000017>

dbVar: <https://www.ncbi.nlm.nih.gov/dbvar>

* A database of human structural genetic variation*

**GWAS Central** [link](http://www.gwascentral.org/)

**GRASP**

[GRASP link](https://grasp.nhlbi.nih.gov/Search.aspx)

GRASP allows you to get records of all SNPs under a particular p-value for any trait of interest for a variety of studies that have been housed at the NCBI.  These can then be loaded into R to be compared with whatever other data you have

**CNV association database**

Just for panic disorder, for Japanese individuals.  [link](https://gwas.biosciencedbc.jp/cgi-bin/cccdb/ccc_top.cgi)

--------------


## GWAS Summary Statistics

[EBI list of available summary statistics](https://www.ebi.ac.uk/gwas/downloads/summary-statistics)

#### Psychiatric genetics consortium

Data are available to download from a few studies including Schizophrenia.  Contains data from a CNV case control study.

https://www.med.unc.edu/pgc/results-and-downloads

#### Height

GIANT consortium data

https://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#GWAS_Anthropometric_2014_Height

This contains information for exomes, and SNPs of the association between height and variants.

The most useful is perhaps GIANT 2-14 from the paper 

[Defining the role of common variation in the genomic and biological architecture of adult human height.](https://www.ncbi.nlm.nih.gov/pubmed/25282103?dopt=Citation) from 253,288 individuals.


#### Eczema

https://data.bris.ac.uk/data/dataset/28uchsdpmub118uex26ylacqm

#### Type 2 diabetes

[DIAGRAM consortium](http://diagram-consortium.org/downloads.html)

#### Inflamatory Bowel Disease

https://www.ibdgenetics.org/

Includes Crohn's disease & Ulcerative colitis meta-analysis.

#### Social Sciences

https://www.thessgac.org/data

Education attainment, cognitive performance, well being, neuroticism and reproductive behaviour.


#### CAD and Myocardial infarction

http://www.cardiogramplusc4d.org/data-downloads/

#### Age at Menarche

http://www.reprogen.org/data_download.html

#### Heart and Ageing

http://www.chargeconsortium.com/main/results



## Reading

**Recent CNV Associations**

**Height** <https://www.nature.com/articles/s41467-017-00556-x>

**Schizophrenia.** <https://www.nature.com/articles/ng.3725.pdf>

and a new unpublished preprint <https://www.biorxiv.org/content/early/2016/08/09/068593.full.pdf+html>

And a GWAS https://www.biorxiv.org/content/early/2016/08/09/068593

**GWAS Catalog** *MacArthur J, Bowler E, et al.*

[The new NHGRI-EBI Catalog of published genome-wide association studies
(GWAS
Catalog).](https://academic.oup.com/nar/article/45/D1/D896/2605751/The-new-NHGRI-EBI-Catalog-of-published-genome-wide)

Nucleic Acids Research, 2017, Vol. 45 (Database issue): D896-D901.




ExAc: https://www.nature.com/articles/nature19057

-------------

A couple more papers that might be useful

http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003709

*Talks about intolerance scores*

https://academic.oup.com/mbe/advance-article/doi/10.1093/molbev/msx305/4675175

*Very interesting paper that shows how a GWAS database can be used to inform hypotheses about the evolutionary genetics of cancer.*

