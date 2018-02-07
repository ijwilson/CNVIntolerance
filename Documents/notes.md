## Notes

* For quite a few of the studies, full summary statistics are available for the data.  This gives p values etc. for each SNP 
might be useful.  This is available from teh EBI and for other databases such as the GIANT consortium, which is particularly useful for BMI and height.

* Further summary statistics are available for other traits, 
    * [GIANT consortium data files](http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files).  
    * Neale LAB and biobank summary statistics [spreadsheet](https://sites.google.com/broadinstitute.org/ukbbgwasresults/home?authuser=0).  Also includes Heritability of >2000 traits from UK biobank
    

* 

Further ideas
==============

Look for SNPs in the promoter regions of genes that are significant
in the CNV GWAS.

```
promoters(mygenes, upstream=2000)
```

Should enable me to get the promoters.