# VizCNV
This is a shiny app for chromosomal copy number variant analysis.
It can parse the vcf file with SV calls, visualize CNV and B-allele frequency 
and genetic phasing information interactively.

## Prerequisites
R version >= 3.6.3
Following R libraries are required:

Shifting level models based segmentation is performed using [SLMSuite](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1734-5). 

Launch app on local:
```r
shiny::runGitHub(repo = "cluhaowie/VizCNV")
```
![fig](/www/screenshot1.png)

Upload the required file from local filesysterm:
![gif1](docs/uploadFile.gif)

If launch the app on cloud or on server, input file need to be upload due to access restriction.

Visualize the CNV calls in table format, read depth plot and B-allele frequency together:
![gif1](docs/view.gif)

The app require read depth file as the input:
A output from mosedepth can be used



