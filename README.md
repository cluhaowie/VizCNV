README
================

<img src="docs/VizCNV_logo.svg" width="150" height="150" />

# VizCNV - An integrated platform for concurrent phased BAF and CNV analysis

------------------------------------------------------------------------

This is a shiny app for chromosomal copy number variant analysis. It can
parse the vcf file with SV calls, visualize CNV and B-allele frequency
and genetic phasing information interactively. The VizCNV is still under
active development, suggestions is appreciated.

[Download current
version](https://github.com/cluhaowie/VizCNV/releases/tag/v4.2.3)

[Pull docker image](https://hub.docker.com/r/duclare123/vizcnv_dev)

[Q&A](https://github.com/cluhaowie/VizCNV/issues)

## Prerequisites

R version \>= 4.2 Following R libraries are required: Shifting level
models based segmentation is performed using
[SLMSuite](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-017-1734-5).

## Citation

------------------------------------------------------------------------

Please cite the following article if you use VizCNV in your research:

- Du, H., Jolly, A., Grochowski, C. M., Yuan, B., Dawood, M.,
  Jhangiani, S. N., Li, H., Muzny, D., Fatih, J. M., Coban-Akdemir, Z.,
  Carlin, M. E., Scheuerle, A. E., Witzl, K., Posey, J. E., Pendleton,
  M., Harrington, E., Juul, S., Hastings, P. J., Bi, W., … Liu, P.
  (2022). The multiple *de novo* copy number variant (M*dn*CNV)
  phenomenon presents with peri-zygotic DNA mutational signatures and
  multilocus pathogenic variation. *Genome Medicine*, *14*(1), 122.

The manuscript of VizCNV is WIP!
