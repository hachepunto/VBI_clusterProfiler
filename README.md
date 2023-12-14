# Introducción a clusterProfiler
## Viernes de Bioinformática
## Instituto Nacional de Medicina Genómica

Este material:
# https://github.com/hachepunto/VBI_clusterProfiler



[clusterProfiler](https://bioconductor.org/packages/release/bioc/html/clusterProfiler.html) in [Bioconductor](https://bioconductor.org/)


![clusterProfiler](https://ars.els-cdn.com/content/image/1-s2.0-S2666675821000667-fx1_lrg.jpg)


### Installation

To install this package, start R (version "4.3") and enter:

```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("clusterProfiler")
```

For older versions of R, please refer to the appropriate [Bioconductor release](https://bioconductor.org/about/release-announcements/).


Citation (from within R, enter `citation("clusterProfiler")`):

Wu T, Hu E, Xu S, Chen M, Guo P, Dai Z, Feng T, Zhou L, Tang W, Zhan L, Fu x, Liu S, Bo X, Yu G (2021). “clusterProfiler 4.0: A universal enrichment tool for interpreting omics data.” The Innovation, 2(3), 100141. [doi:10.1016/j.xinn.2021.100141](https://doi.org/10.1016/j.xinn.2021.100141).

Yu G, Wang L, Han Y, He Q (2012). “clusterProfiler: an R package for comparing biological themes among gene clusters.” OMICS: A Journal of Integrative Biology, 16(5), 284-287. [doi:10.1089/omi.2011.0118](https://doi.org/10.1089/omi.2011.0118).