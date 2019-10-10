# PICS
Probabilistic inference of ChIP-Seq using an empirical Bayes mixture model approach.

This is R implementation of statitical method proposed in my paper
Zhang, Xuekui, Gordon Robertson, Martin Krzywinski, Kaida Ning, Arnaud Droit, Steven Jones, and Raphael Gottardo. "PICS: probabilistic inference for ChIP‐seq." Biometrics 67, no. 1 (2011): 151-163.

Please cite this paper, if you used this R package in your research. Thanks!

To install this R package, please use the following R code
library(devtools)
install_github("ubcxzhang/MTPS")

Alternatively, this package can be installed from Bioconductor using  the following R code
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("PICS")
