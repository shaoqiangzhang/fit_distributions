# fit_distributions
For a gene expression array, fit the array to three distributions (Poisson, NB, ZINB) respectively, and use Chi-Square Goodness of Fit Test to find the best fit of the gene array

## install packages

```
install.packages("pscl")
install.packages("AER")
install.packages("fitdistrplus")
install.packages("osDesign")
library(BiocManager)
BiocManager::install("edgeR")
```
## call libraries 

```
library(edgeR)
library(pscl)
library(AER)   
library(fitdistrplus)
library(osDesign)
```
## source the three R files

