# Two-layer hierarchical LASSO (TLHLASSO) model 

The last version is here: [LINK](https://github.com/YULEITSINGTAO/TLHLASSO)

## Introduction 
Two-layer hierarchical LASSO (2LHLASSO) instrumental variable methodology to integrate multi-omics predictor variables based on the flow of genetic information in the central dogma. Now we mainly apply the 2LHLASSO model on prostate cancer cohort to make the prediction of prostate cancer outcome. The Gleason score (GS), considered the most relevant clinical indicator for biochemical relapse-free survival (BRFS), was treated as the target outcome trait to be predicted, and the DNA methylation profiles were employed as the instrumental variables to estimate the expected transcription profiles. The analysis of the TCGA prostate cancer dataset showed that the new method of 2LHLASSO had higher predictability for outcome traits when compared to the linear model that simply combined multi-omics variables. 

## Installation 
```r
remotes::install_github("YULEITSINGTAO/TLHLASSO")
```
You will need the following R packages installed to run TLHLASSO:
- dplyr
- tidyr
- glmnet
- DT

## Using 2LHLASSO

## Running 2LHLASSO


## License
