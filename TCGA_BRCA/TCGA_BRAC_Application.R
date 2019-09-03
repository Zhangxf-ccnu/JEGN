# Run JEGN on TCGA breast cancer datasets

rm(list=ls())
library(JEGN)

data("TCGA.BRCA")

# run JEGN using JEGN
TCGA.BRCA.JEGN1 = JEGN(TCGA.BRCA$X, lambda =  0.95, alpha = 0.4, model = "nonparanormal")

# run JEGN using JEGN.admm
TCGA.BRCA.JEGN2 = JEGN.admm(TCGA.BRCA$Sigma, lambda =  0.95, alpha = 0.4)


