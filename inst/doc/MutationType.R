## ----echo = FALSE, eval=TRUE--------------------------------------------------

library(BiocStyle)
library(knitr)

## ----echo = FALSE, eval=TRUE,message=FALSE------------------------------------
library(MutationType)
library(VariantAnnotation)
library(BSgenome.Hsapiens.UCSC.hg38)

## ----echo = TRUE, eval=TRUE, results = 'asis'---------------------------------
vcffile <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")
vcffile <- readVcf(vcffile,"hg38")
setOfSNVs <- seekSNVs(vcffile, "hg38")

## ----echo = FALSE, eval=TRUE, results = 'asis'--------------------------------
knitr::kable(setOfSNVs[1:10,])

## ----echo = TRUE, eval=TRUE, results = 'asis'---------------------------------
setOfSNVsInfo <- seekMutAndSeq(setOfSNVs, contextLength = 3, Hsapiens)

## ----echo = FALSE, eval=TRUE, results = 'asis'--------------------------------
knitr::kable(setOfSNVsInfo[1:10,])

## ----echo = TRUE, eval=TRUE, results = 'asis', warning=FALSE------------------
results <- generalStat(setOfSNVsInfo)

## ----echo = FALSE, eval=TRUE, results = 'asis', warning=FALSE-----------------
knitr::kable(results[[1]], align = "c")
results[[2]]

## ----echo = FALSE, eval=TRUE, results = 'asis', warning=FALSE-----------------
knitr::kable(results[[3]][1:10,])
results[[4]]

## ----echo = FALSE, eval=TRUE--------------------------------------------------
sessionInfo()

