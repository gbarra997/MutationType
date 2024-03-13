library(VariantAnnotation)
wrongPath <- paste(system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation"), ".doc", sep="")
rightPath <- system.file("extdata", "chr22.vcf.gz", package="VariantAnnotation")


test_that("wrong type ref genome", {
  expect_error(seekSNVs(rightPath, 3))
})


test_that("wrong vcf file", {
  expect_error(seekSNVs(system.file(wrongPath, "chr22.vcf.gz", package="VariantAnnotation")
                        , "hg18"))
})


