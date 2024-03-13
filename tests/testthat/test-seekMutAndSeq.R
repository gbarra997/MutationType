library(BSgenome.Hsapiens.UCSC.hg38)

df <- data.frame(id = c("rs7410291", "rs1234567"),chr = c("chr22", "chrX"),pos = c("50300078", "100200300"),start = c("0", "0"),end = c("0", "0"),ref = c("A", "C"),alt = c("G", "T"),stra = c("*", "*") )
wrongDf <- data.frame(id = c("rs7410291", "rs1234567"),chr = c("chr22", "chrX"),pos = c("50300078", "100200300"),start = c("0", "0"),end = c("0", "0"),stra = c("*", "*") )

t <- c(1,2)

#Wrong reference type
test_that("reference genome wrong", {
  expect_error(seekMutAndSeq(df, 3, t))
})

#Wrong reference type
test_that("wrong context len", {
  expect_error(seekMutAndSeq(df, 4, Hsapiens))
})



#Wrong df input
test_that("wrong df inpuit", {
  expect_error(seekMutAndSeq(wrongDf, 3, Hsapiens))
})
