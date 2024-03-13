df <- data.frame(id = c("rs7410291", "rs1234567"))
#Wrong input
test_that("wrong input", {
  expect_error(generalStat(df))
})

