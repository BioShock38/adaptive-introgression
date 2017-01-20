library(testthat)
Rcpp::sourceCpp("aiUtils.cpp")

test_that("missing values have been taken into account in median estimation", {
  test <- matrix(c(9, 2, 1 , 9, NA, 2, 2, 2, 0, 2, 1 , 9, NA, 2, 2, 2), nrow = 2)
  res <- rowMedian_cpp(test)
  med <- array(0, dim = nrow(test))
  for (k in 1:nrow(test)){
    rowtest <- test[k, ]
    rowtest <- rowtest[!is.na(rowtest)]
    rowtest <- rowtest[rowtest != 9]
    med[k] <- median(rowtest)
    testthat::expect_equal(res[k], med[k])
  }
})

test_that("the median-based imputation is correct", {
  test <- matrix(c(9, 2, 1 , 9, NA, 2, 2, 2, 0, 2, 1 , 9, NA, 2, 2, 2), nrow = 2)
  xpctd <- matrix(c(1, 2, 1 , 2, 1, 2, 2, 2, 0, 2, 1 , 2, 1, 2, 2, 2), nrow = 2)
  res <- impute_geno(test)
  testthat::expect_equal(res, xpctd)
})