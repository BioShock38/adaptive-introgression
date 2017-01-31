library(testthat)
Rcpp::sourceCpp("imputationUtils.cpp")

test_that("missing values have been taken into account in median estimation", {
  test <- matrix(c(9, 2, 1, 9, NA, 2, 2, 2, 0, 2, 1, 9, NA, 2, 2, 2), nrow = 2)
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

test_that("missing values have been taken into account in median estimation", {
  test <- matrix(c(9, NA, 9, 9, 9, 9, 9, NA, 9, 9, 9, 9, 9, NA, 9, NA), nrow = 2)
  res <- rowMedian_cpp(test)
  for (k in 1:nrow(test)){
    testthat::expect_equal(is.na(res[k]), TRUE)
  }
})

test_that("the median-based imputation is correct", {
  test <- matrix(c(9, 2, 1, 9, NA, 2, 2, 2, 0, 2, 1, 9, NA, 2, 2, 2), nrow = 2)
  xpctd <- matrix(c(1, 2, 1, 2, 1, 2, 2, 2, 0, 2, 1, 2, 1, 2, 2, 2), nrow = 2)
  res <- impute_geno(test)
  testthat::expect_equal(res$x, xpctd)
})

test_that("rows with missing values only are discarded", {
  test <- matrix(c(9, NA, 9, 9, 9, 9, 9, NA, 9, 9, 9, 9, 9, NA, 9, NA), nrow = 2)
  skip <- check_row(test, 0);
  testthat::expect_equal(skip, 1)
})

test_that("rows with only 2's or NA's are discarded", {
  test <- matrix(c(9, 2, 9, 2, NA, 2, 9, 2, 9, 2, 9, 2, NA, NA, NA, NA), nrow = 2)
  skip <- check_row(test, 1);
  testthat::expect_equal(skip, 1)
})

test_that("rows with different finite values are kept", {
  test <- matrix(c(9, 2, 9, 2, NA, 0, 9, NA, 9, 9, 9, 9, NA, NA, NA, NA), nrow = 2)
  skip <- check_row(test, 1);
  testthat::expect_equal(skip, 0)
})

test_that("medians for each SNP and each population are computed correctly", {
  test <- matrix(c(1, 2, 1, 9, 1, 2, 9, 2, 0, 2, 0, 9, NA, 2, 0, 2), nrow = 2)  
  lab <- as.vector(c(1, 1, 1, 1, 2, 2, 2, 2))
  pop <- as.vector(c(1, 2))
  xpctd <- matrix(c(1, 2, 0, 2), nrow = 2)
  res <- median_per_pop(test, lab, pop)
  testthat::expect_equal(res, xpctd)
})