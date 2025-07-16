test_that("impute_knn ignore all na columns", {
  set.seed(123)

  mat <- matrix(rnorm(100), nrow = 10, ncol = 10)
  # NA cols
  mat[, 3] <- NA
  mat[, 7] <- NA
  # NA vals
  mat[2, 2] <- NA
  mat[4, 4] <- NA

  imputed <- impute_knn(obj = mat, k = 5, rowmax = 0.5, colmax = 0.9, rng.seed = 1234)

  testthat::expect_true(all(is.na(imputed[, 3])) && all(is.na(imputed[, 7])))
  testthat::expect_true(!is.na(imputed[2, 2]) && !is.na(imputed[4, 4]))
})
