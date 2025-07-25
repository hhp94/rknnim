test_that("impute_knn ignore all na rows", {
  set.seed(123)

  mat <- matrix(rnorm(100), nrow = 10, ncol = 10, dimnames = list(1:10, 1:10))
  # Introduce all NA rows
  mat[3, ] <- NA
  mat[7, ] <- NA
  # Introduce some specific NA values that should be imputed
  mat[2, 2] <- NA
  mat[4, 4] <- NA

  # Call impute_knn. For all-NA rows to be ignored, their NA proportion (100%)
  # must be >= colmax. Setting colmax to 1 ensures this.
  imputed <- impute_knn(obj = mat, k = 5, rowmax = 0.5, colmax = 1, method = "euclidean", cores = 1)

  # Expect that the all-NA rows remain all NA
  testthat::expect_true(all(is.na(imputed[3, ])) && all(is.na(imputed[7, ])))
  # Expect that the specific NA values are imputed
  testthat::expect_true(!is.na(imputed[2, 2]) && !is.na(imputed[4, 4]))
})

test_that("Exactly replicate impute::impute.knn", {
  data("khanmiss1")
  if (rlang::is_installed("impute")) {
    r1 <- knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn")
    r2 <- t(impute::impute.knn(khanmiss1, k = 3, rowmax = 1, maxp = nrow(khanmiss1))$data)

    expect_equal(r1, r2)
  } else {
    expect_no_error(knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn"))
  }
})
