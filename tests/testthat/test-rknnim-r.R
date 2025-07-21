test_that("impute_knn_r matches impute.knn on small fixed matrix with k=3 and k=nrow(data)-1", {
  sample_matrix <- matrix(
    c(
      28, 90, 40, 40, 43,
      80, 70, NA, NA, NA,
      22, 79, 56, 48, 79,
      9, 78, 67, 3, 54,
      5, 14, 100, NA, 49,
      38, 56, NA, 41, 56,
      16, NA, 93, 100, 51,
      NA, 4, 5, 72, NA,
      98, 4, 66, 32, NA,
      86, 21, 47, 42, NA
    ),
    nrow = 10,
    ncol = 5,
    byrow = TRUE
  )
  rownames(sample_matrix) <- paste0("r", 1:10)
  colnames(sample_matrix) <- paste0("s", 1:5)

  set.seed(1234)
  orig <- impute::impute.knn(
    sample_matrix,
    k = 3,
    rowmax = 0.8,
    colmax = 1.0,
    maxp = nrow(sample_matrix)
  )$data
  custom <- r_impute(
    sample_matrix,
    k = 3,
    rowmax = 0.8,
    colmax = 1.0
  )
  expect_equal(orig, custom)

  k_max <- nrow(sample_matrix) - 1
  orig1 <- impute::impute.knn(
    sample_matrix,
    k = k_max,
    rowmax = 0.8,
    colmax = 1.0,
    maxp = nrow(sample_matrix)
  )$data

  custom1 <- r_impute(
    sample_matrix,
    k = k_max,
    rowmax = 0.8,
    colmax = 1.0
  )
  expect_equal(orig1, custom1)
  expect_error(
    r_impute(sample_matrix, k = nrow(sample_matrix), rowmax = 0.8, colmax = 1.0)
  )
})

test_that("impute_knn_r matches impute.knn on simulated sample_matrix with low missingness", {
  set.seed(123)
  sample_matrix <- sim_mat(n = 20, m = 10, nchr = 2, ngrp = 1, perc_NA = 0.1, perc_col_NA = 0.2)$input
  orig <- impute::impute.knn(
    sample_matrix,
    k = 5,
    rowmax = 0.8,
    colmax = 1.0,
    maxp = nrow(sample_matrix)
  )$data

  custom <- r_impute(
    sample_matrix,
    k = 5,
    rowmax = 0.8,
    colmax = 1.0
  )
  expect_equal(orig, custom)
})

test_that("impute_knn_r matches impute.knn on simulated sample_matrix with higher missingness", {
  set.seed(456)
  sample_matrix <- sim_mat(n = 20, m = 10, nchr = 2, ngrp = 1, perc_NA = 0.3, perc_col_NA = 0.5)$input
  orig <- impute::impute.knn(
    sample_matrix,
    k = 5,
    rowmax = 0.8,
    colmax = 1.0,
    maxp = nrow(sample_matrix)
  )$data

  custom <- r_impute(
    sample_matrix,
    k = 5,
    rowmax = 0.8,
    colmax = 1.0
  )
  expect_equal(orig, custom)
})
