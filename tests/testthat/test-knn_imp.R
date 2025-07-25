test_that("Exactly replicate impute::impute.knn", {
  data("khanmiss1")
  r1 <- knn_imp(t(khanmiss1), k = 3, rowmax = 1, method = "impute.knn")
  r2 <- t(impute::impute.knn(khanmiss1, k = 3, rowmax = 1, maxp = nrow(khanmiss1))$data)

  expect_equal(r1, r2)
})
