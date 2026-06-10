test_that("psmqte returns QTE and QTT estimates", {
  dat <- simulate_psmqte_data(250, seed = 10)
  fit <- psmqte(
    Y ~ X1 + X2,
    treatment = "D",
    data = dat,
    taus = c(0.25, 0.5, 0.75),
    M = c(1, 2),
    se = FALSE
  )
  expect_s3_class(fit, "psmqte")
  expect_equal(nrow(fit$estimates), 6)
  expect_true(all(is.finite(fit$estimates$qte)))
  expect_true(all(is.finite(fit$estimates$qtt)))
})

test_that("psmate returns ATE and ATT estimates", {
  dat <- simulate_psmqte_data(250, seed = 11)
  fit <- psmate(
    Y ~ X1 + X2,
    treatment = "D",
    data = dat,
    M = c(1, 2),
    se = TRUE
  )
  expect_s3_class(fit, "psmate")
  expect_equal(nrow(fit$estimates), 2)
  expect_true(all(is.finite(fit$estimates$ate)))
  expect_true(all(is.finite(fit$estimates$att)))
})

test_that("matching utility reports reuse counts", {
  dat <- simulate_psmqte_data(100, seed = 12)
  mat <- psm_match(dat$p_true, dat$D, M = 2)
  expect_equal(nrow(mat$match_matrix), 100)
  expect_equal(ncol(mat$match_matrix), 2)
  expect_equal(sum(mat$K_all), 200)
})

test_that("safe default example data loads", {
  dat <- psmqte_example_data()
  expect_equal(nrow(dat), 1000)
  expect_true(all(c("Y", "D", "X1", "X2", "X3", "age", "size") %in%
                    names(dat)))
})

test_that("fresh simulated example data loads", {
  dat <- psmqte_example_data("simulated", n = 250, seed = 2)
  expect_equal(nrow(dat), 250)
  expect_true(all(c("Y", "D", "X1", "X2") %in% names(dat)))
})

test_that("optional NHEFS example data loads when causaldata is installed", {
  skip_if_not_installed("causaldata")
  dat <- psmqte_example_data("nhefs")
  expect_true(nrow(dat) > 1000)
  expect_true(all(c("Y", "D", "age", "sex", "race", "school") %in% names(dat)))
})
