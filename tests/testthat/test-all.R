
test_that("ae_sensitivity works", {
  p <- ae_sensitivity(savedir = NULL)
  expect_s3_class(p, "ggplot")
})

test_that("multi_weighting works", {
  p <- multi_weighting(savedir = NULL)
  expect_s3_class(p, "ggplot")
})

test_that("average_monotonicity works", {
  p <- average_monotonicity(savedir = NULL)
  expect_s3_class(p, "ggplot")
})

test_that("tsls weighting matches simulation", {
  TOL <- 1e-3
  set.seed(3534023)
  late <- c(1, -1)
  psc <- c(.05, .2, .35)
  prz <- generate_prz(prtreated = .4, prhigh = .6)
  w <- compute_tsls_weights(psc, prz, ZETA_ID)
  tsls_analytic <- sum(w*late)
  r <- tsls_simulation(psc, prz, late, ZETA_ID, n = 100000, numreps = 2000)
  diff <- abs(r$mean - tsls_analytic)
  expect_gt(TOL - diff, 0)
})
