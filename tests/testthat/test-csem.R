items <- attitude_data[, 1:10]

test_that("all supported models fit and give finite, positive CSEM", {
  for (m in c("graded", "gpcm", "pcm", "grsm", "rsm")) {
    f <- suppressWarnings(csem_estimate(items, model = m, theta_points = 41))
    expect_s3_class(f, "csem_fit")
    a <- csem_analytical(f)
    expect_true(all(is.finite(a$csem) & a$csem > 0))
    expect_equal(a$csem, 1 / sqrt(f$test_info))
  }
})

test_that("aliases map to the same model", {
  f1 <- csem_estimate(items, model = "grm", theta_points = 21)
  expect_identical(f1$model, "graded")
  expect_error(csem_estimate(items, model = "nominal"), "must be one of")
})

test_that("standardization rescales PCM information correctly", {
  f_std <- csem_estimate(items, model = "pcm", theta_points = 21)
  f_raw <- csem_estimate(items, model = "pcm", theta_points = 21,
                         standardize = FALSE)
  s <- sqrt(f_std$latent_var)
  z <- f_std$theta_grid
  raw_at_sz <- as.numeric(mirt::testinfo(f_raw$fit, matrix(s * z)))
  expect_equal(f_std$test_info, s^2 * raw_at_sz, tolerance = 1e-8)
  # GRM: latent variance fixed at 1, so both metrics coincide
  g <- csem_estimate(items, theta_points = 21)
  expect_equal(g$scale_sd, 1)
})

test_that("information-based marginal reliability averages conditional reliability", {
  f <- csem_estimate(items, theta_points = 801, theta_range = c(-8, 8))
  r <- csem_reliability(f)
  dz <- diff(f$theta_grid)[1]
  avg <- sum(r$crel * dnorm(f$theta_grid)) * dz
  expect_equal(avg, unname(f$reliability["information"]), tolerance = 1e-3)
  expect_true(all(r$crel > 0 & r$crel < 1))
  rc <- csem_reliability(f, type = "complement")
  expect_true(all(rc$crel >= 0))
  expect_gt(rc$n_truncated, 0)
})

test_that("delta-method band brackets the estimate", {
  f <- csem_estimate(items, theta_points = 21)
  d <- csem_delta(f)
  expect_true(all(d$lower < d$csem & d$csem < d$upper))
  expect_error(csem_delta(csem_estimate(items, theta_points = 21, SE = FALSE)),
               "covariance")
})

test_that("bootstrap runs for parametric and nonparametric resampling", {
  skip_on_cran()
  f <- csem_estimate(items, model = "gpcm", theta_points = 21)
  for (tp in c("parametric", "nonparametric")) {
    b <- csem_bootstrap(f, n_boot = 20, type = tp, seed = 1)
    expect_s3_class(b, "csem_bootstrap")
    expect_true(all(b$csem_boot_lower <= b$csem_boot_upper))
    expect_equal(b$n_boot + sum(b$failures$n), 20)
  }
  b1 <- csem_bootstrap(f, n_boot = 20, seed = 7)
  b2 <- csem_bootstrap(f, n_boot = 20, seed = 7)
  expect_identical(b1$replicates, b2$replicates)
})

test_that("comparison returns aligned profiles and differences", {
  skip_on_cran()
  cmp <- csem_compare(items, models = c("graded", "pcm"), n_boot = 20,
                      seed = 1, theta_points = 21)
  expect_s3_class(cmp, "csem_compare")
  expect_setequal(cmp$models, c("GRM", "PCM"))
  expect_equal(nrow(cmp$differences), 21)
  expect_equal(min(cmp$fit_table$dAIC), 0)
})

test_that("regions translate reliability targets into CSEM thresholds", {
  f <- csem_estimate(items, theta_points = 81)
  r <- csem_regions(f, reliability = c(.8, .5))
  expect_equal(r$csem_threshold, c(0.5, 1))
  expect_true(r$coverage[2] >= r$coverage[1])
})

test_that("plots build for every type", {
  f <- csem_estimate(items, theta_points = 21)
  for (tp in c("csem", "info", "items", "reliability"))
    expect_s3_class(csem_plot(f, type = tp), "ggplot")
  expect_s3_class(csem_plot(csem_delta(f)), "ggplot")
  expect_s3_class(csem_plot(csem_compare(items, theta_points = 21)), "ggplot")
})
