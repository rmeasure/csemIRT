#' Theta Regions Meeting a Precision Target
#'
#' Identifies the theta intervals in which the CSEM is small enough to
#' attain a target conditional reliability, and the share of the latent
#' population that falls in those intervals. When a confidence band is
#' available, a conservative region is also reported in which even the
#' upper confidence limit of the CSEM meets the target.
#'
#' @param x An object of class \code{"csem_fit"}, \code{"csem_analytical"},
#'   \code{"csem_bootstrap"} or \code{"csem_delta"}.
#' @param reliability Numeric vector of target conditional reliabilities
#'   (ratio definition, see \code{csem_reliability()}).
#'   Default \code{c(0.90, 0.80, 0.70)}.
#'
#' @details
#' Under the ratio definition \eqn{\rho = \sigma^2/(\sigma^2 + e^2)}, a
#' target \eqn{\rho^*} corresponds to the CSEM threshold
#' \eqn{e^* = \sigma\sqrt{1/\rho^* - 1}}; for \eqn{\sigma = 1}, targets of
#' .90, .80 and .70 correspond to CSEM thresholds of 0.333, 0.500 and 0.655.
#' Population coverage is the probability mass of
#' \eqn{N(0, \sigma^2)} inside the region.
#'
#' @return A data frame of class \code{"csem_regions"} with one row per
#'   target: the CSEM threshold, the theta interval(s) meeting it, the
#'   population coverage, and the conservative (band-based) counterparts
#'   when available.
#'
#' @examples
#' data(attitude_data)
#' fit <- csem_estimate(attitude_data[, 1:10])
#' csem_regions(fit)
#'
#' @export
csem_regions <- function(x, reliability = c(0.90, 0.80, 0.70)) {
  if (inherits(x, "csem_fit")) x <- csem_analytical(x)
  if (!inherits(x, c("csem_analytical", "csem_bootstrap", "csem_delta")))
    stop("'x' must be a csem_fit, csem_analytical, csem_bootstrap or ",
         "csem_delta object.", call. = FALSE)
  if (any(reliability <= 0 | reliability >= 1))
    stop("'reliability' targets must lie in (0, 1).", call. = FALSE)

  v  <- x$var_theta %||% 1
  z  <- x$theta_grid
  cs <- if (inherits(x, "csem_bootstrap")) x$csem_analytical else x$csem
  up <- switch(class(x)[1], csem_bootstrap = x$csem_boot_upper,
               csem_delta = x$upper, NULL)
  dz <- diff(z)[1]

  cover <- function(flag) {         # population mass inside flagged grid cells
    if (!any(flag)) return(0)
    lo <- z - dz / 2; hi <- z + dz / 2
    sum((stats::pnorm(hi, 0, sqrt(v)) - stats::pnorm(lo, 0, sqrt(v)))[flag])
  }
  ranges <- function(flag) {
    if (!any(flag)) return("none")
    r <- rle(flag); e <- cumsum(r$lengths); s <- e - r$lengths + 1
    k <- which(r$values)
    paste(sprintf("[%.2f, %.2f]", z[s[k]], z[e[k]]), collapse = " U ")
  }

  out <- do.call(rbind, lapply(reliability, function(r) {
    thr  <- sqrt(v * (1 / r - 1))
    flag <- cs <= thr
    row <- data.frame(target_rel = r, csem_threshold = thr,
                      theta_region = ranges(flag),
                      coverage = cover(flag), stringsAsFactors = FALSE)
    if (!is.null(up)) {
      fc <- up <= thr
      row$region_conservative   <- ranges(fc)
      row$coverage_conservative <- cover(fc)
    }
    row
  }))
  attr(out, "model") <- x$model_label
  attr(out, "range") <- range(z)
  class(out) <- c("csem_regions", "data.frame")
  out
}


#' @export
print.csem_regions <- function(x, ...) {
  cat("=== Precision regions - ", attr(x, "model"), " ===\n", sep = "")
  cat(sprintf("(searched theta in [%.1f, %.1f]; coverage = share of N(0, var) population)\n\n",
              attr(x, "range")[1], attr(x, "range")[2]))
  y <- as.data.frame(unclass(x))
  num <- vapply(y, is.numeric, logical(1))
  y[num] <- lapply(y[num], round, 3)
  print(y, row.names = FALSE)
  invisible(x)
}
