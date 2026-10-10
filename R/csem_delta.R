#' Delta-Method Confidence Band for the CSEM Profile
#'
#' Approximates the sampling uncertainty of the analytical CSEM profile from
#' the asymptotic covariance matrix of the item (and latent variance)
#' parameter estimates. The Jacobian of \eqn{\log \mathrm{CSEM}(\theta)} with
#' respect to the free parameters is obtained by central finite differences,
#' and Wald intervals are formed on the log scale and back-transformed, which
#' keeps the bounds positive and respects the skewness of the CSEM.
#'
#' The delta method is fast (no re-estimation) and serves as an analytic
#' cross-check of the bootstrap band from \code{csem_bootstrap()}; the two
#' typically agree closely in the central theta range and diverge at the
#' extremes, where the linear approximation is least accurate.
#'
#' @param fit An object of class \code{"csem_fit"} estimated with
#'   \code{SE = TRUE}.
#' @param conf_level Numeric. Confidence level. Default 0.95.
#' @param h Numeric. Relative step size for finite differences.
#'   Default 1e-4.
#'
#' @return An object of class \code{"csem_delta"} with the analytical CSEM,
#'   its delta-method standard error and confidence limits on the theta grid
#'   (\code{csem_df}).
#'
#' @references
#' Oehlert, G. W. (1992). A note on the delta method. \emph{The American
#' Statistician, 46}(1), 27-29. \doi{10.1080/00031305.1992.10475842}
#'
#' @examples
#' data(attitude_data)
#' fit <- csem_estimate(attitude_data[, 1:10])
#' dlt <- csem_delta(fit)
#' dlt
#'
#' @export
csem_delta <- function(fit, conf_level = 0.95, h = 1e-4) {
  if (!inherits(fit, "csem_fit"))
    stop("'fit' must be a 'csem_fit' object.", call. = FALSE)
  V <- tryCatch(mirt::extract.mirt(fit$fit, "vcov"), error = function(e) NULL)
  if (is.null(V) || !length(V) || any(!is.finite(diag(V))))
    stop("No usable parameter covariance matrix. Re-run csem_estimate() ",
         "with SE = TRUE.", call. = FALSE)

  sv   <- mirt::mod2values(fit$fit)
  dat  <- fit$data
  z    <- fit$theta_grid
  std  <- fit$standardized

  # Map each covariance-matrix column to the parnum(s) it controls
  # (constrained parameters share a column, e.g. "b1.2.7.12")
  pn <- lapply(strsplit(colnames(V), ".", fixed = TRUE),
               function(p) as.integer(p[-1]))
  p0 <- vapply(pn, function(k) sv$value[match(k[1], sv$parnum)], numeric(1))

  logcsem_at <- function(p) {
    s <- sv
    for (i in seq_along(pn)) s$value[s$parnum %in% pn[[i]]] <- p[i]
    m <- suppressWarnings(suppressMessages(
      mirt::mirt(dat, 1, itemtype = fit$itemtype, pars = s,
                 TOL = NaN, verbose = FALSE)))
    sd <- if (std) sqrt(.latent_var(m)) else 1
    log(.csem_from_info(.info_grid(m, z, sd)))
  }

  f0 <- logcsem_at(p0)
  J  <- vapply(seq_along(p0), function(i) {
    step <- h * max(1, abs(p0[i]))
    up <- dn <- p0
    up[i] <- up[i] + step; dn[i] <- dn[i] - step
    (logcsem_at(up) - logcsem_at(dn)) / (2 * step)
  }, numeric(length(z)))
  J <- matrix(J, nrow = length(z))

  se_log <- sqrt(pmax(rowSums((J %*% V) * J), 0))
  q      <- stats::qnorm(1 - (1 - conf_level) / 2)
  csem   <- exp(f0)

  df <- data.frame(theta = z, csem = csem,
                   se    = csem * se_log,          # delta SE on CSEM scale
                   lower = exp(f0 - q * se_log),
                   upper = exp(f0 + q * se_log))

  structure(list(
    theta_grid  = z,
    csem        = csem,
    se          = df$se,
    lower       = df$lower,
    upper       = df$upper,
    csem_df     = df,
    conf_level  = conf_level,
    model       = fit$model,
    model_label = fit$model_label,
    model_short = fit$model_short,
    var_theta   = fit$var_theta,
    method      = "delta"
  ), class = "csem_delta")
}


#' @export
print.csem_delta <- function(x, ...) {
  cat("=== Delta-method CSEM band - ", x$model_label, " ===\n\n", sep = "")
  cat(sprintf("Confidence level: %.0f%% (log-scale Wald interval)\n\n",
              100 * x$conf_level))
  print(round(x$csem_df[.show_idx(x$theta_grid), ], 4), row.names = FALSE)
  invisible(x)
}
