#' Conditional Reliability Profile
#'
#' Converts the CSEM profile into a conditional reliability coefficient at
#' each theta value.
#'
#' @param x An object of class \code{"csem_fit"}, \code{"csem_analytical"},
#'   \code{"csem_bootstrap"} or \code{"csem_delta"}. When an object with a
#'   confidence band is supplied, the band is carried over to the
#'   reliability scale (the transformation is monotone).
#' @param type Character. \code{"ratio"} (default):
#'   \eqn{\rho(\theta) = \sigma^2 / (\sigma^2 + \mathrm{CSEM}^2(\theta))},
#'   which lies in (0, 1) everywhere and whose average over the latent
#'   density is the information-based marginal reliability reported by
#'   \code{csem_estimate()}. \code{"complement"}:
#'   \eqn{\rho(\theta) = 1 - \mathrm{CSEM}^2(\theta)/\sigma^2}, truncated at
#'   0; this form is negative wherever CSEM exceeds the latent SD, and the
#'   truncation is reported.
#' @param var_theta Numeric or \code{NULL}. Latent variance on the reporting
#'   metric. Default \code{NULL} uses the value stored in \code{x}
#'   (1 on the standardized metric).
#'
#' @return An object of class \code{"csem_reliability"}.
#'
#' @references
#' Green, B. F., Bock, R. D., Humphreys, L. G., Linn, R. L., & Reckase,
#' M. D. (1984). Technical guidelines for assessing computerized adaptive
#' tests. \emph{Journal of Educational Measurement, 21}(4), 347-360.
#' \doi{10.1111/j.1745-3984.1984.tb01039.x}
#'
#' Samejima, F. (1994). Estimation of reliability coefficients using the
#' test information function and its modifications. \emph{Applied
#' Psychological Measurement, 18}(3), 229-244.
#' \doi{10.1177/014662169401800304}
#'
#' @examples
#' data(attitude_data)
#' fit <- csem_estimate(attitude_data[, 1:10])
#' csem_reliability(fit)
#'
#' @export
csem_reliability <- function(x, type = c("ratio", "complement"),
                             var_theta = NULL) {
  type <- match.arg(type)
  if (inherits(x, "csem_fit")) {
    marg <- x$reliability
    x <- csem_analytical(x)
  } else marg <- NULL
  if (!inherits(x, c("csem_analytical", "csem_bootstrap", "csem_delta")))
    stop("'x' must be a csem_fit, csem_analytical, csem_bootstrap or ",
         "csem_delta object.", call. = FALSE)

  v <- var_theta %||% x$var_theta %||% 1
  if (v <= 0) stop("'var_theta' must be positive.", call. = FALSE)

  csem <- if (inherits(x, "csem_bootstrap")) x$csem_analytical else x$csem
  lo <- switch(class(x)[1], csem_bootstrap = x$csem_boot_lower,
               csem_delta = x$lower, NULL)
  up <- switch(class(x)[1], csem_bootstrap = x$csem_boot_upper,
               csem_delta = x$upper, NULL)

  tr <- if (type == "ratio") function(e) v / (v + e^2)
        else function(e) pmax(1 - e^2 / v, 0)

  df <- data.frame(theta = x$theta_grid, csem = csem, crel = tr(csem))
  if (!is.null(lo)) {           # larger CSEM -> lower reliability
    df$crel_lower <- tr(up)
    df$crel_upper <- tr(lo)
  }
  n_trunc <- if (type == "complement") sum(csem^2 >= v) else 0L

  structure(list(
    theta_grid           = x$theta_grid,
    csem                 = csem,
    crel                 = df$crel,
    rel_df               = df,
    type                 = type,
    var_theta            = v,
    n_truncated          = n_trunc,
    marginal_reliability = marg,
    model_label          = x$model_label,
    model_short          = x$model_short
  ), class = "csem_reliability")
}


#' @export
print.csem_reliability <- function(x, ...) {
  cat("=== Conditional reliability - ", x$model_label, " ===\n\n", sep = "")
  cat("Definition: ", if (x$type == "ratio")
    "var / (var + CSEM^2)" else "1 - CSEM^2 / var (truncated at 0)",
    sprintf(" with var = %.3f\n", x$var_theta), sep = "")
  if (!is.null(x$marginal_reliability))
    cat(sprintf("Marginal reliability: %.4f (information-based)\n",
                x$marginal_reliability["information"]))
  if (x$n_truncated > 0)
    cat(sprintf("Note: %d grid points truncated at 0 (CSEM >= latent SD).\n",
                x$n_truncated))
  cat("\n")
  print(round(x$rel_df[.show_idx(x$theta_grid), ], 4), row.names = FALSE)
  invisible(x)
}
