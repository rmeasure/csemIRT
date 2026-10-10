#' Analytical Conditional Standard Error of Measurement
#'
#' Computes the analytical CSEM profile of a fitted polytomous IRT model as
#' the reciprocal square root of the test information function,
#' \eqn{\mathrm{CSEM}(\theta) = 1/\sqrt{I(\theta)}} (Lord, 1980), which is
#' the asymptotic standard error of the maximum likelihood estimate of
#' \eqn{\theta}.
#'
#' @param fit An object of class \code{"csem_fit"} from \code{csem_estimate()}.
#'
#' @return An object of class \code{"csem_analytical"} with the theta grid,
#'   test information, CSEM, a data frame of all three (\code{csem_df}), and
#'   the location and value of minimum CSEM / maximum information.
#'
#' @references
#' Lord, F. M. (1980). \emph{Applications of item response theory to
#' practical testing problems}. Lawrence Erlbaum Associates.
#'
#' @examples
#' data(attitude_data)
#' fit <- csem_estimate(attitude_data[, 1:10], model = "gpcm")
#' csem_analytical(fit)
#'
#' @export
csem_analytical <- function(fit) {
  if (!inherits(fit, "csem_fit"))
    stop("'fit' must be a 'csem_fit' object from csem_estimate().",
         call. = FALSE)

  theta <- fit$theta_grid
  info  <- fit$test_info
  csem  <- .csem_from_info(info)
  i_min <- which.min(csem)

  structure(list(
    theta_grid     = theta,
    test_info      = info,
    csem           = csem,
    csem_df        = data.frame(theta = theta, test_info = info, csem = csem),
    theta_min_csem = theta[i_min],
    min_csem       = csem[i_min],
    theta_max_info = theta[i_min],
    max_info       = info[i_min],
    model          = fit$model,
    model_label    = fit$model_label,
    model_short    = fit$model_short,
    var_theta      = fit$var_theta,
    n_items        = fit$n_items,
    n_persons      = fit$n_persons
  ), class = "csem_analytical")
}


#' @export
print.csem_analytical <- function(x, ...) {
  cat("=== Analytical CSEM - ", x$model_label, " ===\n\n", sep = "")
  cat(sprintf("Peak information: %.4f at theta = %.2f\n",
              x$max_info, x$theta_max_info))
  cat(sprintf("Minimum CSEM:     %.4f at theta = %.2f\n",
              x$min_csem, x$theta_min_csem))
  cat(sprintf("CSEM range:       %.4f to %.4f\n\n",
              min(x$csem), max(x$csem)))
  cat("CSEM at selected theta values:\n")
  df <- x$csem_df[.show_idx(x$theta_grid), ]
  print(round(df, 4), row.names = FALSE)
  invisible(x)
}
