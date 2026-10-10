#' Bootstrap Confidence Band for the CSEM Profile
#'
#' Quantifies the sampling uncertainty of the CSEM profile by resampling.
#' In each replication a new data set is generated, the same polytomous IRT
#' model is re-estimated, and the CSEM profile is recomputed on the common
#' theta grid. The empirical distribution of the replicated profiles yields
#' bootstrap standard errors and percentile confidence intervals at every
#' theta point.
#'
#' @param fit An object of class \code{"csem_fit"} from \code{csem_estimate()}.
#' @param n_boot Integer. Number of bootstrap replications. Default 500;
#'   at least 1000 is recommended for reported 95\% percentile intervals.
#' @param type Character. \code{"parametric"} (default): responses are
#'   simulated from the fitted model with
#'   \eqn{\theta \sim N(0, \hat\sigma^2)}, so the band reflects sampling
#'   error under the model. \code{"nonparametric"}: persons are resampled
#'   with replacement, so the band also reflects model misfit.
#' @param conf_level Numeric. Confidence level. Default 0.95.
#' @param seed Integer or \code{NULL}. Random seed. Default 42.
#' @param cores Integer. Number of cores (forked processes; ignored on
#'   Windows). Default 1.
#' @param verbose Logical. Report progress. Default \code{FALSE}.
#'
#' @details
#' A replication is discarded (and counted in \code{failures}) if the
#' resampled data lose a response category for any item, because the
#' re-estimated model would then have a different parameter structure; if
#' estimation fails; or if the EM algorithm does not converge. Results are
#' based on the successful replications only, and the number and reasons of
#' failures are reported so that users can judge whether the band is
#' trustworthy.
#'
#' @return An object of class \code{"csem_bootstrap"} with the analytical
#'   CSEM, bootstrap mean, standard error and percentile limits
#'   (\code{csem_df}), the matrix of replicated profiles
#'   (\code{replicates}), the number of successful replications
#'   (\code{n_boot}) and a table of failures (\code{failures}).
#'
#' @references
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An introduction to the
#' bootstrap}. Chapman & Hall.
#'
#' @examples
#' \donttest{
#' data(attitude_data)
#' fit  <- csem_estimate(attitude_data[, 1:10], model = "gpcm")
#' boot <- csem_bootstrap(fit, n_boot = 50, seed = 1)
#' boot
#' }
#'
#' @export
csem_bootstrap <- function(fit,
                           n_boot     = 500,
                           type       = c("parametric", "nonparametric"),
                           conf_level = 0.95,
                           seed       = 42,
                           cores      = 1L,
                           verbose    = FALSE) {

  if (!inherits(fit, "csem_fit"))
    stop("'fit' must be a 'csem_fit' object from csem_estimate().",
         call. = FALSE)
  type <- match.arg(type)
  if (n_boot < 20)
    stop("'n_boot' must be at least 20.", call. = FALSE)
  if (!is.null(seed)) set.seed(seed)
  if (cores > 1L) RNGkind("L'Ecuyer-CMRG")

  z     <- fit$theta_grid
  std   <- fit$standardized
  dat   <- fit$data
  ref_k <- fit$n_cats
  N     <- fit$n_persons
  sd_th <- sqrt(fit$latent_var)
  miss  <- is.na(dat)

  one_rep <- function(b) {
    db <- if (type == "parametric") {
      th <- matrix(stats::rnorm(N, 0, sd_th), ncol = 1)
      s  <- tryCatch(suppressMessages(
        mirt::simdata(model = fit$fit, N = N, Theta = th)),
        error = function(e) NULL)
      if (is.null(s)) return("simulation")
      s <- as.data.frame(s)
      # Keep the observed missing-data pattern and category coding
      mins <- vapply(dat, min, numeric(1), na.rm = TRUE)
      s <- as.data.frame(mapply(function(col, m) col + m, s, mins,
                                SIMPLIFY = FALSE))
      s[miss] <- NA
      names(s) <- names(dat)
      s
    } else {
      dat[sample.int(N, N, replace = TRUE), , drop = FALSE]
    }
    if (any(.n_cats(db) != ref_k)) return("category loss")
    m <- .fit_quiet(db, fit$itemtype)
    if (is.null(m)) return("non-convergence")
    sdb <- if (std) sqrt(.latent_var(m)) else 1
    .csem_from_info(.info_grid(m, z, sdb))
  }

  if (verbose) message("Running ", n_boot, " ", type,
                       " bootstrap replications...")
  res <- .lapply_par(seq_len(n_boot), one_rep, cores = cores)

  ok    <- vapply(res, is.numeric, logical(1))
  fails <- unlist(res[!ok])
  if (sum(ok) < 20)
    stop("Only ", sum(ok), " of ", n_boot, " replications succeeded (",
         paste(names(table(fails)), table(fails), sep = ": ",
               collapse = "; "), ").", call. = FALSE)
  if (mean(!ok) > 0.10)
    warning(sprintf("%.0f%% of replications failed; interpret the band ",
                    100 * mean(!ok)), "with caution.", call. = FALSE)

  R     <- do.call(rbind, res[ok])
  alpha <- 1 - conf_level
  csem_obs <- .csem_from_info(fit$test_info)

  df <- data.frame(
    theta           = z,
    csem_analytical = csem_obs,
    csem_boot_mean  = colMeans(R),
    csem_boot_se    = apply(R, 2, stats::sd),
    csem_lower      = apply(R, 2, stats::quantile, probs = alpha / 2),
    csem_upper      = apply(R, 2, stats::quantile, probs = 1 - alpha / 2))
  df$bias <- df$csem_boot_mean - df$csem_analytical

  failures <- as.data.frame(table(reason = factor(
    fails, levels = c("simulation", "category loss", "non-convergence"))),
    responseName = "n")

  structure(list(
    theta_grid      = z,
    csem_analytical = csem_obs,
    csem_boot_mean  = df$csem_boot_mean,
    csem_boot_se    = df$csem_boot_se,
    csem_boot_lower = df$csem_lower,
    csem_boot_upper = df$csem_upper,
    csem_df         = df,
    replicates      = R,
    n_boot          = sum(ok),
    n_requested     = n_boot,
    failures        = failures,
    type            = type,
    conf_level      = conf_level,
    model           = fit$model,
    model_label     = fit$model_label,
    model_short     = fit$model_short,
    var_theta       = fit$var_theta
  ), class = "csem_bootstrap")
}


#' @export
print.csem_bootstrap <- function(x, ...) {
  cat("=== Bootstrap CSEM band - ", x$model_label, " ===\n\n", sep = "")
  cat(sprintf("Type: %s | successful replications: %d of %d | %.0f%% percentile CI\n",
              x$type, x$n_boot, x$n_requested, 100 * x$conf_level))
  if (any(x$failures$n > 0)) {
    f <- x$failures[x$failures$n > 0, ]
    cat("Failed replications: ",
        paste(f$reason, f$n, sep = " = ", collapse = ", "), "\n", sep = "")
  }
  cat("\n")
  print(round(x$csem_df[.show_idx(x$theta_grid), ], 4), row.names = FALSE)
  invisible(x)
}
