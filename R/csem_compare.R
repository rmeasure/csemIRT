#' Compare CSEM Profiles Across Polytomous IRT Models
#'
#' Fits several polytomous IRT models to the same data, expresses their
#' CSEM profiles on a common standardized latent metric, and (optionally)
#' quantifies whether the models imply different measurement precision at
#' each theta value using a paired nonparametric bootstrap.
#'
#' @param data Item response data (see \code{csem_estimate()}).
#' @param models Character vector of models to compare (see
#'   \code{csem_estimate()}). Default \code{c("graded", "gpcm", "pcm")}.
#' @param reference Character. Model against which differences are formed.
#'   Default: the first element of \code{models}.
#' @param n_boot Integer. Number of paired bootstrap replications; 0
#'   (default) skips the bootstrap and returns analytical profiles only.
#' @param conf_level Numeric. Confidence level. Default 0.95.
#' @param min_diff Numeric. Smallest CSEM difference regarded as practically
#'   relevant. Differences are flagged as practically relevant where the
#'   bootstrap interval lies entirely beyond \code{+/- min_diff}.
#'   Default 0.05.
#' @param seed Integer or \code{NULL}. Random seed. Default 42.
#' @param cores Integer. Number of cores. Default 1.
#' @param ... Further arguments passed to \code{csem_estimate()} (e.g.
#'   \code{theta_range}, \code{theta_points}, \code{method}).
#'
#' @details
#' The standardized metric (unit latent variance) is used for every model,
#' so that the PCM and RSM, which estimate the latent variance, are
#' comparable with the GRM, GPCM and GRSM.
#'
#' In each bootstrap replication persons are resampled once and \emph{all}
#' models are refitted to the same resample. The replicated difference
#' \eqn{\mathrm{CSEM}_m(\theta) - \mathrm{CSEM}_{ref}(\theta)} therefore
#' carries the dependence between models induced by the common data, and
#' its percentile interval tests, pointwise, whether model choice changes
#' the conclusion about precision at that theta value. Replications in which
#' any model fails or a category is lost are discarded for all models.
#'
#' @return An object of class \code{"csem_compare"} containing the fitted
#'   \code{csem_fit} objects (\code{fits}), a model-fit table
#'   (\code{fit_table}: log-likelihood, number of parameters, AIC, BIC,
#'   latent variance, marginal reliability, minimum CSEM and its location),
#'   long-format CSEM profiles (\code{profiles}) and, when
#'   \code{n_boot > 0}, pointwise differences from the reference model with
#'   confidence limits (\code{differences}).
#'
#' @examples
#' data(attitude_data)
#' cmp <- csem_compare(attitude_data[, 1:10])
#' cmp
#' \donttest{
#' cmp_b <- csem_compare(attitude_data[, 1:10], n_boot = 50, seed = 1)
#' csem_plot(cmp_b, type = "difference")
#' }
#'
#' @export
csem_compare <- function(data,
                         models     = c("graded", "gpcm", "pcm"),
                         reference  = models[1],
                         n_boot     = 0,
                         conf_level = 0.95,
                         min_diff   = 0.05,
                         seed       = 42,
                         cores      = 1L,
                         ...) {

  models <- vapply(models, function(m) .match_model(m)$model, character(1),
                   USE.NAMES = FALSE)
  if (anyDuplicated(models)) stop("Duplicated models.", call. = FALSE)
  if (length(models) < 2L)
    stop("Provide at least two models to compare.", call. = FALSE)
  reference <- .match_model(reference)$model
  if (!reference %in% models)
    stop("'reference' must be one of 'models'.", call. = FALSE)

  dots <- list(...)
  if (isFALSE(dots$standardize))
    stop("csem_compare() always uses the standardized metric.", call. = FALSE)
  dots$standardize <- NULL

  fits <- lapply(models, function(m)
    do.call(csem_estimate, c(list(data = data, model = m,
                                  standardize = TRUE), dots)))
  names(fits) <- vapply(fits, `[[`, character(1), "model_short")
  z <- fits[[1]]$theta_grid

  fit_table <- do.call(rbind, lapply(fits, function(f) {
    cs <- .csem_from_info(f$test_info)
    data.frame(f$fit_stats,
               latent_var  = f$latent_var,
               rel_marg    = f$reliability[["information"]],
               min_csem    = min(cs),
               theta_min   = z[which.min(cs)],
               converged   = f$converged,
               stringsAsFactors = FALSE)
  }))
  fit_table$dAIC <- fit_table$AIC - min(fit_table$AIC)
  fit_table$dBIC <- fit_table$BIC - min(fit_table$BIC)
  rownames(fit_table) <- NULL

  profiles <- do.call(rbind, lapply(fits, function(f)
    data.frame(model = f$model_short, theta = z,
               test_info = f$test_info,
               csem = .csem_from_info(f$test_info),
               stringsAsFactors = FALSE)))
  rownames(profiles) <- NULL

  ref_short <- .match_model(reference)$short
  out <- list(fits = fits, fit_table = fit_table, profiles = profiles,
              theta_grid = z, reference = ref_short,
              models = names(fits), n_boot = 0, conf_level = conf_level,
              min_diff = min_diff)

  # -- Paired nonparametric bootstrap -------------------------------------
  if (n_boot > 0) {
    if (n_boot < 20) stop("'n_boot' must be 0 or at least 20.", call. = FALSE)
    if (!is.null(seed)) set.seed(seed)
    if (cores > 1L) RNGkind("L'Ecuyer-CMRG")
    dat   <- fits[[1]]$data
    ref_k <- fits[[1]]$n_cats
    N     <- nrow(dat)
    itypes <- vapply(fits, `[[`, character(1), "itemtype")

    one_rep <- function(b) {
      db <- dat[sample.int(N, N, replace = TRUE), , drop = FALSE]
      if (any(.n_cats(db) != ref_k)) return("category loss")
      out <- vector("list", length(itypes))
      for (i in seq_along(itypes)) {
        m <- .fit_quiet(db, itypes[i])
        if (is.null(m)) return("non-convergence")
        out[[i]] <- .csem_from_info(.info_grid(m, z, sqrt(.latent_var(m))))
      }
      do.call(cbind, out)
    }
    res <- .lapply_par(seq_len(n_boot), one_rep, cores = cores)
    ok  <- vapply(res, is.matrix, logical(1))
    if (sum(ok) < 20)
      stop("Only ", sum(ok), " paired replications succeeded.", call. = FALSE)
    fails <- unlist(res[!ok])
    arr <- simplify2array(res[ok])            # theta x model x replicate
    alpha <- 1 - conf_level
    q <- function(x, p) stats::quantile(x, p, names = FALSE)

    # Bands for each model
    band <- do.call(rbind, lapply(seq_along(fits), function(i) {
      R <- arr[, i, ]
      data.frame(model = names(fits)[i], theta = z,
                 csem_boot_se = apply(R, 1, stats::sd),
                 lower = apply(R, 1, q, alpha / 2),
                 upper = apply(R, 1, q, 1 - alpha / 2))
    }))
    out$profiles <- merge(out$profiles, band, by = c("model", "theta"),
                          sort = FALSE)
    out$profiles <- out$profiles[order(match(out$profiles$model, names(fits)),
                                       out$profiles$theta), ]
    rownames(out$profiles) <- NULL

    # Differences from the reference model
    ir <- match(ref_short, names(fits))
    out$differences <- do.call(rbind, lapply(setdiff(seq_along(fits), ir),
      function(i) {
        D   <- arr[, i, ] - arr[, ir, ]
        est <- .csem_from_info(fits[[i]]$test_info) -
               .csem_from_info(fits[[ir]]$test_info)
        lo  <- apply(D, 1, q, alpha / 2)
        hi  <- apply(D, 1, q, 1 - alpha / 2)
        data.frame(model = names(fits)[i], reference = ref_short,
                   theta = z, diff = est,
                   se = apply(D, 1, stats::sd), lower = lo, upper = hi,
                   differs = lo > 0 | hi < 0,
                   relevant = lo > min_diff | hi < -min_diff)
      }))
    rownames(out$differences) <- NULL
    out$n_boot   <- sum(ok)
    out$n_requested <- n_boot
    out$failures <- table(factor(fails,
                                 levels = c("category loss", "non-convergence")))
  }
  class(out) <- "csem_compare"
  out
}


#' @export
print.csem_compare <- function(x, ...) {
  cat("=== CSEM comparison across polytomous IRT models ===\n")
  cat("(standardized latent metric; reference model: ", x$reference, ")\n\n",
      sep = "")
  ft <- x$fit_table
  num <- vapply(ft, is.numeric, logical(1))
  ft[num] <- lapply(ft[num], round, 3)
  print(ft[, c("model", "logLik", "npar", "AIC", "dAIC", "BIC", "dBIC",
               "latent_var", "rel_marg", "min_csem", "theta_min")],
        row.names = FALSE)

  idx <- .show_idx(x$theta_grid)
  wide <- data.frame(theta = x$theta_grid[idx])
  for (m in x$models)
    wide[[m]] <- x$profiles$csem[x$profiles$model == m][idx]
  cat("\nCSEM at selected theta values:\n")
  print(round(wide, 3), row.names = FALSE)

  if (!is.null(x$differences)) {
    cat(sprintf("\nPaired nonparametric bootstrap: %d of %d replications, %.0f%% CI\n",
                x$n_boot, x$n_requested, 100 * x$conf_level))
    for (m in unique(x$differences$model)) {
      d <- x$differences[x$differences$model == m, ]
      cat(sprintf("  %s - %s: CI excludes 0 at %d of %d grid points%s\n",
                  m, x$reference, sum(d$differs), nrow(d),
                  .describe_ranges(d$theta, d$differs)))
      cat(sprintf("  %s - %s: |difference| > %.2f with confidence at %d grid points%s\n",
                  m, x$reference, x$min_diff, sum(d$relevant),
                  .describe_ranges(d$theta, d$relevant)))
    }
  }
  invisible(x)
}

# Summarize a logical vector over a grid as theta intervals
.describe_ranges <- function(theta, flag) {
  if (!any(flag)) return("")
  r <- rle(flag); ends <- cumsum(r$lengths); starts <- ends - r$lengths + 1
  s <- which(r$values)
  paste0(" (theta in ", paste(sprintf("[%.2f, %.2f]", theta[starts[s]],
                                       theta[ends[s]]), collapse = ", "), ")")
}
