#' Fit a Polytomous IRT Model for CSEM Analysis
#'
#' Fits a unidimensional polytomous item response theory model to item
#' response data by marginal maximum likelihood (EM algorithm, via
#' \pkg{mirt}) and returns the item parameters, information functions,
#' person estimates and fit statistics required for conditional standard
#' error of measurement (CSEM) analysis.
#'
#' @param data A \code{data.frame} or \code{matrix} of item responses
#'   (rows = persons, columns = items). Responses must be integer category
#'   codes; non-item columns must be removed beforehand. Missing responses
#'   (\code{NA}) are allowed and handled by full-information estimation.
#' @param model Character. The polytomous IRT model:
#'   \code{"graded"} (Graded Response Model; alias \code{"grm"}),
#'   \code{"gpcm"} (Generalized Partial Credit Model),
#'   \code{"pcm"} (Partial Credit Model; alias \code{"rasch"}),
#'   \code{"grsm"} (Graded Rating Scale Model) or
#'   \code{"rsm"} (Rating Scale Model). Default \code{"graded"}.
#' @param theta_range Numeric vector of length 2: range of the theta grid.
#'   Default \code{c(-4, 4)}.
#' @param theta_points Integer. Number of theta grid points. Default 161
#'   (a step of 0.05 on the default range).
#' @param method Character. Person estimation method passed to
#'   \code{mirt::fscores()}: \code{"EAP"} (default), \code{"MAP"} or
#'   \code{"ML"}.
#' @param standardize Logical. If \code{TRUE} (default), all quantities are
#'   expressed on a standardized latent metric with unit variance. This
#'   matters for the PCM and RSM, whose latent variance is freely estimated
#'   (slopes are fixed to 1): without rescaling, their CSEM values would not
#'   be comparable with those of the GRM or GPCM. See Details.
#' @param SE Logical. Compute the parameter covariance matrix (needed for
#'   \code{csem_delta()}). Default \code{TRUE}.
#' @param verbose Logical. Print estimation progress. Default \code{FALSE}.
#'
#' @details
#' \strong{Common metric.} The GRM, GPCM and GRSM identify the latent scale
#' by fixing \eqn{\mathrm{Var}(\theta) = 1}. The PCM and RSM instead fix all
#' slopes to 1 and estimate \eqn{\sigma^2 = \mathrm{Var}(\theta)}. With
#' \code{standardize = TRUE}, the theta grid is defined on \eqn{z =
#' \theta/\sigma} and information is transformed as \eqn{I_z(z) = \sigma^2
#' I_\theta(\sigma z)}, so that \eqn{\mathrm{CSEM}_z(z) =
#' \mathrm{CSEM}_\theta(\sigma z)/\sigma}. For models with \eqn{\sigma = 1}
#' the transformation is the identity.
#'
#' @return An object of class \code{"csem_fit"}, a list containing the
#'   fitted \code{mirt} model (\code{fit}), the model specification
#'   (\code{model}, \code{itemtype}, \code{model_label}), item parameters
#'   (\code{item_params}, IRT parameterization), the theta grid
#'   (\code{theta_grid}), test and item information on that grid
#'   (\code{test_info}, \code{item_info}), the latent variance
#'   (\code{latent_var}) and scaling constant (\code{scale_sd}), person
#'   estimates and their standard errors on the reporting metric
#'   (\code{theta_persons}, \code{se_persons}), marginal reliability
#'   (\code{reliability}: information-based and empirical), fit statistics
#'   (\code{fit_stats}: log-likelihood, number of parameters, AIC, BIC),
#'   a convergence flag and sample information.
#'
#' @references
#' Samejima, F. (1969). Estimation of latent ability using a response
#' pattern of graded scores. \emph{Psychometrika Monograph Supplement, 17}.
#'
#' Masters, G. N. (1982). A Rasch model for partial credit scoring.
#' \emph{Psychometrika, 47}(2), 149-174. \doi{10.1007/BF02296272}
#'
#' Muraki, E. (1992). A generalized partial credit model: Application of an
#' EM algorithm. \emph{Applied Psychological Measurement, 16}(2), 159-176.
#' \doi{10.1177/014662169201600206}
#'
#' Chalmers, R. P. (2012). mirt: A multidimensional item response theory
#' package for the R environment. \emph{Journal of Statistical Software,
#' 48}(6), 1-29. \doi{10.18637/jss.v048.i06}
#'
#' @examples
#' data(attitude_data)
#' items <- attitude_data[, 1:10]
#' fit_grm  <- csem_estimate(items, model = "graded")
#' fit_pcm  <- csem_estimate(items, model = "pcm")
#' fit_grm
#'
#' @export
csem_estimate <- function(data,
                          model        = "graded",
                          theta_range  = c(-4, 4),
                          theta_points = 161,
                          method       = "EAP",
                          standardize  = TRUE,
                          SE           = TRUE,
                          verbose      = FALSE) {

  spec <- .match_model(model)

  # -- Input validation ----------------------------------------------------
  if (!is.data.frame(data) && !is.matrix(data))
    stop("'data' must be a data.frame or matrix.", call. = FALSE)
  data <- as.data.frame(data)

  not_num <- !vapply(data, is.numeric, logical(1))
  if (any(not_num))
    stop("All columns must be numeric item responses. Non-numeric columns: ",
         paste(names(data)[not_num], collapse = ", "), call. = FALSE)
  if (any(!is.na(unlist(data)) & unlist(data) %% 1 != 0))
    stop("Item responses must be integer category codes.", call. = FALSE)
  if (any(colSums(!is.na(data)) == 0))
    stop("At least one item has no observed responses.", call. = FALSE)

  if (length(theta_range) != 2L || diff(theta_range) <= 0)
    stop("'theta_range' must be c(lower, upper) with lower < upper.",
         call. = FALSE)

  n_persons <- nrow(data)
  n_items   <- ncol(data)
  n_cats    <- .n_cats(data)

  if (any(n_cats < 2L))
    stop("Items with fewer than two observed categories: ",
         paste(names(data)[n_cats < 2L], collapse = ", "), call. = FALSE)
  if (all(n_cats == 2L))
    warning("All items are dichotomous; the polytomous models reduce to ",
            "their binary counterparts.", call. = FALSE)
  if (spec$model %in% c("grsm", "rsm") && length(unique(n_cats)) > 1L)
    stop("The ", spec$short, " requires the same number of categories for ",
         "every item.", call. = FALSE)
  if (n_persons < 200)
    warning("Sample size (N = ", n_persons, ") is below 200; parameter and ",
            "CSEM estimates may be unstable.", call. = FALSE)
  if (n_items < 5)
    warning("Fewer than 5 items (k = ", n_items, ").", call. = FALSE)

  # -- Fit model -----------------------------------------------------------
  if (verbose) message("Fitting ", spec$short, " to ", n_persons,
                       " persons x ", n_items, " items...")
  fit <- tryCatch(
    mirt::mirt(data, model = 1, itemtype = spec$itemtype, SE = SE,
               verbose = verbose),
    error = function(e) stop(spec$short, " estimation failed: ",
                             conditionMessage(e), call. = FALSE))
  converged <- isTRUE(mirt::extract.mirt(fit, "converged"))
  if (!converged)
    warning(spec$short, " estimation did not converge.", call. = FALSE)

  # -- Item parameters (IRT parameterization) -----------------------------
  # (IRT parameterization where mirt supports it; slope-intercept otherwise)
  ip <- tryCatch(mirt::coef(fit, simplify = TRUE, IRTpars = TRUE)$items,
                 error = function(e) mirt::coef(fit, simplify = TRUE)$items)
  ip <- as.data.frame(round(ip, 4))
  ip <- cbind(item = rownames(ip), ip)
  rownames(ip) <- NULL

  # -- Metric -------------------------------------------------------------
  sigma2   <- .latent_var(fit)
  scale_sd <- if (standardize) sqrt(sigma2) else 1
  var_rep  <- if (standardize) 1 else sigma2   # latent variance, reporting metric

  # -- Information on grid ------------------------------------------------
  theta_grid <- seq(theta_range[1], theta_range[2], length.out = theta_points)
  inf <- .info_grid(fit, theta_grid, scale_sd, items = TRUE)

  # -- Person estimates ---------------------------------------------------
  fs <- suppressWarnings(
    mirt::fscores(fit, method = method, full.scores.SE = TRUE))
  theta_p <- as.numeric(fs[, 1]) / scale_sd
  se_p    <- as.numeric(fs[, 2]) / scale_sd

  # -- Marginal reliability -----------------------------------------------
  # Information-based (Green et al., 1984): integrate the conditional
  # reliability var/(var + CSEM^2) over the N(0, var) latent density.
  f_rel <- function(t) {
    I <- .info_grid(fit, t, scale_sd)
    var_rep * I / (var_rep * I + 1) * stats::dnorm(t, 0, sqrt(var_rep))
  }
  rel_info <- stats::integrate(f_rel, -Inf, Inf)$value
  ok <- is.finite(theta_p) & is.finite(se_p)
  rel_emp <- stats::var(theta_p[ok]) /
    (stats::var(theta_p[ok]) + mean(se_p[ok]^2))

  ll   <- mirt::extract.mirt(fit, "logLik")
  npar <- mirt::extract.mirt(fit, "nest")
  fit_stats <- data.frame(
    model  = spec$short,
    logLik = ll,
    npar   = npar,
    AIC    = mirt::extract.mirt(fit, "AIC"),
    BIC    = mirt::extract.mirt(fit, "BIC"),
    stringsAsFactors = FALSE)

  structure(list(
    fit           = fit,
    data          = data,
    model         = spec$model,
    itemtype      = spec$itemtype,
    model_label   = spec$label,
    model_short   = spec$short,
    item_params   = ip,
    theta_grid    = theta_grid,
    test_info     = inf$test,
    item_info     = inf$items,
    latent_var    = sigma2,
    scale_sd      = scale_sd,
    var_theta     = var_rep,
    standardized  = standardize,
    theta_persons = theta_p,
    se_persons    = se_p,
    method        = method,
    reliability   = c(information = rel_info, empirical = rel_emp),
    fit_stats     = fit_stats,
    converged     = converged,
    n_persons     = n_persons,
    n_items       = n_items,
    n_cats        = n_cats,
    theta_range   = theta_range
  ), class = "csem_fit")
}


#' @export
print.csem_fit <- function(x, ...) {
  cat("=== ", x$model_label, " - CSEM analysis ===\n\n", sep = "")
  cat(sprintf("Persons x items:     %d x %d\n", x$n_persons, x$n_items))
  cat(sprintf("Categories per item: %s\n",
              paste(range(x$n_cats), collapse = "-")))
  cat(sprintf("Converged:           %s\n", x$converged))
  cat(sprintf("Latent variance:     %.4f%s\n", x$latent_var,
              if (x$standardized && abs(x$latent_var - 1) > 1e-8)
                " (rescaled to 1 for reporting)" else ""))
  cat(sprintf("logLik / AIC / BIC:  %.2f / %.2f / %.2f\n",
              x$fit_stats$logLik, x$fit_stats$AIC, x$fit_stats$BIC))
  cat(sprintf("Marginal reliability: %.4f (information-based), %.4f (empirical, %s)\n",
              x$reliability["information"], x$reliability["empirical"],
              x$method))
  i <- which.max(x$test_info)
  cat(sprintf("Peak test information: %.3f at theta = %.2f\n\n",
              x$test_info[i], x$theta_grid[i]))
  cat("Item parameters (first 6 items):\n")
  print(utils::head(x$item_params), row.names = FALSE)
  invisible(x)
}
