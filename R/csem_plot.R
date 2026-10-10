#' Plot CSEM Results
#'
#' Produces \pkg{ggplot2} visualizations of CSEM analyses.
#'
#' @param x An object returned by a \pkg{csemIRT} function.
#' @param type Character. One of
#'   \code{"csem"} (CSEM profile with precision benchmarks; any object),
#'   \code{"info"} (test information function; \code{csem_fit},
#'   \code{csem_analytical} or \code{csem_compare}),
#'   \code{"items"} (item information functions; \code{csem_fit}),
#'   \code{"reliability"} (conditional reliability; \code{csem_reliability}
#'   or any object accepted by \code{csem_reliability()}),
#'   \code{"band"} (CSEM with confidence band; \code{csem_bootstrap},
#'   \code{csem_delta} or a \code{csem_compare} with bootstrap),
#'   \code{"difference"} (CSEM differences from the reference model with
#'   bootstrap band; \code{csem_compare} with \code{n_boot > 0}).
#'   The default depends on the class of \code{x}.
#' @param benchmarks Numeric vector of conditional-reliability targets drawn
#'   as horizontal reference lines on CSEM plots. Default
#'   \code{c(0.90, 0.80)} (CSEM = 0.333 and 0.500 on the standardized metric).
#' @param title Optional plot title.
#'
#' @return A \code{ggplot} object.
#'
#' @examples
#' data(attitude_data)
#' fit <- csem_estimate(attitude_data[, 1:10])
#' csem_plot(fit)
#' csem_plot(fit, type = "items")
#' csem_plot(csem_compare(attitude_data[, 1:10]))
#'
#' @export
csem_plot <- function(x, type = NULL, benchmarks = c(0.90, 0.80),
                      title = NULL) {
  if (is.null(type))
    type <- switch(class(x)[1],
                   csem_bootstrap = "band", csem_delta = "band",
                   csem_reliability = "reliability", "csem")
  type <- match.arg(type, c("csem", "info", "items", "reliability", "band",
                            "difference"))
  th_lab <- expression(theta ~ "(standardized latent trait)")
  base <- ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor = ggplot2::element_blank())
  v <- x$var_theta %||% 1
  bench <- data.frame(y = sqrt(v * (1 / benchmarks - 1)),
                      lab = sprintf("rho = %.2f", benchmarks))
  add_bench <- function(p, xmax)
    p + ggplot2::geom_hline(data = bench, ggplot2::aes(yintercept = .data$y),
                            linetype = "dashed", colour = "grey45") +
      ggplot2::annotate("text", x = xmax, y = bench$y, label = bench$lab,
                        hjust = 0, vjust = -0.4, size = 3, colour = "grey30")

  # -- Comparison objects ----------------------------------------------------
  if (inherits(x, "csem_compare")) {
    pr <- x$profiles
    pr$model <- factor(pr$model, levels = x$models)
    pal <- .csem_palette[seq_along(x$models)]
    if (type == "info") {
      p <- ggplot2::ggplot(pr, ggplot2::aes(.data$theta, .data$test_info,
                                            colour = .data$model)) +
        ggplot2::geom_line(linewidth = 1) +
        ggplot2::scale_colour_manual(values = pal, name = NULL) +
        ggplot2::labs(x = th_lab, y = "Test information",
                      title = title %||% "Test information by model") + base
      return(p)
    }
    if (type == "difference") {
      if (is.null(x$differences))
        stop("Run csem_compare() with n_boot > 0 for type = 'difference'.",
             call. = FALSE)
      d <- x$differences
      p <- ggplot2::ggplot(d, ggplot2::aes(.data$theta, .data$diff)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                          ymax = .data$upper),
                             fill = "#0072B2", alpha = 0.25) +
        ggplot2::geom_line(colour = "#0072B2", linewidth = 1) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey30") +
        ggplot2::geom_hline(yintercept = c(-1, 1) * (x$min_diff %||% 0),
                            linetype = "dotted", colour = "grey45") +
        ggplot2::facet_wrap(~ model, labeller = ggplot2::as_labeller(
          function(m) paste(m, "-", x$reference))) +
        ggplot2::labs(x = th_lab, y = "Difference in CSEM",
                      title = title %||% "CSEM difference from reference model",
                      subtitle = sprintf("%.0f%% paired bootstrap band (%d replications); dotted lines: +/- %.2f",
                                         100 * x$conf_level, x$n_boot,
                                         x$min_diff %||% 0)) + base
      return(p)
    }
    p <- ggplot2::ggplot(pr, ggplot2::aes(.data$theta, .data$csem,
                                          colour = .data$model,
                                          fill = .data$model))
    if (type == "band" && "lower" %in% names(pr))
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                                 ymax = .data$upper),
                                    alpha = 0.18, colour = NA)
    p <- p + ggplot2::geom_line(linewidth = 1) +
      ggplot2::scale_colour_manual(values = pal, name = NULL) +
      ggplot2::scale_fill_manual(values = pal, name = NULL) +
      ggplot2::coord_cartesian(ylim = c(0, min(max(pr$csem), 1.5))) +
      ggplot2::labs(x = th_lab, y = "CSEM",
                    title = title %||% "CSEM profiles by model") + base
    return(add_bench(p, min(pr$theta)))
  }

  # -- Single-model objects --------------------------------------------------
  if (type == "items") {
    if (!inherits(x, "csem_fit"))
      stop("type = 'items' requires a csem_fit object.", call. = FALSE)
    ii <- x$item_info
    df <- data.frame(theta = rep(x$theta_grid, ncol(ii)),
                     info = as.vector(ii),
                     item = factor(rep(colnames(ii), each = nrow(ii)),
                                   levels = colnames(ii)))
    return(ggplot2::ggplot(df, ggplot2::aes(.data$theta, .data$info,
                                            group = .data$item)) +
             ggplot2::geom_line(colour = "#0072B2", alpha = 0.7) +
             ggplot2::labs(x = th_lab, y = "Item information",
                           title = title %||% paste("Item information -",
                                                    x$model_short)) + base)
  }

  if (type == "reliability") {
    if (!inherits(x, "csem_reliability")) x <- csem_reliability(x)
    df <- x$rel_df
    p <- ggplot2::ggplot(df, ggplot2::aes(.data$theta, .data$crel))
    if ("crel_lower" %in% names(df))
      p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$crel_lower,
                                                 ymax = .data$crel_upper),
                                    fill = "#0072B2", alpha = 0.2)
    p <- p + ggplot2::geom_line(colour = "#0072B2", linewidth = 1) +
      ggplot2::geom_hline(yintercept = benchmarks, linetype = "dashed",
                          colour = "grey45") +
      ggplot2::scale_y_continuous(limits = c(0, 1)) +
      ggplot2::labs(x = th_lab, y = "Conditional reliability",
                    title = title %||% paste("Conditional reliability -",
                                             x$model_short)) + base
    if (!is.null(x$marginal_reliability))
      p <- p + ggplot2::labs(subtitle = sprintf(
        "Marginal reliability = %.3f", x$marginal_reliability["information"]))
    return(p)
  }

  if (inherits(x, "csem_fit")) x <- csem_analytical(x)

  if (type == "info") {
    if (!inherits(x, "csem_analytical"))
      stop("type = 'info' requires csem_fit or csem_analytical.", call. = FALSE)
    return(ggplot2::ggplot(x$csem_df, ggplot2::aes(.data$theta,
                                                   .data$test_info)) +
             ggplot2::geom_line(colour = "#0072B2", linewidth = 1) +
             ggplot2::geom_vline(xintercept = x$theta_max_info,
                                 linetype = "dotted") +
             ggplot2::labs(x = th_lab, y = "Test information",
                           title = title %||% paste("Test information -",
                                                    x$model_short)) + base)
  }

  # CSEM, with band if available
  if (inherits(x, "csem_bootstrap")) {
    df <- data.frame(theta = x$theta_grid, csem = x$csem_analytical,
                     lower = x$csem_boot_lower, upper = x$csem_boot_upper)
    sub <- sprintf("%.0f%% %s bootstrap band (%d replications)",
                   100 * x$conf_level, x$type, x$n_boot)
  } else if (inherits(x, "csem_delta")) {
    df <- data.frame(theta = x$theta_grid, csem = x$csem,
                     lower = x$lower, upper = x$upper)
    sub <- sprintf("%.0f%% delta-method band", 100 * x$conf_level)
  } else {
    df <- data.frame(theta = x$theta_grid, csem = x$csem)
    sub <- NULL
  }
  if (type == "band" && is.null(df$lower))
    stop("type = 'band' requires a csem_bootstrap or csem_delta object.",
         call. = FALSE)
  p <- ggplot2::ggplot(df, ggplot2::aes(.data$theta, .data$csem))
  if (!is.null(df$lower))
    p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$lower,
                                               ymax = .data$upper),
                                  fill = "#0072B2", alpha = 0.22)
  p <- p + ggplot2::geom_line(colour = "#0072B2", linewidth = 1) +
    ggplot2::coord_cartesian(ylim = c(0, min(max(df$upper %||% df$csem), 1.5))) +
    ggplot2::labs(x = th_lab, y = "CSEM", subtitle = sub,
                  title = title %||% paste("Conditional SEM -",
                                           x$model_short)) + base
  add_bench(p, min(df$theta))
}
