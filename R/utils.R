# -- Internal helpers ----------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

# Supported polytomous models and their mirt itemtype codes
.csem_models <- function() {
  data.frame(
    model    = c("graded", "gpcm", "pcm", "grsm", "rsm"),
    itemtype = c("graded", "gpcm", "Rasch", "grsm", "rsm"),
    label    = c("Graded Response Model (GRM)",
                 "Generalized Partial Credit Model (GPCM)",
                 "Partial Credit Model (PCM)",
                 "Graded Rating Scale Model (GRSM)",
                 "Rating Scale Model (RSM)"),
    short    = c("GRM", "GPCM", "PCM", "GRSM", "RSM"),
    stringsAsFactors = FALSE
  )
}

.match_model <- function(model) {
  if (length(model) != 1L || !is.character(model))
    stop("'model' must be a single character string.", call. = FALSE)
  m <- tolower(model)
  aliases <- c(grm = "graded", samejima = "graded", rasch = "pcm",
               partialcredit = "pcm", ratingscale = "rsm")
  if (m %in% names(aliases)) m <- aliases[[m]]
  tab <- .csem_models()
  if (!m %in% tab$model)
    stop("'model' must be one of: ", paste(tab$model, collapse = ", "),
         " (aliases: 'grm', 'rasch').", call. = FALSE)
  as.list(tab[tab$model == m, ])
}

# Latent variance of a fitted unidimensional mirt model
.latent_var <- function(mod) {
  as.numeric(mirt::coef(mod, simplify = TRUE)$cov[1, 1])
}

# Test and item information on a standardized theta grid z.
# The model metric is theta = s * z, where s = latent SD (s = 1 for models
# that fix the latent variance). Information re-expressed in the z metric is
# I_z(z) = s^2 * I_theta(s * z), so CSEM_z(z) = CSEM_theta(s * z) / s.
.info_grid <- function(mod, z, s, items = FALSE) {
  th   <- matrix(s * z, ncol = 1)
  tinf <- s^2 * as.numeric(mirt::testinfo(mod, th))
  if (!items) return(tinf)
  nit  <- mirt::extract.mirt(mod, "nitems")
  iinf <- vapply(seq_len(nit), function(j)
    s^2 * as.numeric(mirt::iteminfo(mirt::extract.item(mod, j), th)),
    numeric(length(z)))
  colnames(iinf) <- colnames(mirt::extract.mirt(mod, "data"))
  list(test = tinf, items = iinf)
}

.csem_from_info <- function(info) 1 / sqrt(pmax(info, 1e-10))

# Number of observed categories per item
.n_cats <- function(data) {
  vapply(data, function(x) length(unique(x[!is.na(x)])), integer(1))
}

# Fit a unidimensional model quietly; returns NULL on error / non-convergence
.fit_quiet <- function(data, itemtype, SE = FALSE, technical = list()) {
  mod <- tryCatch(
    suppressWarnings(suppressMessages(
      mirt::mirt(data, 1, itemtype = itemtype, SE = SE, verbose = FALSE,
                 technical = technical))),
    error = function(e) NULL)
  if (is.null(mod)) return(NULL)
  if (!isTRUE(mirt::extract.mirt(mod, "converged"))) return(NULL)
  mod
}

# Indices of grid points closest to selected theta values (for printing)
.show_idx <- function(grid, at = c(-3, -2, -1, 0, 1, 2, 3)) {
  at <- at[at >= min(grid) & at <= max(grid)]
  unique(vapply(at, function(a) which.min(abs(grid - a)), integer(1)))
}

# Lapply that optionally runs in parallel (forked; falls back to serial)
.lapply_par <- function(X, FUN, cores = 1L) {
  if (cores > 1L && .Platform$OS.type != "windows")
    parallel::mclapply(X, FUN, mc.cores = cores, mc.preschedule = TRUE)
  else lapply(X, FUN)
}

# Okabe-Ito colour-blind safe palette
.csem_palette <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7",
                   "#E69F00", "#56B4E9", "#000000")
