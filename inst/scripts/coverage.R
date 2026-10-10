# Monte Carlo check: coverage of delta-method and parametric-bootstrap bands
suppressMessages(library(csemIRT)); library(mirt)
dir.create("cov_parts", showWarnings = FALSE)
z <- seq(-3, 3, by = 0.25)
true_csem <- function(dat, itemtype, values_fun, sd_true) {
  sv <- mirt(dat, 1, itemtype = itemtype, pars = "values", verbose = FALSE)
  sv <- values_fun(sv)
  m  <- mirt(dat, 1, itemtype = itemtype, pars = sv, TOL = NaN, verbose = FALSE)
  1 / sqrt(sd_true^2 * testinfo(m, matrix(sd_true * z)))
}
# Condition 1: GRM, attitude-scale generating parameters (N = 500)
a1 <- c(1.8,1.5,2.0,1.6,1.3,1.7,1.9,1.4,1.6,1.8)
b1 <- matrix(c(-2.0,-1.0,0.0,1.0, -1.8,-0.8,0.2,1.2, -2.2,-1.2,-0.2,0.8,
               -1.6,-0.6,0.4,1.4, -1.5,-0.5,0.5,1.5, -1.7,-0.7,0.3,1.3,
               -2.1,-1.1,-0.1,0.9, -1.4,-0.4,0.6,1.6, -1.9,-0.9,0.1,1.1,
               -1.6,-0.6,0.4,1.4), 10, byrow = TRUE)
d1 <- -a1 * b1
# Condition 2: PCM with latent variance 1.5 (N = 500)
set.seed(11)
b2 <- t(sapply(1:10, function(j) sort(c(-1.5,-0.5,0.5,1.5) + rnorm(4, 0, .3))))
d2 <- cbind(0, -t(apply(b2, 1, cumsum)))
conds <- list(
  GRM = list(itemtype = "graded", model = "graded", sd = 1,
             gen = function(N) simdata(matrix(a1), d1, N, "graded") + 1,
             vals = function(sv) { for (j in 1:10) {
               r <- sv$item == sprintf("Item_%d", j)
               sv$value[r & sv$name == "a1"] <- a1[j]
               for (k in 1:4) sv$value[r & sv$name == paste0("d", k)] <- d1[j, k] }
               sv }),
  PCM = list(itemtype = "Rasch", model = "pcm", sd = sqrt(1.5),
             gen = function(N) simdata(matrix(1, 10), d2, N, "gpcm",
                                       Theta = matrix(rnorm(N, 0, sqrt(1.5)))) + 1,
             vals = function(sv) { for (j in 1:10) {
               r <- sv$item == sprintf("Item_%d", j)
               for (k in 0:4) sv$value[r & sv$name == paste0("d", k)] <- d2[j, k + 1] }
               sv$value[sv$name == "COV_11"] <- 1.5
               sv }))
R <- 60; B <- 200
for (cn in names(conds)) {
  cd <- conds[[cn]]
  for (r in seq_len(R)) {
    fn <- sprintf("cov_parts/%s_%03d.rds", cn, r)
    if (file.exists(fn)) next
    set.seed(10000 * match(cn, names(conds)) + r)
    dat <- as.data.frame(cd$gen(500)); names(dat) <- sprintf("Item_%d", 1:10)
    if (any(sapply(dat, function(x) length(unique(x))) < 5)) { saveRDS(NULL, fn); next }
    truth <- true_csem(dat, cd$itemtype, cd$vals, cd$sd)
    f  <- suppressWarnings(csem_estimate(dat, cd$model, theta_range = c(-3, 3),
                                         theta_points = length(z)))
    dl <- csem_delta(f)
    bt <- csem_bootstrap(f, n_boot = B, seed = r, cores = 2)
    saveRDS(data.frame(cond = cn, rep = r, theta = z, truth = truth, est = dl$csem,
               d_lo = dl$lower, d_hi = dl$upper,
               b_lo = bt$csem_boot_lower, b_hi = bt$csem_boot_upper,
               d_se = dl$se, b_se = bt$csem_boot_se), fn)
    cat(cn, r, format(Sys.time()), "\n")
  }
}
