suppressMessages(library(csemIRT)); library(ggplot2)
t0 <- Sys.time()
it <- attitude_data[, 1:10]
f   <- csem_estimate(it, "graded")
an  <- csem_analytical(f)
dl  <- csem_delta(f)
bp  <- csem_bootstrap(f, n_boot = 1000, type = "parametric",    seed = 2026, cores = 2)
bn  <- csem_bootstrap(f, n_boot = 1000, type = "nonparametric", seed = 2026, cores = 2)
saveRDS(list(f=f, an=an, dl=dl, bp=bp, bn=bn), "ex1.rds")
cat("ex1 done", format(Sys.time()-t0), "\n")
cmp <- csem_compare(wellbeing_data[, 1:12], models = c("graded","gpcm","pcm"),
                    n_boot = 1000, seed = 2026, cores = 2)
saveRDS(cmp, "ex2.rds")
cat("ex2 done", format(Sys.time()-t0), "\n")
