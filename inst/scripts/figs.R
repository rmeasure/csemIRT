suppressMessages(library(csemIRT)); library(ggplot2)
x   <- readRDS("ex1.rds"); cmp <- readRDS("ex2.rds")
cmp$min_diff <- 0.05
d <- cmp$differences; d$relevant <- d$lower > 0.05 | d$upper < -0.05; cmp$differences <- d
th <- theme(plot.title = element_text(size = 11, face = "bold"),
            plot.subtitle = element_text(size = 9))
# Figure 1: GRM CSEM with parametric bootstrap band, delta band overlaid
bp <- x$bp$csem_df; dl <- x$dl$csem_df
df <- rbind(data.frame(theta = bp$theta, csem = bp$csem_analytical, lower = bp$csem_lower,
                       upper = bp$csem_upper, band = "Parametric bootstrap (B = 1000)"),
            data.frame(theta = dl$theta, csem = dl$csem, lower = dl$lower, upper = dl$upper,
                       band = "Delta method"))
bench <- data.frame(y = c(1/3, 0.5), lab = c("rho = .90", "rho = .80"))
p1 <- ggplot(df, aes(theta, csem)) +
  geom_ribbon(data = df[df$band != "Delta method", ], aes(ymin = lower, ymax = upper),
              fill = "#0072B2", alpha = .25) +
  geom_line(data = df[df$band == "Delta method", ], aes(y = lower), linetype = "dashed", colour = "#D55E00") +
  geom_line(data = df[df$band == "Delta method", ], aes(y = upper), linetype = "dashed", colour = "#D55E00") +
  geom_line(data = df[df$band != "Delta method", ], colour = "#0072B2", linewidth = 1) +
  geom_hline(data = bench, aes(yintercept = y), linetype = "dotted", colour = "grey40") +
  annotate("text", x = -4, y = bench$y, label = bench$lab, hjust = 0, vjust = -0.4, size = 3, colour = "grey30") +
  coord_cartesian(ylim = c(0, 1.2), xlim = c(-4, 4)) +
  labs(x = expression(theta), y = "CSEM",
       title = "GRM, attitude_data: analytical CSEM with 95% bands",
       subtitle = "Shaded: parametric bootstrap percentile band (B = 1000); dashed: delta-method band") +
  theme_minimal(base_size = 11) + theme(panel.grid.minor = element_blank()) + th
ggsave("Figure1.png", p1, width = 6.5, height = 4, dpi = 300, bg = "white")
# Figure 2: model comparison (two panels)
p2a <- csem_plot(cmp, type = "band", title = "(a) CSEM profiles with 95% bootstrap bands") +
  coord_cartesian(ylim = c(0, 1.2), xlim = c(-4, 3)) + th
p2b <- csem_plot(cmp, type = "difference", title = "(b) Difference from GRM") +
  coord_cartesian(ylim = c(-0.45, 0.25), xlim = c(-4, 3)) + th
g <- gridExtra::arrangeGrob(p2a, p2b, ncol = 1, heights = c(1, 0.85))
ggsave("Figure2.png", g, width = 6.5, height = 7.5, dpi = 300, bg = "white")
print(cmp)
