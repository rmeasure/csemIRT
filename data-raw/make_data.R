# Regenerates the built-in datasets of csemIRT.
# attitude_data and character_data are carried over unchanged from csemGRM
# 0.1.0 (simulated from the GRM; see data.R), except that the 'group'
# variable is now assigned at random rather than by the sign of theta.
library(mirt)
load("data/attitude_data.rda"); load("data/character_data.rda")
set.seed(99)
attitude_data$group  <- sample(c("Group_1", "Group_2"), nrow(attitude_data), TRUE)
character_data$group <- sample(c("Group_1", "Group_2"), nrow(character_data), TRUE)

# wellbeing_data: 12 items, 4 categories, N = 600, generated from the GPCM
# with heterogeneous slopes and thresholds concentrated below the mean
# (a positively worded well-being scale with a ceiling effect).
set.seed(2026)
k  <- 12
a  <- c(0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 2.4, 0.7, 2.6)
b  <- t(sapply(seq_len(k), function(j) sort(c(-2.4, -1.4, -0.3) + rnorm(3, 0, 0.25))))
# mirt gpcm intercept parameterization: d_k = -a * sum_{v<=k} b_v (d0 = 0)
d  <- cbind(0, t(apply(b, 1, cumsum)) * -a)
th <- matrix(rnorm(600), ncol = 1)
w  <- simdata(a = matrix(a), d = d, N = 600, itemtype = "gpcm", Theta = th)
wellbeing_data <- as.data.frame(w + 1L)
names(wellbeing_data) <- sprintf("W%02d", 1:k)
wellbeing_data$group <- sample(c("Group_1", "Group_2"), 600, TRUE)

save(attitude_data,  file = "data/attitude_data.rda",  compress = "xz")
save(character_data, file = "data/character_data.rda", compress = "xz")
save(wellbeing_data, file = "data/wellbeing_data.rda", compress = "xz")
