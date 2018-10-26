### inc
library(dplyr)
library(magrittr)

library(ggplot2)
library(cowplot)

### par
num_ind <- 100
num_rep <- 5
num_obs <- num_ind * num_rep

num_gr <- 5

### scenario 1: y ~ time + (1|id) + (0 + time||id_slope)
# simulate data
id_sim <- seq(1, num_ind) %>% as.character
gr_sim <- seq(1, num_gr) %>% as.character
time_sim <- seq(1, num_rep)

dat_sim <- data_frame(
  id = rep(id_sim, each = num_rep), 
  id_slope = id,
  gr = sample(gr_sim, num_ind, replace = TRUE) %>% rep(each = num_rep),
  time = rep(time_sim, num_ind),
  slope_fixed = 1,
  slope_rand = rep(rnorm(num_ind), each = num_rep),
  intercept_rand = rep(rnorm(num_ind), each = num_rep))

dat_sim_time1 <- filter(dat_sim, time == time_sim[1]) %>% select(-time)

dat_sim <- within(dat_sim, {
  y <- intercept_rand + (slope_fixed + slope_rand) * time + rnorm(num_obs)
})

# fit model
m <- lmer(y ~ time + (1|id) + (0 + time||id_slope), dat_sim)

# get estmiates
slope_rand_pred <- ranef(m, whichel = "id_slope", drop = TRUE) %$% id_slope

dat_pred <- data_frame(id = names(slope_rand_pred), slope_rand_pred = slope_rand_pred) %>%
  left_join(dat_sim_time1, by = "id")

p1 <- ggplot(dat_pred, aes(gr, slope_rand)) + geom_boxplot() + geom_hline(yintercept = 0, linetype = 3)

### scenario 2: y ~ time + (0 + time||gr)
# simulate data
id_sim <- seq(1, num_ind) %>% as.character
gr_sim <- seq(1, num_gr) %>% as.character
time_sim <- seq(1, num_rep)

intercept_rand_sim <- rnorm(num_gr)
names(intercept_rand_sim) <- gr_sim

dat_sim2 <- data_frame(
  id = rep(id_sim, each = num_rep), 
  gr = sample(gr_sim, num_ind, replace = TRUE) %>% rep(each = num_rep),
  time = rep(time_sim, num_ind),
  slope_fixed = 1,
  slope_rand = intercept_rand_sim[gr])

dat_sim2_time1 <- filter(dat_sim2, time == time_sim[1]) %>% select(-time)

dat_sim2 <- within(dat_sim2, {
  y <- (slope_fixed + slope_rand) * time + rnorm(num_obs)
})

# fit model
m2 <- lmer(y ~ time + (0 + time||gr), dat_sim2)

# get estmiates
slope_rand_pred <- ranef(m2, whichel = "gr", drop = TRUE) %$% gr

dat_pred2 <- data_frame(gr = names(slope_rand_pred), slope_rand_pred = slope_rand_pred) %>%
  right_join(dat_sim2_time1, by = "gr")

p2 <- ggplot(dat_pred2, aes(gr, slope_rand_pred)) + geom_boxplot() + geom_hline(yintercept = 0, linetype = 3)

### scenario 3: y ~ time + (0 + time||gr) 
# - the same data simulation (`dat_sim2` is used again)
# - fitted by relamtLmer
# simulate data
mat_id <- model.matrix(~ -1 + gr, dat_sim2_time1) %>% Matrix %>% tcrossprod
rownames(mat_id) <- dat_sim2_time1$id
colnames(mat_id) <- dat_sim2_time1$id

# fit model
m3 <- relmatLmer(y ~ time + (0 + time||id), dat_sim2, relmat = list(id = mat_id))

# get estmiates
slope_rand_pred <- ranef(m3, whichel = "id", drop = TRUE) %$% id
names_slope <- names(slope_rand_pred)

relfac <- m3@optinfo$relmat$relfac %$% id 
relfac <- relfac[names_slope, names_slope] # to make sure names are the same

# see https://github.com/cran/pedigreemm/blob/master/R/pedigree.R#L367 
slope_rand_pred <- crossprod(slope_rand_pred, relfac) %>% as.numeric
names(slope_rand_pred) <- names_slope

dat_pred3 <- data_frame(id = names(slope_rand_pred), slope_rand_pred = slope_rand_pred) %>%
  left_join(dat_sim2_time1, by = "id")

p3 <- ggplot(dat_pred3, aes(gr, slope_rand_pred)) + geom_boxplot() + geom_hline(yintercept = 0, linetype = 3)

### final plot: models `m2` and `m3` produce the same results
#plot_grid(p2, p3)
