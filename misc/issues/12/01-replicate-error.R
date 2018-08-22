### inc
library(lme4)
library(magrittr)

### (0) simulated data (grouping)
set.seed(1)
n_total <- 750
n_gr <- 10
n_miss <- 10
var_explained <- 0.5

ids <- seq(1, n_total) %>% as.character
x_gr <- sample(1:n_gr, n_total, replace = TRUE) # grouping variable (random effect)
x_rand <- rnorm(n_total) # random variable (fixed effect), not linked to outcome

x_rand_miss <- x_rand # add missing values
ind_miss <- sample(seq(1, n_total), n_miss, replace = FALSE)
x_rand_miss[ind_miss] <- NA

u_gr <- rnorm(n_gr)[x_gr] # simulate `n_gr` random values and use `x_gr` as group indicator
u_res <- rnorm(n_total)

y <- sqrt(var_explained) * u_gr + sqrt(1 - var_explained) * u_res

dat <- data.frame(y = y, id = ids, x_gr = x_gr, x_rand = x_rand, x_rand_miss = x_rand_miss)

# var-cov matrix for grouping variable `x_gr` (needed to test `lme4qtl` models)
R_gr <- sapply(x_gr, function(x) as.numeric(x == x_gr)) 
rownames(R_gr) <- ids
colnames(R_gr) <- ids

### (1) fit data (no missing) with `lme4`
m0 <- lmer(y ~ (1|x_gr), dat)
summary(m0)

m1 <- lmer(y ~ x_rand + (1|x_gr), dat)

anova(m0, m1)

### (2) fit data (no missing) with `lme4qtl`
m0 <- relmatLmer(y ~ (1|id), dat, relmat = list(id = R_gr))
m1 <- relmatLmer(y ~ x_rand + (1|id), dat, relmat = list(id = R_gr))
anova(m0, m1)

### (3) fit data (missing) with `lme4`
m0 <- lmer(y ~ (1|x_gr), dat)
m1 <- lmer(y ~ x_rand_miss + (1|x_gr), dat)

try(anova(m0, m1))

# the reason of the error
nobs(m0)
nobs(m1)

### (4) fit data (missing) with `lme4qtl`
m0 <- relmatLmer(y ~ (1|id), dat, relmat = list(id = R_gr))
m1 <- relmatLmer(y ~ x_rand_miss + (1|id), dat, relmat = list(id = R_gr))
try(anova(m0, m1))

### (5) fit data (missing) with `lme4qtl` + use `subset` 
# to filter out missing values of a variable under testing
dat_subset <- subset(dat, !is.na(x_rand_miss))

m0 <- relmatLmer(y ~ (1|id), dat_subset, relmat = list(id = R_gr))
m1 <- relmatLmer(y ~ x_rand_miss + (1|id), dat_subset, relmat = list(id = R_gr))
anova(m0, m1)



