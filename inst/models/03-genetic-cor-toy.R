### inc
library(lme4)

library(dplyr)
library(magrittr)

### par
p <- 3 # pedigree size
k <- 1 # n. pedigrees
t <- 2 # n. traits
nt <- p*k # n. samples per trait
n <- p*k*t # n. samples overall

### data
set.seed(1)
dat <- data_frame(y = runif(n), t = rep(paste0("t", 1:t), each = nt))
