# Addressing issue #5 using data & models
# from http://www.animalsciencepublications.org/publications/jas/abstracts/88/2/497

### inc
#library(lme4qtl)
library(devtools)
load_all("~/git/variani/lme4qtl/")

### download data
urls <- c("https://github.com/variani/lme4qtl/files/1632181/data.count.repeat.txt",
  "https://github.com/variani/lme4qtl/files/1632182/data.count.txt")

#ret <- lapply(urls, function(x) download.file(x, basename(x)))

### load data
data(dat40, package = "lme4qtl")

data.count <- read.table("data.count.txt", sep="\t")
data.count.repeat <- read.table("data.count.repeat.txt", sep="\t")

### models
m1 <- relmatGlmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), data.count.repeat, relmat = list(ID = kin2), family = poisson)

# drop `(1|FAMID)` (captures little variance)
m2 <- relmatGlmer(trait1 ~ AGE + SEX + (1|ID), data.count.repeat, relmat = list(ID = kin2), family = poisson)

# drop `AGE + SEX` (not significant)
m3 <- relmatGlmer(trait1 ~ (1|ID), data.count.repeat, relmat = list(ID = kin2), family = poisson)

### print
#> packageVersion("lme4qtl")
#[1] ‘0.1.10’

### print warnings
#> m1@optinfo$conv$lme4$messages
#[[1]]
#[1] "Model failed to converge with max|grad| = 0.00340346 (tol = 0.001, component 1)"
#
#[[2]]
#[1] "Model is nearly unidentifiable: very large eigenvalue\n - Rescale variables?"

#> m2@optinfo$conv$lme4$messages
#[[1]]
#[1] "Model is nearly unidentifiable: very large eigenvalue\n - Rescale variables?"

#> m3@optinfo$conv$lme4$messages
# NULL
