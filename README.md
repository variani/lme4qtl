## About lme4qtl

lme4qtl extends the [lme4](https://github.com/lme4/lme4) R package for quantitative trait locus (qtl) mapping. It is all about the covariance structure of random effects. `lme4qtl` allows to define this structure via matrices,
e.g. the kinship or IBD matrices.

See slides [bit.ly/1UiTZvQ](http://bit.ly/1UiTZvQ) introducing the lme4qtl R package.

The release of lme4qtl is scheduled according to manuscript submission. 

|  Package | Continuous response |
|----------|---------------------|
| stats   | `lm(myTrait ~ myCovariate, myData)` |
| lme4    | `lmer(myTrait ~ myCovariate + (1\|myID), myData)` |
| lme4qtl | `relmatLmer(myTrait ~ myCovariate + (1\|myID), myData, relmat = list(myID = myMatrix))` |

|  Package | Binary response |
|----------|---------------------|
| stats    | `glm(myStatus ~ 1, myData, family = binomial)` |
| lme4    | `glmer(myStatus ~ (1\|myID), myData, family = binomial)` |
| lme4qtl | `relmatGlmer(myStatus ~ (1\|myID), myData, relmat = list(myID = myMatrix), family = binomial)` |


## Quick start

Source: [inst/examples/package-lme4qtl.R](inst/examples/package-lme4qtl.R)

```
library(lme4) # needed for `VarCorr` function
library(lme4qtl)

### load the synthetic data: 
# - table of phenotypes `dat40`
# - the double kinship matrix `kin2`
data(dat40)

### model the continiuous trait `trait1`
mod <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))

# get the estimation of h2
(vf <- as.data.frame(VarCorr(mod))[, c("grp", "vcov")])
#       grp      vcov
#1       ID 5.2845001
#2    FAMID 0.0000000
#3 Residual 0.6172059

prop <- with(vf, vcov / sum(vcov))

(h2 <- prop[1]) 
#[1] `0.895419`

### model the binary trait `trait1bin`
# - Model is nearly unidentifiable, when `(1|FAMID)` effect is added
gmod <- relmatGlmer(trait1bin ~ AGE + SEX + (1|ID), dat40, relmat = list(ID = kin2), family = binomial)
```

## Install

### Install from a local file

A source file is something like `lme4qtl_0.1.5.zip` or `lme4qtl_0.1.5_R_x86_64-pc-linux-gnu.tar.gz`.

```
# install dependencies, if necessary
p <- c("Matrix", "lme4", "kinship2", "plyr", "ggplot2")
install.packages(p)

# install lme4qtl from a local file
# - the first argument is path to the installation file
# - the second argument (`repos = NULL`) says to install locally rather than from a repository
install.packages("lme4qtl_0.1.5.zip", repos = NULL)
```
