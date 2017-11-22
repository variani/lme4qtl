## About lme4qtl

[![travis-ci build status](https://travis-ci.org/variani/lme4qtl.svg?branch=master)](https://travis-ci.org/variani/lme4qtl)


lme4qtl extends the [lme4](https://github.com/lme4/lme4) R package for quantitative trait locus (qtl) mapping. It is all about the covariance structure of random effects. `lme4qtl` supports user-defined matrices for that,
e.g. kinship or IBDs.

See slides [bit.ly/1UiTZvQ](http://bit.ly/1UiTZvQ) introducing the `lme4qtl` R package or read our preprint at [biorxiv](http://www.biorxiv.org/content/early/2017/08/31/139816).

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

Note that rownames/colnames of `myMatrix` have to be values of `myID` variable, so matching between relationship matrix and grouping variable is possible. The order doesn't matter.

## Quick start

```
library(lme4) # needed for `VarCorr` function
library(lme4qtl)

# load synthetic data set `dat40` distributed within `lme4qtl`
# - table of phenotypes `dat40`
# - the double kinship matrix `kin2`
data(dat40)

# (1) model continiuous trait `trait1`
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

# (2) model binary trait `trait1bin`
gmod <- relmatGlmer(trait1bin ~ (1|ID), dat40, relmat = list(ID = kin2), family = binomial)
```

## Installation

```
# install.packages("devtools")
library(devtools)
install_github("variani/lme4qtl")
```

## Contact

You are welcome to submit suggestions and bug-reports at https://github.com/variani/lme4qtl/issues.
