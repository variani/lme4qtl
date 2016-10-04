## About lme4qtl

lme4qtl extends the [lme4](https://github.com/lme4/lme4) R package for quantitative trait locus (qtl) mapping.

Please see the slides [bit.ly/1UiTZvQ](http://bit.ly/1UiTZvQ) introducing the lme4qtl R package.

The release of lme4qtl is scheduled according to the manuscript submission. 

## Quick start

Source: [inst/examples/package-lme4qtl.R](inst/examples/package-lme4qtl.R)

```
library(lme4) # needed for `VarCorr` function
library(lme4qtl)

# load the synthetic data: 
# - table of phenotypes `dat40`
# - the double kinship matrix `kin2`
data(dat40)

# run the model
mod <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))

# get the estimation of h2
(vf <- as.data.frame(VarCorr(mod))[, c("grp", "vcov")])
#       grp      vcov
#1       ID 5.2845001
#2    FAMID 0.0000000
#3 Residual 0.6172059

prop <- with(vf, vcov / sum(vcov))

(h2 <- prop[1]) 
#[1] `0.8954191`
```

## Install

### Install from a local file

Instalation from a source file, e.g. `lme4qtl_0.1.4.tar.gz`, includes two parts:

* install all dependencies of the `lme4qtl` package;
* install the `lme4qtl` package itself from the source.

The standard `install.packages` R function will not help. The right function to do the job is `devtools::install`, which in turn needs a source directory istead of a compressed file. Thus, there are two steps:

Step 1. Unpack the source file to a a directory, e.g. from the Linux terminal:

```
tar -xf lme4qtl_0.1.4.tar.gz
```

The output directory will be `lme4qtl`.

Step 2. Install the package from the R console:

```
library(devtools)
install("lme4qtl/", dependencies = TRUE)
```

If the fresh packages-dependencies are needed (suggested for the first-time installation), then type:

```
install("lme4qtl/", dependencies = TRUE, upgrade_dependencies = TRUE)
```

If `devtools` is not installed, type `install.packages("devtools")` before.
