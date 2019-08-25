
<!-- README.md is generated from README.Rmd. Please edit that file -->
lme4qtl
=======

<!-- badges: start -->
[![travis-ci build status](https://travis-ci.org/variani/lme4qtl.svg?branch=master)](https://travis-ci.org/variani/lme4qtl) <!-- badges: end -->

`lme4qtl` extends the [lme4](https://github.com/lme4/lme4) R package for quantitative trait locus (qtl) mapping. It is all about the covariance structure of random effects. `lme4qtl` supports user-defined matrices for that, e.g. kinship or IBDs.

See slides [bit.ly/1UiTZvQ](http://bit.ly/1UiTZvQ) introducing the `lme4qtl` R package or read our [article](http://dx.doi.org/10.1186/s12859-018-2057-x) / [preprint](http://www.biorxiv.org/content/early/2017/08/31/139816).

<table style="width:46%;">
<colgroup>
<col width="15%" />
<col width="30%" />
</colgroup>
<thead>
<tr class="header">
<th>Package</th>
<th>Continuous response</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>stats</td>
<td><code>lm(myTrait ~ myCovariate, myData)</code></td>
</tr>
<tr class="even">
<td>lme4</td>
<td><code>lmer(myTrait ~ myCovariate + (1\|myID), myData)</code></td>
</tr>
<tr class="odd">
<td>lme4qtl</td>
<td><code>relmatLmer(myTrait ~ myCovariate + (1\|myID), myData, relmat = list(myID = myMatrix))</code></td>
</tr>
</tbody>
</table>

<table style="width:46%;">
<colgroup>
<col width="15%" />
<col width="30%" />
</colgroup>
<thead>
<tr class="header">
<th>Package</th>
<th>Binary response</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>stats</td>
<td><code>glm(myStatus ~ 1, myData, family = binomial)</code></td>
</tr>
<tr class="even">
<td>lme4</td>
<td><code>glmer(myStatus ~ (1\|myID), myData, family = binomial)</code></td>
</tr>
<tr class="odd">
<td>lme4qtl</td>
<td><code>relmatGlmer(myStatus ~ (1\|myID), myData, relmat = list(myID = myMatrix), family = binomial)</code></td>
</tr>
</tbody>
</table>

Note that rownames/colnames of `myMatrix` have to be values of `myID` variable, so matching between relationship matrix and grouping variable is possible. The order doesn't matter.

Installation
------------

You can install the development version from [GitHub](https://github.com/variani/lme4qtl) with:

``` r
# install.packages("devtools")
devtools::install_github("variani/lme4qtl")
```

The official release on [CRAN](https://CRAN.R-project.org) is [pending](https://github.com/variani/lme4qtl/issues/9).

Example
-------

``` r
library(lme4qtl)

# load simulated data: phenotypes in `dat40` and kinship in `kin2`
data(dat40, package = "lme4qtl")

# fit a model for continuous trait 
m1 <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))
#> boundary (singular) fit: see ?isSingular
#> boundary (singular) fit: see ?isSingular
m1
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: trait1 ~ AGE + SEX + (1 | FAMID) + (1 | ID)
#>    Data: dat40
#> REML criterion at convergence: 963.3853
#> Random effects:
#>  Groups   Name        Std.Dev.
#>  ID       (Intercept) 2.2988  
#>  FAMID    (Intercept) 0.0000  
#>  Residual             0.7856  
#> Number of obs: 224, groups:  ID, 224; FAMID, 39
#> Fixed Effects:
#> (Intercept)          AGE         SEX2  
#>    7.563248     0.008314    -0.364197  
#> convergence code 0; 1 optimizer warnings; 0 lme4 warnings

# get an estimate of h2, the proportion of variance explained by (1|ID)
lme4::VarCorr(m1)
#>  Groups   Name        Std.Dev.
#>  ID       (Intercept) 2.29880 
#>  FAMID    (Intercept) 0.00000 
#>  Residual             0.78562
lme4qtl::VarProp(m1)
#>        grp        var1 var2      vcov     sdcor      prop
#> 1       ID (Intercept) <NA> 5.2845002 2.2988041 0.8954191
#> 2    FAMID (Intercept) <NA> 0.0000000 0.0000000 0.0000000
#> 3 Residual        <NA> <NA> 0.6172059 0.7856245 0.1045809

# fir a model model for binary trait
m2 <- relmatGlmer(trait1bin ~ (1|ID), dat40, relmat = list(ID = kin2), family = binomial)
m2
#> Generalized linear mixed model fit by maximum likelihood (Laplace
#>   Approximation) [glmerMod]
#>  Family: binomial  ( logit )
#> Formula: trait1bin ~ (1 | ID)
#>    Data: dat40
#>      AIC      BIC   logLik deviance df.resid 
#> 175.5894 182.5000 -85.7947 171.5894      232 
#> Random effects:
#>  Groups Name        Std.Dev.
#>  ID     (Intercept) 51.09   
#> Number of obs: 234, groups:  ID, 234
#> Fixed Effects:
#> (Intercept)  
#>      -27.88
```

Citation
--------

To cite the `lme4qtl` package in publications use:

      Ziyatdinov et al., lme4qtl: linear mixed models with flexible
      covariance structure for genetic studies of related individuals, 
      BMC Bioinformatics (2018)

Contact
-------

You are welcome to submit suggestions and bug-reports at <https://github.com/variani/lme4qtl/issues>.
