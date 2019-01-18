# A version of fast GWAS, when 
# (i) mixed model is fitted once (without any SNP as covariate);
# (ii) each SNP is tested using variance components etimated at step (i).

### install
# library(devtools)
# install_github("variani/matlm")
#
# install_github("variani/wlm")
# install_github("variani/qq")
# 
# re-install lme4qtl if necessary (e.g. for `varcov` fun. was not exported prev.)
# install_github("variani/lme4qtl")

### inc
library(Matrix)
library(MASS)
library(ggplot2)
# for Step 1: linear mixed model (no SNPs)
library(lme4qtl) 
# for Step 2: association tests
library(matlm) 
library(wlm) 
# for Step 3: explore GWAS results
library(qq)

### par
M <- 1500 # the number of SNPs

### data
data(dat40, package = "lme4qtl")
N <- nrow(dat40) # the number of ind.
ids <- dat40$ID # individual IDs

# simulate (null) random genotypes (0/1 with 50% prob.; not linked to `kin2`)
gdat40 <- matrix(rbinom(N*M, 1, 0.5), nrow = N, ncol = M)
rownames(gdat40) <- ids
colnames(gdat40) <- paste0("snp", seq(1, ncol(gdat40)))

# simulate (null) random cont. predictors (linked to `kin2`) by MVN ~ (0, kin2)
pdat40 <- t(mvrnorm(M, rep(0, N), kin2))
rownames(gdat40) <- ids
colnames(gdat40) <- paste0("pred", seq(1, ncol(gdat40)))

#-----------------------------------------------
# Step 1: 
# - fit base model (no SNPs) 
# - extract variance-covariance matrix `V`
#-----------------------------------------------
mod <- lme4qtl::relmatLmer(trait1 ~ AGE + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))
# in case of convergence problems, one might try/experiment with:
# - relmatLmer(trait1 ~ AGE + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2), calc.derivs = FALSE)
# - relmatLmer(trait1 ~ AGE + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2), calc.derivs = FALSE, control = lmerControl(optimizer = "bobyqa", check.conv.grad = list(action = "warning", tol = 0.005, relTol = NULL)))

V <- lme4qtl::varcov(mod, idvar = "ID")

Matrix::image(V[1:20, 1:20], main = "Estimated V (with artifacts)") # some artifacts close to zero due to limited numeric precision

# get rid of the artifacts and see the expected matrix of family blocks
V_thr <- V
V_thr[abs(V) < 1e-10] <- 0
Matrix::image(V_thr[1:20, 1:20], main = "Estimated V (with artifacts removed)")

#-------
# Step 2:
# - perform association tests on M predictors
# - examimed several combinations:
#   - linear models: least squares (LS) vs. generalized least squares (GLS) that takes V as input
#   - binary genotypes (simulated with no structure) vs. cont. predictors (linked to kin2)
#-------
# transformation on data (due to structure in V) needs to be computed once 
# (note: EVD (not Cholesky) is required; otherwise, missing data would produce messy results)
decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
W <- decomp$transform

gassoc_lm <- matlm::matlm(trait1 ~ AGE, dat40, pred = gdat40, ids = ids,
  batch_size = 100, verbose = 2)
gassoc_gls <- matlm::matlm(trait1 ~ AGE, dat40, pred = gdat40, ids = ids, transform = W, 
  batch_size = 100, verbose = 2)

passoc_lm <- matlm::matlm(trait1 ~ AGE, dat40, pred = pdat40, ids = ids,
  batch_size = 100, verbose = 2)
passoc_gls <- matlm::matlm(trait1 ~ AGE, dat40, pred = pdat40, ids = ids, transform = W, 
  batch_size = 100, verbose = 2)
  
#-------
# Step 3:
# - QQ plots
#   - the first two plots show that both LS & GLS approches give valid results,
#     because binary genotype predictors (gdat40) were simulated without any link to 
#     data structure in kin2
#   - the last plots show that GLS produces a valid distribution of p-values,
#     while LS shows an inflated Type I error rate
#-------
qq::qq_plot(gassoc_lm$tab$pval) + ggtitle("LS: (null) random binary genotypes (not linked to kin2)")
qq::qq_plot(gassoc_gls$tab$pval) + ggtitle("LMM/WLS: (null) random binary genotypes (not linked to kin2)")

qq::qq_plot(passoc_lm$tab$pval) + ggtitle("LS: (null) random cont. predictors (linked to kin2)")
qq::qq_plot(passoc_gls$tab$pval) + ggtitle("LMM/WLS: (null) random cont. predictors (linked to kin2)")



