### install
# library(devtools)
# install_github("variani/matlm")
# install_github("variani/wlm")
# install_github("variani/qq")

### inc
library(Matrix)
library(MASS)
library(ggplot2)
# for Step 1: linear mixed model
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

# simulate (null) random genotypes (0/1 with 50% prob.)
gdat40 <- matrix(rbinom(N*M, 1, 0.5), nrow = N, ncol = M)
rownames(gdat40) <- ids
colnames(gdat40) <- paste0("snp", seq(1, ncol(gdat40)))

pdat40 <- t(mvrnorm(M, rep(0, N), kin2))
rownames(gdat40) <- ids
colnames(gdat40) <- paste0("pred", seq(1, ncol(gdat40)))

#-----------------------------------------------
# Step 1: 
# - fit base model (not SNPs) 
# - extract variance-covariance matrix `V`
#-----------------------------------------------
mod <- lme4qtl::relmatLmer(trait1 ~ AGE + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))

V <- lme4qtl::varcov(mod, idvar = "ID")

Matrix::image(V[1:20, 1:20], main = "Estimated V (with artifacts)") # some artifacts close to zero due to limited numeric precision

# get rid of the artifacts and see the expected matrix of family blocks
V_thr <- V
V_thr[abs(V) < 1e-10] <- 0
Matrix::image(V_thr[1:20, 1:20], main = "Estimated V (with artifacts removed)")

#-------
#
#-------
decomp <- wlm::decompose_varcov(V, method = "evd", output = "all")
W <- decomp$transform

gassoc_lm <- matlm::matlm(trait1 ~ AGE, dat40, pred = gdat40, ids = ids,
  batch_size = 100, verbose = 2)
gassoc_wlm <- matlm::matlm(trait1 ~ AGE, dat40, pred = gdat40, ids = ids, transform = W, 
  batch_size = 100, verbose = 2)

passoc_lm <- matlm::matlm(trait1 ~ AGE, dat40, pred = pdat40, ids = ids,
  batch_size = 100, verbose = 2)
passoc_wlm <- matlm::matlm(trait1 ~ AGE, dat40, pred = pdat40, ids = ids, transform = W, 
  batch_size = 100, verbose = 2)
  
#-------
#
#-------
qq_plot(gassoc_lm$tab$pval) + ggtitle("LS: (null) random binary genotypes (not linked to kin2)")
qq_plot(gassoc_wlm$tab$pval) + ggtitle("LMM/WLS: (null) random binary genotypes (not linked to kin2)")

qq_plot(passoc_lm$tab$pval) + ggtitle("LS: (null) random cont. predictors (linked to kin2)")
qq_plot(passoc_wlm$tab$pval) + ggtitle("LMM/WLS: (null) random cont. predictors (linked to kin2)")
