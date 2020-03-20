data(dat40, package = "lme4qtl")
K <- as.matrix(kin2[1:5, 1:5])

# Covariance matrix K and its decompoisiton K = L L',
#   where L is lower triangular  matrix
# Also, K = R'R, where R is upper triangular matrix, i.e. R = L'  
#
# See Section "Implementation of lme4qtl" in
# https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-018-2057-x
#
# In lme4/lme4qtl many matrix producs use crossprod, 
# so we compute R rather than L

## 1. Cholesky decomposition, K = R'R
# chol function returns R
# see also lme4qtl:::relfac.chol
R <- Matrix::chol(K)
stopifnot(all(crossprod(R) == K))

## 2. Eigen value decomposition (EVD), K = V D V'
# chol function returns R
# see also lme4qtl:::relfac.evd
evd <- eigen(K, symmetric = TRUE)

D <- diag(evd$values)
V <- evd$vectors
R <- tcrossprod(sqrt(D), V) # D^{0.5} V'
try(stopifnot(all(crossprod(R) == K)))
stopifnot(all(round(crossprod(R), 10) == K))

## 3. Example of low-rank cov. matrix Klr
# see also lme4qtl:::relfac.evd
Klr <- K
Klr[3:5, 3:5] <- 1

# dimension vs. rank (5 vs. 3)
dim(Klr)
qr(Klr)$rank

# try Cholesky on Klr
try(chol(Klr))

# EVD on Klr + filter out zero eigen-values
evd <- eigen(Klr, symmetric = TRUE)

# filtering
ind_zero <- (abs(evd$values) < 1e-10) # near-zero eigen values can be negative
ev <- evd$values
ev[ind_zero] <- 0

# decomposition
D <- diag(ev)
V <- evd$vectors

R <- tcrossprod(sqrt(D), V) # D^{0.5} V'

# check
try(stopifnot(all(crossprod(R) == Klr)))
stopifnot(all(round(crossprod(R), 10) == Klr))
