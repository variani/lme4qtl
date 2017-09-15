### inc
library(Matrix)
library(lme4)
library(lme4qtl)

### data
load("Sample_file_for_testing.RData") # -> `FH2`, `maatriks2`

### check
isSymmetric(maatriks2, tol = 0)
#[1] FALSE

isSymmetric(maatriks2, tol = 0.01)
#[1] FALSE

isSymmetric(maatriks2, tol = 0.02)
#[1] TRUE

### model
#mudel2 <- relmatLmer(CTG_LDL_statin ~ Age+ Sugu + (1|Vcode1), FH2, relmat = list(Vcode1 = maatriks2))

### solution
mat <- forceSymmetric(maatriks2)
# caution: See `?forceSymmetric`, e.g. "result can be pretty surprising".

mod <- relmatLmer(CTG_LDL_statin ~ Age+ Sugu + (1|Vcode1), FH2, relmat = list(Vcode1 = mat), verbose = 2)

VarCorr(mod)
# Groups   Name        Std.Dev.
# Vcode1   (Intercept) 0.00000
# Residual             0.96904
