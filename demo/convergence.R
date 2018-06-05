# help("troubleshooting", package = "lme4")

### inc
library(lme4) # ?lme4::lmerControl for controls; `args(lmerControl)`
library(minqa) # `?minqa::bobyqa` for controls

### data
data(dat40, package = "lme4qtl")

### model
ctr1 <- lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 30))
mod1 <- relmatLmer(trait1 ~ (1|ID), dat40, relmat = list(ID = kin2), control = ctr1)

ctr2 <- ctr1
mod2 <- relmatLmer(trait1 ~ (1|ID), dat40, relmat = list(ID = kin2), control = ctr1, calc.derivs = FALSE)

ctr3 <- lmerControl(optimizer = "bobyqa", boundary.tol = 0.5, optCtrl = list(maxfun = 30, iprint = 5), check.conv.grad = list(action = "warning", tol = 0.5, relTol = NULL))
mod3 <- relmatLmer(trait1 ~ (1|ID), dat40, relmat = list(ID = kin2), control = ctr3, calc.derivs = FALSE)
