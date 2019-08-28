# [R-sig-ME] extracting gradient and hessian matrix from lme4: https://stat.ethz.ch/pipermail/r-sig-mixed-models/2012q2/018094.html
#   - a similar example: https://stackoverflow.com/a/31704646
# Use case of profile: https://stackoverflow.com/a/37243460
# Delta method in R package msm: https://rdrr.io/rforge/msm/man/deltamethod.html
# Use case of msm::deltamethod: https://github.com/sashagusev/SKK-REML-sim/blob/master/func_reml.R#L197
# [R-sig-ME] Error in Profile likelihood based confidence intervals in glmer(): https://stat.ethz.ch/pipermail/r-sig-mixed-models/2014q3/022394.html

# data
data(dat40, package = "lme4qtl")

# down-sample to check of se is getting larger
#dat40 <- head(dat40, 100)

# fit a model
m <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))

# try to profile and get CI for theta...
prof <- profile(m1, which = "theta_", signames = FALSE)
head(prof) # many NAs
try(confint(prof))

# point estimate of variance components
vf <- VarProp(m)
vc <- vf[["prop"]]
vc <- head(vc, -1) # the last one is `sigma(m)`

# hessian
f <- update(m, devFunOnly = TRUE)

library(numDeriv)
th <- getME(m, "theta")
h <- hessian(f, th)

library(MASS)
cov <- ginv(h)

# delta method
nth <- length(th)
terms <- paste0("x", seq(nth), "*", "x", seq(nth)) # "x1*x1" "x2*x2"
terms <- c(terms, "1") # # "x1*x1" "x2*x2", "1"

total <- paste(terms, collapse = " + ") # "x1*x1 + x2*x2 + 1"

forms <- lapply(seq(nth), function(i) 
 paste("~", terms[i], "/ (", total, ")"))
  
ses <- sapply(forms, function(f)
  deltamethod(as.formula(f), th, cov, ses = TRUE))

# show results
tab <- data.frame(grp = head(vf$grp, -1), est = vc, se = ses)
print(tab)

