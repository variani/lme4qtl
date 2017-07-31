context("logLik")

test_that("start in lmer (modular)", {
  library(lme4)
  
  data(sleepstudy)
  f <- Reaction ~ Days + (1 | Subject)
  
  lmod <- lFormula(f, sleepstudy)
  devfun <- do.call(mkLmerDevfun, lmod)
  
  opt <- optimizeLmer(devfun)
  mod1 <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)

  start <- list(theta = getME(mod1, "theta"))
  opt2 <- list(par = as.numeric(start$theta), fval = NA,conv = 1000, message="start copied")
  mod2 <- mkMerMod(environment(devfun), opt2, lmod$reTrms, fr = lmod$fr)

  expect_true(getME(mod2, "theta") == start)
  expect_true(is.na(as.numeric(logLik(mod2))))
  expect_equal(as.numeric(logLik(mod1)), logLikNum(mod1))
  expect_equal(as.numeric(logLik(mod1)), logLikNum(mod2))
})

