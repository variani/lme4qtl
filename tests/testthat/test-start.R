context("start")

test_that("start in lmer (modular)", {
  data(sleepstudy, package = "lme4")
  
  f <- Reaction ~ Days + (1 | Subject)
  lmod <- lFormula(f, sleepstudy)
  devfun <- do.call(mkLmerDevfun, lmod)
  opt <- optimizeLmer(devfun)

  mod <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)
  theta <- as.numeric(getME(mod, "theta"))
    
  lmod2 <- lFormula(f, sleepstudy)
  devfun2 <- do.call(mkLmerDevfun, lmod2)
  opt2 <- list(par = 3, fval = NA,conv = 1000, message="start copied")
  devfun2(opt2$par <- 2)
  mod2 <- mkMerMod(environment(devfun2), opt2, lmod2$reTrms, fr = lmod2$fr)

  lmod3 <- lFormula(f, sleepstudy)
  devfun3 <- do.call(mkLmerDevfun, lmod3)
  opt3 <- list(par = 3, fval = NA,conv = 1000, message="start copied")
  devfun3(opt3$par <- 3)
  mod3 <- mkMerMod(environment(devfun3), opt3, lmod3$reTrms, fr = lmod3$fr)

  #print(getME(mod, "theta"))  
  #print(getME(mod2, "theta"))
  #print(getME(mod3, "theta"))

  expect_true(getME(mod, "theta") == theta)
  expect_true(getME(mod2, "theta") == 2)
  expect_true(getME(mod3, "theta") == 3)  
})

test_that("start in lmer", {
  if(packageVersion("lme4") >= "1.1.14") {
    data(dat30, package = "solarius")

    # model 1
    mod1 <- lmer(trait1 ~ 1 + (1|famid), dat30)

    # model 2
    start <- list(theta = as.numeric(getME(mod1, "theta")))
    control <- lmerControl(optimizer = NULL)

    mod2 <- lmer(trait1 ~ 1 + (1|famid), dat30, control = control, start = start)

    expect_equal(getME(mod2, "theta"), getME(mod2, "theta"))
  }
})

test_that("start in lmer (modular) on dat30 with (1|famid)", {
  data(dat30, package = "solarius")
    
  f <- trait1 ~ 1 + (1|famid)
  lmod <- lFormula(f, dat30)
  devfun <- do.call(mkLmerDevfun, lmod)
  theta1 <- with(environment(devfun), pp$theta)
  
  opt <- optimizeLmer(devfun)
  theta2 <- with(environment(devfun), pp$theta)
  
  mod1 <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)

  start <- list(theta = getME(mod1, "theta"))
  opt2 <- list(par = as.numeric(start$theta), fval = NA,conv = 1000, message="start copied")
  mod2 <- mkMerMod(environment(devfun), opt2, lmod$reTrms, fr = lmod$fr)

  expect_true(getME(mod2, "theta") == start)
})

test_that("start in relmatLmer", {
  stopifnot(require(kinship2))
  
  data(dat30, package = "solarius")
  kin2 <- 2 * with(dat30, kinship(id, fa, mo, sex))

  # model 1
  mod1 <- relmatLmer(trait1 ~ 1 + (1|id), dat30, relmat = list(id = kin2))

  # model 2
  start <- list(theta = as.numeric(getME(mod1, "theta")))
  control <- lmerControl(optimizer = "none")

  mod2 <- relmatLmer(trait1 ~ 1 + (1|id), dat30, relmat = list(id = kin2), control = control, start = start)

  expect_true(getME(mod2, "theta") == start)
})

test_that("start in relmatLmer: REML = FALSE", {
  stopifnot(require(kinship2))
  
  data(dat30, package = "solarius")
  kin2 <- 2 * with(dat30, kinship(id, fa, mo, sex))

  # model 1
  mod1 <- relmatLmer(trait1 ~ 1 + (1|id), dat30, REML = FALSE, relmat = list(id = kin2))

  # model 2
  start <- list(theta = as.numeric(getME(mod1, "theta")))
  control <- lmerControl(optimizer = "none")

  mod2 <- relmatLmer(trait1 ~ 1 + (1|id), dat30, REML = FALSE, relmat = list(id = kin2), control = control, start = start)

  expect_true(getME(mod2, "theta") == start)
})

test_that("start in relmatLmer: 2 random effects", {
  stopifnot(require(kinship2))
  stopifnot(require(solarius))
  
  data(dat30, package = "solarius")
    
  data(dat30)
  kin2 <- 2 * with(dat30, kinship(id, fa, mo, sex))

  # model 1
  mod1 <- relmatLmer(trait1 ~ 1 + (1|famid) + (1|id), dat30, relmat = list(id = kin2))

  # model 2
  start <- list(theta = as.numeric(getME(mod1, "theta")))
  control <- lmerControl(optimizer = "none")

  mod2 <- relmatLmer(trait1 ~ 1 + (1|famid) + (1|id), dat30, relmat = list(id = kin2), control = control, start = start)

  expect_true(all(as.numeric(getME(mod2, "theta")) == as.numeric(start$theta)))
})

test_that("start in hacked lmer (modular) on dat30 with (1|id)", {
  stopifnot(require(kinship2))
  
  data(dat30, package = "solarius")
  kin2 <- 2 * with(dat30, kinship(id, fa, mo, sex))
  relmat <- list(id = kin2)
  
  # formula
  f <- trait1 ~ 1 + (1|id)
  control <- lmerControl(check.nobs.vs.rankZ = "ignore", 
    check.nobs.vs.nlev = "ignore", check.nobs.vs.nRE = "ignore")
  
  lmod <- lFormula(f, dat30, control = control)
  
  #-------------------------------
  # start of solaris-specific code
  #-------------------------------
  stopifnot(is.list(relmat), length(names(relmat)) == length(relmat))
  relnms <- names(relmat)
  relfac <- relmat
  flist <- lmod$reTrms[["flist"]]   ## list of factors

  ind <- (relnms %in% names(flist))
  if(any(ind)) {
    ## random-effects design matrix components
    Ztlist <- lmod$reTrms[["Ztlist"]]
    
    asgn <- attr(flist, "assign")
    for(i in seq_along(relnms[ind])) {
      
      relmati <- relnms[ind][i]
      if(!(relmati %in% names(flist))) {
        stop("a relationship matrix must be (", relmati, ")",
          " associated with only one random effects term (", paste(names(flist), collapse = ", "), ")")
      }
      tn <- which(relmati == names(flist))
      fn <- names(flist)[tn]
      
      ### option 1
      #zn <- rownames(Ztlist[[tn]]) 
      # here was the correction: index value is `tn`, rather than `i`
      ### option 2
      # > packageVersion("lme4")
      # [1] ‘1.1.8’
      zn <- lmod$fr[, fn]
      if(class(rownames(relmat[[i]])) == "character") {
        zn <- as.character(zn)
      }
      
      relmat[[i]] <- Matrix::Matrix(relmat[[i]][zn,zn], sparse = TRUE)
      stopifnot(nrow(relmat[[i]]) == nrow(lmod$fr))
    
      relfac[[i]] <- Matrix::chol(relmat[[i]])
      # ?Matrix::chol
      # Returned value: a matrix of class ‘Cholesky’, i.e., upper triangular: R such that R'R = x.
      # Note that another notation is equivalent x = L L', where L is a lower triangular 
      # @ http://en.wikipedia.org/wiki/Cholesky_decomposition
      Ztlist[[i]] <-  relfac[[i]] %*% Ztlist[[i]]
      # If the substitution is Z* = Z L, then Z*' = L' Z' = R Z'
      # @ http://www.journalofanimalscience.org/content/88/2/497.long%5BVazquez%20et%20al.,%202010%5D
    }
  
    lmod$reTrms[["Ztlist"]] <- Ztlist
    lmod$reTrms[["Zt"]] <- do.call(rBind, Ztlist)
  }
    
  # mod1
  devfun <- do.call(mkLmerDevfun, lmod)
  opt <- optimizeLmer(devfun)
  mod1 <- mkMerMod(environment(devfun), opt, lmod$reTrms, fr = lmod$fr)

  # mod2
  start <- list(theta = getME(mod1, "theta"))
  devfun <- do.call(mkLmerDevfun, c(lmod, start = start))
  
  opt2 <- list(par = as.numeric(start$theta), fval = NA,conv = 1000, message="start copied")
  mod2 <- mkMerMod(environment(devfun), opt2, lmod$reTrms, fr = lmod$fr)

  expect_true(getME(mod2, "theta") == start)
})
