context("update")

test_that("start + optimizer.none", {
  data(sleepstudy)
  
  f1 <- Reaction ~ 1 + (Days | Subject)
  f2 <- Reaction ~ Days + (Days | Subject)
  
  t1 <- system.time(mod1 <- relmatLmer(f1, sleepstudy))
  
  start <- list(theta = getME(mod1, "theta"))
  t2 <- system.time(mod2 <- relmatLmer(f2, sleepstudy, control = lmerControl(optimizer = "none"), start = start))
  
  expect_equivalent(getME(mod1, "theta"), getME(mod2, "theta"))
  expect_true(is.na(mod2@optinfo$feval))
  #expect_lt(t2[3], t1[3])
})

test_that("update + start + optimizer.none", {
  data(sleepstudy)
  
  f1 <- Reaction ~ 1 + (Days | Subject)
  f2 <- Reaction ~ Days + (Days | Subject)
  
  t1 <- system.time(mod1 <- relmatLmer(f1, sleepstudy))
  
  start <- list(theta = getME(mod1, "theta"))
  t2 <- system.time(mod2 <- update(mod1, f2, control = lmerControl(optimizer = "none"), start = start))
  
  expect_equivalent(getME(mod1, "theta"), getME(mod2, "theta"))
  expect_true(is.na(mod2@optinfo$feval))
  #expect_lt(t2[3], t1[3])
})
