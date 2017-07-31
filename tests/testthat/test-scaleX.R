context("scaleX")

test_that("check.scaleX", {
  if(FALSE) 
  {
  
  stopifnot(require(kinship2))
  
  data(dat30)
  kin2 <- 2 * with(dat30, kinship(id, fa, mo, sex))

  control1 <- lmerControl(check.scaleX  = "warning")
  control2 <- lmerControl(check.scaleX = "stop")
  control3 <- lmerControl(check.scaleX  = "ignore")

  dat30$x <- 0
  dat30$x[1] <- 0.001
  
  # models
  mod3 <- relmatLmer(trait1 ~ x + (1|id), dat30, relmat = list(id = kin2), control = control3)

  # tests  
  expect_warning({
    mod1 <- relmatLmer(trait1 ~ x + (1|id), dat30, relmat = list(id = kin2), control = control1)
  }, "Some predictor variables are on very different scales: consider rescaling")
  expect_error({
    mod <- relmatLmer(trait1 ~ x + (1|id), dat30, relmat = list(id = kin2), control = control2)
  })
  expect_true("lmerMod" %in% class(mod3))
  
  }
})
