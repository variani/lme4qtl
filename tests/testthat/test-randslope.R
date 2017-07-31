context("Random slope")

test_that("GxE in dat30", {
  stopifnot(require(kinship2))
  
  data(dat30, package = "solarius")
  dat30 <- mutate(dat30,
    sex = as.factor(sex))

  kin2 <- 2 * with(dat30, kinship(id, fa, mo, sex))

  mod1 <- relmatLmer(trait1 ~ 1 + (1|id), dat30, relmat = list(id = kin2))

  mod2 <- relmatLmer(trait1 ~ 1 + (0+sex|id), dat30, relmat = list(id = kin2))
  
  # h2
  vf1 <- as.data.frame(VarCorr(mod1))
  vc <- vf1$vcov / sum(vf1$vcov)
  h2 <- vc[1]

  # rho
  vf2 <- as.data.frame(VarCorr(mod2))
  rho <- vf2$sdcor[3]
  
  expect_equal(round(h2, 2), 0.84)  
  expect_equal(round(rho, 2), 0.75)
})
