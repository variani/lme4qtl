context("Function relmatGlmer")

test_that("no relmat structure", {
  data(cbpp, package = "lme4")
  
  f <- cbind(incidence, size - incidence) ~ period + (1 | herd)
  
  m1 <- glmer(f, cbpp, family = binomial)
  m2 <- relmatGlmer(f, cbpp, family = binomial)
  
  expect_equal(logLikNum(m1), logLikNum(m2))
})

test_that("nAGQ > 1", {
  stopifnot(require(HSAUR2))

  data(toenail, package = "HSAUR2")
  
  f <- outcome ~ treatment*visit + (1|patientID)
  #m1 <- glmer(f, toenail, family = binomial, nAGQ = 5)
  #m2 <- relmatGlmer(f, toenail, family = binomial, nAGQ = 5)
  
  #expect_equal(logLikNum(m1), logLikNum(m2))
})

test_that("compare with `pedigreemm` pkg", {
  stopifnot(require(pedigreemm))
  
  # data
  data(milk, package = "pedigreemm")
  milk <- within(milk, {
    sdMilk <- milk / sd(milk)
  })
  milk <- within(milk, {
    highMilk <- as.numeric(sdMilk > 6)
  })

  A <- getA(ped = pedCowsR)

  # subset
  #dat <- subset(milk, sire %in% c("330", "331", "332", "333", "334", "335") & lact == 3)
  dat <- subset(milk, lact == 2)
  
  # models
  f <- highMilk ~ log(dim) + (1|id)

  #m1 <- pedigreemm(f, dat, family = binomial, pedigree = list(id = pedCowsR))
  #m2 <- relmatGlmer(f, dat, family = binomial, relmat = list(id = A)) 
})
