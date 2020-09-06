context("input data")

data(dat40)
dat <- head(dat40, 70)
dat <- within(dat, {
  ID1 <- ID
  ID2 <- ID
})


Z1 <- model.matrix(~ SEX - 1, dat) 
K1 <- tcrossprod(Z1)
rownames(K1) <- colnames(K1) <- dat$ID1

Z2 <- model.matrix(~ FAMID - 1, dat) 
K2 <- tcrossprod(Z2)
rownames(K2) <- colnames(K2) <- dat$ID2

test_that("two custom covariance matrices", {
  m <- relmatLmer(trait1 ~ (1|ID1) + (1|ID2), dat, relmat = list(ID1 = K1, ID2 = K2))

  expect_equal(length(ranef(m)), 2) # two random effects ID1 & ID2
})

test_that("error when some IDs in K are missing", {
  K1 <- K1[1:2, 1:2]

  expect_error({
    m <- relmatLmer(trait1 ~ (1|ID1), dat, relmat = list(ID1 = K1))
  })
})

test_that("error when rownames(K) is missing", {
  rownames(K1) <- NULL

  expect_error({
    m <- relmatLmer(trait1 ~ (1|ID1), dat, relmat = list(ID1 = K1))
  }, "rownames")
})

test_that("error when colnames(K) is missing", {
  colnames(K1) <- NULL

  expect_error({
    m <- relmatLmer(trait1 ~ (1|ID1), dat, relmat = list(ID1 = K1))
  }, "subscript out of bounds")
})

test_that("extra IDs in K", {
  m <- relmatLmer(trait1 ~ (1|ID1), head(dat, 10), relmat = list(ID1 = K1))

  expect_equal(length(ranef(m)), 1) # one random effect ID1
})

test_that("ID variable is a factor with #levels > #observations", {
  dat <- within(dat, ID1 <- factor(ID1)) 

  dat <- head(dat, 10)
  K1 <- K1[1:10, 1:10]
  
  m <- relmatLmer(trait1 ~ (1|ID1), dat, relmat = list(ID1 = K1))

  expect_true(nlevels(dat$ID1) > 10)
  expect_equal(length(ranef(m)), 1) # one random effect ID1
})
