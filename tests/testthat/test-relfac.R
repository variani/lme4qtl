context("relfac")

test_that("relfac methods", {
  K <- matrix(0.9, 2, 2)
  diag(K) <- 1
  K <- Matrix(K)
  
  K.chol <- crossprod(relfac.chol(K))
  K.svd <- crossprod(relfac.svd(K, 1e-10))
  K.evd <- crossprod(relfac.evd(K, 1e-10))
  
  expect_equal(as.vector(K), as.vector(K.chol), as.vector(K.svd), as.vector(K.evd))
})

test_that("relfac methods for rank deficient matrices", {
  stopifnot(require(Matrix))

  gr <- factor(rep(1:2, each = 2))
  Z <- model.matrix(~ gr - 1)
  K <- tcrossprod(Matrix(Z))
  
  rf <- relfac(K)
  R <- crossprod(rf)

  expect_equal(as.vector(K), as.vector(R))  
})

test_that("equivalence on (1|plate) and K", {
  stopifnot(require(Matrix))

  data(Penicillin)
  Penicillin <- within(Penicillin, {
    id <- as.character(seq(1, nrow(Penicillin)))
  })
  
  # lme4 model
  m1 <- lmer(diameter ~ (1|plate), Penicillin)

  # lme4qtl model
  Z <- with(Penicillin, model.matrix(~ plate - 1)) 
  Z <- Matrix(Z)

  K <- tcrossprod(Z)
  rownames(K) <- Penicillin$id
  colnames(K) <- Penicillin$id

  m2 <- relmatLmer(diameter ~ (1|id), Penicillin, relmat = list(id = K))

  expect_true(abs(getME(m1, "theta") - getME(m2, "theta")) < 1e-6) 
})

test_that("relfac methods for h2 estimation", {
  stopifnot(require(kinship2))
  
  data(dat30, package = "solarius")

  dat30 <- mutate(dat30,
    id = factor(id))
  
  kin2 <- 2 * with(dat30, kinship(id, fa, mo, sex))
  relmat <- list(id = kin2)
  
  mod1 <- relmatLmer(trait1 ~ 1 + (1|id), dat30, relmat = relmat, method.relfac = "chol")
  mod2 <- relmatLmer(trait1 ~ 1 + (1|id), dat30, relmat = relmat, method.relfac = "svd")

  vcf <- as.data.frame(VarCorr(mod1))[, c("grp", "vcov")]
  vc <- vcf$vcov / sum(vcf$vcov)
  h2.mod1 <- vc[1]

  vcf <- as.data.frame(VarCorr(mod2))[, c("grp", "vcov")]
  vc <- vcf$vcov / sum(vcf$vcov)
  h2.mod2 <- vc[1]

  expect_true(h2.mod1 > 0.8)
  expect_true(abs(h2.mod1 - h2.mod2) < 1e-5)
})
