context("Random slope")

test_that("cont. slope; one ranef", {
  stopifnot(require(Matrix))
    
  data(dat40)
  dat40 <- subset(dat40, FAMID %in% 1:3) # a smaller subset
  
  # model: random slope on FAMID; fitted by lme4  
  m1 <- lmer(trait1 ~ AGE + (0 + AGE|FAMID), dat40)

  # model: random slope on FAMID; fitted by lme4qtk using ID grouping and rel. matrix
  # two models are equvalent (m1 and m2)  
  mat_famid <- model.matrix(~ -1 + FAMID, dat40) %>% Matrix %>% tcrossprod
  rownames(mat_famid) <- dat40$ID
  colnames(mat_famid) <- dat40$ID

  m2 <- relmatLmer(trait1 ~ AGE + (0 + AGE|ID), dat40, relmat = list(ID = mat_famid))
  
  # get BLUP estimates 
  u1_hat <- round(relmatRanef(m1, "FAMID"), 5)
  u2_hat <- round(relmatRanef(m2, "ID"), 5)
  u2_hat_fam <- unique(u2_hat)

  expect_equal(u1_hat, u2_hat_fam)  
})

test_that("cont. slope; 2 ranef + 2 formula terms", {
  stopifnot(require(Matrix))
    
  data(dat40)
  dat40 <- subset(dat40, FAMID %in% 1:5) # a smaller subset
  dat40 <- within(dat40, ID2 <- ID)
  
  # model: random slope on FAMID; fitted by lme4  
  m1 <- lmer(trait1 ~ AGE + (1 + AGE||FAMID), dat40)

  # model: random slope on FAMID; fitted by lme4qtk using ID grouping and rel. matrix
  # two models are equvalent (m1 and m2)  
  mat_famid <- model.matrix(~ -1 + FAMID, dat40) %>% Matrix %>% tcrossprod
  rownames(mat_famid) <- dat40$ID
  colnames(mat_famid) <- dat40$ID

  m2 <- relmatLmer(trait1 ~ AGE + (1|ID) + (0 + AGE|ID2), dat40, relmat = list(ID = mat_famid, ID2 = mat_famid))
  
  # get BLUP estimates 
  u1_hat <- ranef(m1, which = "FAMID")[[1]][["AGE"]] %>% round(5)
  u2_hat <- relmatRanef(m2, "ID2") %>% round(5)
  u2_hat_fam <- unique(u2_hat)

  expect_equal(u1_hat, u2_hat_fam)  
})

test_that("cont. slope; 2 ranef + 1 formula term", {
  stopifnot(require(Matrix))
    
  data(dat40)
  dat40 <- subset(dat40, FAMID %in% 1:5) # a smaller subset
  dat40 <- within(dat40, ID2 <- ID)
  
  # model: random slope on FAMID; fitted by lme4  
  m1 <- lmer(trait1 ~ AGE + (1 + AGE||FAMID), dat40)

  # model: random slope on FAMID; fitted by lme4qtk using ID grouping and rel. matrix
  # two models are equvalent (m1 and m2)  
  mat_famid <- model.matrix(~ -1 + FAMID, dat40) %>% Matrix %>% tcrossprod
  rownames(mat_famid) <- dat40$ID
  colnames(mat_famid) <- dat40$ID

  m2 <- relmatLmer(trait1 ~ AGE + (1|ID) + (0 + AGE|ID2), dat40, relmat = list(ID = mat_famid, ID2 = mat_famid))
  
  m3 <- relmatLmer(trait1 ~ AGE + (1 + AGE||ID), dat40, relmat = list(ID = mat_famid))
  m4 <- relmatLmer(trait1 ~ AGE + (1|ID) + (0 + AGE|ID2), dat40, relmat = list(ID = mat_famid))

  # get BLUP estimates 
  u1_hat <- ranef(m1, which = "FAMID")[[1]][["AGE"]] %>% round(5)
  u2_hat <- relmatRanef(m2, "ID2") %>% round(5)
  u2_hat_fam <- unique(u2_hat)

  u3_hat <- ranef(m3, which = "ID")[[1]][["AGE"]] %>% round(5)
  u3_hat_fam <- unique(u3_hat)
  u4_hat <- relmatRanef(m4, "ID2") %>% round(5)
  u4_hat_fam <- unique(u4_hat)
  
  expect_equal(u3_hat, u4_hat_fam)  
})

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
