context("modular implementation of relmat*`")

test_that("ids order by relmatLmer", {
  if(suppressWarnings(require(lmmlite))) {
    data(recla, package = "lmmlite")
    dat <- data.frame(y = recla$pheno[,1], int = recla$covar[, 1], sex = recla$covar[, 2],
      id = rownames(recla$kinship))
    
    mod <- relmatLmer(y ~ sex + (1|id), dat, relmat = list(id = recla$kinship))
    h2 <- with(as.data.frame(VarCorr(mod)), vcov[1] / sum(vcov))

    expect_true(h2 > 0.7)
  }
})

test_that("ids order by relmatGlmer", {
  if(suppressWarnings(require(lmmlite))) {
    data(recla, package = "lmmlite")
    
    dat <- data.frame(y = recla$pheno[,1], int = recla$covar[, 1], sex = recla$covar[, 2],
      id = rownames(recla$kinship))

    dat <- mutate(dat,
      ycnt = round(y))
    
    mod <- relmatGlmer(ycnt ~ (1|id), dat[1:70, ], relmat = list(id = recla$kinship), family = poisson)
    vcov <- with(as.data.frame(VarCorr(mod)), vcov)
    
    expect_true(vcov < 0.2)
  }
})
