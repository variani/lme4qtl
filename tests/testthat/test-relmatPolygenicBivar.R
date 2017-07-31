context("Function lmerPolygenicBivar")

test_that("basic polygenic bivar. model", {
  #data(dat40)
  #relmat <- list(ID = kin2)
  #
  #f1 <- trait1 ~ AGE + (1 | ID)
  #f2 <- trait2 ~ SEX + (1 | ID)
  
  #mod1 <- relmatLmer(f1, dat40, relmat = relmat)
  #mod2 <- relmatLmer(f2, dat40, relmat = relmat)

  #poly <- lmerPolygenicBivar(c(f1, f2), dat40, relmat = relmat)
  
  #vcf <- as.data.frame(VarCorr(mod))[, c("grp", "vcov")]
  #vc <- vcf$vcov / sum(vcf$vcov)
  #h2 <- vc[1]

  #expect_true(h2 > 0.8)
})
