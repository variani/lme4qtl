context("Function relmatLmer")

test_that("h2 of trait1 in dat30", {
  stopifnot(require(kinship2))
  
  if(require("solarius")) {
    data(dat30, package = "solarius")
    kin2 <- 2 * with(dat30, kinship(id, fa, mo, sex))

    mod <- relmatLmer(trait1 ~ 1 + (1|id), dat30, relmat = list(id = kin2), REML = FALSE)

    vcf <- as.data.frame(VarCorr(mod))[, c("grp", "vcov")]
    vc <- vcf$vcov / sum(vcf$vcov)
    h2 <- vc[1]

    expect_true(round(h2, 4) == 0.8343)
  }
})
