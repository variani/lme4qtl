context("Polygenic random slope")

test_that("polygenic random slope #18", {
  library(lme4)
  library(lme4qtl)

  data(dat40)
  dat40 <- within(dat40, { 
    ID2 <- ID
    AGE <- AGE
  })

  # create a random kinship matrix or the FAMIDs
  fams <- unique(dat40$FAMID)
  set.seed(1)
  kin_fam <- tcrossprod(matrix(rnorm(length(fams)^2), length(fams))) + diag(1, length(fams))
  rownames(kin_fam) <- colnames(kin_fam) <- fams
  
  polyslope_id <- relmatLmer(trait1 ~ AGE + (1 + AGE|ID), dat40, relmat = list(ID = kin2))
  polyslope_famid <- relmatLmer(trait1 ~ AGE + (1 + AGE|FAMID), dat40, relmat = list(FAMID = kin_fam)) 

  expect_true(class(polyslope_id) == "lmerMod")
  expect_true(class(polyslope_famid) == "lmerMod")
})
