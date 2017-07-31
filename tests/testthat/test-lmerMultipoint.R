context("lmerMultipoint")

test_that("linkage on dat30", {
  #stopifnot(require(kinship2))
  #stopifnot(require(solarius))  
  
  #data(dat30, package = "solarius")
  #mibddir <- system.file('extdata', 'solarOutput',
  #  'solarMibdsCsv', package = 'solarius')
  
  #kin2 <- 2 * with(dat30, kinship(id, fa, mo, sex))
  
  #dat30 <- mutate(dat30,
  #  ID = id)
  #relmat <- list(ID = kin2)
  
  #mibd <- solarius:::solarMIBD(mibddir, chr = 5, nmibd = 1)
  
  ### model
  #link <- lmerMultipoint(trait1 ~ 1 + (1|ID), dat30, relmat = relmat, mibd = mibd, verbose = 0) 
    
  #expect_true(link$lodf$LOD > 1.4)
})
