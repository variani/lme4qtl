### data
data(dat40)

### models with two RE
m1 <- relmatLmer(trait1 ~ (1|ID) + (1|FAMID), dat40, relmat = list(ID = kin2))
# boundary (singular) fit: see ?isSingular
# (1|FAMID) ~ zero variance

m2 <- relmatLmer(trait1 ~ (1|ID) + (1|FAMID), dat40, relmat = list(MyID = kin2))
# Model failed to converge: degenerate  Hessian with 1 negative eigenvalues

### models with one RE
m3 <- relmatLmer(trait1 ~ (1|ID), dat40, relmat = list(MyID = kin2))
# no errors 
# long computation time
