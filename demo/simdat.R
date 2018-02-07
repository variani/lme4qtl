### par
num_ind <- 500

size_fam <- 5
num_fam <- num_ind / size_fam

### kinship matrix for sib-pairs
kin_fam <- matrix(0.5, size_fam, size_fam) 
diag(kin_fam)	<- 1  

kin <- kronecker(diag(num_fam), kin_fam)
	
### simulate data: y ~ 0.1 pred1 + 0.1 pred2 + sqrt(0.3) rand_fam + sqrt(0.3) rand_kin + sqrt(0.4) rand_res
set.seed(1)
pred1 <- rnorm(num_ind)
pred2 <- rnorm(num_ind)

rand_fam <- rep(rnorm(num_fam), each = size_fam)

ch <- chol(kin)
Z <- matrix(rnorm(num_ind), ncol = 1)
rand_kin <- as.numeric(crossprod(ch, Z))

rand_res <- rnorm(num_ind)

rand_total <- sqrt(0.3) * rand_fam + sqrt(0.3) * rand_kin + sqrt(0.4) * rand_res

y <- 0.1 * pred1 + 0.1 * pred2 + rand_total

### fit the model
simdat <- data.frame(id = seq(1, num_ind), fam = rep(seq(1, num_fam), each = size_fam),
  y = y, pred1 = pred1, pred2 = pred2)

rownames(kin) <- seq(1, num_ind)
colnames(kin) <- seq(1, num_ind)

mod <- relmatLmer(y ~ pred1 + pred2 + (1|fam) + (1|id), simdat, relmat = list(id = kin))

VarProp(mod)

