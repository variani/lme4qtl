# code example for https://github.com/variani/lme4qtl/issues/16
library(lme4)
library(lme4qtl)

data(dat40)

# fitting a model: y = mu + u1 + u2 + e
mod <- relmatLmer(trait1 ~ (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))

# prediction, see ?predict.merMod
# (1) BLUPs
u12_hat <- predict(mod, random.only = TRUE) # length = #individuals
blup_12 <- ranef(mod) # length (FAM) = # families; length(ID) = #individuals

# check if blup_12[["ID"]] = u12_hat, as blup_12[["FAMID"]] are all zeros
u12_hat %>% unique %>% sort %>% head
blup_12[["ID"]][, 1] %>% sort %>% head

# (2) predict with/without random effects
y_hat <- predict(mod, re.form = NULL)
y_hat_fixed <- predict(mod, re.form = NA) # u12_hat/blup_12 are not used

# y_hat_fixed = mu_hat = contstant
plot(y_hat, y_hat_fixed)

