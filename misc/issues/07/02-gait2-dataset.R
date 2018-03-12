###
# Reproducing the results reported in the GAIT2 vignette
# - it requires access to the GAIT2 data (not publicaly available)

### inc
library(plyr)
library(dplyr)

library(gait)
library(solarius)

### GAIT2 data
dir_phen <- "~/Data/GAIT2/phen/"
dir_snp <- "~/Data/GAIT2/ncdf/"

phen <- gait2.phen(dir_phen, transforms = "tr1", id.alert = TRUE, traits = "tr1_APTT")

phen <- rename(phen,
  aptt = tr1_APTT,
  gender = SEXf, age = AGEc,
  id = ID, famid = FAMID, hhid = HHID)
	
phen <- mutate(phen, rid = id, id7 = id)

dkin <- Matrix(solarKinship2(phen))

### fit 
m1 <- relmatLmer(aptt ~ age + gender + (1|id) + (1|hhid), phen, relmat = list(id = dkin))

### profile deviance
prof <- profile(m1, which = "theta_", prof.scale = "varcov") 
prof_prop <- varpropProf(prof) 

ci <- confint(prof_prop, level = 0.95)

### results
#> ci
#               2.5 %     97.5 %
#.sigprop01 0.4450731 0.84293219
#.sigprop02 0.0000000 0.06472766
#.sigmaprop 0.1461354 0.50813557

