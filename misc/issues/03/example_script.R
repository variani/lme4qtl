dataX <- read.table("lme4qtl_data.csv", sep=",", dec=".", header=TRUE, skipNul = TRUE)
library(ape)
tree <- read.nexus("tree_all.tre")
phylo_all <- vcv.phylo(tree, cor = T) 

library(coxme)
m_lmekin <- lmekin(variable ~ mass*habitat + (1|species)+(1|site), data=dataX,
                      varlist=phylo_all, method = "REML")
print(m_lmekin)

# lmekin produces a low SE for mass:
#
#                   Value  Std Error     z       p
# mass           0.1682963 0.05734047  2.94 3.3e-03
# mass:habitatU -0.1192794 0.02719739 -4.39 1.2e-05

library(lme4qtl)
m_lme4qtl <- relmatLmer(variable ~ mass*habitat + (1|species)+(1|site), data=dataX, 
                     relmat = list(species = phylo_all), REML = T)
summary(m_lme4qtl)

# lme4qtl produces a much higher estimate and SE for mass:

#                Estimate Std. Error t value
# mass           0.38195    0.15330   2.491
# mass:habitatU -0.12362    0.02711  -4.561

library(visreg)
# CI bands for mass are very wide even though the model indicates a signficant interaction:
visreg(m_lme4qtl, "mass", "habitat", overlay=F, partial=T, band=T, xlab="centered and scaled log10(mass)")

# testing without species phylogeny
m_nophylo <- relmatLmer(variable ~ mass*habitat + (1|species)+(1|site), data=dataX, REML = T)
summary(m_nophylo)
visreg(m_nophylo, "mass", "habitat", overlay=F, partial=T, band=T, xlab="centered and scaled log10(mass)")

# without phylogeny the estimate and SE are similar to lmekin
#
#               Estimate Std. Error t value
# mass           0.11821    0.06049   1.954
# mass:habitatU -0.12054    0.02711  -4.446

# Why do lmekin and relmatLmer produce such different outputs for the parameters of body mass?
 