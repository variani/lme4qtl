### inc
library(ape)
library(lme4qtl)
library(coxme)

### data
dataX <- read.table("lme4qtl_data.csv", sep=",", dec=".", header=TRUE, skipNul = TRUE)

dataX <- within(dataX,
  species_chr <- as.character(species))
  
tree <- read.nexus("tree_all.tre")
phylo_all <- vcv.phylo(tree, cor = T) 

ids_phylo <- rownames(phylo_all)
ids_dat <- as.character(dataX$species)

### `phylo_all_ordered`: make the job which `lmekin` is supposed to do
# - order rownames in `phylo_all` matrix according the the order in `species` column in data frame `dataX`
species_vals <- levels(dataX$species)

phylo_all_ordered <- phylo_all[species_vals, species_vals]

### models
m01 <- lmekin(variable ~ mass*habitat + (1|species), data=dataX, varlist=phylo_all, method = "REML")
m02 <- lmekin(variable ~ mass*habitat + (1|species), data=dataX, varlist=phylo_all_ordered, method = "REML")

m1 <- relmatLmer(variable ~ mass*habitat + (1|species), data=dataX, relmat = list(species = phylo_all), REML = T)
m2 <- relmatLmer(variable ~ mass*habitat + (1|species_chr), data=dataX, relmat = list(species_chr = phylo_all), REML = T)

### compare outputs of the models
coef(m01)[["fixed"]]
#  (Intercept)          mass      habitatU mass:habitatU 
#    1.0926508     0.1483706    -0.2663249    -0.1065064 
coef(m02)[["fixed"]]
#  (Intercept)          mass      habitatU mass:habitatU 
#    0.9814955     0.3642342    -0.2655203    -0.1103906 

fixef(m1)
#  (Intercept)          mass      habitatU mass:habitatU 
#    0.9814953     0.3642359    -0.2655202    -0.1103905 
fixef(m2)
#  (Intercept)          mass      habitatU mass:habitatU 
#    0.9814953     0.3642359    -0.2655202    -0.1103905

### conclusion
# - lme4qtl works properly whenever the order of IDs in covariance matrix is 
# - lmekin (coxme package) can produce a mess 
