############################################
# analysis of nesting data for Bates et al.
# Written by M. Fidino
############################################

package_load<-function(packages = NULL, quiet=TRUE, verbose=FALSE, warn.conflicts=FALSE){
  
  # download required packages if they're not already
  
  pkgsToDownload<- packages[!(packages  %in% installed.packages()[,"Package"])]
  if(length(pkgsToDownload)>0)
    install.packages(pkgsToDownload, repos="http://cran.us.r-project.org", quiet=quiet, verbose=verbose)
  
  # then load them
  for(i in 1:length(packages))
    require(packages[i], character.only=T, quietly=quiet, warn.conflicts=warn.conflicts)
}

# required packages
packs <- c("dplyr", "runjags", "MCMCpack", "mcmcplots",
           "runjags", 'parallel', "phytools", "ape")

# load the packages
package_load(packs)




# read in the data
# migt = species migratory status info
# fdt = species nesting data
migt <- read.csv("./data/bates_2017_migratory_status.csv", header = TRUE)
migt$species <- factor(migt$species)
fdt <- read.csv("./data/bates_2017_bird_lay_dates.csv", header = TRUE)
fdt$species <- factor(fdt$species)

# length of data
n <- nrow(fdt)

# get unique species
span <- unique(fdt$species)

# get number of species
nspec <- length(span)

# get number of coefs in group lvl
# intercept, time, climate residuals
ncof <- 3

# create X matrix to include in linear predictor
# adding residuals afterwards
# Setting up year so that the intercept represents
# nesting on the first year data was collected
X <- matrix(c(rep(1, n), fdt$year - (min(fdt$year) )), ncol = 2, nrow = n)

# read in residual data
ress <- read.csv("./data/climate_residuals.csv", header = TRUE)

# add 1 to index ress because X[,2] ranges from 0 to 143
ress_rep <- ress[X[,2]+1,]

# bind climate residuals to covariates
X <- cbind(X,ress_rep)

# make species a numberic vector
species <- as.numeric(factor(fdt$species))

# keep the species info though so we can relate
# the numeric vector to the actual species
spf <- levels(factor(fdt$species))

# make U matrix
# this is the group level matrix
# for the mvn distribution

# make sp code a factor
spcode <- factor(migt$species)
migt$species <- factor(migt$species)
migt <- migt[order(migt$species),]
# make a vector so we can use model.matrix
migt$y <- 1
# reorder the migratory status factor
migt$migstat <- factor(migt$migstat, levels = c("resident", "short", "long"))
# make the U matrix
U <- matrix(model.matrix(migt$y~migt$migstat), nrow = nspec, ncol = 3)

# make W matrix, which is the 
# initial value 
W <- diag(ncof)

# make julian date a numeric vector
y <- as.numeric(fdt$jdate)

# make a matrix and other data for predictions
yrs <- c(seq(0, 140, 5), 143)
nyr <- length(yrs)

x_pred <- matrix(c(rep(1, nspec*nyr), rep(yrs, nspec)), ncol = 2,
                 nrow = nspec*nyr)
pr <- ress[x_pred[,2]+1,]

# add predicted to x_pred
x_pred <- cbind(x_pred, pr)

at_start <- which(x_pred[,2]==0)
at_end <- which(x_pred[,2]==143)

# read in phylogeny
my_tree <- ape::read.nexus(
  "./data/bird_consensus_tree.nex"
)

# generate an ultrametric tree
tree_ultra=ape::chronos(my_tree, lambda=0) 
tree_unmatched <- multi2di(tree_ultra, random=TRUE)

# Get the variance covariance matrix, and invert it
ivcv <- solve(
  ape::vcv(
    tree_unmatched
  )
)

data_list <- list(
  y = y, X = X, U = U, W = W,
  ncof = ncof, species = species,
  nspec = nspec, n = n, nmig = 3, 
  ivcv = ivcv,
  species_pred = rep(1:nspec, each = nyr),
                  x_pred = x_pred, n_pred = nrow(x_pred),
                  mx = X[,2] - mean(X[,2])
)

# make intial values function
inits <- function(chain){
  gen_list <- function(chain = chain){
  list(
    B.raw=array(
      rnorm(
        nspec*ncof
      ),
      dim = c(nspec, ncof)
    ),
    Tau.B.raw=rwish(ncof+1, diag(ncof)),
    xi = runif(ncof),
    TauPhylo = runif(ncof, 0.01, 1),
       .RNG.name = switch(chain,
                          "1" = "base::Wichmann-Hill",
                          "2" = "base::Marsaglia-Multicarry",
                          "3" = "base::Super-Duper",
                          "4" = "base::Mersenne-Twister",
                          "5" = "base::Wichmann-Hill",
                          "6" = "base::Marsaglia-Multicarry",
                          "7" = "base::Super-Duper",
                          "8" = "base::Mersenne-Twister"),
       .RNG.seed = sample(1:1e+06, 1))
  }
  return(switch(chain,           
                "1" = gen_list(chain),
                "2" = gen_list(chain),
                "3" = gen_list(chain),
                "4" = gen_list(chain),
                "5" = gen_list(chain),
                "6" = gen_list(chain),
                "7" = gen_list(chain),
                "8" = gen_list(chain)
  )
  )
  }
# parameters to track
params <- c("B", "G", "sigma.y", "sigma.B", "rho.B", "y_pred", "nu",
            "t0_log", "t1_log", "TauPhylo", "B_pre")

# set up arguments for JAGS
n_chains = detectCores()-2
adapt_steps = 100000
burn_in = 100000
sample_steps = 10000
thin_steps = 20

mod_mcmc_more_sp <- as.mcmc.list(run.jags( model= "./jags_models/bates_2017_robust_t_phylo_model.R" , 
                                           monitor=params , 
                                           data=data_list ,  
                                           inits=inits , 
                                           n.chains=n_chains ,
                                           adapt=adapt_steps ,
                                           burnin=burn_in , 
                                           sample=ceiling(sample_steps / n_chains) ,
                                           thin=thin_steps ,
                                           summarise=FALSE ,
                                           plots=FALSE,
                                           method = "parallel"))

saveRDS(mod_mcmc_more_sp, "./mcmc_output/bates_2017_model_output_phylo.RDS")
mm <- as.matrix(mod_mcmc_more_sp, chains = TRUE)
# saving the output of the model
write.csv(mm, "C:/Users/mfidino/Documents/bates_2017_model_output_2.csv")

ph <- mm[,grep("Phylo", colnames(mm))]


# Trying it via brms

install.packages("brms")
library(brms)


?brm



#make the data.frame
my_df <- data.frame(
  y = data_list$y,
  co = data_list$X[,2],
  yr = data_list$X[,3],
  migstat = migt$migstat[data_list$species],
  species = factor(fdt$species),
  phylo = factor(fdt$species),
  mx = data_list$mx
)
phylo <- tree_unmatched
A <- ape::vcv.phylo(phylo)

# the formula
# y ~ yr + co + yr:migstat + co:migstat + (1 + yr + co | gr(phylo, cov = A)) + 
# (1 + yr + co | species)
#
# sigma ~ mx
#
# nu ~ offset(1)

longshot <- brm(
  bf(
    y ~ yr*migstat + co*migstat + 
    (1 + yr + co|gr(phylo, cov=A)) + (1 + yr + co |species),
    sigma ~ mx,
    nu ~ offset(1)
  ),
  data = my_df,
  family = student(),
  data2 = list(A = A),
  prior = c(
    prior(normal(0, 1000), "b"),
    prior(normal(0, 1000), "Intercept"),
    prior(exponential(1/29), "Intercept", dpar = "nu"),
    prior(uniform(-10,10), "Intercept", dpar = "sigma"),
    prior(uniform(-10,10), "b", dpar = "sigma")
  ),
  iter = 5000
)


Sys.getenv("PATH")
Sys.getenv("BINPREF")
readLines("~/.R/Makevars.win")
readLines("~/.Rprofile")
readLines("~/.Renviron")
devtools::session_info("rstan")

  hm <- get_prior(
    bf(
      y ~ yr*migstat + co*migstat + 
        (1 + yr + co|gr(phylo, cov=A)) + (1 + yr + co |species),
      sigma ~ mx,
      nu ~ offset(1)
    ),
    data = my_df,
    family = student(),
    data2 = list(A = A)
  )
