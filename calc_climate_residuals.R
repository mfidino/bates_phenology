#########################
# Calculate climate residuals for bates et al.
# Written by M. Fidino 2017
#

# this function loads packages and downloads them if you do not have them

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
packs <- c("dplyr", "reshape2", "MCMCpack", "mcmcplots",
           "runjags", 'parallel')

# load the packages
package_load(packs)

# read in the co2 data
co <- read.csv("./data/bates_2017_co2.csv", header = TRUE)

# remove any data before 1868
co <- co[co$yr>1868,]

# scale year so that 1869 = zero
co$yr <- co$yr - min(co$yr)

# interpolate values that we are missing
mlow <- approx(co$yr, co$co2, xout = c(0:144))

# get 1872 to 2015 data
just_for_birds <- tail(mlow$y, length(1872:2015))

# exploratory plot
plot(just_for_birds~ c(1872:2015), type = "l", xlim = c(1872,2020),
     ylim = c(275, 400),
     lwd = 3, bty = "n", xlab = "Year", ylab = "co2 concentration (ppmv)",
     cex.lab = 1.5, xaxt = "n", yaxt = "n", main = "co2 levels over time",
     col = "red")


axis(2, at = seq(275, 400, 25), labels=F, tck=-.035)
axis(2, at = seq(275, 400, 12.5), labels=F, tck=-.025)
axis(2, at = seq(275, 400, 12.5/2), labels=F, tck=-.015)
mtext(text = seq(275, 400, 25), 2, line = .75, at = seq(275, 400, 25), las = 1, cex = 1)

axis(1, at = seq(1870,2020, 25), labels=F, tck=-.035)
axis(1, at = seq(1870,2020, 12.5), labels=F, tck=-.025)
axis(1, at = seq(1870,2020, 6.25), labels=F, tck=-.015)
mtext(text = seq(1870, 2020, 25), 1, 
      line = 0.75, at = seq(1870, 2020, 25), 
      las = 1, cex = 1)

# prepare data for analysis in jags
dl <- list(N = length(just_for_birds), Y = just_for_birds, x = 0:143)

# parameters to track
params <- c( "Yp", "Y")

# set up arguments for JAGS
n_chains = detectCores() - 1
adapt_steps = 5000
burn_in = 20000
sample_steps = 10000
thin_steps = 2
mod_mcmc <- as.mcmc.list(run.jags( model= "./jags_models/bates_2017_climate_resid_model.R" , 
                                           monitor=params , 
                                           data=dl ,  
                                           n.chains=n_chains ,
                                           adapt=adapt_steps ,
                                           burnin=burn_in , 
                                           sample=ceiling(sample_steps / n_chains) ,
                                           thin=thin_steps ,
                                           summarise=FALSE ,
                                           plots=FALSE,
                                           method = "parallel"))

# convert to a matrix and get median estimate
mmat <- as.matrix(mod_mcmc, chains = TRUE)
mmat <- apply(mmat, 2, median)

# get only the predictions
y2 <- mmat[grep("Yp\\[", names(mmat))]
# calculate residual
ress <- just_for_birds - y2
# save the residuals for analysis
write.csv(ress, "./data/climate_residuals.csv", row.names = FALSE)

# make a plot of the residuals
pdf("bates_2017_climate_resid_plot.pdf", height = 5, width = 5)

plot(ress~ c(1872:2015), type = "l", xlim = c(1872,2020),
     ylim = c(-30, 40),
     lwd = 3, bty = "n", xlab = "Year", ylab = "Difference from linear trend",
     cex.lab = 1.5, xaxt = "n", yaxt = "n", main = "Non-linear component\n of rising co2 levels",
     col = "blue")


axis(2, at = seq(-30, 40, 10), labels=F, tck=-.035)
axis(2, at = seq(-30, 40, 5), labels=F, tck=-.025)
axis(2, at = seq(-30, 40, 2.5), labels=F, tck=-.015)
mtext(text = seq(-30, 40, 10), 2, line = .75, at = seq(-30, 40, 10), las = 1, cex = 1)

axis(1, at = seq(1870,2020, 25), labels=F, tck=-.035)
axis(1, at = seq(1870,2020, 12.5), labels=F, tck=-.025)
axis(1, at = seq(1870,2020, 6.25), labels=F, tck=-.015)
mtext(text = seq(1870, 2020, 25), 1, 
      line = 0.75, at = seq(1870, 2020, 25), 
      las = 1, cex = 1)


dev.off()

