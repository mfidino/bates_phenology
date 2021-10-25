######################
#
# bates
#
#



mm <- data.table::fread("C:/Users/mfidino/Documents/bates_2017_model_output_2.csv", 
                        header = TRUE, data.table = FALSE)
# remove sigma.y
mm <- mm[,-grep("sigma.y", colnames(mm))]
m2 <- apply(mm, 2, median)

# remove y_pred from the covariate samples
covs <- mm[,-grep("y_pred|sigma.y", colnames(mm))]
# get just the predictions
preds <- mm[,grep("y_pred", colnames(mm))]

# calculate Bayesian p-value
bayes_p <-  preds > data_list$y 
mean(bayes_p) # this is the bayesian p-value

# parameters from changing standard deviation
t0_log <- mm[,grep("t0_log", colnames(mm))]
t1_log <- mm[,grep("t1_log", colnames(mm))]
aa <- 0:143
aa <- aa - mean(aa)
est_mat <- matrix(0, ncol = length(0:143), nrow = 10000)
# calculate standard deviation through time
for(i in 0:143){
  est_mat[,i+1] <- 1/sqrt(exp(t0_log + t1_log * aa[i+1]))
}
my_med <- apply(est_mat,2, quantile, probs = c(0.025, 0.5, 0.975))

# plot out the standard deviation through time
tiff("sd_plot_time.tiff", height = 4, width = 4, units = "in",
     res = 600, compression = "lzw")
plot(my_med[2,], type = "l", ylab = "standard deviation",
     xlab ="year", xlim = c(0, 150))
lines(my_med[1,], lty = 2)
lines(my_med[3,], lty = 2)
dev.off()

# get quantiles of the intercept
meds <- apply(covs[,-c(1:2)], 2, quantile, probs = c(0.025, 0.5, 0.975))
# get just the beta values
betas <- meds[,grep("^B", colnames(meds))]
# make it into a matrix
b1 <- matrix(betas[2,], ncol = 3, nrow = nrow(betas)/3)
betas <- round(betas, 2)

# for use in a supplementary table
beta_print <- paste0(betas[2,], " (", betas[1,], " - ", betas[3,], ")")

beta_print <- data.frame(matrix(beta_print, ncol = 3, nrow = nrow(betas)/3))
colnames(beta_print) <- c("Intercept", "Linear Time", "C02")


# calculate the number of days each species has changed through
# time with our predictions
y_diff <- round(apply(preds[, at_start] - preds[,at_end],
                      2, quantile, probs = c(0.025, 0.5, 0.975)), 2)



# get signs from meds (used for calculating best fit model for each species)
ms <- apply(t(meds), 2, sign)
# if they are all the same then the signs will by all negative or positive
msum <- rowSums(ms)

negs <- which(msum==-3) # negative effects
pos <- which(msum ==3) # positive effects

# significant year species
sigyr <- c(colnames(meds[,negs])[grep("^B\\[\\d\\d?\\,2",colnames(meds[,negs]))],
           colnames(meds[,pos])[grep("^B\\[\\d\\d?\\,2",colnames(meds[,pos]))])
# significant CO2 species
sigco <- c(colnames(meds[,negs])[grep("^B\\[\\d\\d?\\,3",colnames(meds[,negs]))],
           colnames(meds[,pos])[grep("^B\\[\\d\\d?\\,3",colnames(meds[,pos]))])
# get the species numeric codes
yrnum <- as.numeric(gsub("B\\[(\\d\\d?)..*", "\\1", sigyr))
# get the species numeric codes for CO2
conum <- as.numeric(gsub("B\\[(\\d\\d?)..*", "\\1", sigco))
# a few species had significant effects on both
bonum <- yrnum[yrnum %in% conum]
# a vector of species who have at least one significant effect
anum <- c(yrnum, conum, bonum)
# species who have no significant effects
none <- c(1:73)[-which(c(1:73) %in% anum)]
# remove bonum species from yrnum and conum
yrnum <- yrnum[-which(yrnum %in% bonum)]
conum <- conum[-which(conum %in% bonum)]
# get the 4 character species codes
yrspc <- levels(fdt$species)[yrnum]
cospc <- levels(fdt$species)[conum]
bothspc <- levels(fdt$species)[bonum]
nonespc <- levels(fdt$species)[none]
# make a table for species betas

ints <- matrix(meds[2,1:216], ncol = 3, nrow = ncol(betas)/nrow(betas)) # intercepts
lows <- matrix(meds[1,1:216], ncol = 3, nrow = ncol(betas)/nrow(betas)) # lower quantile
highs <- matrix(meds[3,1:216], ncol = 3, nrow = ncol(betas)/nrow(betas)) # upper quantile

# starting to make the species data.frame for suppliementary material
sp_frame <- data.frame(species = migt$cmn, group = migt$migstat,
                       code = migt$species)
# round to 2
sp_frame$int <- sprintf("%.2f", round(ints[,1],2))
b <- paste(sprintf("%.2f", round(lows[,1],2)), " - ", 
           sprintf("%.2f", round(highs[,1],2)), sep = "")
sp_frame$int <- paste(sp_frame$int, " [", b, "]", sep = "")
# do for year effect
sp_frame$time <- sprintf("%.2f",round(ints[,2],2))
b <- paste(sprintf("%.2f", round(lows[,2],2)), " - ", 
           sprintf("%.2f", round(highs[,2],2)), sep = "")

sp_frame$time <- paste(sp_frame$time, " [", b, "]", sep = "")

# creating a new column that will have a star if there is a significant
# effect
sp_frame$ts <- " "
sp_frame$ts[yrnum] <- "*"
# do for CO2
sp_frame$c02 <- sprintf("%.2f",round(ints[,3],2))
b <- paste(sprintf("%.2f", round(lows[,3],2)), " - ", 
           sprintf("%.2f", round(highs[,3],2)), sep = "")
sp_frame$c02 <- paste(sp_frame$c02, " [", b, "]", sep = "")
# CO2 column with star
sp_frame$cs <- " "
sp_frame$cs[conum] <- "*"
# calculate advancement for yrnum species, holding CO2 at its median value
# of zero.
sp_frame$day_increase <- " "
sp_frame$day_increase[yrnum] <- sprintf("%.2f", round(143 * ints[yrnum,2],2))
b <- paste(sprintf("%.2f", round(highs[yrnum,2]*143,2) ), " - ", 
           sprintf("%.2f", round(lows[yrnum,2]*143,2) ), sep = "")
sp_frame$day_increase[yrnum] <- paste(sp_frame$day_increase[yrnum],
                                      " [", b, "]", sep = "")

bc <- 40.2785 # CO2 at last year
sp_frame$day_increase[bonum] <- 
  sprintf("%.2f", round(143 * ints[bonum,2] + bc * ints[bonum,3],2))
b <- paste(sprintf("%.2f", round((lows[,2]*143) + (bc * lows[,3]),2) ), " - ", 
           sprintf("%.2f", round(highs[,2]*143 + bc * highs[,3],2) ), sep = "")
sp_frame$day_increase[bonum] <- paste(sp_frame$day_increase[bonum],
                                      " [", b[bonum], "]", sep = "")

ll <- lows[,2]*143 + bc * lows[,3]
hh <- highs[,2]*143 + bc * highs[,3]
bb <- cbind(ll, hh)
bb <- rowSums(sign(bb))


getsig <- function(noest){
  for(i in 1:length(noest)){
    a <- noest[[i]][[1]][,144] - noest[[i]][[1]][,1]
  }
}

sp_frame <- sp_frame[order(sp_frame$group),]
write.csv(sp_frame, "species_betas_no_eabl.csv", row.names = FALSE)

# get just the species number
# sig spec yr


# type 1 = time, 2 = c02, 3 = both, 4 = none


make_proj <- function(nv = NULL, mm = NULL, type = NULL, ress = NULL,
                      s = NULL){
  
  qt_ls <- function(prob, df, mu, a) qt(prob, df)*a + mu
  res2 <- as.numeric(ress$x)
  yr <- 0:143
  pm <- pm_low <- pm_high <-  array(0, dim = c(nrow(mm), length(nv), length(yr)))
  pm_list <-pmlow_list <- pmhigh_list<-pm_prob <-  vector("list", length = length(nv))
  nu <- median(mm[,grep("nu", colnames(mm))])
  for(i in 1:length(nv)){
    tk <- paste0("^B\\[", nv[i],",1\\]|^B\\[", nv[i],",2\\]|^B\\[", nv[i],",3\\]")
    sp_beta <- as.matrix(mm[,grep(tk, colnames(mm))])
    
    
    
    for(k in 1:length(yr)){
      if(type == 1){
        pm[,i,k] <- sp_beta %*% c(1, yr[k], median(res2))
        pm_high[,i,k] <- qt_ls(0.975, nu,pm[,i,k], median(s[,k]))
        pm_low[,i,k] <- qt_ls(0.025, nu, pm[,i,k], median(s[,k]))
      }
      if(type == 2){
        pm[,i,k] <- sp_beta %*% c(1, median(yr), res2[k] )
        pm_high[,i,k] <- qt_ls(0.975, nu, pm[,i,k], median(s[,k]))
        pm_low[,i,k] <- qt_ls(0.025,  nu, pm[,i,k],median(s[,k]))
      }
      if(type ==3){
        pm[,i,k] <- sp_beta %*% c(1, yr[k], res2[k])
        pm_high[,i,k] <- qt_ls(0.975,median(nu), pm[,i,k], median(s[,k]))
        pm_low[,i,k] <- qt_ls(0.025, median(nu), pm[,i,k],  median(s[,k]))
      }
      if(type ==4){
        pm[,i,k] <- sp_beta %*% c(1, median(yr), median(res2))
        pm_high[,i,k] <- qt_ls(0.975,nu, pm[,i,k], median(s[,k]))
        pm_low[,i,k] <- qt_ls(0.025,nu, pm[,i,k],  median(s[,k]))
      }
    }
    pm_list[[i]] <- apply(pm[,i,], 2, quantile, probs = c(0.025, 0.5, 0.975))
    pmlow_list[[i]] <- apply(pm_low[,i,], 2, median)
    pmhigh_list[[i]] <- apply(pm_high[,i,], 2, median)
    pm_prob[[i]] <- median(pm[,i,144]/pm[,i,1])
  }
  return(list(pm_list, pmlow_list, pmhigh_list, pm_prob))
}

# to do. Make the inits figure with 4, 7, 72, 15 to illustrate 
# varying responses.

yrest <- make_proj(nv = yrnum, mm = covs, type = 1, ress = ress,
                   s = est_mat)
coest <- make_proj(nv = conum, mm = covs, type = 2, ress = ress, s = my_med)
boest <- make_proj(nv = bonum, mm = covs, type = 3, ress = ress, my_med)
nochange <- 1:72
nochange <- nochange[-which(nochange %in% c(yrnum, conum, bonum))]

noest <- make_proj(nv = nochange, mm = covs, type = 4, ress = ress, my_med)

changes <- matrix(0, ncol = 3, nrow = 72)

for(i in 1:length(yrnum)) {
  b1 <- c(yrest[[1]][[i]][2,144] - yrest[[1]][[i]][2,1],
          yrest[[2]][[i]][144] - yrest[[2]][[i]][1],
          yrest[[3]][[i]][144] - yrest[[3]][[i]][1])
  
  changes[yrnum[i],] <- b1
}
for(i in 1:length(conum)) {
  b1 <- c(coest[[1]][[i]][2,144] - coest[[1]][[i]][2,1],
          coest[[2]][[i]][144] - coest[[2]][[i]][1],
          coest[[3]][[i]][144] - coest[[3]][[i]][1])
  
  changes[conum[i],] <- b1
}
for(i in 1:length(bonum)) {
  b1 <- c(boest[[1]][[i]][2,144] - boest[[1]][[i]][2,1],
          boest[[2]][[i]][144] -   boest[[2]][[i]][1],
          boest[[3]][[i]][144] -   boest[[3]][[i]][1])
  
  changes[bonum[i],] <- b1
}

for(i in 1:length(nochange)) {
  b1 <- c(noest[[1]][[i]][2,144] - noest[[1]][[i]][2,1],
          noest[[2]][[i]][144] -   noest[[2]][[i]][1],
          noest[[3]][[i]][144] -   noest[[3]][[i]][1])
  
  changes[nochange[i],] <- b1
}

changes <- round(changes, 2)

tch <- which(changes[,1] == 0)

cha <- paste0(changes[,1], " [", changes[,3], " - ", changes[,2], "]")


jj2 %>% group_by(group) %>% summarise(me = mean(test), mi = min(test),
                                      ma = max(test))

sp_frame$day_increase <- 0
sp_frame$day_increase[-nochange] <- cha[-nochange]
sp_frame$rule <- 0
#sp_frame$rule[sp_frame$code %in% first_inclu] <- 1
#sp_frame$rule <- 2 - sp_frame$rule
yrspc <- spf[yrnum]

new_aou <- foreign::read.dbf("LIST17.DBF")

new_aou <- new_aou[new_aou$SPEC %in% sp_frame$code,]
new_levs <- new_aou$SPEC
colnames(new_aou)[2] <- "code"


sp_frame <- dplyr::left_join(sp_frame, new_aou[,c(2,5)], by = "code" )
sp_frame$code <- factor(sp_frame$code, levels = new_levs)
sp_frame <- sp_frame[with(sp_frame, order(group, code)),]

write.csv(sp_frame, "field_table1.csv", row.names = FALSE)
pdf("yr_effect.pdf", width = 6, height = 10)
windows(6,10)
par(mfrow = c(2,1))
for(i in 1:length(yrnum)){
  dat <- fdt[fdt$species== yrspc[i],]
  plot(1~1, xlim = c(0,145), ylim = c(50, 250) ,
       xlab = "Year", ylab = "Initial lay date", bty = "n", 
       xaxt = "n", yaxt = "n", main = migt$cmn[migt$code== dat$species[1]] )
  axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
  axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
  axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
  mtext(text = seq(50, 250, 50), 2, line = .75, at = seq(50, 250, 50), las = 1, cex = 1)
  
  axis(1, at = seq(0, 150, 25), labels=F, tck=-.035)
  axis(1, at = seq(0, 150, 12.5), labels=F, tck=-.025)
  axis(1, at = seq(-0, 150, 6.25), labels=F, tck=-.015)
  mtext(text = c(floor(seq(1872, 2015, length.out = 7))), 1, line = .75, at = seq(0, 150, 25), las = 1, cex = 1)
  
  lines(yrest[[1]][[i]][2,] , lwd = 2)
  lines(yrest[[1]][[i]][1,] , lwd = 2, lty = 2)
  lines(yrest[[1]][[i]][3,] , lwd = 2, lty = 2)
  
  lines(yrest[[2]][[i]] , lwd = 2, lty = 3)
  lines(yrest[[3]][[i]] , lwd = 2, lty = 3)
  
  
  
  #lines(c(yrest[[i]][2,] + abs(median(res2))), lwd = 2)
  #lines(c(yrest[[i]][1,] + abs(median(res2))), lwd = 2, lty = 2)
  #lines(c(yrest[[i]][3,] + abs(median(res2))), lwd = 2, lty = 2)
  points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
         pch = 16)
  
}
dev.off()
# include warbler
my_plots <- 72

pdf("all_sp.pdf", width = 6, height = 10)
par(mfrow = c(2,1))
for(i in 1:72){
  dat <- fdt[fdt$species== migt$code[i],]
  plot(1~1, xlim = c(0,145), ylim = c(50, 250) ,
       xlab = "Year", ylab = "Initial lay date", bty = "n", 
       xaxt = "n", yaxt = "n", main = migt$cmn[i])
  axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
  axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
  axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
  mtext(text = seq(50, 250, 50), 2, line = .75, at = seq(50, 250, 50), las = 1, cex = 1)
  
  axis(1, at = seq(0, 150, 25), labels=F, tck=-.035)
  axis(1, at = seq(0, 150, 12.5), labels=F, tck=-.025)
  axis(1, at = seq(-0, 150, 6.25), labels=F, tck=-.015)
  mtext(text = c(floor(seq(1872, 2015, length.out = 7))), 1, line = .75, at = seq(0, 150, 25), las = 1, cex = 1)
  
  lines(noest[[1]][[i]][2,] , lwd = 2)
  lines(noest[[1]][[i]][1,] , lwd = 2, lty = 2)
  lines(noest[[1]][[i]][3,] , lwd = 2, lty = 2)
  
  lines(noest[[2]][[i]] , lwd = 2, lty = 3)
  lines(noest[[3]][[i]] , lwd = 2, lty = 3)
  
  
  
  #lines(c(yrest[[i]][2,] + abs(median(res2))), lwd = 2)
  #lines(c(yrest[[i]][1,] + abs(median(res2))), lwd = 2, lty = 2)
  #lines(c(yrest[[i]][3,] + abs(median(res2))), lwd = 2, lty = 2)
  points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
         pch = 16)
  
}
dev.off()

x <- 
  wireframe()

pdf("both_est.pdf", width = 6, height = 10)
par(mfrow = c(2,1))
for(i in 1:length(bonum)){
  dat <- fdt[fdt$species== bothspc[i],]
  plot(1~1, xlim = c(0,145), ylim = c(50, 250) ,
       xlab = "Year", ylab = "Initial lay date", bty = "n", 
       xaxt = "n", yaxt = "n", main = dat$species[i])
  axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
  axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
  axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
  mtext(text = seq(50, 250, 50), 2, line = .75, at = seq(50, 250, 50), las = 1, cex = 1)
  
  axis(1, at = seq(0, 150, 25), labels=F, tck=-.035)
  axis(1, at = seq(0, 150, 12.5), labels=F, tck=-.025)
  axis(1, at = seq(-0, 150, 6.25), labels=F, tck=-.015)
  mtext(text = c(floor(seq(1872, 2015, length.out = 7))), 1, line = .75, at = seq(0, 150, 25), las = 1, cex = 1)
  
  lines(boest[[1]][[i]][2,] , lwd = 2)
  lines(boest[[1]][[i]][1,] , lwd = 2, lty = 2)
  lines(boest[[1]][[i]][3,] , lwd = 2, lty = 2)
  
  lines(boest[[2]][[i]] , lwd = 2, lty = 3)
  lines(boest[[3]][[i]] , lwd = 2, lty = 3)
  
  
  
  #lines(c(yrest[[i]][2,] + abs(median(res2))), lwd = 2)
  #lines(c(yrest[[i]][1,] + abs(median(res2))), lwd = 2, lty = 2)
  #lines(c(yrest[[i]][3,] + abs(median(res2))), lwd = 2, lty = 2)
  points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
         pch = 16)
  
}
dev.off()


pdf("co_est.pdf", width = 6, height = 10)
par(mfrow = c(2,1))
for(i in 1:length(conum)){
  dat <- fdt[fdt$species== cospc[i],]
  plot(1~1, xlim = c(0,145), ylim = c(50, 250) ,
       xlab = "Year", ylab = "Initial lay date", bty = "n", 
       xaxt = "n", yaxt = "n", main = dat$species[i])
  axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
  axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
  axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
  mtext(text = seq(50, 250, 50), 2, line = .75, at = seq(50, 250, 50), las = 1, cex = 1)
  
  axis(1, at = seq(0, 150, 25), labels=F, tck=-.035)
  axis(1, at = seq(0, 150, 12.5), labels=F, tck=-.025)
  axis(1, at = seq(-0, 150, 6.25), labels=F, tck=-.015)
  mtext(text = c(floor(seq(1872, 2015, length.out = 7))), 1, line = .75, at = seq(0, 150, 25), las = 1, cex = 1)
  
  lines(coest[[1]][[i]][2,] , lwd = 2)
  lines(coest[[1]][[i]][1,] , lwd = 2, lty = 2)
  lines(coest[[1]][[i]][3,] , lwd = 2, lty = 2)
  
  lines(coest[[2]][[i]] , lwd = 2, lty = 3)
  lines(coest[[3]][[i]] , lwd = 2, lty = 3)
  
  
  
  #lines(c(yrest[[i]][2,] + abs(median(res2))), lwd = 2)
  #lines(c(yrest[[i]][1,] + abs(median(res2))), lwd = 2, lty = 2)
  #lines(c(yrest[[i]][3,] + abs(median(res2))), lwd = 2, lty = 2)
  points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
         pch = 16)
  
}
dev.off()


pdf("no_est.pdf", width = 6, height = 10)
par(mfrow = c(2,1))
for(i in 1:length(none)){
  sp <- levels(fdt$species)[none[i]]
  dat <- fdt[fdt$species== sp,]
  plot(1~1, xlim = c(0,145), ylim = c(50, 250) ,
       xlab = "Year", ylab = "Initial lay date", bty = "n", 
       xaxt = "n", yaxt = "n", main = dat$species[i])
  axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
  axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
  axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
  mtext(text = seq(50, 250, 50), 2, line = .75, at = seq(50, 250, 50), las = 1, cex = 1)
  
  axis(1, at = seq(0, 150, 25), labels=F, tck=-.035)
  axis(1, at = seq(0, 150, 12.5), labels=F, tck=-.025)
  axis(1, at = seq(-0, 150, 6.25), labels=F, tck=-.015)
  mtext(text = c(floor(seq(1872, 2015, length.out = 7))), 1, line = .75, at = seq(0, 150, 25), las = 1, cex = 1)
  
  lines(noest[[1]][[i]][2,] , lwd = 2)
  lines(noest[[1]][[i]][1,] , lwd = 2, lty = 2)
  lines(noest[[1]][[i]][3,] , lwd = 2, lty = 2)
  
  lines(noest[[2]][[i]] , lwd = 2, lty = 3)
  lines(noest[[3]][[i]] , lwd = 2, lty = 3)
  
  
  
  #lines(c(yrest[[i]][2,] + abs(median(res2))), lwd = 2)
  #lines(c(yrest[[i]][1,] + abs(median(res2))), lwd = 2, lty = 2)
  #lines(c(yrest[[i]][3,] + abs(median(res2))), lwd = 2, lty = 2)
  points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
         pch = 16)
  
}
dev.off()


# x, y, z variables
x <- mtcars$wt
y <- mtcars$disp
z <- mtcars$mpg
# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane


# final figure for ms
windows(8,8)

pdf("nesting_figure1_noeabl.pdf", height = 8, width = 8)
m1 <- matrix(c(5,5,rep(1,8), rep(2, 8)), ncol = 18, nrow = 8,
             byrow = TRUE)
#m4 <- rep(8, 18)
m2 <- matrix(c(5,5,rep(3,8), rep(4, 8)), ncol = 18, nrow = 8,
             byrow = TRUE)
m3 <- matrix(c(7,7, rep(6, 16)), ncol = 18, nrow = 2, byrow = TRUE)

m <- rbind(m1, m2, m3)
layout(m)
par( mar = c(1.75,1.75, 1.25,0.5),
     oma = c(1.75,1.75,1.25,0.5) + 0.1)

sp <- levels(fdt$species)[71]
dat <- fdt[fdt$species== sp,]
plot(1~1, xlim = c(0,145), ylim = c(50, 250) ,
     xlab = "Year", ylab = "Initial lay date", bty = "n", 
     xaxt = "n", yaxt = "n")
axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
mtext(text = seq(50, 250, 50), 2, line = 1, 
      at = seq(50, 250, 50), las = 1, cex = 1.5)
mtext(text = "a)", 3, line = 0.1, at = 0, las = 1, cex = 1.5)
text(x = -50, y = 35, labels = "Initial lay date (Julian day)", 
     cex = 4, xpd = NA, srt = 90)

lines(yrest[[1]][[13]][2,] , lwd = 2)
lines(yrest[[1]][[13]][1,] , lwd = 2, lty = 2)
lines(yrest[[1]][[13]][3,] , lwd = 2, lty = 2)

lines(yrest[[2]][[13]] , lwd = 2, lty = 3)
lines(yrest[[3]][[13]] , lwd = 2, lty = 3)

points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
       pch = 16, cex = 1.5)

axis(1, at = seq(0, 150, 25), labels=F, tck=-.035)
axis(1, at = seq(0, 150, 12.5), labels=F, tck=-.025)
axis(1, at = seq(-0, 150, 6.25), labels=F, tck=-.015)
#text(x = seq(0, 150, 25), y = rep(30, 7), 
#     labels = c(floor(seq(1872, 2015, length.out = 7))),
#     adj = 1, las = 1, cex = 1.5, srt = 45, xpd = NA)

sp <- levels(fdt$species)[33]
dat <- fdt[fdt$species== sp,]

plot(1~1, xlim = c(0,145), ylim = c(50, 250) ,
     xlab = "Year", ylab = "Initial lay date", bty = "n", 
     xaxt = "n", yaxt = "n")
axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
#mtext(text = seq(50, 250, 50), 2, line = .75, 
#      at = seq(50, 250, 50), las = 1, cex = 1.5)
mtext(text = "b)", 3, line = 0.1, at = 0, las = 1, cex = 1.5)

axis(1, at = seq(0, 150, 25), labels=F, tck=-.035)
axis(1, at = seq(0, 150, 12.5), labels=F, tck=-.025)
axis(1, at = seq(-0, 150, 6.25), labels=F, tck=-.015)
#text(x = seq(0, 150, 25), y = rep(30, 7), 
#     labels = c(floor(seq(1872, 2015, length.out = 7))),
#     adj = 1, las = 1, cex = 1.5, srt = 45, xpd = NA)

lines(coest[[1]][[4]][2,] , lwd = 2)
lines(coest[[1]][[4]][1,] , lwd = 2, lty = 2)
lines(coest[[1]][[4]][3,] , lwd = 2, lty = 2)

lines(coest[[2]][[4]] , lwd = 2, lty = 3)
lines(coest[[3]][[4]] , lwd = 2, lty = 3)

points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
       pch = 16, cex = 1.5)

sp <- levels(fdt$species)[15]
dat <- fdt[fdt$species== sp,]

plot(1~1, xlim = c(0,145), ylim = c(50, 250) ,
     xlab = "Year", ylab = "Initial lay date", bty = "n", 
     xaxt = "n", yaxt = "n")
axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
mtext(text = seq(50, 250, 50), 2, line = 1, 
      at = seq(50, 250, 50), las = 1, cex = 1.5)
mtext(text = "c)", 3, line = 0.1, at = 0, las = 1, cex = 1.5)

axis(1, at = seq(0, 150, 25), labels=F, tck=-.035)
axis(1, at = seq(0, 150, 12.5), labels=F, tck=-.025)
axis(1, at = seq(-0, 150, 6.25), labels=F, tck=-.015)
text(x = seq(0, 150, 25), y = rep(30, 7), 
     labels = c(floor(seq(1872, 2015, length.out = 7))),
     adj = 1, las = 1, cex = 2, srt = 45, xpd = NA)


lines(noest[[1]][[11]][2,] , lwd = 2)
lines(noest[[1]][[11]][1,] , lwd = 2, lty = 2)
lines(noest[[1]][[11]][3,] , lwd = 2, lty = 2)

lines(noest[[2]][[11]] , lwd = 2, lty = 3)
lines(noest[[3]][[11]] , lwd = 2, lty = 3)

points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
       pch = 16, cex = 1.5)

sp <- levels(fdt$species)[7]
dat <- fdt[fdt$species== sp,]
plot(1~1, xlim = c(0,145), ylim = c(50, 250) ,
     xlab = "Year", ylab = "Initial lay date (Julian day)", bty = "n", 
     xaxt = "n", yaxt = "n")
axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
#mtext(text = seq(50, 250, 50), 2, line = .75, 
#     at = seq(50, 250, 50), las = 1, cex = 1.5)
mtext(text = "d)", 3, line = 0.1, at = 0, las = 1, cex = 1.5)

axis(1, at = seq(0, 150, 25), labels=F, tck=-.035)
axis(1, at = seq(0, 150, 12.5), labels=F, tck=-.025)
axis(1, at = seq(-0, 150, 6.25), labels=F, tck=-.015)
text(x = seq(0, 150, 25), y = rep(30, 7), 
     labels = c(floor(seq(1872, 2015, length.out = 7))),
     adj = 1, las = 1, cex = 2, srt = 45, xpd = NA)
text(x = -15, y = -25, labels = "Year", cex = 4, xpd = NA)

lines(boest[[1]][[1]][2,] , lwd = 2)
lines(boest[[1]][[1]][1,] , lwd = 2, lty = 2)
lines(boest[[1]][[1]][3,] , lwd = 2, lty = 2)

lines(boest[[2]][[1]] , lwd = 2, lty = 3)
lines(boest[[3]][[1]] , lwd = 2, lty = 3)

points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
       pch = 16, cex = 1.5)

dev.off()
###



#lines(c(yrest[[i]][2,] + abs(median(res2))), lwd = 2)
#lines(c(yrest[[i]][1,] + abs(median(res2))), lwd = 2, lty = 2)
#lines(c(yrest[[i]][3,] + abs(median(res2))), lwd = 2, lty = 2)
points(dat$jdate ~ c(dat$year-1872), col = alpha("black", 0.5),
       pch = 16, cex = 1.5)

cept <- t(meds[,1:72])
hm <- t(meds[,73:216])
yer <- hm[1:72,2]
co <- hm[-c(1:72),2]
windows(4,4)
pdf("bates_figure2.pdf", height = 4, width = 4)
par( mar = c(5,5.5, 1.25,1))
plot(1~1, xlim = c(-0.4, 0.35), ylim = c(-0.75, 0.75) ,
     xlab = "Species specific change in mean lay date\nto year (Julian day)",
     ylab = "Species specific change in mean lay date\nto changes in global CO2 (Julian day)", bty = "n", 
     xaxt = "n", yaxt = "n")
axis(2, at = seq(-0.75, 0.75, .25), labels=F, tck=-.035)
axis(2, at = seq(-0.75, 0.75, .125), labels=F, tck=-.025)
mtext(seq(-0.75, 0.75, 0.25), 2, 0.8, at = seq(-0.75, 0.75, 0.25),
      las = 2)

axis(1, at = seq(-0.4, 0.3, .1), labels=F, tck=-.035)
axis(1, at = seq(-0.4, 0.3, .05), labels=F, tck=-.025)
mytext <- seq(-1, 1, 0.1)
mytext <- mytext[mytext>-0.4 & mytext <0.4]
mtext(mytext, 1, 0.8, at = seq(-0.4, 0.3, 0.1) )
#abline(coef(lm(cept~slo)))
abline(v = 0, lty = 2)
abline(h = 0, lty =2 )
lines(x=c(.3,3.5), y = c(0,0), col = "white", lwd = 2)
spl <- which(levels(spcode)=="YEWA")
lines(y = c(co[spl], co[spl]-0.05), x = c(yer[spl], yer[spl]-0.05)) 
spl <- which(levels(spcode)=="AMRO")
lines(y = c(co[spl], co[spl]+0.05), x = c(yer[spl], yer[spl]+0.05))
spl <- which(levels(spcode)=="BLJA")
lines(y = c(co[spl], co[spl]+0.05), x = c(yer[spl], yer[spl]+0.05))
spl <- which(levels(spcode)=="FISP")
lines(y = c(co[spl], co[spl]-0.05), x = c(yer[spl], yer[spl]-0.05))
points(co[none]~yer[none], pch = 21, bg = "#717e9e",   cex = 1.3)
points(co[yrnum]~yer[yrnum], pch = 24, bg = "#f3bf2b", cex = 1.3)
points(co[conum]~yer[conum], pch = 22, bg = "#9d3a35", cex = 1.3)
points(co[bonum]~yer[bonum], pch = 23, bg = "black",   cex = 1.3)
legend(x=-0.015, y = 0.75, pch = c(24,22,23,21), bty = "n",
       legend = c("Response to year", "Response to CO2", 
                  "Response to both", "No response"), cex = 0.75,
       pt.bg = c("#f3bf2b", "#9d3a35", "black", "#717e9e"),
       pt.cex = 1.3)



spl <- which(levels(spcode)=="YEWA")
#points(co[spl]~yer[spl], pch = 19)
text(x = yer[spl]-0.05, y = co[spl]-0.05, labels = "YEWA", cex = 0.75, pos = 1,
     offset = 0)


#abline(a=0,b = -1.875, lty = 2)
#abline(a=0, b=1.875, lty = 2)
  # amgo
spl <- which(levels(spcode)=="FISP")
#points(co[spl]~yer[spl], pch = 19)
text(x = yer[spl]-0.085, y = co[spl]-0.04, labels = "FISP", cex = 0.75, pos = 1,
     offset = 0)
#amro
spl <- which(levels(spcode)=="AMRO")
#points(co[spl]~yer[spl], pch = 19)
text(x = yer[spl]+0.05, y = co[spl]+0.05, labels = "AMRO", cex = 0.75, pos = 3,
     offset = 0)
#text(x = -0.4, y = 0.70, labels = as.roman(1), family = "serif")
#text(x = 0.03, y = 0.70, labels = as.roman(2), family = "serif")
#text(x = -0.4, y = -0.05, labels = as.roman(3), family = "serif")
#text(x = , y = 0.70, labels = as.roman(4), family = "serif")
#yewa
spl <- which(levels(spcode)=="BLJA")
#points(co[spl]~yer[spl], pch = 19)
text(x = yer[spl]+0.1, y = co[spl]+0.05, labels = "BLJA", cex = 0.75,
     pos = 3,
     offset = 0)
dev.off()
#BLJA
spl <- which(levels(spcode)=="BLJA")
points(co[spl]~yer[spl], pch = 19)
text(x = cept[spl]+0.02, y = slo[spl]+0.05, labels = "BLJA", cex = 1, pos = 4,
     offset = 0.3)
spl <- which(levels(spcode)=="BHCO")
points(co[spl]~yer[spl], pch = 19)
text(x = cept[spl]+0.02, y = slo[spl]+0.05, labels = "BLJA", cex = 1, pos = 4,
     offset = 0.3)


plot(1~1, xlim = c(-0.75, 0.5), ylim = c(-0.75, 0.5) ,
     xlab = "Change in lay date per year (Julain day)",
     ylab = "Initial lay date (Julian day)", bty = "n", 
     xaxt = "n", yaxt = "n")

axis(2, at = seq(-0.75, 0.5, 0.25), labels=F, tck=-.035)
axis(2, at = seq(-0.75, 0.5, .125), labels=F, tck=-.025)
axis(2, at = seq(-0.75, 0.5, 0.0625), labels=F, tck=-.015)
mtext(text = seq(-0.75, 0.5, 0.25), 2, line = .75, 
      at = seq(-0.75, 0.5, 0.25), las = 1, cex = 1.5)

axis(1, at = seq(-0.75, 0.5, 0.25), labels=F, tck=-.035)
axis(1, at = seq(-0.75, 0.5, .125), labels=F, tck=-.025)
axis(1, at = seq(-0.75, 0.5, 0.0625), labels=F, tck=-.015)
mtext(text = seq(-0.75, 0.5, 0.25), 1, line = .75, 
      at = seq(-0.75, 0.5, 0.25), las = 1, cex = 1.5)

long <- which(migt$migstat == "long")
shor <-which(migt$migstat == "short")
res <- which(migt$migstat == "resident")


points(x = slo[shor,2], y = cept[shor,2], 
       pch = 21, bg = "#d95f02", cex = 1.5)
points(x = slo[res,2], y = cept[res,2], 
       pch = 23, bg = "#7570b3", cex = 1.5)
points(x = slo[long,2], y = cept[long,2], 
       pch = 24, bg = "#1b9e77", cex = 1.5)

plot(1~1, xlim = c(-0.75, 0.5), ylim = c(50, 250) ,
     xlab = "Change in lay date per year (Julain day)",
     ylab = "Initial lay date (Julian day)", bty = "n", 
     xaxt = "n", yaxt = "n")

axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
axis(2, at = seq(50, 250, 12.5), labels=F, tck=-.015)
mtext(text = seq(50, 250, 50), 2, line = .75, 
      at = seq(50, 250, 50), las = 1, cex = 1.5)

axis(1, at = seq(-0.75, 0.5, 0.25), labels=F, tck=-.035)
axis(1, at = seq(-0.75, 0.5, .125), labels=F, tck=-.025)
axis(1, at = seq(-0.75, 0.5, 0.0625), labels=F, tck=-.015)
mtext(text = seq(-0.75, 0.5, 0.25), 1, line = .75, 
      at = seq(-0.75, 0.5, 0.25), las = 1, cex = 1.5)

long <- which(migt$migstat == "long")
shor <-which(migt$migstat == "short")
res <- which(migt$migstat == "resident")


points(x = slo[shor,2], y = cept[shor,2], 
       pch = 21, bg = "#d95f02", cex = 1.5)
points(x = slo[res,2], y = cept[res,2], 
       pch = 23, bg = "#7570b3", cex = 1.5)
points(x = slo[long,2], y = cept[long,2], 
       pch = 24, bg = "#1b9e77", cex = 1.5)
plot()

x <- ress
y <- 
  library(plot3D)

x <- dat$year
y <- ress[c(x-1872+1),]
z <- dat$jdate
scatter3D(x, y, z, xlim = range(x), ylim = range(y),
          zlim = range(z), pch = 20, cex = 2, theta = 45,
          phi = 0, xlab = "Year", ylab = "c02", zlab = "Julain date",
          col = NULL), surf = list(x = c(1872:2015), y = ress[,1],
                                   z = coest[[1]][[1]][2,],
                                   fitpoints = coest[[1]][[1]][2,]))


scatter3D(x, y, z, pch = 18, cex = 2, 
          theta = 20, phi = 20, ticktype = "detailed",
          xlab = "wt", ylab = "disp", zlab = "mpg",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints), main = "mtcars")


upq <- apply(mm, 2, quantile, probs = c(0.975))

quants <- apply(mm, 2, quantile, probs = c(0.025, 0.5, 0.975))

bs <- grep("^B\\[", colnames(mm))
# put the B's in a matrix
b <- data.frame(matrix(meds[bs], nrow = nspec, ncol = 3))
bq <- data.frame(matrix(upq[bs], nrow = nspec, ncol = 3))
bwm <- cbind(b, migt[,-3], bq[,2])
colnames(bwm) <- c("X1", "X2", "X3", "cmn", "migstat", "remove", "code", "upq")

to_c <- bwm[bwm$upq<0,]

bwm$cmn[bwm$upq>0]
plot(1,1, ylim = c(100,250), xlim = c(-0.3,0.2), type = "n", xlab = 
       "rate of increase (days per year)", ylab = "initial lay date", bty = "n", xaxt = "n", yaxt = "n")

axis(2, at = seq(100, 250, 25), labels=F, tck=-.035)
axis(2, at = seq(100, 250, 12.5), labels=F, tck=-.025)
axis(2, at = seq(100, 250, 12.5/2), labels=F, tck=-.015)
mtext(text = seq(100, 250, 25), 2, line = .75, at = seq(100, 250, 25), las = 1, cex = 1)

axis(1, at = seq(-0.3, 0.2, 0.1), labels=F, tck=-.035)
axis(1, at = seq(-0.3, 0.2, 0.05), labels=F, tck=-.025)
axis(1, at = seq(-0.3, 0.2, 0.025), labels=F, tck=-.015)
mtext(text = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2), 1, line = .75, at = signif(seq(-0.3, 0.2, 0.1), 2), las = 1, cex = 1)

# plot residents
with(bwm, points(X1[migstat=="resident"]~X2[migstat=="resident"], pch = 21, bg = "gray50", cex = 1.5))
with(bwm, points(X1[migstat=="long"]~X2[migstat=="long"], pch = 21, bg = "green", cex = 1.5))
with(bwm, points(X1[migstat=="short"]~X2[migstat=="short"], pch = 21, bg = "purple", cex = 1.5))
legend("topright", c("resident", "long", "short"), pch = c(16,16,16), col = c("gray50", "green", "purple") )


boxplot(bwm$X1~bwm$migstat)
boxplot(bwm$X2~bwm$migstat)
boxplot(bwm$X3~bwm$migstat)
with(to_c, points(X1~X2, pch = 1, col = "black", cex = 3))

to_c$cmn

# time to make species specific plots.

# we have the meds quantile so lets get y_preds from it and bind it
# with x pred

yp <- data.frame(cbind(t(quants[,grep("y_pred", colnames(quants))]),x_pred))
colnames(yp) <- c("low", "med", "high", "spn", "yr")
yp$sp <- rep(levels(span), each = nyr)

sp_cof <- data.frame(cbind(t(quants[,grep("^B\\[\\w\\w?,1", colnames(quants))]),
                           t(quants[,grep("^B\\[\\w\\w?,2", colnames(quants))])))
sp_cof$sp <- levels(span)
colnames(sp_cof) <- c("il", "im", "ih", "sl", "sm", "sh", "sp")    

# get order to which we want to plot them

spor <- levels(span)
yv <- c(1815:2015)

pdf("species_plots.pdf")
par(mfrow = c(2,2))
for(i in 1:nspec){
  sp2 <- sp_cof[i,]
  plot(1,1, ylim = c(50,250), xlim = c(1815, 2015), type = "n", xlab = 
         "rate of increase (days per year)", ylab = "initial lay date", bty = "n", xaxt = "n", yaxt = "n",
       main = spor[i])
  
  axis(2, at = seq(50, 250, 50), labels=F, tck=-.035)
  axis(2, at = seq(50, 250, 25), labels=F, tck=-.025)
  axis(2, at = seq(50, 250, 25/2), labels=F, tck=-.015)
  mtext(text = seq(50, 250, 50), 2, line = .75, at = seq(50, 250, 50), las = 1, cex = 1)
  
  axis(1, at = seq(1815, 2015, 10), labels=F, tck=-.035)
  axis(1, at = seq(1815, 2015, 5), labels=F, tck=-.025)
  axis(1, at = seq(1815, 2015, 2.5), labels=F, tck=-.015)
  mtext(text = seq(1815, 2015, 30), 1, line = .75, at = seq(1815, 2015, 30), las = 1, cex = 0.5)
  fh <- sp2$ih + sp2$sh * 0:200
  fm <- sp2$im + sp2$sm * 0:200
  fl <- sp2$il + sp2$sl * 0:200
  lines(fm ~ yv, lwd = 2)
  lines(fh~ yv)
  lines(fl ~ yv)
  points(fdt[grep(spor[i], fdt$species),2]~
           fdt[grep(spor[i], fdt$species),3])
  
}
dev.off()


gs <- mm[,grep("^G\\[", colnames(mm))]
pdf("group_lvl_effects.pdf")
par(mfrow = c(2, 3))

my_95 <- function(x){x[which(x>quantile(x, probs = c(0.025)) & x<quantile(x, probs = 0.975))]}
hist(my_95(gs[,1]), xlab = c("estimate"), main = "resident mean")
hist(my_95(gs[,3]), xlab = c("difference from resident"), main = "short distance intercept")
hist(my_95(gs[,5]), xlab = c("difference from resident"), main = "long distance intercept")

hist(my_95(gs[,2]), xlab = "estimate", main = "resident slope")
hist(my_95(gs[,4]), xlab = "difference from resident", main = "short distance slope")
hist(my_95(gs[,6]), xlab = "difference from resident", main = "long distance slope")
dev.off()

pdf("correlation.pdf")
hist(my_95(mm[,grep("rho.B\\[1,2", colnames(mm))]), xlab = "correlation coefficient",
     main = "correlation between intercept and slope")
dev.off()