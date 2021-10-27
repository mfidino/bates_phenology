######################
#
# Note: You need to run most of 'analysis_script.R'
#  such that you have the data_list in your global
#  environment.
#

mm <- data.table::fread("./mcmc_output/model_output.csv", 
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

# Get G parameters
G <- meds[,grep("^G", colnames(meds))]

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
  pm_list <-pmlow_list <- pmhigh_list<-pm_prob <-pm_diff <-  vector("list", length = length(nv))
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
    pm_prob[[i]] <- median(pm[,i,144]<pm[,i,1])
    pm_diff[[i]] <- quantile(pm[,i,144] - pm[,i,1],probs = c(0.025, 0.5, 0.975))
  }
  return(list(pm_list, pmlow_list, pmhigh_list, pm_prob, pm_diff))
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



# final figure for ms
windows(8,8)

pdf("./figures/figure1.pdf", height = 8, width = 8)
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
text(x = -50, y = 35, labels = "Initial lay date (Ordinal day)", 
     cex = 4, xpd = NA, srt = 90)

lines(yrest[[1]][[10]][2,] , lwd = 2)
lines(yrest[[1]][[10]][1,] , lwd = 2, lty = 2)
lines(yrest[[1]][[10]][3,] , lwd = 2, lty = 2)

lines(yrest[[2]][[10]] , lwd = 2, lty = 3)
lines(yrest[[3]][[10]] , lwd = 2, lty = 3)

points(dat$jdate ~ c(dat$year-1872), col = scales::alpha("black", 0.5),
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

lines(coest[[1]][[3]][2,] , lwd = 2)
lines(coest[[1]][[3]][1,] , lwd = 2, lty = 2)
lines(coest[[1]][[3]][3,] , lwd = 2, lty = 2)

lines(coest[[2]][[3]] , lwd = 2, lty = 3)
lines(coest[[3]][[3]] , lwd = 2, lty = 3)

points(dat$jdate ~ c(dat$year-1872), col = scales::alpha("black", 0.5),
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


lines(noest[[1]][[12]][2,] , lwd = 2)
lines(noest[[1]][[12]][1,] , lwd = 2, lty = 2)
lines(noest[[1]][[12]][3,] , lwd = 2, lty = 2)

lines(noest[[2]][[12]] , lwd = 2, lty = 3)
lines(noest[[3]][[12]] , lwd = 2, lty = 3)

points(dat$jdate ~ c(dat$year-1872), col = scales::alpha("black", 0.5),
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

points(dat$jdate ~ c(dat$year-1872), col = scales::alpha("black", 0.5),
       pch = 16, cex = 1.5)

dev.off()
###



#lines(c(yrest[[i]][2,] + abs(median(res2))), lwd = 2)
#lines(c(yrest[[i]][1,] + abs(median(res2))), lwd = 2, lty = 2)
#lines(c(yrest[[i]][3,] + abs(median(res2))), lwd = 2, lty = 2)

yer <- betas[2, grep("^B\\[\\d\\d?\\,2", colnames(betas))]
co <- betas[2, grep("^B\\[\\d\\d?\\,3", colnames(betas))]

windows(4,4)
pdf("./figures/figure2.pdf", height = 4, width = 4)
par( mar = c(5,5.5, 1.25,1))
plot(1~1, xlim = c(-0.4, 0.35), ylim = c(-0.75, 0.75) ,
     xlab = "Species specific change in mean lay date\nto year (Ordinal day)",
     ylab = "Species specific change in mean lay date\nto changes in global CO2 (Ordinal day)", bty = "n", 
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

text(x = yer[spl]-0.05, y = co[spl]-0.05, labels = "YEWA", cex = 0.75, pos = 1,
     offset = 0)

# fisp
spl <- which(levels(spcode)=="FISP")

text(x = yer[spl]-0.085, y = co[spl]-0.04, labels = "FISP", cex = 0.75, pos = 1,
     offset = 0)
#amro
spl <- which(levels(spcode)=="AMRO")
text(x = yer[spl]+0.05, y = co[spl]+0.05, labels = "AMRO", cex = 0.75, pos = 3,
     offset = 0)
#blga
spl <- which(levels(spcode)=="BLJA")
text(x = yer[spl]+0.1, y = co[spl]+0.05, labels = "BLJA", cex = 0.75,
     pos = 3,
     offset = 0)
dev.off()




samp_size <- data.frame(table(data_list$species))

# 
all_change <- make_proj(nv = 1:72, mm = covs, type = 3, ress = ress, my_med)


eff_size <- matrix(unlist(all_change[[5]]), ncol = 3, byrow = TRUE)

sigs <- apply(eff_size, 1 ,sign)
sigs <- apply(sigs, 2, function(x) abs(sum(x)))
sigs <- which(sigs == 3)

my_cols <- rep("gray80", 72)
my_cols[sigs] <- "gray40"
my_pch <- rep(21, 72)
my_pch[sigs] <- 23
set.seed(20)
to_move <- rnorm(72, 0, 0.125)

#axis(2, at= seq(1,9), labels = F, tck = -.025)

# This is a function to add axis ticks on log2 scale
minor.ticks.axis <- function(ax,n,t.ratio=0.5,mn,mx,...){
  
  lims <- par("usr")
  if(ax %in%c(1,3)) lims <- lims[1:2] else lims[3:4]
  
  major.ticks <- pretty(lims,n=5)
  if(missing(mn)) mn <- min(major.ticks)
  if(missing(mx)) mx <- max(major.ticks)
  
  major.ticks <- major.ticks[major.ticks >= mn & major.ticks <= mx]
  
  axis(ax,at=major.ticks,tck = -0.025, labels = FALSE)
  mtext(text = c("8", sprintf("%.0f", 2^seq(4,10, 2))), 
        1, line = 0.55, at = c(3, seq(4,10, 2)), cex = 1.2)
  
  
  axis(ax, at = seq(3,10,1), labels = FALSE, tck = -0.0125)
  
}
windows(5,5)
tiff("./figures/figure3.tiff",
     height = 5, width = 5, units = "in", res = 600, compression= "lzw")
par(mar = c(4,6,0.5,0.5))
plot(1~1, xlim = c(3,10), ylim = c(-60, 20) ,
     xlab = "", ylab = "", bty = "l", 
     xaxt = "n", yaxt = "n", type = "n")

minor.ticks.axis(1,1,mn=3, mx=10 )

axis(2, seq(-60,20, 10), labels = FALSE, tck = -0.025)
axis(2, seq(-60,20, 10/2), labels = FALSE, tck = -0.025/2)

mtext(sprintf("%.0f", seq(-60,20, 10)), 2, line = 0.7, las = 1,
      at = seq(-60,20, 10), cex = 1.2)

for(i in 1:nrow(eff_size)){
  lines(x = c(log(rep(samp_size[i, 2], 2), 2) + to_move[i]),
        y = eff_size[i, -2],
        col = scales::alpha("black", 0.5),
        lwd = 1.5
  )
}

points(eff_size[,2]~ c(log(samp_size[,2], base = 2) + to_move),
       pch = my_pch,
       bg = my_cols,
       cex = 1.2)
abline(h = 0, lty = 2, lwd = 2)

mtext(expression("Sample size (log"[2]*" scale)"), 1, 
      line = 2.5, at = 6.5, cex = 1.5  )
mtext("Difference in initial lay date\n between 1872 and 2015",
      2, at = -20, line = 2.5, cex = 1.5)

legend("bottomright", legend = c("Yes", "No"), pch = c(21,23), 
       pt.bg = c("gray80", "gray40"), title = "95% CI includes zero",
       bty = "n")

dev.off()
# do the same thing with regression slopes
betas <- t(betas)

tmp <- betas[73:144,]
set.seed(20)
to_move <- rnorm(72, 0, 0.125)

sigs <- apply(tmp, 1 ,sign)
sigs <- apply(sigs, 2, function(x) abs(sum(x)))
sigs <- which(sigs == 3)

my_cols <- rep("gray80", 72)
my_cols[sigs] <- "black"
my_pch <- rep(21, 72)
my_pch[sigs] <- 23


tiff("./figures/supp_fig1.tiff",
     height = 5, width = 5, units = "in", res = 600, compression= "lzw")
par(mar = c(4,6,0.5,0.5))
plot(1~1, xlim = c(3,10), ylim = c(-0.6, 0.6) ,
     xlab = "", ylab = "", bty = "l", 
     xaxt = "n", yaxt = "n", type = "n")

minor.ticks.axis(1,1,mn=3, mx=10 )

axis(2, seq(-0.6,0.5, 0.1), labels = FALSE, tck = -0.025)
axis(2, seq(-0.6,0.6, 0.1/2), labels = FALSE, tck = -0.025/2)

mtext(sprintf("%.1f", seq(-0.6,0.6, 0.2)), 2, line = 0.7, las = 1,
      at = seq(-0.6,0.6, 0.2), cex = 1.2)

for(i in 1:nrow(eff_size)){
  lines(x = c(log(rep(samp_size[i, 2], 2), 2) + to_move[i]),
        y = tmp[i, -2],
        col = scales::alpha("black", 0.5),
        lwd = 1.5
  )
}

points(tmp[,2]~ c(log(samp_size[,2], base = 2) + to_move),
       pch = my_pch,
       bg = my_cols,
       cex = 1.2)
abline(h = 0, lty = 2, lwd = 2)

mtext(expression("Sample size (log"[2]*" scale)"), 1, 
      line = 2.5, at = 6.5, cex = 1.3  )
mtext("Species specific change in\nmean lay date to year",
      2, at = -0, line = 2.9, cex = 1.3)

legend("bottomright", legend = c("Yes", "No"), pch = c(21,23), 
       pt.bg = c("gray80", "black"), title = "95% CI includes zero",
       bty = "n")

dev.off()


tmp <- betas[145:nrow(betas),]
set.seed(20)
to_move <- rnorm(72, 0, 0.125)

sigs <- apply(tmp, 1 ,sign)
sigs <- apply(sigs, 2, function(x) abs(sum(x)))
sigs <- which(sigs == 3)

my_cols <- rep("gray80", 72)
my_cols[sigs] <- "black"
my_pch <- rep(21, 72)
my_pch[sigs] <- 23




tiff("./figures/supp_fig2.tiff",
     height = 5, width = 5, units = "in", res = 600, compression= "lzw")
par(mar = c(4,6,0.5,0.5))
plot(1~1, xlim = c(3,10), ylim = c(-1.6, 1.6) ,
     xlab = "", ylab = "", bty = "l", 
     xaxt = "n", yaxt = "n", type = "n")

minor.ticks.axis(1,1,mn=3, mx=10 )

axis(2, seq(-1.6,1.6, 0.4), labels = FALSE, tck = -0.025)
axis(2, seq(-1.6,1.6, 0.4/2), labels = FALSE, tck = -0.025/2)

mtext(sprintf("%.1f", seq(-1.6,1.6, 0.4)), 2, line = 0.7, las = 1,
      at = seq(-1.6,1.6, 0.4), cex = 1.2)

for(i in 1:nrow(eff_size)){
  lines(x = c(log(rep(samp_size[i, 2], 2), 2) + to_move[i]),
        y = tmp[i, -2],
        col = scales::alpha("black", 0.5),
        lwd = 1.5
  )
}

points(tmp[,2]~ c(log(samp_size[,2], base = 2) + to_move),
       pch = my_pch,
       bg = my_cols,
       cex = 1.2)
abline(h = 0, lty = 2, lwd = 2)

mtext(expression("Sample size (log"[2]*" scale)"), 1, 
      line = 2.5, at = 6.5, cex = 1.3  )
mtext("Species specific change in\nmean lay date to change in global CO2",
      2, at = -0, line = 2.9, cex = 1.3)

legend("bottomright", legend = c("Yes", "No"), pch = c(21,23), 
       pt.bg = c("gray80", "black"), title = "95% CI includes zero",
       bty = "n")

dev.off()


plot(tmp[,2] ~ c(log(samp_size[,2], 10) + to_move), ylim = c(-0.6, 0.5),
     pch = 21, bg = my_cols,
     bty = "l",
     xlab = "log10 sample size",
     ylab = "Year effect")
for(i in 1:nrow(tmp)){
  lines(x = c(log(rep(samp_size[i, 2], 2), 10) + to_move[i]),
        y = tmp[i, -2],
        lwd = 1.2
  )
}
points(tmp[,2] ~ c(log(samp_size[,2], 10) + to_move),
       pch = 21, bg = my_cols,cex = 1.15)
abline(h = 0, lty = 2)


tmp <- betas[145:nrow(betas),]
set.seed(20)
to_move <- rnorm(72, 0, 0.125)

sigs <- apply(tmp, 1 ,sign)
sigs <- apply(sigs, 2, function(x) abs(sum(x)))
sigs <- which(sigs == 3)

my_cols <- rep("gray80", 72)
my_cols[sigs] <- "black"

plot(tmp[,2] ~ c(log(samp_size[,2], 10) + to_move), ylim = c(-2, 2),
     pch = 21, bg = my_cols,
     bty = "l",
     xlab = "log10 sample size",
     ylab = "CO2 effect")
for(i in 1:nrow(tmp)){
  lines(x = c(log(rep(samp_size[i, 2], 2), 10) + to_move[i]),
        y = tmp[i, -2],
        lwd = 1.2
  )
}
points(tmp[,2] ~ c(log(samp_size[,2], 10) + to_move),
       pch = 21, bg = my_cols,cex = 1.15)
abline(h = 0, lty = 2)

