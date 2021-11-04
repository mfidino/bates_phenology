
samp_size <- data.frame(table(data_list$species))

# Use make_pro
longshot <- make_proj(nv = 1:72, mm = covs, type = 3, ress = ress, my_med)

eff_size <- matrix(unlist(all_change[[5]]), ncol = 3, byrow = TRUE)

sigs <- apply(eff_size, 1 ,sign)
sigs <- apply(sigs, 2, function(x) abs(sum(x)))
sigs <- which(sigs == 3)

my_cols <- rep("gray80", 72)
my_cols[sigs] <- "gray40"
my_pch <- rep(21, 72)
my_pch[sigs] <- 23
set.seed(20)

plot(eff_size[,2] ~ c(log(samp_size[,2], base = 2) + to_move), bty = 'l',
     xlab = "log 10 sample size",
     ylab = "Number of days intial lay day has changed relative to historical baseline",
     ylim = c(-60,20),
     las = 1)


plot(1~1, type = "n", xlim = c(-2,8), ylim = c(1,9), xlab = "",
     ylab = "", xaxt = "n", yaxt="n", bty = "n", bty = "n")





# sorting the species by their baseline detection probability.
mu_t <- order(colMeans(model_array[,2,]))

fancy_sp <- fancy_sp[mu_t]
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

tiff("eff_to_sample.tiff",
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
betas <- t(betas)

tmp <- betas[73:144,]
set.seed(20)
to_move <- rnorm(72, 0, 0.125)

sigs <- apply(tmp, 1 ,sign)
sigs <- apply(sigs, 2, function(x) abs(sum(x)))
sigs <- which(sigs == 3)

my_cols <- rep("gray80", 72)
my_cols[sigs] <- "black"

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
