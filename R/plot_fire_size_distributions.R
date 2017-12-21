# plot_fire_size_distributions.R

library(sfsmisc) # For engineering axis 

# TODO: subset by region and ecoregion

# TODO: Monte Carlo method to show that the environmental conditions of igition
# TODO: types are the same when you account for elevation and month.

# TODO: Further demonstrate that distributions are the same when looking at fires
# TODO: that end up larger than 10 acres. 

# TODO: Update these data to include many more years of FPA-FOD! No reason not 
# TODO: to be looking at 1992_2013 version
load("Data/FPA_FOD/FPA_FOD_ecmwf_1992_2015.RData")

burn_area <- FPA_FOD$FIRE_SIZE
fire_cause <- FPA_FOD$STAT_CAUSE_DESCR
lightningMask <- fire_cause == "Lightning"

total_burn_area  <- sum(burn_area) 
burn_area_sorted <- sort(burn_area)
cumulative_sum   <- cumsum(burn_area_sorted)

maxFireSize <- max(burn_area)
logLim <- log10(maxFireSize)

# TODO: Make each of these lines for all fires, human started fires, and lightning
# TODO: started fires seperately.
acre_bins <- 10^(seq(0, logLim, length.out=80))
nBins <- length(acre_bins)
percent_burn_area <- rep(NA, length(acre_bins))
percent_lightning <- rep(NA, length(acre_bins))
percent_human     <- rep(NA, length(acre_bins))

for (i in 1:nBins){
  
  sizeMask <- burn_area <= acre_bins[i]
  
  percent_burn_area[i] <- sum(burn_area[sizeMask]) / total_burn_area * 100
  
  percent_lightning[i] <- sum(burn_area[sizeMask & lightningMask]) / total_burn_area * 100
  percent_human[i]     <- sum(burn_area[sizeMask & !lightningMask]) / total_burn_area * 100

}


png(filename="Figures/burn_area_comparison/cumulativeBurnAreaBySize_all_US.png",
    width=2200, height=1800, res=250)

par(mar=c(4,12,4,4))

plot(acre_bins, percent_burn_area,
     type="l", bty="n", 
     log="x", lwd=4,
     xaxt='n',
     ylab="",
     xlab="",
     las=1,
     cex.axis=1.6)



mtext("% of total\n burn area",2 , line=3, cex=2, las=1)

lines(acre_bins, percent_lightning, lwd=4, col="gray")
lines(acre_bins, percent_human, lwd=4, col="orange")

# label percent totals of each part
lightningPercent <- round(percent_lightning[nBins], 1)
humanPercent     <- round(percent_human[nBins], 1)


text(acre_bins[nBins], percent_lightning[nBins]+3, labels=paste(lightningPercent, "%"), col="gray", cex=2, xpd=T)
text(acre_bins[nBins], percent_human[nBins]+3, labels=paste(humanPercent, "%"), col="orange", cex=2, xpd=T)

eaxis(1, cex.axis=1.6)
mtext("Fire Size (acres)", 1 , line=3, cex=2)

dev.off()

abline(v=257280)


# TODO: Make the same plot but instead of % of total burn area on x-axis make it
# TODO: total burn area in acres.


