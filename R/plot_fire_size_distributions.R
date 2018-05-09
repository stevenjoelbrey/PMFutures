# plot_fire_size_distributions.R

# ------------------------- Description ---------------------------------------
# This script is used to create plots that show the difference in size of burn
# area by ignition type. This plot will use burn area from FPA_FOD dataframe
# for all available years. 

lightningCol <- adjustcolor("cyan4", alpha.f = 1)
humanCol     <- adjustcolor("orange", alpha.f = 1)

# ----------------- User selected data subset parameters  ----------------------

# TODO: Figure out why colorado subset does not add up to 100% for both types

regionName <- "southeast" # "west" "southeast"

if (regionName == "west"){
  
  minLat <- 31.2
  maxLat <- 49
  minLon <- -125
  maxLon <- -104
  
} else if(regionName == "southeast"){
  
  minLat <- 25 
  maxLat <- 39 
  minLon <- -91   
  maxLon <- -75.5  
  
}

library(sfsmisc) # For engineering axis 

# TODO: subset by region and ecoregion but only do that if using ggplot2

# Path is project relative
dataFile <- "Data/FPA_FOD/FPA_FOD_1992_2015_eco.RData"
load(dataFile)

# Select components needed for subsetting the FPA-FOD data
fireLon <- FPA_FOD$LONGITUDE
fireLat <- FPA_FOD$LATITUDE

# Spatial subset the data 
latMask <- (fireLat <= maxLat) & (fireLat >= minLat)
lonMask <- (fireLon <= maxLon) & (fireLon >= minLon)
spatialMask <- latMask & lonMask

# Also, mask out fires with "Unknown/Undefined" cause. 
hasStartInfo <- FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"

# Subset the FOD
FPA_FOD <- FPA_FOD[spatialMask & hasStartInfo,]

# Extract key pieces of subset FPA-FOD data 
burn_area     <- FPA_FOD$FIRE_SIZE
fire_cause    <- FPA_FOD$STAT_CAUSE_DESCR

# re-select components needed for subsetting the FPA-FOD data
lightningMask <- fire_cause == "Lightning"
fireLon <- FPA_FOD$LONGITUDE
fireLat <- FPA_FOD$LATITUDE

################################################################################
# Create burn area time series for western U.S. by ignition type and show
# frequency histogram of fire size by type underneath
################################################################################
fireYear <- FPA_FOD$FIRE_YEAR
years    <- min(fireYear):max(fireYear)
nYears   <- length(years)
humanTotal     <- rep(NA, nYears) # Array to store human-ignited burn area
lightningTotal <- rep(NA, nYears)

# Sum the annual burn area totals 
for (i in 1:nYears){
  
  y <- years[i]
  yearMask <- fireYear == y
    
  lightningTotal[i] <- sum( FPA_FOD$FIRE_SIZE[yearMask & lightningMask] )
  humanTotal[i]     <- sum( FPA_FOD$FIRE_SIZE[yearMask & !lightningMask] )
  
}

################################################################################
fileName <- paste0("Figures/burn_area_comparison/ignition_size_histogram_",
                   regionName,".png")
png(filename=fileName, height=2000, width=2000, res=250)

par(mfrow=c(2,1), mar=c(0, 6, 4, 2))

yMax <- max(c(lightningTotal,humanTotal))

plot(years, lightningTotal,
     bty="n",
     yaxt="n",
     ylab="",
     xlab="",
     pch="",
     ylim=c(0, yMax))

points(years, humanTotal, pch=19, col=humanCol)
lines(years, humanTotal, pch=19, col=humanCol, lty=2)

points(years, lightningTotal, pch=19, col=lightningCol)
lines(years, lightningTotal, pch=19, col=lightningCol, lty=2)

#axis(1, at=years)
eaxis(2)
mtext("Annual Burn Area [Acres]", side=2,  line=4)

legend("topleft",
       bty="n",
       legend=c("Human", "Lightning"),
       fill=c(humanCol, lightningCol),
       cex=1.5
)

# Make the assumption that if a fire is less than 0.01 (10^-2) that it is 0.01. 
fireSize <- FPA_FOD$FIRE_SIZE
fireSize[FPA_FOD$FIRE_SIZE < 0.01] <- 0.01

# Make custom log breaks
logBreaks <- -2:6
  
# Now show the distribution of fire size
h_human <- hist(log10(fireSize[!lightningMask]), 
                breaks=logBreaks,
                plot=FALSE)

h_lightning <- hist(log10(fireSize[lightningMask]), 
                    breaks=logBreaks,
                    plot=FALSE)

par(mar=c(4, 6, 2, 2))

plot(h_human$mids, h_human$counts, 
     col=humanCol, pch=19,
     bty="n", xaxt="n", yaxt="n",
     ylab="", xlab="")
lines(h_human$mids, h_human$counts, col=humanCol, lwd=2)

points(h_lightning$mids, h_lightning$counts, col=lightningCol, pch=19)
lines(h_lightning$mids, h_lightning$counts, col=lightningCol, lwd=2)

labels <- 10^(logBreaks) 
axis(side=1, at=logBreaks, labels=labels)
eaxis(side=2)

mtext("Fire Size [Acres]", side=1, line=2)
mtext("Frequency", side=2,  line=4)

legend("topright",
       bty="n",
       legend=c(paste("n = ", sum(!lightningMask)), paste("n = ", sum(lightningMask)) ),
       fill=c(humanCol, lightningCol),
       cex=1.5
)


dev.off()


################################################################################
# Create cumulative sum distribution figure 
################################################################################

# Sort the burn area arrays and place into perspective of total burn area
total_burn_area  <- sum(burn_area) 
burn_area_sorted <- sort(burn_area)            # small to large
cumulative_sum   <- cumsum(burn_area_sorted)   

# For plotting, get the min and max fire size for bounds 
maxFireSize <- max(burn_area)
logLim      <- log10(maxFireSize) # consider cieling 

# TODO: Make each of these lines for all fires, human started fires, and lightning
# TODO: started fires seperately.
acre_bins <- 10^(seq(0, logLim, length.out=80))
nBins <- length(acre_bins)
percent_burn_area <- rep(NA, length(acre_bins))
percent_lightning <- rep(NA, length(acre_bins)) # area implied
percent_human     <- rep(NA, length(acre_bins))

# Count what percent of burn area is less than a given acre_bin
for (i in 1:nBins){
  
  # What fires are less than acre_bin[i], starts asking for < 1 acres. 
  sizeMask <- burn_area <= acre_bins[i]
  
  # When i is equal to nBins (largest bin), make sure all values of size mask
  # are TRUE
  if(i == nBins){
    TEST <- sum(sizeMask)/length(sizeMask)
    print(paste("Burn area account for proportion is:", TEST))
  }
  
  # Fires of either kind of ignition
  percent_burn_area[i] <- (sum(burn_area[sizeMask]) / total_burn_area) * 100
  
  # Now for each ignition type, so burn_area masked with ignition type too
  percent_lightning[i] <- (sum(burn_area[sizeMask & lightningMask]) / total_burn_area) * 100
  percent_human[i]     <- (sum(burn_area[sizeMask & !lightningMask]) / total_burn_area) * 100

}


# Plot and save the figure 

fileName <- paste0("Figures/burn_area_comparison/",
                   "cumulativeBurnAreaBySize_", regionName, ".png")
png(filename=fileName,
    width=2200, height=1800, res=250)

par(mar=c(4,12,4,4))

# Set up an empty plot 
plot(acre_bins, percent_burn_area,
     type="l", bty="n", 
     log="x", lwd=4,
     xaxt='n',
     ylab="",
     xlab="",
     las=1,
     cex.axis=1.6)

mtext("% of total\n burn area",2 , line=3, cex=2, las=1)

lines(acre_bins, percent_lightning, lwd=4, col=lightningCol)
lines(acre_bins, percent_human, lwd=4, col=humanCol)

# label percent totals of each part
lightningPercent <- round(percent_lightning[nBins], 1)
humanPercent     <- round(percent_human[nBins], 1)

legend("topleft",
       bty="n",
       legend=c("Lightning-ignition",
                "Human-ignition",
                "Combined"),
       lty=1,
       col=c(lightningCol,humanCol,"black"),
       cex=2,
       lwd=4
       )

text(acre_bins[nBins], percent_lightning[nBins]+3, labels=paste(lightningPercent, "%"), col=lightningCol, cex=2, xpd=T)
text(acre_bins[nBins], percent_human[nBins]+3, labels=paste(humanPercent, "%"), col=humanCol, cex=2, xpd=T)

eaxis(1, cex.axis=1.6)
mtext("Fire Size (acres)", 1 , line=3, cex=2)


y1 <- sum(burn_area[(burn_area<=1000)]) / total_burn_area * 100


# Vertical line at 1000 acres fire size
segments(x0=1000, y0=-5, x1 = 1000, y1 = y1, col="black", lty=2, lwd=3)
# horizonal line to show the percent of total burn area these fires account for 
segments(x0=0.1, y0=y1, x1 = 1000, y1 = y1, col="black", lty=2, lwd=3)
text(50, y1+4, labels=paste(round(y1,2), "%"), col="black", cex=2, xpd=T)

# y2 <- sum(burn_area[(burn_area<=10^5)]) / total_burn_area * 100
# # Vertical line at 10^5 acres fire size
# segments(x0=10^5, y0=-5, x1 = 10^5, y1 = y2, col="black", lty=2, lwd=3)
# # horizonal line to show the percent of total burn area these fires account for 
# segments(x0=0.1, y0=y2, x1 = 10^5, y1 = y2, col="black", lty=2, lwd=3)

dev.off()

