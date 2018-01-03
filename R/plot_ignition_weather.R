# assign_ignition_weather.R 

# ------------------------- Description ---------------------------------------
# This script will be used to determine if the weather on ignition days is 
# different for different ignition types. 

# Load the required libraries
library(maps)
library(sfsmisc)

# ----------------------- Subset arguments -------------------------------------

# Lat lon extent
minLat <- 30
maxLat <- 50
minLon <- -125
maxLon <- -100

# acres
minFireSize <- 1000 

# Season
includeMonths <- c(5:10)

# ecoregion 
ecoregion <- 6.2

# Leave this line uncommented if you want to retain all ecoregions
#ecoregion <- unique(FPA_FOD$NA_L2CODE)

dataFile <- paste0("Data/FPA_FOD/FPA_FOD_ecmwf_1992_2015.RData")
load(dataFile)

# C is much nicer than K. What am I a chemist?
FPA_FOD$t2m <- FPA_FOD$t2m - 273.15  

# Get the spatial coordinates of the fires 
lat <- FPA_FOD$LATITUDE
lon <- FPA_FOD$LONGITUDE
startMonth <- FPA_FOD$START_MONTH

# spatial subset 
latMask <- lat > minLat & lat < maxLat
lonMask <- lon > minLon & lon < maxLon

# time subset 
timeMask <- FPA_FOD$START_MONTH %in% includeMonths

# Ecoregion ecoregion
ecoregionMask <- FPA_FOD$NA_L2CODE %in% ecoregion

# # Fire size mask 
# sizeMask <- FPA_FOD$FIRE_SIZE >= minFireSize

# Subset the data, one mask to rule them all. 
m <- timeMask & latMask & lonMask & ecoregionMask 
FPA_FOD <- FPA_FOD[m,]

# Get the cuase so we can create comparisons 
fireCause     <- FPA_FOD$STAT_CAUSE_DESCR
lightningMask <- fireCause == "Lightning"

################################################################################
# Sample mean differences vs. Julain day
################################################################################

# Time difference tolerance in JDays for fires to be comparable
JDayTol       <- 10
minSampleSize <- 50

# TODO: preset the elevation bins and store differences for these bins to see
# TOOD: if differences in ignition temp occur at different elevations

JDay <- FPA_FOD$DISCOVERY_DOY
JDaysArray <- sort(unique(JDay))
nDays <- length(JDaysArray)

# Find the edges of the julain day experiment using a histogram of JDay counts
# TODO: Do this by ignition type 
hist(JDay, breaks = JDaysArray)
hist(JDay[!lightningMask], breaks = JDaysArray, col="orange", add=T)
span <- as.numeric(quantile(JDay, c(0.10, 0.90) ) ) 
abline(v=span, col="red", lwd=3)

minJDay <- span[1]
maxJDay <- span[2]

# Create expiment data storage arrays
loopJDays <- minJDay:maxJDay
JDayFireCount <- rep(NA, length(loopJDays))

# Loop through each day and perform experiment
# TODO: Should I loop over every JDay or by the tolerance? Probably Tol.
for (i in 1:length(loopJDays)){ 
  
  # Subset fires to a given julain day (or maybe range of +-10 or something?)
  JDayMask <- abs(JDay - loopJDays[i]) < JDayTol
  JDayFireCount[i] <- sum(JDayMask)
  
  # Subset the dataframe by those that are within the JDay tolerance
  df <- FPA_FOD[JDayMask,]
  lMask <- df$STAT_CAUSE_DESCR == "Lightning"
  
  COL <- rep("orange", sum(JDayMask))
  COL[lMask] <- "gray"
  
  # TODO: boxplot of these data by elevation bin! 
  # TODO: Also show the 
  plot(df$elevation, df$t2m, col=COL, pch=19, cex=0.5)
  
  # All elevations difference of means 
  meanLT <- mean(df$t2m[lMask])
  meanHT <- mean(df$t2m[!lMask])
  
  abline(h=meanLT, col="gray", lwd=3)
  abline(h=meanHT, col="orange", lwd=3)
  
  allDT <- meanLT - meanHT
  
  # Now account for elevation 
  # IDEA: Adjust temperature based on elevation using standard lapse rate?
  sampleElavations <- df$elevation
  
  elevationBins <- seq(min(sampleElavations), 
                       max(sampleElavations), 
                       by = 300)
  # Need to make sure that these bins span the range of eleavations in the data
  # TODO: Maybe better to change existing limit to maximum elevation value? 
  if(max(elevationBins) < max(sampleElavations)){
    elevationBins <- append(elevationBins, max(elevationBins) + 300)
  }
  
  abline(v=elevationBins)
  
  # Count total in each elevation bin for each ignition type. We will exclude
  # the bins where there are not enough of each type
  lightningBinCount <- hist(sampleElavations[lMask], 
                            breaks=elevationBins, plot=F)$counts
  
  humanBinCount <- hist(sampleElavations[!lMask], 
                            breaks=elevationBins, plot=F)$counts
  
  # Mask out the elevation bins that do not have enough data
  fullBins <- humanBinCount > minSampleSize & lightningBinCount > minSampleSize
  
  # Take the difference of the means in each elevation bin
  differenceOfMeans <- rep(NA, length(fullBins))
  for (m in 1:(length(elevationBins)-1) ){
    
    elevMask <- sampleElavations >= elevationBins[m] & 
      sampleElavations < elevationBins[m+1]
    
    # mean human Temperature sample  
    meanHTs <- mean(df$t2m[elevMask & !lMask])
    # mean lightning Temperature sample 
    meanLTs <- mean(df$t2m[elevMask & lMask])
    
    differenceOfMeans[m] <- meanLTs - meanHTs
    
    # is SD interesting for anything but a significance test?
    
  }
  
  mean(differenceOfMeans[fullBins])
  
}



# Function for plotting two distributions
plot_ignition_weather_distribution <- function(VAR="t2m", varLabel, humanMask){
  
  # Create the density estimates for the two types. 
  dHuman     <- density(FPA_FOD_subset[[VAR]][humanMask])
  dLightning <- density(FPA_FOD_subset[[VAR]][!humanMask])
  
  par(las=1)
  
  # Get the needed max value for the y-axis
  maxVal <- max(c(dHuman$y,  dLightning$y))
  
  plot(dHuman$x, dHuman$y, 
       type="l", col='orange',
       xlab=varLabel,
       ylab="", 
       ylim=c(0, maxVal),
       lwd=4,
       bty="n", 
       cex.axis=2,
       cex.lab=2,
       yaxt="n")
  eaxis(2, cex.axis=2)
  
  mtext("Probability Density", 2,  line=5, cex=2)
  
  lines(dLightning$x, dLightning$y, col="darkgray", lwd=4)
  
  title(paste("Mean", varLabel, "at Ignition Locations"), cex.main=1.4)
  
  legend_human <- paste("Human, n=", sum(humanMask))
  legend_lightning <- paste0("Lightning, n=", sum(!humanMask))
  
  legend("topleft",
         legend=c(legend_human, legend_lightning),
         col=c("orange", "darkgray"),
         lty=1,
         cex=1,
         lwd=4,
         bty="n")


}

quartz(width=15, height=5)

png(filename=paste0("Figures/summary/ignitionTemperature_6_2_Fires>1000.png"),
    width=3000, height=1000, res=200)
par(mfrow=c(1,3), mar=c(4,5,4,4))
plot_ignition_weather_distribution("t2m", "Temperature [K]", humanMask)

# Show the locations of the fires
# TODO: Dot size is size of fire!
map("state", xlim=c(-125,-102), ylim=c(33,50))
points(FPA_FOD_subset$LONGITUDE[humanMask], FPA_FOD_subset$LATITUDE[humanMask], 
       col="orange", pch=".")
points(FPA_FOD_subset$LONGITUDE[!humanMask], FPA_FOD_subset$LATITUDE[!humanMask], 
       col="darkgray", pch=".")
title(paste("Fires locations"), cex.main=2)


# TODO: histogram of the JDay that these fires are occuring on. 
minDay <- min(FPA_FOD_subset$DISCOVERY_DOY)
maxDay <- max(FPA_FOD_subset$DISCOVERY_DOY)

hist(FPA_FOD_subset$DISCOVERY_DOY[humanMask], 
     col=adjustcolor("orange",0.5),
     breaks=c(minDay:maxDay), las=1, xlab="Day of Year",
     ylim=c(0,30),
     border="transparent",
     main="Day of year ignition counts",
     cex.main=2)
hist(FPA_FOD_subset$DISCOVERY_DOY[!humanMask], 
     col=adjustcolor("darkgray", 0.5), breaks=c(minDay:maxDay), add=T,
     border="transparent")
dev.off()

plot_ignition_weather_distribution("tp_lastMonth", "Total Precip", humanMask)



quartz()
plot(FPA_FOD_subset$windSpeed_MonthAfter, FPA_FOD_subset$FIRE_SIZE)

#windSpeed_MonthAfter



