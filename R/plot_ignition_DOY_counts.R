# plot_ignition_DOY_counts.R 

# ------------------------- Description ---------------------------------------
# This script plots ignition count by julaion day and shows how much burn area 
# occurs up to that DOY for the entire FPA FOD time period. 
# Code in this script was inhereted from the deprecated plot_ignition_weather.R
# script. 

# TODO: Create a dataframe that calculates the IQR of wildfire season, but use
# TODO: the dates that account for middle 50% of burn area, not fire counts!. 
# TODO: Do this by ecoregion. 

# Load the required libraries
library(maps)
library(sfsmisc)

# ----------------------- Subset arguments -------------------------------------

regionName <- "southeast"

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


# # Bounds take care of this subestting specification 
# ecoregion <- c(6.2,  9.2,  9.3, 10.1, 9.4, 13.1, 12.1, 10.2, 11.1,  7.1,  8.4,  8.3,  8.5,  9.5,  5.2,  8.1, 5.3,  3.1,  6.1, 15.4,  0.0 , 2.2,  3.2,  2.3,  8.2,  9.6)
# 

# Season
includeMonths <- c(1:12)

# Min fire size. Recal, even though fire managers view every small fire as a 
# potentially large fire, in terms on when real emissions occur, there is a 
# good reason to ignore lots of really small fires. 
minFireSize <- 0

# ----------------------- Setup Figure Dirs ------------------------------------
# THESE ARGUMENTS SET UP A UNIQUE EXPERIMENT. MAKE A DIRECTORUY TO STORE EACH
# EXISTING EXPERIMENT. 
experimentVARS <- paste0("region=", regionName,"_",
                         "months=", min(includeMonths), "-", max(includeMonths), "_",
                         "sizeMin=", minFireSize
)

experimentDir <- paste0("Figures/ignitionCompare_", experimentVARS, "/")
if(!dir.exists(experimentDir)){
  
  dir.create(experimentDir)

} else{
  
  print("The Direcory exists. Overwriting files within")
  
}

# Load the fire occurence data, "FPA_FOD" dataframe
dataFile <- paste0("Data/FPA_FOD/FPA_FOD_gridMet_1992-2015.RData")
load(dataFile)

# Get the spatial coordinates of the fires 
lat <- FPA_FOD$LATITUDE
lon <- FPA_FOD$LONGITUDE
startMonth <- FPA_FOD$START_MONTH

# spatial subset 
latMask <- (lat > minLat) & (lat < maxLat)
lonMask <- (lon > minLon) & (lon < maxLon)

# time subset 
timeMask <- FPA_FOD$START_MONTH %in% includeMonths

# Fire size mask 
sizeMask <- FPA_FOD$FIRE_SIZE >= minFireSize

# Mask out the fires that do not have a known cause (NEW)
causeMask <- FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"

# Elevation Mask
# NOTE: Because the data are not always precise for location of fire, I have at
# NOTE: at least one fire that is just off coast of San Bern CA in the ocean and 
# NOTE: was assigned an elevation of -557m (FOD_ID: 1276313). For fires less than
# NOTE: 100m I will assume something similar happenned, the data are suspect,
# NOTE: and throw them out in this analysis where elevation are used. I know
# NOTE: that Death Valley is lowest point, at about -85m.
ELEV_MASK <- FPA_FOD$elevation <= -100 # meters

# Subset the data, one mask to rule them all. 
m <- timeMask & latMask & lonMask & sizeMask & causeMask & !ELEV_MASK
FPA_FOD <- FPA_FOD[m,]

# Get the cuase so we can create comparisons 
fireCause     <- FPA_FOD$STAT_CAUSE_DESCR
lightningMask <- fireCause == "Lightning"
humanMask     <- fireCause != "Lightning" 

################################################################################
# Set DOY plot 
################################################################################

# Get the fire occurence day of year for ignitions information 
JDay       <- FPA_FOD$DISCOVERY_DOY
JDaysArray <- sort(unique(JDay)) # limited by selected months
nDays      <- length(JDaysArray)

# Do the same for burn area of fires started by Julian day. 
# TODO: Do this with dplr or aggreate 
# TODO: https://stackoverflow.com/questions/7560671/aggregate-data-in-one-column-based-on-values-in-another-column
JDayBurnArea_H <- rep(NA, nDays)
JDayBurnArea_L <- rep(NA, nDays)
for (j in 1:nDays){
  # Mask the data by each JDay and count the burn area that occurs from fires
  # that start on that JDay. 
  jDayMask        <- FPA_FOD$DISCOVERY_DOY == JDaysArray[j] 
  JDayBurnArea_H[j] <- sum(FPA_FOD$FIRE_SIZE[jDayMask & humanMask])
  JDayBurnArea_L[j] <- sum(FPA_FOD$FIRE_SIZE[jDayMask & lightningMask])
  
}

# Make arrays that show the cumulative burn area by igntion type
cumulative_BA_H <- cumsum(JDayBurnArea_H)
cumulative_BA_L <- cumsum(JDayBurnArea_L)
cumulative_BA <- cumulative_BA_H + cumulative_BA_L

# Estimate the 10th and 90th percent accumulation of burn area JDays
get_span <- function(BAC){
  percentOfBurn <- BAC/max(BAC)
  x10 <- which.min(abs(percentOfBurn-.10))
  x90 <- which.min(abs(percentOfBurn-.90))
  SPAN <- c(JDaysArray[x10], JDaysArray[x90])
  return(SPAN)
}

# Find the edges of the julain day experiment using a histogram of ignition counts
# by JDay 
png(filename=paste0(experimentDir, "JDayIgnitionCounts.png"),
    res=250, height=2000, width=3700)

par(mar=c(4,5.5,4,8))

# Histrogram of all ignitions by Julian day
h <- hist(JDay, breaks = JDaysArray, 
          main="", las=1, xaxt="n", 
          xlab="", ylab="Ignition Count",
          cex.lab=1.8,
          col="gray",
          ylim=c(0,7200))
title(paste(regionName, "ignition counts: 1992-2015"))

# Non-lightning only 
hist(JDay[!lightningMask], breaks = JDaysArray, 
     col="orange", add=T, xlim=c(0, 370), xaxt="n")
mtext("Day of year", 1, line=2.4, cex=1.5)
# TODO: Make this by burn area? That way we are choosing the area based on 
# TODO: significance? 
span <- as.numeric(quantile(JDay, c(0.10, 0.90) ) ) 
#abline(v=span, col="black", lwd=3, lty=3)

# Add a date axis
t0 <- as.POSIXct("2012-01-01", tz="UTC")
DOY <- t0 + JDaysArray*24*60^2

axis(1, at=JDaysArray, labels=format(DOY, "%m/%d"))
# TODO: Also consider axis.Date()
# TODO: burn area in background on second axis would make compelling argument..
par(new=TRUE)

# Set up the ylim, because sometimes there is more burn area from humanm ignited
yMax <- max(c(cumulative_BA_L, cumulative_BA_H))

plot(JDaysArray, cumulative_BA_L, 
     ylim=c(0, yMax),
     pch=19, col="gray", bty="n", yaxt="n", xaxt="n",
     ylab="", xlab="", type="l", lwd=5)
lines(JDaysArray, cumulative_BA_H, col="orange", pch=19, lwd=5)

eaxis(4)
mtext("Cumulative Burn Area (Acres) -", side=4, col="black", line=-2, cex=1.5)

# Show date of 10% and 90th percent of burn area met for both ignitions
abline(v=get_span(cumulative_BA_L), col="gray",lty=3, lwd=3)
abline(v=get_span(cumulative_BA_H), col="orange",lty=3, lwd=3)

# Give this a white background so vertical lines do not go through this text
legend("topleft", 
       legend=c("Lightning Ignitions", "Human Ignitions"),
       fill=c("gray", "orange"), 
       #bty="n", 
       box.col="white",
       cex=2,
       bg="white"
)

dev.off()
