# assign_ignition_weather.R 

# This script is designed to assignb ecmwf reanalysis weather and other
# environmental factors to fires documented in the FPA-FOD database. 

# Load the required libraries
library(maps)
library(sfsmisc)

load("Data/FPA_FOD/FPA_FOD_ecmwf_2003_2013.RData")
FPA_FOD$t2m <- FPA_FOD$t2m - 273.15 # C is much nicer


##################
# Subset arguments
##################

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


lat <- FPA_FOD$LATITUDE
lon <- FPA_FOD$LONGITUDE
startMonth <- FPA_FOD$START_MONTH

#################
# Subset the data 
#################

# spatial subset 
latMask <- lat > minLat & lat < maxLat
lonMask <- lon > minLon & lon < maxLon

# time subset 
timeMask <- FPA_FOD$START_MONTH %in% includeMonths

# Ecoregion ecoregion
ecoregionMask <- FPA_FOD$NA_L2CODE %in% ecoregion

# Fire size mask 
sizeMask <- FPA_FOD$FIRE_SIZE >= minFireSize

# Subset the data, one mask to rule them all. 
m <- timeMask & latMask & lonMask & ecoregionMask & sizeMask
FPA_FOD_subset <- FPA_FOD[m,]

fireCause <- FPA_FOD_subset$STAT_CAUSE_DESCR
humanMask <- fireCause != "Lightning"

# t2m
# rh2m
# days_since_rain
# tp 
# mean temperature last 10 days
# total rain last 10 days. 




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



