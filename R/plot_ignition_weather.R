# assign_ignition_weather.R 

# This script is designed to assignb ecmwf reanalysis weather and other
# environmental factors to fires documented in the FPA-FOD database. 

# Load the required libraries
library(maps)

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
minFireSize <- 0 

# Season
includeMonths <- c(5:10)

# ecoregion 
ecoregion <- 6.2
# Leave this line uncommented if you want to retain all ecoregions
ecoregion <- unique(FPA_FOD$NA_L2CODE)


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

# Seperate ignition location variables by the humanMask vs. !humanMask
t2m_human     <- density(FPA_FOD_subset$t2m[humanMask])
t2m_lightning <- density(FPA_FOD_subset$t2m[!humanMask])

rh2m_human     <- density(FPA_FOD_subset$rh2m[humanMask]) 
rh2m_lightning <- density(FPA_FOD_subset$rh2m[!humanMask]) 

dsr_human     <- density(FPA_FOD_subset$days_since_rain[humanMask]) 
dsr_lightning <- density(FPA_FOD_subset$days_since_rain[!humanMask]) 

tp_human     <- density(FPA_FOD_subset$tp[humanMask]) 
tp_lightning <- density(FPA_FOD_subset$tp[!humanMask]) 

elev_human <- density(FPA_FOD_subset$elevation[humanMask]) 
elev_lightning <- density(FPA_FOD_subset$elevation[!humanMask]) 


# t2m
# rh2m
# days_since_rain
# tp 
# mean temperature last 10 days
# total rain last 10 days. 

# Function for plotting two distributions
plot_ignition_weatther_distribution <- function(VAR="t2m", varLabel, humanMask){
  
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
       cex.lab=2)
  
  mtext("Probability Density", 2,  line=5, cex=2)
  
  lines(dLightning$x, dLightning$y, col="darkgray", lwd=4)
  
  title(paste("Mean Daily", varLabel, "at Ignition Locations"), cex.main=2)
  
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

plot_ignition_weatther_distribution(t2m_human, t2m_lightning, "Temperature [K]", humanMask)
plot_ignition_weatther_distribution(rh2m_human, rh2m_lightning, "RH %", humanMask)
plot_ignition_weatther_distribution(dsr_human, dsr_lightning, "days since rain > 0 inch", humanMask)
plot_ignition_weatther_distribution(tp_human, tp_lightning, "Total Precip [m]", humanMask)

plot_ignition_weatther_distribution(elev_human, elev_lightning, "[m]", humanMask)


# TODO: Dot size is size of fire!
map("state", xlim=c(-125,-102), ylim=c(33,50))
points(FPA_FOD_subset$LONGITUDE[humanMask], FPA_FOD_subset$LATITUDE[humanMask], 
       col="red", pch=".")

points(FPA_FOD_subset$LONGITUDE[!humanMask], FPA_FOD_subset$LATITUDE[!humanMask], 
       col="black", pch=".")
title(paste("Fires in ecoregion", ecoregion))

# TODO: histogram of the JDay that these fires are occuring on. 
minDay <- min(FPA_FOD_subset$DISCOVERY_DOY)
maxDay <- max(FPA_FOD_subset$DISCOVERY_DOY)

hist(FPA_FOD_subset$DISCOVERY_DOY[humanMask], 
     col=adjustcolor("red",0.5),
     breaks=c(minDay:maxDay), las=1, xlab="Day of Year",
     ylim=c(0,1200),
     border="transparent",
     main="Day of year ignition counts")
hist(FPA_FOD_subset$DISCOVERY_DOY[!humanMask], 
     col=adjustcolor("black", 0.5), breaks=c(minDay:maxDay), add=T,
     border="transparent")






