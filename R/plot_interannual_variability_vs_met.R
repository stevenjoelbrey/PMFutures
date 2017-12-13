# plot_interannual_variability_vs_met.R

# Description: 
# This script will be used to plot regional interannual variability in burn area
# segregated by ecoregion and other grid attributes. 

# goal, show if weather variables explain variance in interannual variability 
# and if relationships differ between different emission inventories. 

# This is done using ecmwf reanalysis fields and FPA-FOD fire occurance data 

library(stringr)
library(maps)
library(ncdf4)
library(fields)
library(sfsmisc)

year1 <- 1992
year2 <- 2013
ecoregion_select <- 9.4 
month_select     <- 5:10

years  <- year1:year2
nYears <- length(years)

# Load FPA-FOD data, the one with all attributes (ecoregions) assgined.
# TODO: load the one with weather and make sure they match! 
load(paste0("Data/FPA_FOD/FPA_FOD_", year1,"_", year2,".RData"))
fireDate      <- FPA_FOD$DISCOVERY_DATE
fireLat       <- FPA_FOD$LATITUDE
fireLon       <- FPA_FOD$LONGITUDE
fireYear      <- FPA_FOD$FIRE_YEAR
fireEcoregion <- FPA_FOD$NA_L2CODE
fireCause     <- FPA_FOD$STAT_CAUSE_DESCR
fireMonth     <- FPA_FOD$START_MONTH      
fireSize      <- FPA_FOD$FIRE_SIZE

# Handy arrays
humanStart <- fireCause != "Lightning"
nFires     <- length(fireSize)

# Load ecoregion borders
load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData")
ecoregion_polygon <- SPDF[SPDF@data$NA_L2CODE==ecoregion_select,]
#################################################################
# create a fire size attribute bin that will be used for plotting 
#################################################################
sizeBins <- c(0, 10^c(1:6))
cexBins  <- c(0, 0.1, 0.4, 1, 2, 3.8)
nBins <- length(sizeBins)
sizeClass <- rep(NA, nFires)
# assign the fires to size bins
for (b in 1:(nBins-1)){
  
  m <- fireSize >= sizeBins[b] & fireSize < sizeBins[b+1]
  sizeClass[m] <- cexBins[b]
  
}

# TODO: Make this less terrible... Use scientific notation
legendText <- c("(10 - 99]", "(100 - 1,000]", 
                "(1,000 - 10,000]", "(10,000 - 100,000]", "(100,000 - 1,100,000]")

# Create annual totals for specified ecoregion and months 

FPA_BA_lightning <- rep(NA, nYears)
FPA_BA_human     <- rep(NA, nYears)

for (i in 1:nYears){
  
  yearMask <- years[i] == fireYear
  monthMask <- fireMonth %in% month_select
  ecoRegionMask <- fireEcoregion == ecoregion_select
  m <- yearMask & monthMask & ecoRegionMask
  
  FPA_BA_human[i]     <- sum(fireSize[m & humanStart])
  FPA_BA_lightning[i] <- sum(fireSize[m & !humanStart])
  
}



pdf(file=paste0("Figures/summary/FPA_FOD_interannual_variability_mapped_",ecoregion_select,".pdf"),
    height=10, width=24)

par(mfrow=c(1,2), xpd=T, mar=c(4,10,4,0))

plot(years, FPA_BA_lightning/10^6, col="purple", 
     pch=19, bty="n", yaxt="n", xaxt="n",
     ylab="", xlab="", cex=2) 
lines(years, FPA_BA_lightning/10^6, col="purple", lty=2)

# Add axis to the plot 
eaxis(1, cex.axis=2)
eaxis(2, sub10=5, cex.axis=2)
mtext("Acres Burned (millions)", side=2, line=4, cex=2)

points(years, FPA_BA_human/10^6,pch=19, col="orange", cex=2)
lines(years, FPA_BA_human/10^6, col="orange", lty=2)


legend("topleft",
       legend=c("Lightning", "Human"),
       pch=19,
       pt.cex=2,
       col=c("purple", "orange"),
       cex=2,
       bty = "n"
       )

# Map the fires in the time series
# TODO: make sure to plot the fires in ascending size! That way we can see small
# TODO: fires next to big fires. 

par(xpd=T, mar=c(4,4,4,4))

monthMask      <- fireMonth %in% month_select
ecoRegionMask  <- fireEcoregion == ecoregion_select

humanLocations <- monthMask & ecoRegionMask & humanStart

plot(FPA_FOD$LONGITUDE[humanLocations], FPA_FOD$LATITUDE[humanLocations], 
     cex=sizeClass[humanLocations],
     col=adjustcolor("orange", 0.5), pch=1,
     bty="n", xaxt="n", yaxt="n",
     ylab="",
     xlab="")

map("state", add=T)

lightningLocations <- monthMask & ecoRegionMask & !humanStart

points(FPA_FOD$LONGITUDE[lightningLocations], FPA_FOD$LATITUDE[lightningLocations], 
     cex=sizeClass[lightningLocations],
     col=adjustcolor("purple",0.5), pch=1,
     bty="n", xaxt="n", yaxt="n",
     ylab="",
     xlab="")

# Show the ecoregion border 
plot(ecoregion_polygon, add=T, lty=1, border="lightgray")


legend("bottomleft", bty="n",
       #title="Acres",
       inset=c(0,0),
       xpd=TRUE,
       legendText,
       pt.cex=cexBins[2:(nBins-1)],
       pch=1,
       col="black",
       horiz=F)

dev.off()

# # NOTE: This only works since all fires are treated as points, and only exist
# # NOTE: in the western hemisphere. 
# fireLonAdjusted <- fireLon + 360
# 
# 
# # Get temperature
# ncDir <- "/Volumes/Brey_external/era_interim_nc_daily_merged/"
# 
# print("Loading lots of nc data")
# 
# nc_file <- paste0(ncDir,"t2m_2003_2016.nc")
# nc <- nc_open(nc_file)
# 
# # Handle ecmwf time with origin 1900-01-01 00:00:0.0
# ecmwf_hours <- ncvar_get(nc, "time")
# ecmwf_seconds <- ecmwf_hours * 60^2
# 
# # make time useful unit
# t0 <- as.POSIXct("1900-01-01 00:00:0.0", tz="UTC")
# ecmwfDate <- t0 + ecmwf_seconds
# 
# # We only want to load through 2013
# tf <- which(ecmwfDate == as.POSIXct("2013-12-31", tz="UTC"))
# 
# # Now actually load the data
# ecmwf_latitude <- ncvar_get(nc, "latitude")
# ecmwf_longitude <- ncvar_get(nc, "longitude")
# nLat <- length(ecmwf_latitude)
# nLon <- length(ecmwf_longitude)
# 
# t2m <- ncvar_get(nc, "t2m", start=c(1,1,1), count=c(nLon, nLat, tf))
# nc_close(nc)
# 
# 
# timeLT <- as.POSIXlt(ecmwfDate)[1:tf]
# mon <- timeLT$mon + 1
# yr  <- timeLT$year + 1900
# 
# 
# 
# # TODO: make function to apply to any met variable!
# latMask <- ecmwf_latitude >= minLat & ecmwf_latitude <= maxLat
# lonMask <- ecmwf_longitude >= minLon & ecmwf_longitude <= maxLon
# 
# mean_temperature <- rep(NA, nYears)
# for (i in 1:nYears){
#  
#   yearMask <- years[i] == yr
#   monthMask <- mon %in% 5:10
#   
#   # Spatial First
#   spatialSubset <- t2m[lonMask, latMask, ]
#   
#   mean_temperature[i] <- mean(spatialSubset[,, yearMask & monthMask])
#   
# }
# 
# 
# plot(mean_temperature, FPA_summer_BA, yaxt="n", bty="n", col="black", pch=19, 
#      ylab="", xlab="mean summer temperature")
# points(mean_temperature, FPA_human_summer_BA, yaxt="n", bty="n", col="orange", pch=19)
# 
# points(mean_temperature, GFED_summer_BA, col="green", pch=19)
# 
# cor(mean_temperature, FPA_summer_BA)
# cor(mean_temperature, FPA_human_summer_BA)
# 
# 
# eaxis(2)
# 





