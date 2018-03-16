#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  year1 <- 1992
  year2 <- 2015
}else{
  year1 <- args[1]
  year2 <- args[2]
}

#------------------------------- Description ---------------------------------_
# This script takes the merged FPA FOD fire data created by 
# R/merge_ecoregion_year_Data.R and assigns reanalysis and elevation data. 

print("-----------------------------------------------------------------------")
print("This script assigns ecmwf reanalysis and elevation data to FPA-FOD fires.")
print(args)
print("-----------------------------------------------------------------------")

print("Spanning:")
print( c(year1, year2) )

# execute via command line, this takes about 3 hours, e.g. below  
# Rscript --vanilla R/assign_ecmwf_reanalysis_to_FPA_FOD.R 1992 2015
################################################################################

# This script is assigns ecmwf daily weather attributes to FPA-FOD, fire 
# ingition conditions. 

# Assgin:
# RH
# T
# precip, TODO: how much precip in season so far
# days_since_rain_1992_2016.nc
# dew point 

library(stringr)
library(maps)
library(ncdf4)
library(fields)
library(geosphere)

# Load FPA-FOD data, the one with all attributes (ecoregions) assgined.
FPAFOD_file <- paste0("Data/FPA_FOD/FPA_FOD_", year1, "_", year2, "_eco.RData")
load(FPAFOD_file)
nRow     <- dim(FPA_FOD)[1]
fireDate <- FPA_FOD$DISCOVERY_DATE 
fireLat  <- FPA_FOD$LATITUDE
fireLon  <- FPA_FOD$LONGITUDE

# This whole chunk of code inside this if statement is for assigning dates, 
# which for some reason were not available last time all the assign_ecoregion_to_FPAFOD.R
# codes were run. 
# Also check START_MONTH field in here...
if(sum(is.na(fireDate)) > 0){
  print("Working with ecoregion version of file where DISCOVERY_DATE all NA")
  
  df_old <- get(load(paste0("Data/FPA_FOD/FPA_FOD_", year1, "_", year2, ".RData")))
  rm(FPA_FOD)
  # reload version with ecoregion
  load(FPAFOD_file)
  
  i <- base::match(df_old$FOD_ID, FPA_FOD$FOD_ID)
  
  # Make sure these matching indicies (i) are set up correct
  if (sum(df_old$FOD_ID == FPA_FOD$FOD_ID[i]) == nRow){
    print("Dates correctly assigned")
    FPA_FOD$DISCOVERY_DATE[i] <- df_old$DISCOVERY_DATE
  } else{
    stop("No way of recovering dates for the merged ecoregion file.")
  }
  
  # Old no longer needed
  rm(df_old)
  
  # Now redefine what was above
  fireDate <- FPA_FOD$DISCOVERY_DATE # I do not know why these are NA. I can't find why. 
  fireLat  <- FPA_FOD$LATITUDE
  fireLon  <- FPA_FOD$LONGITUDE
  
}

print("unique Months in the data START_MONTH column")
print(sort(unique(FPA_FOD$START_MONTH)))

print("unique years in the data FIRE_YEAR column")
print(sort(unique(FPA_FOD$FIRE_YEAR)))

print("The total number of NA dates in these data are")
print(sum(is.na(FPA_FOD$DISCOVERY_DATE)))

# NOTE: This only works since all fires are treated as points, and only exist
# NOTE: in the western hemisphere. This adjustment is so 0:360 nc data can 
# NOTE: interact with -180:180 wildfire data
fireLonAdjusted <- fireLon + 360

# Set path to ecmwf era-interim data 
ncDir <- "/Volumes/Brey_external/era_interim_nc_daily_merged/"

print("Loading lots of large nc data")

# NOTE: All the daily_merged nc data end in the year 2016 so that number will
# NOTE: be hard coded in the file names throughout
nc_file <- paste0(ncDir,"t2m_", year1, "_", 2016, ".nc")
nc <- nc_open(nc_file)

# Handle ecmwf time with origin 1900-01-01 00:00:0.0
ecmwf_hours <- ncvar_get(nc, "time") # unites are hours from origin
ecmwf_seconds <- ecmwf_hours * 60^2 # want secinds for easy POSIXct conversions

# make time useful unit
t0 <- as.POSIXct("1900-01-01 00:00:0.0", tz="UTC")
ecmwfDate <- t0 + ecmwf_seconds

# We only want to load through year2 (usualy 2015). So figure out the index of 
# the last date in the year2 is. 
tf <- which(ecmwfDate == as.POSIXct(paste0(year2, "-12-31"), tz="UTC"))

# Now actually load the data fields from the nc file connection
ecmwf_latitude <- ncvar_get(nc, "latitude")
ecmwf_longitude <- ncvar_get(nc, "longitude")
nLat <- length(ecmwf_latitude)
nLon <- length(ecmwf_longitude)

# Figure out the lon and lat subsets. These hard coded numbers are known limits
# to the spatial extent of North America when in lat lon [0 360] space. 
lon1 <- which(ecmwf_longitude == 180)
lon2 <- which(ecmwf_longitude == 310.5)
lonCount <- (lon2 - lon1) + 1

lat1 <- which(ecmwf_latitude == 87)
lat2 <- which(ecmwf_latitude == 15)
latCount <- (lat2 - lat1) + 1

# Get the spatial dims again using the new lonCount and latCount variables 
ecmwf_latitude  <- ncvar_get(nc, "latitude", start=lat1, count=latCount)
ecmwf_longitude <- ncvar_get(nc, "longitude", start=lon1, count=lonCount)

# Finally, get temperature, but only for specified spatial and temporal limits.
# If we tried to load these nc files as a whole we would run out of memory. 
t2m <- ncvar_get(nc, "t2m", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))

nc_close(nc)

# To keep things as clear as possible, subset the time array so that they ALL
# match in terms of dimensions. 
ecmwfDate <- ecmwfDate[1:tf]
if(length(ecmwfDate) != dim(t2m)[3]){
  stop("The ecmwfDate date array and variable array do not match in length")
}else{
  print(paste("The max time in the ecmwfDate array is:", ecmwfDate[tf]))
}

################################################################################
# Plot a single time slice of t2m to make sure that you are getting the land
# area you think you are...
quartz()
f <- length(ecmwf_latitude):1
image.plot(ecmwf_longitude, ecmwf_latitude[f], t2m[,f,400])
title(paste("This map should show T over North America for", ecmwfDate[400]))

# Create array of starts and finishes to make loading nc code cleaner
starts <- c(lon1,lat1, 1)
counts <- c(lonCount, latCount, tf)

# RH% 
nc_file <- paste0(ncDir,"rh2m_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
rh2m <- ncvar_get(nc, "rh2m", start=starts, count=counts)
nc_close(nc)

# days since rain 
nc_file <- paste0(ncDir,"days_since_rain_",year1,"_", 2016,".nc")
nc <- nc_open(nc_file)
daysSinceRain <- ncvar_get(nc, "days_since_rain", start=starts, count=counts)
nc_close(nc)

# total precip
nc_file <- paste0(ncDir,"tp_",year1,"_", 2016,".nc")
nc <- nc_open(nc_file)
tp <- ncvar_get(nc, "tp", start=starts, count=counts)
nc_close(nc)

# dew point
nc_file <- paste0(ncDir,"d2m_",year1,"_", 2016,".nc")
nc <- nc_open(nc_file)
d2m <- ncvar_get(nc, "d2m", start=starts, count=counts)
nc_close(nc)

# Evaporation
nc_file <- paste0(ncDir,"e_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
evap <- ncvar_get(nc, "e", start=starts, count=counts)
nc_close(nc)

# Wind speed
nc_file <- paste0(ncDir,"u10_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
u10 <- ncvar_get(nc, "u10", start=starts, count=counts)
nc_close(nc)

nc_file <- paste0(ncDir,"v10_",year1,"_",2016,".nc")
nc <- nc_open(nc_file)
v10 <- ncvar_get(nc, "v10", start=starts, count=counts)
nc_close(nc)

windSpeed <- sqrt(u10^2 + v10^2)
# save workspace memory 
rm(v10, u10)

################################################################################
# Plot a time slice of wind. At this scale, winds are high over ocean, not so 
# much anywhere else. 
quartz()
f <- length(ecmwf_latitude):1
image.plot(ecmwf_longitude, ecmwf_latitude[f], windSpeed[,f,400])
title(paste("This map should show windspeed over North America for", ecmwfDate[400]))

# Load the elevation grid. Grid spacing of about ~1.8 km
# NOTE: The elevation grid being loaded is on a lon -180:180 grid!!! 
nc_file <- "Data/GIS/ETOPO1.nc"
nc <- nc_open(nc_file)
elev  <- ncvar_get(nc, "elevation")
elev_lat <- ncvar_get(nc, "latitude")
elev_lon <- ncvar_get(nc, "longitude")
nc_close(nc)

# # Too big to plot the whole thing. Can only afford to plot a subset
# quartz()
# f <- length(elev_lat):1
# image.plot(elev_lon[1:200], elev_lat[f][1:200], elev[,f])
# title(paste("This is the elevation to be assigned. lon = -180:180"))

# Get the nc grid attributes. 

# # Make sure the fireLon and elev_lon can be paired by ploting together. 
# # Plot the grids and data together to make sure the grids are the same
# # Sanity check the placement of everything by visualized ecmwf and fire locations
# jpeg(filename="testElevation.jpg") #png(filename="testElevation.png", height=1000, width=2000, res=200)
# #quartz(width=8, height=5)
# elev_flipper <- length(elev_lat):1
# image.plot(elev_lon, elev_lat[elev_flipper], elev[,elev_flipper])
# points(fireLon, fireLat, pch=".", col="black")
# title("The fire locations (black dots, should be over thge U.S. )")
# dev.off()

# Make sure the fireLonAdjusted and ecmwf_longitude can be paired by ploting 
# together. Ensure these imags make sense.
quartz(width=8, height=5)
flipper                <- length(ecmwf_latitude):1
ecmwf_latitude_flipped <- ecmwf_latitude[flipper]
t2m_flipped            <- t2m[,flipper, 180] # and time sliced 

image.plot(ecmwf_longitude, ecmwf_latitude_flipped, t2m_flipped)
points(fireLonAdjusted, fireLat, pch=".", col="black")
title("The fire locations (black dots, should be over the U.S. )")

# NOTE: If the previous two images show fires and data agreeing of spatial coords
# NOTE: then we are good-to-go with these long coord pairs. 

# Create arrays the length of the number of wildfires for storing fire attributes
# to later be assigned to the FPA_FOD dataframe. 
v <- rep(NA, nRow)
t2m_assigned           <- v
rh2m_assigned          <- v
daysSinceRain_assigned <- v
tp_assigned            <- v
d2m_assigned           <- v
elev_assigned          <- v

# Long term predictive measures
tm2_lastMonth  <- v
rh2m_lastMonth <- v
tp_lastMonth   <- v
evap_lastMonth <- v

# After fire start days metric. Namely I am going to look at precip and wind
windSpeed_MonthAfter <- v
tp_MonthAfter        <- v
t2m_monthAfter       <- v

# Loop through every single fire and assign its past current and future weather
nECMWFTime <- length(ecmwfDate) # This is for checking for end of all dates for month after assignments 
for (i in 1:nRow){ 
  
  # Find the fire match in space and time in the reanalysis data 
  # TODO: Woud havesine of close boxes help this calculation? 
  xi <- which.min(abs(fireLonAdjusted[i] - ecmwf_longitude))
  yi <- which.min(abs(fireLat[i] - ecmwf_latitude))
  ti <- which(fireDate[i] == ecmwfDate)
  
  if(fireDate[i] != ecmwfDate[ti]){
    stop("Something is wrong with ti (date matching)")
  }
  
  # Check the distance on the assigned grid points. The minus 180 on the fire
  # coords are because haversine only works in -180:180 coordinates
  dist_meters <- geosphere::distHaversine(c( (fireLonAdjusted[i]-180), fireLat[i]), 
                               c((ecmwf_longitude[xi]-180), ecmwf_latitude[yi]))
  dist_km <- dist_meters/1000.0
  
  distLon <- 0.75*111*cos(fireLat[i]*pi/180) # use fireLat as reference 
  distLat <- 0.75*111
  distHypotenuse <- sqrt(distLon^2 + distLat^2)
  maxDistAllowed <- distHypotenuse/2 + 0.18 # less than 1% estimated tolerance using fire lat
  
  if(dist_km > maxDistAllowed){
    
    # Plot the error
    quartz()
    f <- length(ecmwf_latitude):1
    image.plot(ecmwf_longitude, ecmwf_latitude[f], t2m[,f,ti])
    points(fireLonAdjusted[i], fireLat[i], pch=3, cex=2, col="black")
    points(ecmwf_longitude[xi], ecmwf_latitude[yi])
    title("Fire Point and lat point are too far apart.")
    # Stop the code
    stop(paste("At row:", i, "distance=", dist_km, "km. Exceeds estimated tolerance."))
    
  }
  
  # Assign each environmental variable that lives on the adjusted lon grid.
  # This should include all ecmwf variables to be assigned.
  # "i" is the wildfire index, xi, yi, ti belong to ecmwf data 
  t2m_assigned[i]           <- t2m[xi, yi, ti]
  rh2m_assigned[i]          <- rh2m[xi, yi, ti]
  daysSinceRain_assigned[i] <- daysSinceRain[xi, yi, ti]
  tp_assigned[i]            <- tp[xi, yi, ti]
  d2m_assigned[i]           <- d2m[xi, yi, ti]
  
  # Assign environmental variables with a timescale greater than say  30
  # Recal can't make these assignments near the edges of data, near beginning
  # 1992 and end of year 2015
  pastIndicies <- (ti-29):ti # makes length == 30
  if(pastIndicies[1] > 0){
    
    tm2_lastMonth[i] <- mean(t2m[xi, yi, pastIndicies])
    tp_lastMonth[i]  <- sum(tp[xi, yi, pastIndicies])
    evap_lastMonth[i]<- sum(evap[xi, yi, pastIndicies])
    rh2m_lastMonth[i]<- mean(rh2m[xi, yi, pastIndicies])
    
  }
  
  # After fire start date indicies 
  futureIndicies <- ti:(ti+29) # length == 30
  if(futureIndicies[length(futureIndicies)] <= nECMWFTime){
    
    windSpeed_MonthAfter[i] <- mean(windSpeed[xi, yi, futureIndicies])
    tp_MonthAfter[i]        <- sum(tp[xi, yi, futureIndicies])
    t2m_monthAfter[i]       <- mean(t2m[xi, yi, futureIndicies])
    
  }
  
  # Not ecmwf, but elevation too, no time dimension
  # NOTE: non adjusted lon as these elevation data live on -180:180 grid
  xxi <- which.min(abs(fireLon[i] - elev_lon))
  yyi <- which.min(abs(fireLat[i] - elev_lat))
  
  # Check the distance on the assigned grid points
  dist_meters <- distHaversine(c(fireLon[i], fireLat[i]), 
                               c(elev_lon[xxi], elev_lat[yyi]))
  dist_km <- dist_meters/1000.0
  
  # Figure out what max distance should be
  distLon <- 0.01666667*111*cos(fireLat[i]*pi/180)
  distLat <- 0.01666667*111
  distHypotenuse <- sqrt(distLon^2 + distLat^2)
  maxDistAllowed <- distHypotenuse/2 + .1 # Very small nudge, ~100m less than 1% estimated tolerance using fire lat
  
  # Different grid different tolerance (27.82987 between centers). length of 
  # hypotenuse is then 19.67869 km 
  if(dist_km > maxDistAllowed){
    stop(paste("Elevation assignent for i =", i, "is not within error tolerance."))
  }
  
  # Assign these elevation data 
  elev_assigned[i] <- elev[xxi, yyi]
  
  # Output progress to the screen
  if(i %% 1000 == 0){
    print(paste("Percent Complete: ", i/nRow*100))
  }
  
}

# Assign each of the new environmental arrays to FPA_FOD
FPA_FOD$t2m             <- t2m_assigned
FPA_FOD$rh2m            <- rh2m_assigned
FPA_FOD$days_since_rain <- daysSinceRain_assigned
FPA_FOD$tp              <- tp_assigned
FPA_FOD$d2m             <- d2m_assigned

# Past environment
FPA_FOD$tm2_lastMonth  <- tm2_lastMonth
FPA_FOD$tp_lastMonth   <- tp_lastMonth
FPA_FOD$rh2m_lastMonth <- rh2m_lastMonth
FPA_FOD$evap_lastMonth <- evap_lastMonth

# Future environment
FPA_FOD$windSpeed_MonthAfter <- windSpeed_MonthAfter
FPA_FOD$tp_MonthAfter        <- tp_MonthAfter
FPA_FOD$t2m_monthAfter       <- t2m_monthAfter

# static in time 
FPA_FOD$elevation <- elev_assigned

save(FPA_FOD, file = paste0("Data/FPA_FOD/FPA_FOD_ecmwf_",year1,"_",year2,".RData"))
# The end of assignments. Plot the locations of different elevation bins of 
# FPA FOD fires. Make sure high is high low is low. 


library(ggplot2)
library("ggmap")
library(dplyr)

pdf(file="fire_elevation_mapped.pdf", width=15, height=10)
# Make color depend on yval
df <- data.frame(lon=FPA_FOD$LONGITUDE, lat=FPA_FOD$LATITUDE, z=FPA_FOD$elevation)
ggplot(df, aes(x=lon, y=lat, color=z)) + geom_point(shape=".") + 
  theme_bw()
# TODO: Add map of US to this. 

dev.off()


