# plot_interannual_variability_vs_met.R


# This script will be used to plot regional interannual variability in burn area
# segregated by ecoregion and other grid attributes. 

# goal, show if weather variables explain variance in interannual variability 
# and if relationships differ between different emission inventories. 

library(stringr)
library(maps)
library(ncdf4)
library(fields)

# Met data spatial domain 
minLat     <- 30.
maxLat     <- 49.5
minLon     <- 234.0
maxLon     <- 258.75




ncDir <- "/Volumes/Brey_external/GFED4s/"

# Load the GFED4s burn area
nc_file <- paste0(ncDir ,"GFED4.1s_ecmwf_burned_area_2003_2016.nc")
nc <- nc_open(nc_file)
GFED_time <- ncvar_get(nc, "time")

nc_close(nc)


years <- 2003:2013
nYears <- length(years)
year  <- sort(rep(years, 12)) # 12 months per year

GFED_burn_area <- array(NA, dim=c(480, 241, length(year)))
gfed_time <- rep(NA, length(year))

for (i in 1:nYears){

  nc_file <- paste0(ncDir ,"GFED4.1s_ecmwf_monthly_burned_area_",years[i],".nc")
  nc <- nc_open(nc_file)
  
  GFED_burn_area[,, years[i]==year] <- ncvar_get(nc, "monthly_burned_area")
  gfed_time[years[i]==year] <- ncvar_get(nc, "time")
  
  nc_close(nc)
  
}

sum(is.na(GFED_burn_area))

# Convert to acres 
GFED_burn_area <- GFED_burn_area / 4046.86

# Get other attributes
nc_file <- paste0(ncDir ,"GFED4.1s_ecmwf_monthly_burned_area_",years[i],".nc")
nc <- nc_open(nc_file)

gfed_lon  <- ncvar_get(nc, "longitude")
gfed_lat <- ncvar_get(nc, "latitude") 

nc_close(nc)

gfedMonth <- as.numeric(str_sub(as.character(gfed_time),5,6))



# Load the grid attributes
nc_file <- paste0("Data/grid_attributes/grid_attributes_75x75.nc")
nc <- nc_open(nc_file)
ecoregion <- round(ncvar_get(nc, "ecoregion"),2)
attribute_lat <- ncvar_get(nc, "latitude")
attribute_lon <- ncvar_get(nc, "longitude")
nc_close(nc)

# Make western U.S., ecoregion 6.2 summer sums. 
GFED_summer_BA <- rep(NA, nYears)
latMask <- gfed_lat <= 49
eco_subset <- ecoregion[, latMask]

for(i in 1:nYears){
  
  yearMask  <- years[i] == year
  monthMask <- gfedMonth %in% 5:10
  tMask <- yearMask & monthMask
  
  # Sum over the required years 
  yearSum <- apply(GFED_burn_area[,,tMask], 1:2, sum, na.rm=T)
  
  # lat subset (get rid of canada)
  yearSum <- yearSum[, latMask]
  
  # ecoregion sum
  ecoRegionMask <- eco_subset==6.2
  ecoRegionMask[is.na(ecoRegionMask)] <- FALSE
  
  ecoregionSubset <- yearSum[ecoRegionMask]
  
  # Take the sum of this very specific subset of the data 
  GFED_summer_BA[i] <- sum(ecoregionSubset)
  
}


# Load FPA-FOD data, the one with all attributes (ecoregions) assgined.
# TODO: load the one with weather and make sure they match! 
load("Data/FPA_FOD/FPA_FOD_2003_2013.RData")
fireDate      <- FPA_FOD$DISCOVERY_DATE
fireLat       <- FPA_FOD$LATITUDE
fireLon       <- FPA_FOD$LONGITUDE
fireYear     <- FPA_FOD$FIRE_YEAR
fireEcoregion <- FPA_FOD$NA_L2CODE
fireCause     <- FPA_FOD$STAT_CAUSE_DESCR
fireMonth     <- FPA_FOD$START_MONTH      
fireSize      <- FPA_FOD$FIRE_SIZE

# NOTE: This only works since all fires are treated as points, and only exist
# NOTE: in the western hemisphere. 
fireLonAdjusted <- fireLon + 360

# Create summer totals for ecoregions and then plot interannual variabilitity
# against temperature 

FPA_summer_BA <- rep(NA, nYears)
FPA_human_summer_BA <- rep(NA, nYears)

for (i in 1:nYears){
  
  yearMask <- years[i] == fireYear
  monthMask <- fireMonth %in% 5:10
  ecoRegionMask <- fireEcoregion == 6.2
  fireLat <- fireLat <= 49
  m <- yearMask & monthMask & ecoRegionMask & fireLat
  
  FPA_summer_BA[i] <- sum(fireSize[m])
  
  # Now human only 
  FPA_human_summer_BA[i] <- sum(fireSize[m & (fireCause != "Lightning")])
  
}

# Get temperature
ncDir <- "/Volumes/Brey_external/era_interim_nc_daily_merged/"

print("Loading lots of nc data")

nc_file <- paste0(ncDir,"t2m_2003_2016.nc")
nc <- nc_open(nc_file)

# Handle ecmwf time with origin 1900-01-01 00:00:0.0
ecmwf_hours <- ncvar_get(nc, "time")
ecmwf_seconds <- ecmwf_hours * 60^2

# make time useful unit
t0 <- as.POSIXct("1900-01-01 00:00:0.0", tz="UTC")
ecmwfDate <- t0 + ecmwf_seconds

# We only want to load through 2013
tf <- which(ecmwfDate == as.POSIXct("2013-12-31", tz="UTC"))

# Now actually load the data
ecmwf_latitude <- ncvar_get(nc, "latitude")
ecmwf_longitude <- ncvar_get(nc, "longitude")
nLat <- length(ecmwf_latitude)
nLon <- length(ecmwf_longitude)

t2m <- ncvar_get(nc, "t2m", start=c(1,1,1), count=c(nLon, nLat, tf))
nc_close(nc)


timeLT <- as.POSIXlt(ecmwfDate)[1:tf]
mon <- timeLT$mon + 1
yr  <- timeLT$year + 1900



# TODO: make function to apply to any met variable!
latMask <- ecmwf_latitude >= minLat & ecmwf_latitude <= maxLat
lonMask <- ecmwf_longitude >= minLon & ecmwf_longitude <= maxLon

mean_temperature <- rep(NA, nYears)
for (i in 1:nYears){
 
  
  yearMask <- years[i] == yr
  monthMask <- mon %in% 5:10
  
  # Spatial First
  spatialSubset <- t2m[lonMask, latMask, ]
  
  mean_temperature[i] <- mean(spatialSubset[,, yearMask & monthMask])
  
}


plot(mean_temperature, FPA_summer_BA, yaxt="n", bty="n", col="black", pch=19, 
     ylab="", xlab="mean summer temperature")
points(mean_temperature, FPA_human_summer_BA, yaxt="n", bty="n", col="orange", pch=19)

points(mean_temperature, GFED_summer_BA, col="green", pch=19)

cor(mean_temperature, FPA_summer_BA)
cor(mean_temperature, FPA_human_summer_BA)


eaxis(2)






