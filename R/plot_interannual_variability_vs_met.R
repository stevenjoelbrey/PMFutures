# plot_interannual_variability_vs_met.R

# ------------------------- Description ---------------------------------------
# This script will be used to plot regional interannual variability in burn area
# segregated by ecoregion and other grid attributes. 

# The main goal is to show if weather variables explain variance in interannual 
# variability and if relationships differ between ignition types. This is done 
# using ecmwf reanalysis fields and FPA-FOD fire occurance data.
# Timescales of interest include season (e.g. 5-10), monthly, and lagged monthly. 

# TODO: Download ecmwf reanalysis on monthly scale, or average what I have already

library(stringr)
library(maps)
library(ncdf4)
library(fields)
library(sfsmisc)

year1 <- 1992
year2 <- 2015 # extent of current FPA-FOD, will be 2015 soon. 
ecoregion_select <- 6.2
month_select     <- 5:10

years  <- year1:year2
nYears <- length(years)

# Load FPA-FOD data, the one with all attributes (ecoregions) assgined.
# TODO: load the one with weather and make sure they match! 
load(paste0("Data/FPA_FOD/FPA_FOD_ecmwf_",year1,"_",year2,".RData"))

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

# Load ecoregion borders, for plotting, analysis etc. 
load("Data/GIS/na_cec_eco_l2/na_cec_eco_level_2.RData") # <- SPDF
ecoregion_polygon <- SPDF[SPDF@data$NA_L2CODE==ecoregion_select,]

#################################################################
# create a fire size attribute bin that will be used for plotting 
#################################################################
sizeBins <- c(0, 10^c(1:6)) # a clean span of 6 orders of magnitude 
cexBins  <- c(0, 0.1, 0.4, 1, 2, 3.8)
nBins <- length(sizeBins)
sizeClass <- rep(NA, nFires)

# Assign the fires to size bins
for (b in 1:(nBins-1)){
  
  # Mask them 
  m <- fireSize >= sizeBins[b] & fireSize < sizeBins[b+1]
  sizeClass[m] <- cexBins[b]
  
}

# TODO: Make this less terrible... Use scientific notation
legendText <- c("(10 - 99]", "(100 - 1,000]", 
                "(1,000 - 10,000]", "(10,000 - 100,000]", 
                "(100,000 - 1,100,000]")

# Create annual totals for specified ecoregion and months 
FPA_BA_lightning <- rep(NA, nYears)
FPA_BA_human     <- rep(NA, nYears)

# Keep track of monthly totals too. ( "_m" indicates monthly variable )
nMonths <- nYears * length(month_select)
monthlyTimeArray   <- rep(as.POSIXct("1750-04-26", tz="UTC"), nMonths)
FPA_BA_lightning_m <- rep(NA, nMonths)
FPA_BA_human_m     <- rep(NA, nMonths)

monthCount <- 0 # Advance the month counter 
for (i in 1:nYears){
  
  yearMask <- years[i] == fireYear
  monthMask <- fireMonth %in% month_select
  ecoRegionMask <- fireEcoregion == ecoregion_select
  m <- yearMask & monthMask & ecoRegionMask
  
  # Sum to burn area for the mask made above! 
  FPA_BA_human[i]     <- sum(fireSize[m & humanStart])
  FPA_BA_lightning[i] <- sum(fireSize[m & !humanStart])
  
  # Keep track of monthly totals too. Go within the year
  for (m in month_select){

    # Advance the month index 
    monthCount <- monthCount + 1

    # Make a nice time array for plotting monthly data
    if(m < 10){
      mm <- paste0("0", m)
    } else{
      mm <- as.character(m)
    }
    dateString <- paste0(years[i], "-", mm, "-15")
    monthlyTimeArray[monthCount] <- as.POSIXct(dateString, tz="UTC")
    
    # Mask isloated month in year timeframe    
    oneMonthMask <- fireMonth == m
    m2 <- oneMonthMask & yearMask & ecoRegionMask
    
    FPA_BA_human_m[monthCount]     <- sum(fireSize[m2 & humanStart])
    FPA_BA_lightning_m[monthCount] <- sum(fireSize[m2 & !humanStart])
    
  }
  
}

################################################################################
# Plot Monthly time series of burn area
################################################################################
if(FALSE){
png(filename=paste0("Figures/summary/FPA_FOD_monthly_timeSeries_",
                ecoregion_select,".png"),
    height=2000, width=3000, res=250)

par(mar=c(4,8,4,4), lty=1, cex=2)

# Set the limits on the y-axis for the plot so we can create a correct sized
# blank. 
yMin <- min( c(FPA_BA_human_m, FPA_BA_lightning_m) )
yMax <- max( c(FPA_BA_human_m, FPA_BA_lightning_m) )

# Create the blank space (baby, I'll right your name)
plot(monthlyTimeArray, FPA_BA_lightning_m, col="gray", 
     bty="n", yaxt="n", xaxt="n",
     ylab="", xlab="", cex=1,
     ylim = c(yMin, yMax),
     xlim = c(min(monthlyTimeArray), max(monthlyTimeArray)),
     pch="")
lines(monthlyTimeArray, FPA_BA_lightning_m, col="gray")

# Label the x-axis
BY <- seq(1, length(monthlyTimeArray), by = 12)
axis(1, at = monthlyTimeArray[BY], labels = years)


# vertical axis
eaxis(2)
mtext("Acres Burned", side=2, line=5, cex=2)

# Human
#points(monthlyTimeArray, FPA_BA_human_m, col="orange", pch=1)
lines(monthlyTimeArray, FPA_BA_human_m, col="orange")

title(paste("Monthly Burn area in ecoregion", ecoregion_select),
      cex.main=2)


legend("topleft",
       bty="n",
       legend=c("Human", "Lightning"),
       pch= 19,
       col=c("Orange", "gray"),
       cex=1.3
       )

dev.off()
}

################################################################################
# Plot Interannual variability over time for selected months 
################################################################################
png(filename=paste0("Figures/summary/FPA_FOD_interannual_variability_mapped_",
                ecoregion_select,
                ".png"),
    height=2000, width=4800, res=250)

par(mfrow=c(1,2), xpd=T, mar=c(4,10,4,0))

maxValue <- max(c(FPA_BA_lightning/10^6, FPA_BA_human/10^6))


plot(years, FPA_BA_lightning/10^6, col="gray", 
     pch=19, bty="n", yaxt="n", xaxt="n",
     ylab="", xlab="", cex=2,
     ylim=c(0, maxValue)) 
lines(years, FPA_BA_lightning/10^6, col="gray", lty=2)

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
       col=c("gray", "orange"),
       cex=2,
       bty = "n"
       )

# Map the fires in the time series
# TODO: make sure to plot the fires in ascending size! That way we can see small
# TODO: fires next to big fires. 

par(xpd=F, mar=c(4,4,4,4))

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
     col=adjustcolor("gray",0.5), pch=1,
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


################################################################################
# Get corrosponding environmental data to make seasonal comparisons for
# specific variables.
# ERA-interim ecmwf averaging definitions can be found at the link below:
# https://software.ecmwf.int/wiki/display/CKB/ERA-Interim%3A+monthly+means
################################################################################

# Apply ecmwf gridded ecoregion mask. Want to only look at meteorology in 
# locations where we are looking at fire! 

gridAttributesFile <- "Data/grid_attributes/grid_attributes_75x75.nc"
grid_attributes    <- nc_open(gridAttributesFile)

# Spatial first
grid_latitude <- grid_attributes$dim$latitude$vals
grid_longitude <- grid_attributes$dim$longitude$vals

# Actual attributes next 
grid_ecoregion     <- round(ncvar_get(grid_attributes, "ecoregion"),2)
grid_state         <- ncvar_get(grid_attributes, "state")
nc_close(grid_attributes)

# Get rid of pesky NA values, they make it hard to make masks 
grid_ecoregion[is.na(grid_ecoregion)] <- -1
grid_ecoRegionMask <- grid_ecoregion == ecoregion_select
grid_latMask <- t(replicate(length(grid_longitude), grid_latitude < 50))

# ecoregion mask where Canada is excluded
grid_mask <- grid_ecoRegionMask & grid_latMask

# NOTE: This only works since all fires are treated as points, and only exist
# NOTE: in the western hemisphere.
fireLonAdjusted <- fireLon + 360

# Location of local nc data batch
ncDir <- "/Volumes/Brey_external/era_interim_nc_daily_merged/"

# Get temperature
# TODO: Consider making this a function, as it is getting used all over the place
nc_file <- paste0(ncDir, "t2m_", year1, "_", 2016, ".nc")
nc <- nc_open(nc_file)

# Handle ecmwf time with origin 1900-01-01 00:00:0.0
ecmwf_hours <- ncvar_get(nc, "time")
ecmwf_seconds <- ecmwf_hours * 60^2

# make time useful unit
t0 <- as.POSIXct("1900-01-01 00:00:0.0", tz="UTC")
ecmwfDate <- t0 + ecmwf_seconds

# We only want to load through 2013
tf <- which(ecmwfDate == as.POSIXct(paste0(year2, "-12-31"), tz="UTC"))

# Now actually load the data
ecmwf_latitude <- ncvar_get(nc, "latitude")
ecmwf_longitude <- ncvar_get(nc, "longitude")
nLat <- length(ecmwf_latitude)
nLon <- length(ecmwf_longitude)

# Figure out the lon and lat subsets that cover North America, but not more.
lon1 <- which(ecmwf_longitude == 180)
lon2 <- which(ecmwf_longitude == 310.5)
lonCount <- (lon2 - lon1) + 1

lat1 <- which(ecmwf_latitude == 87)
lat2 <- which(ecmwf_latitude == 15)
latCount <- (lat2 - lat1) + 1

# Get the spatial dims again using the new lonCount and latCount variables
ecmwf_latitude  <- ncvar_get(nc, "latitude", start=lat1, count=latCount)
ecmwf_longitude <- ncvar_get(nc, "longitude", start=lon1, count=lonCount)

t2m <- ncvar_get(nc, "t2m", start=c(lon1,lat1, 1), count=c(lonCount, latCount, tf))

nc_close(nc)

# To keep things as clear as possible, subset the time array so that they ALL
# match in terms of dimensions.
ecmwfDate <- ecmwfDate[1:tf]

# Get time is useful subsets and masks
timeLT <- as.POSIXlt(ecmwfDate)[1:tf]
mon <- timeLT$mon + 1
yr  <- timeLT$year + 1900

# Make sure these coordinates match!
flipper <- length(ecmwf_latitude):1
quartz(width=8, height=5)
image.plot(ecmwf_longitude, ecmwf_latitude[flipper], t2m[,flipper, 180])
# Add fire
points(fireLonAdjusted, fireLat, pch=".")
map("state", add=T)

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





