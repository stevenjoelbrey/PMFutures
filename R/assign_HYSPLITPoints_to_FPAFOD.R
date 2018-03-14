#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args=2007
}

# TODO: Assign Unique FPA ID to the HMS hysplit points that are assigned to 
# TODO: wildfires. This tells us how many are associated with all lands wildfires
# TODO: and allows easy sanity checking of pairs via plotting. 

# assign_HYSPLITPoints_to_FPAFOD.R 

# TO RUN THIS SCRIPT
# Rscript --vanilla R/assign_HYSPLITPoints_to_FPAFOD.R 2007

year <- args[1]
print(paste("WORKING ON ASSIGNING HMS TO FPA FOD YEAR:", year))

# ------------------------- Description ---------------------------------------
# This script is used to assigh HMS HYPLIT points to FPA FOD fires. The goal of
# these assignments is to get an idea of what FPA FOD fires were associated with
# air quality forecasts. Later, we will learn about the variability associated 
# with these assignments. Make sure to indicate where no comparison is available. 
#
# Description of FPA FOD is from Short 2015
# https://www.fs.usda.gov/rds/archive/Product/RDS-2013-0009
# 
# Description of the HMS HYPSLIT points can be found in Brey et al. 2017
# https://www.atmos-chem-phys.net/18/1745/2018/acp-18-1745-2018-discussion.html
#
#
# The FPA FOD data used by this script was processes by:
# R/assign_ecoregion_to_FPAFOD.R
# R/assign_ecmwf_reanalysis_to_FPA_FOD.R
# R/assign_metGrid_to_FPAFOD.R 

library(geosphere)
library(maps)
library(lubridate)

# -------------------------- Set tolerance values  -----------------------------
distanceTol <- 10  # haversinse kilometers required for HP to be assigned to FPAFOD
timeBefireTol <- 1 # days before FPAFOD discovery_date HP allowed to occur for match
timeAfterTol <- 7  # days after FPAFOD discovery_date HP allowed to be assigned


# Load the HMS Hysplit Points dataframe. To be sure nothing strange happens in
# assignments, allow assignments to be made using only this years data. This 
# will impact assignments at the end and beginning of year. 
load("Data/HMS/hysplitPoints_2007_2017.RData")
yearMask <- lubridate::year(hysplitPoints$DATE) == year
hysplitPoints <- hysplitPoints[yearMask,]

# Keep track of how many wildfires this is assocaited with
hysplitPoints$nFPAFODPaired <- rep(0, dim(hysplitPoints)[1])

# HEADERS <- c("Fire1", "Fire2", "Fire3", "Fire4", "Fire5", "Fire6")
# df <-  data.frame(matrix(ncol = length(HEADERS), nrow = dim(hysplitPoints)[1]))
# names(df) <- HEADERS
# hysplitPoints <- cbind(hysplitPoints, df)

# Load the FPA-FOD wildfires that will be assigned Hysplit Point Stats
load("Data/FPA_FOD/FPA_FOD_gridMet_1992-2015.RData")

# Loop through wildfires, making assignments one at a time. 
yearMask  <- FPA_FOD$FIRE_YEAR == year
FPA_FOD   <- FPA_FOD[yearMask, ]
nWildfire <- dim(FPA_FOD)[1]

# Create columns in FPA FOD that will store Hysplit point information associated 
# with that wildfire (row). 
#   n_HP: The number of hysplit points associated with the wildfire
#   totalDurationHours: Sum of the durationHours of associated wildfires
#   minDate: minDate of Hysplit point associated with wildfire
#   maxDate: max date of hysplit point associated with wildfire 
blank <- rep(NA, nWildfire)
blankTime <- rep(as.POSIXct("1900-01-01", tz="UTC"), nWildfire)
FPA_FOD$n_HP <- blank
FPA_FOD$totalDurationHours <- blank
FPA_FOD$minDate <- blankTime
FPA_FOD$maxDate <- blankTime

for (i in 1:nWildfire){ # nWildfire
  
  # Get fire location
  fireLat <- FPA_FOD$LATITUDE[i]
  fireLon <- FPA_FOD$LONGITUDE[i]
  # Get fire date
  fireDate <- FPA_FOD$DISCOVERY_DATE[i]
  
  # Now, make the first round of masks for hysplit points that could potentially
  # be associated with this wildfire. Calculate difference in date in days.
  # First make POSIXct as.numeric for seconds since epoch, then seconds to days.
  
  # Positive DT values mean hysplit point occurs that nmany days AFTER FPA FOD
  DT <- ( as.numeric(hysplitPoints$DATE - as.numeric(fireDate) )) / (24*60*60)
  tMask <- (DT <= timeAfterTol) & (DT >= -1*timeBefireTol)
  
  # What points are close enough? First subset down to a degree
  DX <- abs(hysplitPoints$Lon - fireLon)
  DY <- abs(hysplitPoints$Lat - fireLat)
  distMask <- (DX <= 0.2) & (DY <= 0.2)
  distMaskClone <- distMask
  
  # haversine the points that passed the 0.3 deg restriction)
  # Only make haversine calculation when needed
  if( sum(distMask & tMask) > 0 ){
    
    # print(paste("index:", i, "has data"))
    
    # Calculate hypotenuse distance accounting for changing km with lat accross
    # longitude
    dx <- cos(fireLat*pi/180) * 111 * DX[distMask]
    dy <- DY[distMask] * 111
    dist_km <- (dx^2 + dy^2)^(0.5)
    
    # Haversine is slow so using the trig shown above instead. 
    # p1 <- c(fireLon, fireLat)
    # p2 <- cbind(hysplitPoints$Lon[distMask], hysplitPoints$Lat[distMask])
    # dist_km <- geosphere::distHaversine(p1, p2, r=6378137)/1000

    # Mask projected onto already subset data 
    dist_km_mask <- dist_km <= distanceTol 
    
    # Update the distMask, this maps back to all hysplit points using the km
    # distance masking requirements 
    distMask[distMask] <- dist_km_mask
    
  }
  
  # See if there are any overlap of the time and space mask. This needs to be
  # done again now that the km distance restriction has been applied
  timeAndSpace <- distMask & tMask
  n_HP <- sum(timeAndSpace)
  
  
  if(n_HP > 0){
    
    # print(paste("index:", i, "has data that meets dist cuttoff"))
    
    # There are hysplit points to be assigned to this wildfire
    FPA_FOD$n_HP[i] <- n_HP
    
    # Assign other column data 
    FPA_FOD$totalDurationHours[i] <- sum(hysplitPoints$DurHours[timeAndSpace])
    FPA_FOD$minDate[i] <- min(hysplitPoints$DATE[timeAndSpace])
    FPA_FOD$maxDate[i] <- max(hysplitPoints$DATE[timeAndSpace])
    
    # We also want to know what wildfire the HYSPLIT point was paired to. 
    # TODO: Get those fires ID, that would be epic!
    hysplitPoints$nFPAFODPaired[timeAndSpace] <- hysplitPoints$nFPAFODPaired[timeAndSpace] + 1
    
    # plot(hysplitPoints$Lon[distMaskClone], hysplitPoints$Lat[distMaskClone],
    #      pch=19, col="blue") # Plots hysplit points that are close enough
    # map("state", add=T)
    # points(fireLon, fireLat, col="red") # plot the wildfire of interest
    # points(hysplitPoints$Lon[timeAndSpace], hysplitPoints$Lat[timeAndSpace],
    #        pch=19, cex=0.3, col="red") # plots the wildfires that are close enough in space and time. 
    # 
    
  }else{
    # All other rows left NA, there are no hysplit points associated with this
    # wildfire 
    FPA_FOD$n_HP[i] <- 0
  }
  
  # Update the user with how far along we are
  if( (i %% 1000) == 0){
    print(paste("Percent Complete:", (i/nWildfire)*100))
  }
  
} # End of wildfire loop

# Where no meaningful date was assigned, turn back to NA
FPA_FOD$minDate[FPA_FOD$minDate < as.POSIXct("1901-01-01", tz="UTC")] <- NA
FPA_FOD$maxDate[FPA_FOD$maxDate < as.POSIXct("1901-01-01", tz="UTC")] <- NA

# save the appended data 
saveName1 <- paste0("Data/HMS/FPA_FOD_with_HP_", year, ".RData")
save(FPA_FOD, file = saveName1)

# Save the hysplitpoints dataframe with the number of FOD fires it was paired to
saveName2 <- paste0("Data/HMS/hysplitPoints_with_FPAFOD_", year, ".RData")
save(hysplitPoints, file = saveName2)
