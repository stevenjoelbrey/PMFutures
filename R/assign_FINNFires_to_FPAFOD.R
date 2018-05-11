#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args=2002
}


# assign_FINNFires_to_FPAFOD.R 

# TO RUN THIS SCRIPT
# Rscript --vanilla R/assign_FINNFires_to_FPAFOD.R  YEAR 

year <- args[1]
print(paste("WORKING ON ASSIGNING FINN FIRES TO FPA FOD YEAR:", year))

# ------------------------- Description ---------------------------------------
# This script is used to assigh FPA FIRE detections to FPA FOD fires. The goal of
# these assignments is to get an idea of what FPA FOD fires were associated with
# these satellite detected emissions. 
#
# Description of FPA FOD is from Short 2015
# https://www.fs.usda.gov/rds/archive/Product/RDS-2013-0009
# 
# Description of the FINN can be found in:
# https://www2.acom.ucar.edu/modeling/finn-fire-inventory-ncar
#
#
# The FPA FOD data used by this script was processes by:
# R/assign_ecoregion_to_FPAFOD.R
# R/assign_ecmwf_reanalysis_to_FPA_FOD.R
# R/assign_metGrid_to_FPAFOD.R 
#
# After this script is run for each year, combine the yearly data, adding 
# desired attributes, using R/merge_paired_fires.R 

library(geosphere)
library(maps)
library(lubridate) # for month() year() etc. 

# -------------------------- Set tolerance values  -----------------------------
distanceTol   <- 10  # haversinse kilometers required for match
timeBefireTol <- 1   # days before FPAFOD discovery_date FINN allowed to occur for match
timeAfterTol  <- 7   # days after FPAFOD discovery_date allowed to be assigned

# Load the FPA-FOD wildfires that will be assigned FINN Stats.
# TODO: Load the FPA-FOD wildfires that already have HMS assignments made? 
load("Data/HMS/FPA_FOD_gridMet_1992-2015.RData")

# Loop through wildfires, making assignments one at a time. 
yearMask  <- FPA_FOD$FIRE_YEAR == year
FPA_FOD   <- FPA_FOD[yearMask, ]
nWildfire <- dim(FPA_FOD)[1]

# Load the FINN fires. 
load("Data/FINN/FINN_CONUS_2002-2016.RData")

# Keep track of how many wildfires this is assocaited with
FINN$nFPAFODPaired <- rep(0, dim(FINN)[1])

# Create columns in FPA FOD that will store Hysplit point information associated 
# with that wildfire (row). 
#   TOTAL_PM25: The total PM25 emissions associated with this file
#   FINN_AREA: The total FINN area associated with these wildfires 
#   minDate: minDate of FINN fires associated with wildfire
#   maxDate: max date of FINN fires associated with wildfire 
blank <- rep(NA, nWildfire)
blankTime <- rep(as.POSIXct("1900-01-01", tz="UTC"), nWildfire)
FPA_FOD$n_FINN     <- blank
FPA_FOD$TOTAL_PM25 <- blank
FPA_FOD$FINN_AREA  <- blank
FPA_FOD$minDate    <- blankTime
FPA_FOD$maxDate    <- blankTime

# Handy conversion for time comparisons 
secondsPerDay <- (24*60*60)

for (i in 1:nWildfire){ # nWildfire
  
  # Get FPA FOD fire location
  fireLat <- FPA_FOD$LATITUDE[i]
  fireLon <- FPA_FOD$LONGITUDE[i]
  # Get fire date
  fireDate <- FPA_FOD$DISCOVERY_DATE[i]
  
  # Now, make the first round of masks for hysplit points that could potentially
  # be associated with this wildfire. Calculate difference in date in days.
  # First make POSIXct as.numeric for seconds since epoch, then seconds to days.
  
  # Positive DT values mean FINN fire occurs that DT days AFTER FPA FOD
  DT <- ( as.numeric(FINN$DATE) - as.numeric(fireDate)  ) / secondsPerDay
  
  # The DT we want are DT positive up to timeAFterTol and DT -1 and 0. 
  tMask <- (DT <= timeAfterTol) & (DT >= -1*timeBefireTol)
  
  # What points are close enough? First subset down to a degree of lat/lon
  DX <- abs(FINN$LONGI - fireLon)
  DY <- abs(FINN$LATI - fireLat)
  distMask <- (DX <= 0.2) & (DY <= 0.2)
  distMaskClone <- distMask
  
  # haversine the points that passed the 0.2 deg restriction)
  # Only make haversine calculation when needed
  if( sum(distMask & tMask) > 0 ){
    
    # print(paste("index:", i, "has data"))
    
    # Calculate hypotenuse distance accounting for changing km with lat accross
    # longitude
    dx <- cos(fireLat*pi/180) * 111 * DX[distMask]
    dy <- DY[distMask] * 111
    dist_km <- (dx^2 + dy^2)^(0.5)

    # Mask projected onto already subset data 
    dist_km_mask <- dist_km <= distanceTol 
    
    # Update the distMask, this maps back to all hysplit points using the km
    # distance masking requirements 
    distMask[distMask] <- dist_km_mask
    
  }
  
  # See if there are any overlap of the time and space mask. This needs to be
  # done again now that the km distance restriction has been applied
  timeAndSpace <- distMask & tMask
  n_paired     <- sum(timeAndSpace)
  
  if(n_paired > 0){
    
    # print(paste("index:", i, "has data that meets dist cuttoff"))
    
    # There are FINN fires to be assigned to this wildfire
    FPA_FOD$n_FINN[i] <- n_paired
    
    # Assign other column data 
    FPA_FOD$TOTAL_PM25[i] <- sum(FINN$PM25[timeAndSpace])
    FPA_FOD$FINN_AREA[i]  <- sum(FINN$AREA[timeAndSpace])
    FPA_FOD$minDate[i]    <- min(FINN$DATE[timeAndSpace])
    FPA_FOD$maxDate[i]    <- max(FINN$DATE[timeAndSpace])
    
    # We also want to know what wildfire the HYSPLIT point was paired to. 
    # TODO: Get those fires ID, that would be epic!
    FINN$nFPAFODPaired[timeAndSpace] <- FINN$nFPAFODPaired[timeAndSpace] + 1
    
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
    FPA_FOD$n_FINN[i] <- 0
  }
  
  # Update the user with how far along we are
  if( (i %% 1000) == 0){
    print(paste("Percent Complete:", (i/nWildfire)*100))
  }
  
} # End of wildfire loop

# Where no meaningful date was assigned, turn back to NA
FPA_FOD$minDate[FPA_FOD$minDate < as.POSIXct("1901-01-01", tz="UTC")] <- NA
FPA_FOD$maxDate[FPA_FOD$maxDate < as.POSIXct("1901-01-01", tz="UTC")] <- NA

# TODO: include the selected tolerance settings for the save names so that they are descriptive 

# save the appended data 
appendedFPA_FOD <- paste0("Data/FINN/FPA_FOD_with_FINN_dxdy=",distanceTol,"_", 
                          "DT=", timeBefireTol, "_", timeAfterTol,"_",
                          year, ".RData")
save(FPA_FOD, file = appendedFPA_FOD)

# Save the FINN dataframe with the number of FOD fires it was paired to
saveName2 <- paste0("Data/FINN/FINN_with_FPAFOD_dxdy=",distanceTol,"_", 
                    "DT=", timeBefireTol, "_", timeAfterTol,"_",
                    year, ".RData")
FINN <- FINN[year(FINN$DATE)==year,]
save(FINN, file = saveName2)
