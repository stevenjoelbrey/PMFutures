#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
if(length(args)==0){
  args=c(2005, "no")
}
print(args)

# -----------------------assign_FINNFires_to_FPAFOD.R--------------------------- 
# TO RUN THIS SCRIPT:
# Rscript --vanilla R/assign_FINNFires_to_FPAFOD.R  YEAR 

year <- as.numeric(args[1])
conservePM25 = as.character(args[2])
print(paste("WORKING ON ASSIGNING FINN FIRES TO FPA FOD FOR YEAR:", year))
print(paste("conservePM25:", conservePM25))

# ------------------------- Description ---------------------------------------
# This script is used to assign FINN FIRE detections to FPA FOD fires. The goal 
# of these assignments is to get an idea of what FPA FOD fires were associated 
# with these satellite detected emissions (primaraly PM25). 
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
# TODO: Consider handling very small wildfires differently, e.g. same day
# TODO: required for time match pairing. 
# TODO: consider making these arguments.
distanceTol   <- 10  # haversinse kilometers required for match
timeBefireTol <- 1   # days before FPAFOD discovery_date FINN allowed to occur for match
timeAfterTol  <- 7   # days after FPAFOD discovery_date allowed to be assigned

# Load the FPA-FOD wildfires that will be assigned FINN stats where paired.
# TODO: Load the FPA-FOD wildfires that already have HMS assignments made? 
load("Data/FPA_FOD/FPA_FOD_gridMet_1992-2015.RData")

# Loop through wildfires, making assignments one at a time. 
yearMask  <- FPA_FOD$FIRE_YEAR == year
FPA_FOD   <- FPA_FOD[yearMask, ]
nWildfire <- dim(FPA_FOD)[1]
fire_cause <- rep("-", nWildfire)
fire_cause[FPA_FOD$STAT_CAUSE_DESCR=="Lightning"] <- "Lightning"
fire_cause[FPA_FOD$STAT_CAUSE_DESCR!="Lightning" & FPA_FOD$STAT_CAUSE_DESCR!="Missing/Undefined"] <- "Human"
fire_cause[FPA_FOD$STAT_CAUSE_DESCR=="Missing/Undefined"] <- "Missing/Undefined"


# Load the FINN fires. They come in one big batch.
load("Data/FINN/FINN_CONUS_2002-2016.RData")

# Subset FINN fires by the dates that will be allowed for pairing. 
tMask <- FINN$DATE >= as.POSIXct(paste0(year-1, "-12-01"), tz="UTC") &
         FINN$DATE <= as.POSIXct(paste0(year+1, "-01-30"), tz="UTC")
FINN <- FINN[tMask,]

# Keep track of how many FPA FOD wildfires each FINN detectection are assocaited 
# with
FINN$nFPAFODPaired       <- rep(0, dim(FINN)[1])
FINN$nFPAHumanPaired     <- rep(0, dim(FINN)[1])
FINN$nFPALightningPaired <- rep(0, dim(FINN)[1])
FINN$nMissingPaired      <- rep(0, dim(FINN)[1])
FINN$FPA_FOD_ID          <- ""
# Preserve these values because they are cleared out as we loop through wildfires
# when conservePM=yes
PM25                     <- FINN$PM25

# Create columns in FPA FOD that will store paired FINN information  
#   TOTAL_PM25: The total PM25 emissions associated with this file [units?]
#   FINN_AREA: The total FINN area associated with these wildfires [units?]
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
kmPerDeg <- 111
rad_earth <- 6378137/1000 # radius of earth in km (from distHaversine function)

for (i in 1:nWildfire){ # nWildfire
  
  # Get FPA FOD fire location and discovery date
  fireLat <- FPA_FOD$LATITUDE[i]
  fireLon <- FPA_FOD$LONGITUDE[i]
  fireDate <- FPA_FOD$DISCOVERY_DATE[i]

  # Now, make the first round of masks for FINN detections that could 
  # potentially be associated with this wildfire. Calculate difference in date 
  # in days. First make POSIXct as.numeric for seconds since epoch, then seconds 
  # to days.
  
  # Positive DT values mean FINN fire occurs that DT days AFTER FPA FOD
  DT <- ( as.numeric(FINN$DATE) - as.numeric(fireDate)  ) / secondsPerDay
  
  # The DT we want are DT positive up to timeAFterTol and DT -1 and 0. 
  tMask <- (DT <= timeAfterTol) & (DT >= -1*timeBefireTol)
  
  # What points are close enough? First subset down to 0.2 degree of lat/lon
  DX <- abs(FINN$LONGI - fireLon)
  DY <- abs(FINN$LATI - fireLat)
  distMask <- (DX <= 0.20) & (DY <= 0.20)
  distMaskClone <- distMask
  
  # Account for differences in longitude as a function of latitude but only 
  # bother with these calculations when you have condidates that pass the 
  # existing tMask and distMask
  if( sum(distMask & tMask) > 0 ){
    
    # Calculate hypotenuse distance accounting for changing km with lat accross
    # longitude.
    # NOTE: These distances will only be done where distMask is TRUE, so it
    # NOTE: is being done on a subset of the distances calculated above
    # dx <- cos(fireLat*pi/180) * kmPerDeg * DX[distMask]
    # dy <- DY[distMask] * kmPerDeg
    # dist_km <- (dx^2 + dy^2)^(0.5)

    # distHaversine
    p1 <- c(fireLon, fireLat)
    p2 <- cbind(FINN$LONGI[distMask], FINN$LATI[distMask])
    dist_km <- geosphere::distHaversine(p1, p2, r=rad_earth)
    
    # Mask projected onto already subset data 
    dist_km_mask <- (dist_km <= distanceTol) 
    
    # Update the distMask[distMask], aka where distMask was TRUE, which are the 
    # points the calculations were done for. This maps back to all FINN points 
    # using the km distance masking requirements. 
    distMask[distMask] <- dist_km_mask
    
  }
  
  # See if there are any overlap of the time and space mask. This needs to be
  # done again now that the km distance restriction has been applied (more strict)
  timeAndSpace <- distMask & tMask
  n_paired     <- sum(timeAndSpace)
  
  if(n_paired > 0){
    
    # There are FINN fires to be assigned to this wildfire
    FPA_FOD$n_FINN[i] <- n_paired
    
    # Assign other column data 
    # NOTE: nothing to stop FINN detections (rows) from being assigned to multiple FPA
    # NOTE: FOD wildfires. This could result in double, triple, etc, counting of FINN PM25. 
    FPA_FOD$TOTAL_PM25[i] <- sum(FINN$PM25[timeAndSpace]) # Could be used on more than one wildfire 
    FPA_FOD$FINN_AREA[i]  <- sum(FINN$AREA[timeAndSpace]) # ""
    FPA_FOD$minDate[i]    <- min(FINN$DATE[timeAndSpace])
    FPA_FOD$maxDate[i]    <- max(FINN$DATE[timeAndSpace])

    if (conservePM25 == "yes"){
      # remove the PM25 and AREA that have already been assigned
      FINN$PM25[timeAndSpace] <- 0
      FINN$AREA[timeAndSpace] <- 0
    }
    
    #-------------------------------------------------------------------------#
    # We also want to know what wildfire the FINN detection was paired to. 
    # TODO: Get those fires ID, that would be epic!
    FINN$nFPAFODPaired[timeAndSpace] <- FINN$nFPAFODPaired[timeAndSpace] + 1
    
    # Keep track of the FPA_FOD IDs of the paired
    FINN$FPA_FOD_ID[timeAndSpace] <- paste(FINN$FPA_FOD_ID[timeAndSpace], FPA_FOD$FOD_ID[i], sep="_") 
    
    # I want to know the ignition types of FPA wildfires that were paired to these
    # FINN fires
    if(fire_cause[i]=="Lightning"){
      FINN$nFPALightningPaired[timeAndSpace] <- FINN$nFPALightningPaired[timeAndSpace] + 1
    } else if(fire_cause[i]=="Human"){
      FINN$nFPAHumanPaired[timeAndSpace] <- FINN$nFPAHumanPaired[timeAndSpace] + 1
    } else{
      FINN$nMissingPaired[timeAndSpace] <- FINN$nMissingPaired[timeAndSpace] + 1
    }
    
    # Plotting code for testing and development
    # plot(hysplitPoints$Lon[distMaskClone], hysplitPoints$Lat[distMaskClone],
    #      pch=19, col="blue") # Plots hysplit points that are close enough
    # map("state", add=T)
    # points(fireLon, fireLat, col="red") # plot the wildfire of interest
    # points(hysplitPoints$Lon[timeAndSpace], hysplitPoints$Lat[timeAndSpace],
    #        pch=19, cex=0.3, col="red") # plots the wildfires that are close enough in space and time. 
    
  }else{
    # All other rows left NA, there are no FINN points associated with this
    # wildfire. There are also no emissions.  
    FPA_FOD$n_FINN[i]     <- 0
    FPA_FOD$TOTAL_PM25[i] <- 0
    FPA_FOD$FINN_AREA[i]  <- 0
  }
  
  # Update the user with how far along we are
  if( (i %% 1000) == 0){
    print(paste("Percent Complete:", (i/nWildfire)*100))
  }
  
} # End of wildfire loop

# When conservePM25=yes, any FINN fire paired with an FPA_FOD wildfire will be
# made 0, so we want to put those PM25 values back before saving. 
FINN$PM25 <- PM2

# Where no meaningful date was assigned, turn back to NA
FPA_FOD$minDate[FPA_FOD$minDate < as.POSIXct("1901-01-01", tz="UTC")] <- NA
FPA_FOD$maxDate[FPA_FOD$maxDate < as.POSIXct("1901-01-01", tz="UTC")] <- NA


# -------------------------- save the appended data ----------------------------
appendedFPA_FOD <- paste0("Data/FINN/FPA_FOD_with_FINN_dxdy=",
                          distanceTol,"_", 
                          "DT=", timeBefireTol, "_", timeAfterTol,"_",
                          year,"_",
                          "conservePM25=", conservePM25,".RData")
save(FPA_FOD, file = appendedFPA_FOD)

# Save the FINN dataframe with the number of FOD fires each was paired to.
saveName2 <- paste0("Data/FINN/FINN_with_FPAFOD_dxdy=",
                    distanceTol,"_", 
                    "DT=", timeBefireTol, "_", timeAfterTol,"_",
                    year,"_", 
                    "conservePM25=", conservePM25, ".RData")

# We want to save these one year at a time to be appended later. 
FINNYearMask <- year(FINN$DATE)==year # This now only gets rid of 2 months
FINN <- FINN[FINNYearMask,]
save(FINN, file = saveName2)
