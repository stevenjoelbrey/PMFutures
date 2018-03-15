# plot_interannual_variability.R

# ------------------------- Description ---------------------------------------
# This script will be used to plot regional interannual variability in burn area
# segregated by ecoregion and other grid attributes. 

# The main goal is to show if weather variables explain variance in interannual 
# variability and if relationships differ between ignition types. R values for
# burn are between ignition types are compared. 

library(stringr)
library(maps)
library(ncdf4)
library(fields)
library(sfsmisc)
library(lubridate) # for month()

# What region are you investigating? 
regionName <- "west"
ecoregion_select <- 11.1
year1 <- 1992
year2 <- 2015
month_select  <- 1:12 # THIS MAY BE VERY WRONG FOR HUMAN OR NON WEST REGIONS!
# Load the region lat and lon bounds
load(paste0("Data/GIS/", regionName,"_bounds.RData"))

years  <- year1:year2
nYears <- length(years)

################################################################################
# Save figures in directory specific to how the data are subset
################################################################################
experimentDir <- paste0("eco=", ecoregion_select, "_", 
                        "months=", min(month_select),"_", max(month_select), 
                        "/" )
figureDir     <- paste0("Figures/regional_interannual_variability/", experimentDir)

if(!dir.exists(figureDir)){
  dir.create(figureDir)
} else{
  print("Saving figures into an already existing directory")
}


# Load FPA-FOD data, the one with all attributes (ecoregions) assgined.
# TODO: load the one with weather and make sure they match! 
load(paste0("Data/FPA_FOD/FPA_FOD_ecmwf_", year1,"_", year2,".RData"))

# Get rid of fires that have "Missing/Undefined" start types, cleaner analysis
# with fewer assumptions. 
HasStartInfo <- FPA_FOD$STAT_CAUSE_DESCR != "Missing/Undefined"

# There are parts of regions in the SE we do not want to include in the analysis
# so the lat lon bounds need to be applied
latMask <- (FPA_FOD$LATITUDE > minLat) & (FPA_FOD$LATITUDE < maxLat)
lonMask <- (FPA_FOD$LONGITUDE > minLon) & (FPA_FOD$LONGITUDE < maxLon)

m <- latMask & lonMask & HasStartInfo
   
FPA_FOD <- FPA_FOD[m, ]

# Get fire parameters too be used in subsetting the data 
fireDate      <- FPA_FOD$DISCOVERY_DATE
fireLat       <- FPA_FOD$LATITUDE
fireLon       <- FPA_FOD$LONGITUDE
fireYear      <- lubridate::year(fireDate)
fireEcoregion <- FPA_FOD$NA_L2CODE
fireCause     <- FPA_FOD$STAT_CAUSE_DESCR
fireMonth     <- lubridate::month(fireDate)      
fireSize      <- FPA_FOD$FIRE_SIZE

# Handy arrays
humanStart <- fireCause != "Lightning" # because "Missing/Undefined" have been removed
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
  ecoRegionMask <-  fireEcoregion %in% ecoregion_select
  #stateMask <- fire_state
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
# Plot Interannual variability over time for selected months 
################################################################################
png(filename=paste0(figureDir,"FPA_FOD_interannual_variability_",
                    ecoregion_select,
                    ".png"),
    height=2000, width=4800, res=250)

par(mfrow=c(1,1), xpd=T, mar=c(4,10,4,0))

maxValue <- max(c(FPA_BA_lightning/10^6, FPA_BA_human/10^6))

# What is the correlation between the two time series? (this is the point of the
# figure)

r_annual <- cor(FPA_BA_lightning, FPA_BA_human, method="spearman")
r_annual_pretty <- round(r_annual, 2)


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
       legend=c("Lightning", "Human", paste("r =", r_annual_pretty) ),
       pch=19,
       pt.cex=2,
       col=c("gray", "orange","transparent"),
       cex=2.7,
       bty = "n"
       )

title(paste("Ecoregion", ecoregion_select), line=0, 
      cex.main=2.5)

dev.off()

