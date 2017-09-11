# met_distributions_by_emissions.R

# The goal of this script is to determine if metereology is different on days
# with varying levels of emissions. 

# TODO: All "_NA_" in met path strings need to be changed to "_west_", or better
# TODO: yet, a dynamic region argument. 
# TODO: This currently only gives results for months 6 - 9, consider making season an 
# TODO: argument. 

library(stringr)
library(ncdf4)
library(fields)
library(maps)


# Select which fire emission inventory to investigate
emissions   <- 'GFED4s' # GFED4s | HMS
emissionVAR <- 'C' # This is the emission variable name read out of the emission nc file
region    <- '_west_' 
startYear <- 2003
endYear   <- 2016

# Figure out which hard drive to read data from 
# TODO: Make emissions file creation more dynamic 
# TODO: Make work for GFED again and run an update to make sure latest version shown
if ( str_detect( getwd(), "Google") ){

  eraDir <- '/Volumes/Brey_external/era_interim_nc_daily_merged/'
  emissionsDir <- paste0("/Volumes/Brey_external/", emissions,"/")
  emissionFile <- paste0(emissionsDir, 'HMS_ecmwf_SPDH', region, 
                         startYear, '_', endYear, '.nc')
  
} else {

  eraDir <- "/barnes-scratch/sbrey/era_interim_nc_daily_merged/"
  emissionsDir <- paste0("/barnes-scratch/sbrey/", emissions, "/")
  emissionFile <- paste0(emissionsDir, 'HMS_ecmwf_SPDH', region, 
                         startYear, '_', endYear, '.nc')
  
}

figureDir <- paste0('../Figures/', emissions, '_era_interm_analysis/')

################################################################################
# Read in the nc data of interest
# TODO: These files will eventually contain "_west_" instead of NA
# TODO: Time subset the met files to match the emissions files! 
################################################################################
v_nc <- nc_open(paste0(eraDir, 'v10_NA_2003_2016.nc'))
v   <- ncvar_get(v_nc, "v10")
longitude <- ncvar_get(v_nc, "longitude")
latitude  <- ncvar_get(v_nc, "latitude")
time      <- ncvar_get(v_nc, "time")
nc_close(v_nc)

# Handle masking the time of all of these arrays based on startYear and endYear

# First, get time into time R format 
t0 <- as.POSIXct("1900-01-01", tz="UTC")
timeCT <- t0 + time*60^2
timeLT <- as.POSIXlt(timeCT)
year <- timeLT$year + 1900
timeMask <- year >= startYear & year <= endYear

# Apply this new mask across the board to the time variables created just above
timeCT <- timeCT[timeMask]
timeLT <- timeLT[timeMask]
year   <- year[timeMask]
time   <- time[timeMask]
 
# v is already in the workspace, subset the time dimension.  
v <- v[,,timeMask] 
  
u_nc <- nc_open(paste0(eraDir, 'u10_NA_2003_2016.nc'))
u   <- ncvar_get(u_nc, "u10")[,, timeMask]
nc_close(u_nc)

srf_wind <- sqrt(u^2 + v^2)
rm(u,v) # save memory

tp_nc <- nc_open(paste0(eraDir, 'tp_NA_2003_2016.nc'))
tp   <- ncvar_get(tp_nc, "tp")[,, timeMask]
nc_close(tp_nc)

RH_nc <- nc_open(paste0(eraDir, 'RH_NA_2003_2016.nc'))
RH   <- ncvar_get(RH_nc, "RH")[,, timeMask]
nc_close(RH_nc)

t2m_nc <- nc_open(paste0(eraDir, 't2m_NA_2003_2016.nc'))
t2m   <- ncvar_get(t2m_nc, "t2m")[,, timeMask]
nc_close(t2m_nc)

# Read in z 500 mb level
z_nc <- nc_open(paste0(eraDir, 'z_NA_2003_2016.nc'))
z  <- ncvar_get(z_nc, "z")
z  <- z[,, 2, timeMask] # retain 500 mb level only, also cut in time
nc_close(z_nc)

################################################################################
# Load the fire emissions
# TODO: Make the 'C' more generic, like a 'emissions' variable instead since
# TODO: not always carbon 
################################################################################
C_nc <- nc_open(emissionFile)
C  <- ncvar_get(C_nc, emissionVAR)
emission_longitude <- ncvar_get(C_nc, "longitude")
emission_latitude  <- ncvar_get(C_nc, "latitude")
emission_time      <- ncvar_get(C_nc, "time")
nc_close(C_nc)

# Make sure the emission and met grids match in space and time
dy <- unique(emission_latitude - latitude)
dx <- unique(emission_longitude - longitude)
dt <- unique(emission_time - time)
if (dy != 0 | dx != 0 | dt != 0){
	stop('The fire emission grid and met grids do not match!')
}

############################################
# We only want summer for this trial work. 
############################################
month <- timeLT$mon + 1 # NOTE: LT time objects define January as month zero.  
monthMask = month >= 6 & month <= 9

# TODO: insert dynamic spatial masks here

srf_wind <- srf_wind[,,monthMask]
tp       <- tp[,,monthMask]
RH       <- RH[,,monthMask]
t2m      <- t2m[,,monthMask]
z        <- z[,,monthMask]
C        <- C[,,monthMask]
timeCT   <- timeCT[monthMask]
month    <- month[monthMask]


################################################################################
# Now that the emissions are loaded, we need create masks for varying levels
# of emitting grid box days. 
################################################################################
zeroEmissionMask <- C == 0

largeCutValue <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.95))
largestEmissionMask <- C >= largeCutValue

cut_40 <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.40))
cut_60 <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.60))
midEmissionMask <- (C >= cut_40) & (C <= cut_60)

# Get the histogram counts using the built in hist() command
noEmissionProbs <- hist(t2m[zeroEmissionMask], plot=FALSE)
midEmissionProbs <- hist(t2m[midEmissionMask], plot=FALSE)
largestEmissionProbs <- hist(t2m[largestEmissionMask],plot=FALSE)

par(las=1, lwd=3, bty="n")
plot(noEmissionProbs$mids, noEmissionProbs$density, type="l",
     ylim=c(0,0.1),
     xlab="Temperature [K]")
lines(midEmissionProbs$mids, midEmissionProbs$density, col="purple")
lines(largestEmissionProbs$mids, largestEmissionProbs$density, col="red")

# TODO: convery how different these emissions are


plotMetByEmission <- function(C, metVar, metVarName, metVarUnits, region){
  
  zeroEmissionMask <- C == 0
  
  largeCutValue <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.95))
  largestEmissionMask <- C >= largeCutValue
  
  cut_40 <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.40))
  cut_60 <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.60))
  midEmissionMask <- C >= cut_40 & C <= cut_60
  
  noEmissionProbs      <- hist(metVar[zeroEmissionMask], plot=FALSE)
  midEmissionProbs     <- hist(metVar[midEmissionMask], plot=FALSE)
  largestEmissionProbs <- hist(metVar[largestEmissionMask], plot=FALSE)
  
  # Get largest density from all possible values
  yMax <- max(c(noEmissionProbs$density,
                midEmissionProbs$density,
                largestEmissionProbs$density))
  
  png(filename=paste0(figureDir, metVarName, region, "distribution.png"),
      width=700, height=400
      )
  
  par(las=1, lwd=3, bty="n", mar=c(4,6,6,4))
  plot(noEmissionProbs$mids, noEmissionProbs$density, type="l",
       ylim=c(0,yMax),
       xlab=metVarUnits,
       ylab="Probability Density",
       cex.lab=2)
  lines(midEmissionProbs$mids, midEmissionProbs$density, col="purple")
  lines(largestEmissionProbs$mids, largestEmissionProbs$density, col="red")
  
  # TODO: display N for each distribution 
  legend(
    "topright",
    legend = c("no emissions", "[40,60] emission percentile", "top 5%"),
    col = c("black", "purple", "red"),
    lty=1,
    bty='n',
    cex=1.2
  )
  title(paste(metVarName, "distribution by emissions group") ,cex.main=2)
  
  dev.off()
  
  
}

plotMetByEmission(C, t2m, "t2m", "K", region)
plotMetByEmission(C, z, "z", "500 mb geopotential m**2/s**2", region)
plotMetByEmission(C, tp, "tp", "inches/day", region) # TODO: Make sure this is inches not meters
plotMetByEmission(C, RH, "RH", "RH%", region)
plotMetByEmission(C, srf_wind, "srf_wind", "surface wind speed m/s", region)



