# met_distributions_by_emissions.R

# The goal of this script is to determine if metereology is different on days
# with varying levels of emissions. 

library(stringr)
library(ncdf4)
library(fields)
library(maps)


# Figure out which hard drive to read data from 
if ( str_detect( getwd(), "Google") ){
  eraDir <- '/Volumes/Brey_external/era_interim_nc_daily_merged/'
  emissionsDir <- "/Volumes/Brey_external/GFED4s/"
} else {
  eraDir <- "/barnes-scratch/sbrey/era_interim_nc_daily_merged/"
  emissionsDir <- "/barnes-scratch/sbrey/GFED4s/"
}

figureDir <- '../Figures/GFED_era_interm_analysis/'

################################################################################
# Read in the nc data of interest
################################################################################
v_nc <- nc_open(paste0(eraDir, 'v10_NA_2003_2016.nc'))
v   <- ncvar_get(v_nc, "v10")
longitude <- ncvar_get(v_nc, "longitude")
latitude  <- ncvar_get(v_nc, "latitude")
time      <- ncvar_get(v_nc, "time")
nc_close(v_nc)

# Get time into time R format 
t0 <- as.POSIXct("1900-01-01", tz="UTC")
timeCT <- t0 + time*60^2
timeLT <- as.POSIXlt(timeCT)
  
u_nc <- nc_open(paste0(eraDir, 'u10_NA_2003_2016.nc'))
u   <- ncvar_get(u_nc, "u10")
nc_close(u_nc)

srf_wind <- sqrt(u^2 + v^2)
rm(u,v) # save memory

tp_nc <- nc_open(paste0(eraDir, 'tp_NA_2003_2016.nc'))
tp   <- ncvar_get(tp_nc, "tp")
nc_close(tp_nc)

RH_nc <- nc_open(paste0(eraDir, 'RH_NA_2003_2016.nc'))
RH   <- ncvar_get(RH_nc, "RH")
nc_close(RH_nc)

t2m_nc <- nc_open(paste0(eraDir, 't2m_NA_2003_2016.nc'))
t2m   <- ncvar_get(t2m_nc, "t2m")
nc_close(t2m_nc)

# Read in z 500 mb level
z_nc <- nc_open(paste0(eraDir, 'z_NA_2003_2016.nc'))
z  <- ncvar_get(z_nc, "z")
z  <- z[,,2,] # retain 500 mb level only 
nc_close(z_nc)

# Load the emissions! 
C_nc <- nc_open(paste0(emissionsDir, 'GFED4.1s_METGrid_C_NA_2003_2016.nc'))
C  <- ncvar_get(C_nc, "C")
nc_close(C_nc)

# We only want july for this trail work. 
month <- timeLT$mon + 1
timeMask = month >= 6 & month <= 9

# TODO: insert dynamic spatial masks here

srf_wind <- srf_wind[,,timeMask]
tp       <- tp[,,timeMask]
RH       <- RH[,,timeMask]
t2m      <- t2m[,,timeMask]
z        <- z[,,timeMask]
C        <- C[,,timeMask]
timeCT   <- timeCT[timeMask]
month    <- month[timeMask]


################################################################################
# Now that the emissions are loaded, we need create masks for varying levels
# of emitting grid box days. 
################################################################################
zeroEmissionMask <- C == 0

largeCutValue <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.95))
largestEmissionMask <- C >= largeCutValue

cut_40 <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.40))
cut_60 <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.60))
midEmissionMask <- C >= cut_40 & C <= cut_60


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


plotMetByEmission <- function(C, metVar, metVarName, metVarUnits){
  
  zeroEmissionMask <- C == 0
  
  largeCutValue <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.95))
  largestEmissionMask <- C >= largeCutValue
  
  cut_40 <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.40))
  cut_60 <- as.numeric(quantile(C[!zeroEmissionMask], probs=0.60))
  midEmissionMask <- C >= cut_40 & C <= cut_60
  
  
  noEmissionProbs      <- hist(metVar[zeroEmissionMask],  plot=FALSE)
  midEmissionProbs     <- hist(metVar[midEmissionMask], plot=FALSE)
  largestEmissionProbs <- hist(metVar[largestEmissionMask], plot=FALSE)
  
  # Get largest density from all possible values
  yMax <- max(c(noEmissionProbs$density,midEmissionProbs$density,
                largestEmissionProbs$density))
  
  png(filename=paste0(figureDir, metVarName, "_distribution.png"),
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
    "topleft",
    legend = c("no emissions", "[40,60] emission percentile", "top 5%"),
    col = c("black", "purple", "red"),
    lty=1,
    bty='n',
    cex=1.2
  )
  title(paste(metVarName, "distribution by emissions group") ,cex.main=2)
  
  dev.off()
  
  
}

plotMetByEmission(C, t2m, "t2m", "K")
plotMetByEmission(C, z, "z", "500 mb geopotential m**2/s**2")
plotMetByEmission(C, tp, "tp", "inches/day")
plotMetByEmission(C, RH, "RH", "RH%")
plotMetByEmission(C, srf_wind, "srf_wind", "surface wind speed m/s")



