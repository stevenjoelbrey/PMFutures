# plotParameterSpace.R 

# This script will plot and enable exploring the relationships between
# meteorology events (numeric values) and fire emissions. I want the 
# functionality of this script to be the forerunner to a shiny app that will
# enable interactive map exploring. 

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

png(filename='test.png')
plot(RH, C)
dev.off()


