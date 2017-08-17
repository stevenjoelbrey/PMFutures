# FPA_FOD_analysis.R

# The goal of this script is to determine if metereology is different on days
# with varying levels of emissions. I want to see the conclusions of Balch et
# al. for myself, that human started fires are able to start in more wet and 
# cooler locations than their lightning started counterparts. 

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

