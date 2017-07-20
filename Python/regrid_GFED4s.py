#!/usr/bin/env python2


###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# The purpose of this script is to regrid GFED data to ecmwf era-interim grid
# (or any lat lon passed) and save as a single daily nc file. 

dataDir  = '/barnes-scratch/sbrey/'
year     = 2003
species  = 'C' # 'C' , 'DM', 'small_fire_fraction' (These have daily fraction est.)

import numpy as np
import sys
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
#import scipy.interpolate 
from mpl_toolkits import basemap


# Get the era-interim (or other MET grid). They ALL live on the same grid, except
# time
# TODO: Handle time matching in this script as well!
metFile = dataDir + 'era_interim_nc_6_hourly/d2m_2003.nc' 
met_nc = Dataset(metFile, 'r')
met_lat = met_nc.variables['latitude'][:]
met_lon = met_nc.variables['longitude'][:]
met_time= met_nc.variables['time']

# Grid the coordinate arrays
met_lon_g, met_lat_g = np.meshgrid(met_lon, met_lat)

# Get emissions file
fire_file = dataDir + 'GFED4s/GFED4.1s_C_2003_2016.nc'
fire_nc   = Dataset(fire_file, 'r')
fire_lat  = fire_nc.variables['latitude'][:] 
fire_lon  = fire_nc.variables['longitude'][:] 
fire_time = fire_nc.variables['time']

# Adjust fire lon from -180:180 lon to 0:360 lon to match met
fire_lon = fire_lon + 180.
fire_lon_g, fire_lat_g = np.meshgrid(fire_lon, fire_lat)

# Get the first instance of fire data to regrid
data = fire_nc.variables['C'][0,:,:]



lons_sub, lats_sub = np.meshgrid(lons[::4], lats[::4])

fire_data_coarse = basemap.interp(data, fire_lon_g, fire_lat_g, met_lon_g, met_lat_g, order=1)


# https://stackoverflow.com/questions/25544110/regridding-regular-netcdf-data


# https://docs.scipy.org/doc/scipy/reference/generated/
# scipy.interpolate.interp2d.html#scipy.interpolate.interp2d
#scipy.interpolate.interp2d()


