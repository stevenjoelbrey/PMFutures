#!/usr/bin/env python2


###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# The purpose of this script is to regrid GFED data to ecmwf era-interim grid
# (or any lat lon passed) and save as a single daily nc file. 

# Example using scipy 
# http://christopherbull.com.au/python/scipy-interpolate-griddata/

dataDir  = '/barnes-scratch/sbrey/'
year     = 2003
species  = 'C' # 'C' , 'DM', 'small_fire_fraction' (These have daily fraction est.)


import sys
from netCDF4 import Dataset 
import matplotlib.pyplot as plt

from mpl_toolkits import basemap


import numpy as np
import scipy.interpolate


old_grid_data=np.random.rand(4,3)

#old grid dim
loni = np.array([109.94999695, 110.05000305, 110.15000153])
depi=np.array([3.04677272, 9.45404911, 16.36396599, 23.89871025])

#new grid dim
lon=np.arange(110.,110.3,.1) #NB: 110.2 outside of convex hull of old so will produce nan
depth=np.array([3.1,9,16,23])

#create mesh
X, Y = np.meshgrid(loni, depi)
XI, YI = np.meshgrid(lon,depth)

#interp
new_grid=scipy.interpolate.griddata((X.flatten(),Y.flatten()),old_grid_data.flatten() ,\ 
                                    (XI,YI),method='linear')

print "this is original"
print old_grid_data.reshape(4,3)
print ""
print "this is interp' by cubic"
print new_grid

print
print "this is diff"
print new_grid-old_grid_data.reshape(4,3)

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

######################
#def regridOneDay():
######################

old_grid = data

#old grid dim
loni = fire_lon
lati = fire_lat

#new grid dim
lon = met_lon #NB: 110.2 outside of convex hull of old so will produce nan
lat = met_lat

#create mesh
X, Y   = np.meshgrid(loni, lati)
XI, YI = np.meshgrid(lon, lat)

#interp
new_grid=scipy.interpolate.griddata((X.flatten(),Y.flatten()), old_grid.flatten() ,\
                                    (XI,YI), method='linear')










# https://docs.scipy.org/doc/scipy/reference/generated/
# scipy.interpolate.interp2d.html#scipy.interpolate.interp2d
#scipy.interpolate.interp2d()


