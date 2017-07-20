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
import time as timer

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

print 'Working on loading the really big emissions file' 
timeStart = timer.time()
emissions = fire_nc.variables['C'][:]
dt = (timer.time() - timeStart) / 60.
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to load emissions array into workspace'
print '----------------------------------------------------------------------'

# Adjust fire lon from -180:180 lon to 0:360 lon to match met
fire_lon = fire_lon + 180.
fire_lon_g, fire_lat_g = np.meshgrid(fire_lon, fire_lat)

# flatten grids that do not change
fire_lon_f = fire_lon_g.flatten()
fire_lat_f = fire_lat_g.flatten()
met_lon_f = met_lon_g.flatten()
met_lat_f = met_lat_g.flatten()

# Now make these into tuple for griddata function
fire_points = (fire_lon_f, fire_lat_f)
met_points = (met_lon_f, met_lat_f)

nTStep = len(fire_time)

fire_new_grid = np.zeros((nTStep, len(met_lat), len(met_lon)))
deltaMass = np.zeros(nTStep)

# For each time step, find the best match in terms of min lon lat
# distance 
# TODO: Only loop over non zero locations? Over subregion?
timeStart = timer.time()
for t in range(4):
	print t
	# Loop over fire lon values
	for x in range(len(fire_lon)):
		# to make it a point, loop over each lat for each lon
		for y in range(len(fire_lat)):
		
			# Difference in regular grid arrays, units of deg
			dx = np.abs(fire_lon[x] - met_lon)
			dy = np.abs(fire_lat[y] - met_lat)
	
			# Index of best match
			#xi = np.argmin(dx)
			#yi = np.argmin(dy)
			xi = np.where(dx == dx.min())[0][0]
			yi = np.where(dy == dy.min())[0][0]

			# assign the data
			fire_new_grid[t, yi, xi] = fire_new_grid[t, yi, xi] + emissions[t, y, x]

	# How much mass is lost or gained?
	deltaMass[t] = np.sum(fire_new_grid[t, :, :]) - np.sum(emissions[t, :, :])

dt = (timer.time() - timeStart) / 60.
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to run the re-grid loop'
print '----------------------------------------------------------------------'


# Before looping handle time matching of fire and met data

# Get the first instance of fire data to regrid
#nTStep = len(fire_time)

#fire_new_grid = np.zeros((nTStep, len(met_lat), len(met_lon)))

#timeStart = timer.time()
#for i in range(10):
	
#	print i 

#	old_grid = emissions[i,:,:]
	
#	# NOTE: This does NOT conserve mtotal emission mass
#	new_grid = scipy.interpolate.griddata((fire_lon_f, fire_lat_f), old_grid.flatten(),
#	                                      (met_lon_g, met_lat_g), method='linear')

#	fire_new_grid[i,:,:] = new_grid

#dt = (timer.time() - timeStart) / 60.
#print '----------------------------------------------------------------------'
#print 'It took ' + str(dt) + ' minutes to run the loop'
#print '----------------------------------------------------------------------'
	




# TODO: see what percent change in total grams of burning for the whole days grid
# TODO: occur using this method. I do not want to loose mass or gain mass when
# TODO: moving to a new grid. 







# https://docs.scipy.org/doc/scipy/reference/generated/
# scipy.interpolate.interp2d.html#scipy.interpolate.interp2d
#scipy.interpolate.interp2d()


