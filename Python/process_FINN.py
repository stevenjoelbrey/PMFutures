#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# The purpose of this script is to make gridded FINN data made by the GEOS-Chem
# community suitable for this projects research goals. I want to make FINN
# format as close to GFED as possible, for easy comparisons. 

# TODO: 
# 1) Handle the different time array
# 2) Make sure it makes sense to add the different eco regions together. 
# 3) re-grid final desired product, global
# 3) figure out if you can save carbon out. 


import sys # for reading command line arguments
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
import numpy as np
import scipy.interpolate
import time as timer
from datetime import date
from datetime import timedelta
from datetime import datetime

# Consider making this directory a command line argument
dataDir  = '/barnes-scratch/sbrey/FINN/'
f = dataDir + 'FINN_daily_2012_0.25x0.25.compressed.nc'

# Explore FINN data before actually re-gridding. 
nc = Dataset(f, 'r')

# I want all burning as an easy to grab variable. But re-gridding all individual sounds 
# good too. 
time = nc.variables['time']
lon  = nc.variables['lon'][:]
lat  = nc.variables['lat'][:]

##########################################################################################
# Make time into same units as ecmwf and GFED4s
##########################################################################################
if time.units == 'hours since 1985-01-01 00:00:00':
	
	t0 = datetime(year=1985, month=1, day=1, hour=0, minute=0, second=0)
	
	# Desired origin              
	t0_new = datetime(year=1900, month=1, day=1, hour=0, minute=0, second=0)
	dt = t0 - t0_new 
	dt_hours = dt.total_seconds() / 60**2
	
	# This array is hours from t0_new
	time_new = time[:] + dt_hours

	# Make a nice month and time array for masking 
	nTime = len(time_new)
	t = []
	month = []
	year = []
	for i in range(nTime):
		dt = timedelta(hours=int(time_new[i]))
		t_new  = t0_new + dt 
		t.append(t_new)
		month.append(t_new.month)
		year.append(t_new.year)
	t = np.array(t)
	month = np.array(month)
	year = np.array(year)
	# replace the hours from origin values of time
	time = time_new # units of hours from 1990-01-01 00:00:00


##########################################################################################
# Calculate grid cell area in meters squared use website below for math reference
# https://badc.nerc.ac.uk/help/coordinates/cell-surf-area.html
##########################################################################################

nLon = len(lon)
nLat = len(lat)
grid_area = np.zeros( (nLat, nLon) )

R = 6371.* 1000 # radius of the earth in meters
# Math needs to be in radians
lat_rad = lat * np.pi / 180.
lon_rad = lon * np.pi / 180.
 
# NOTE: difference in longitude never changes, always 0.25, but will leave in the loop
# NOTE: for consistency  
for i in range(nLat-1):
	for j in range(nLon-1):
		dx = (lon_rad[j+1] - lon_rad[j])
		dy = ( np.sin(lat_rad[i+1]) - np.sin(lat_rad[i]) )
		S = R**2 * dx * dy
		grid_area[i, j] = S

# These are the veg types
# 1 SavannaGrasslands
# 2 WoodySavannah
# 3 TropicalForest
# 4 TemperateForest
# 5 Boreal
# 9 Crops
    
SavannaGrasslands = nc.variables['fire_vegtype1'][:]
WoodySavannah     = nc.variables['fire_vegtype2'][:]
TropicalForest    = nc.variables['fire_vegtype3'][:]
TemperateForest   = nc.variables['fire_vegtype4'][:] 
Boreal            = nc.variables['fire_vegtype5'][:]
Crops             = nc.variables['fire_vegtype9'][:]

# Add the emissions, units of kg/m2/s together for total emissions of CO2
CO2 = SavannaGrasslands + WoodySavannah +\
      TropicalForest + TemperateForest +\
      Boreal + Crops

def convertUnits(var, grid_area):
	"""Desired units are g C / grid / day, from kg/m2/s"""
	#TODO: look into np.tensordot() function
	for i in range(var.shape[0]):
		var[i,:,:] = var[i,:,:] * grid_area
	
	var = var * 1000. * 86400.0 # g/kg * seconds/day
	
	return var 

SavannaGrasslands = convertUnits(SavannaGrasslands, grid_area)
WoodySavannah     = convertUnits(WoodySavannah, grid_area)
TropicalForest    = convertUnits(TropicalForest, grid_area)
TemperateForest   = convertUnits(TemperateForest, grid_area)
Boreal            = convertUnits(Boreal, grid_area)
Crops             = convertUnits(Crops, grid_area)

CO2 = convertUnits(CO2, grid_area)

# Plot the total emissions as a sanity check on the grid 
all_emissions = np.sum(CO2, axis=0)
m = Basemap(projection='robin',lon_0=0,resolution='c')
lons, lats = np.meshgrid(lon, lat)
x,y=m(lons,lats)
m.drawcoastlines()
#m.fillcontinents(color='coral',lake_color='aqua')
# draw parallels and meridians.
#plt.pcolormesh(x,y, sp_slice)
c = plt.pcolormesh(x, y, all_emissions )
#m.drawmapboundary(fill_color='aqua')
plt.title("2012 FINN emissions")
plt.show(block=False)

plt.close()

##########################################################################################
# Save with names matching GFED as closely as possible, e.g., latitude not lat for dim 
# name 
##########################################################################################

outputFile = 'FINN_CO2_2012.nc'

ncFile = Dataset(outputFile, 'w', format='NETCDF4')
ncFile.description = 'Gridded FINN data created for GEOS-Chem. No Wiki.'
ncFile.location = 'Global'
ncFile.createDimension('time',  len(time) )
ncFile.createDimension('latitude', nLat )
ncFile.createDimension('longitude', nLon )

CO2_ = ncFile.createVariable('CO2','f4',('time','latitude','longitude'))
CO2_.units = 'g C02 / grid cell / day'	            

grid_area_ = ncFile.createVariable("grid_area", 'f4', ('latitude', 'longitude')) 	            
grid_area_.units = 'm**2'
	            
# Create time variables	
time_ = ncFile.createVariable('time', 'f4', ('time',))
time_.units = 'hours since 1900-01-01 00:00:00'

# create lat variable
latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
latitude_.units = 'degrees north'

# create longitude variable
longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
longitude_.units = 'degrees east'

# Write the actual data to these dimensions
VAR_[:]       = CO2
grid_area_[:] = grid_area
latitude_[:]  = lat
longitude_[:] = lon
time_[:]      = time

ncFile.close()		

##########################################################################################
# Save the same thing only where each species is a variable! 
##########################################################################################




