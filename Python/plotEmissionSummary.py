#!/usr/bin/env python2

# plotEmissionSummary.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to summarize and plot emissions from CESM. 
# A Multipanel plot will sumarize the present day, 2050, 2100 emissions for a 
# selected variable. Later on, this script will be able to make the standard 
# plots for desired time and season subsets. Much later, this code will be
# incorperated into a shiny app that will allow easy exploration of the model 
# output. 

# Load resources
import cesm_nc_manager as cnm
import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr


###############################################################################
# ------------------------- USER INPUT ---------------------------------------- 
###############################################################################
# TODO: Make sure these can work as agruments passed by the command line. 

variable = 'BC'    # Anything listed in outputFromYellowstone/FireEmissions
scenario = 'RCP85' # Will be a loop
year     = '2100'  # Will be a loop

ncFile   = cnm.makeEmissionNCFile(variable, scenario, year)
#dims     = cnm.getNCVarDims(ncFile, variable)

nc       = Dataset(ncFile, 'r')
bb       = nc.variables['bb'][:]
bbUnits  = nc.variables['bb'].units

# Handle the massive fill_value for bb so math does not go crazy later
bbMask = bb.mask
bb.data[bbMask==1] = 0

###############################################################################
# Handle time dimension. NOTE: Leap years do not exist in the model. 
###############################################################################
         
daysSincet0 = nc.variables['time'][:]   # only used for length of time dim
 
months      = [1,2,3,4,5,6,7,8,9,10,11,12]
daysInMonth = [31,28,31,30,31,30,31,31,30,31,30,31] # Never a 29 day February
nMonth      = len(daysInMonth)
daysInYear  = np.sum(daysInMonth)
nYears      = len(daysSincet0) / daysInYear 

# start year is always 10 years before year selection string. 
# start month and day are always Jan 1 for emissions
startYear = int(year) - 10
t0        = date(startYear, 1, 1) 

t = [] # where time will be appended
for deltaYear in range(nYears):
	# Advance the year
	YEAR = t0.year + deltaYear
	for m in range(nMonth):
		dim = daysInMonth[m]
		for d in np.linspace(1, dim, dim, dtype='int'):
			newDate = date(YEAR, months[m], d)
			t.append(newDate)
			
# Make useful array
t = np.array(t)

###############################################################################
# Get spacial dimensions
###############################################################################
lon = nc.variables['lon'][:]
lat = nc.variables['lat'][:]

nlon = len(lon)
nlat = len(lat)

nc.close()

###############################################################################
# Load model grid cell atributes. These do not depend on time. Must be subset 
# to match bb grid that is already loaded into workspace 
###############################################################################
dataDir = '/fischer-scratch/sbrey/outputFromYellowstone/FireEmissions/'
moduleLayers = ['area','landfrac','landmask']
layer = {}
layerUnits = {}
for f in moduleLayers:
	fileName = dataDir + 'cesm130_clm5_firemodule_'+ f + '_f09x125.nc'
	nc = Dataset(fileName, 'r')
	layer[f] = nc.variables[f][:]
	allLons = nc.variables['lon'][:]
	allLats = nc.variables['lat'][:]
	nc.close()

# Subset each based on the dimensions where lat and lon match

# Get area from km**2 to m*2
m2Perkm2 = 1e6 * 1.
area = layer['area']
area = area * m2Perkm2
landMask = layer['landmask']

# Get rid of the insanely large fill value that messes up multiplication. 
area.data[landMask==0] = 0

lonI = np.where(np.in1d(allLons, lon))[0]
latI = np.where(np.in1d(allLats, lat))[0]

# TODO: Handle this masking in a non horrible way. 
areaSubset = area[latI.min():(latI.max()+1) , lonI.min():(lonI.max()+1)]
#fill_value = area.fill_value
area = areaSubset
#area.fill_value = 0

# create boolean mask that is size of area but true at these indicies
#areaMask = np.zeros(area.shape, dtype=bool)
#areaMask[latIndices, lonIndices] = True

#latSubset  = area[latIndices,:]
#areaSubset = latSubset[:,lonIndices]
#area = areaSubset

# Get emissions from kg/m2/sec to kg in time period of interest
secondsPerDay = 24. * 60**2
kgPerDay = np.zeros(bb.shape)
for i in range(len(t)):
	kgPerDay[i, :, :] = bb[i,:,:] * (area.data * secondsPerDay)

# Now sum over time to come up with maximum values
kgTotal = np.sum(kgPerDay,0) * len(t)



















