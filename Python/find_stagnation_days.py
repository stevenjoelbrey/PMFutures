#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script is designed to create a mask of stagnation days as defined by 
# Wang and Angell [1999]. The three variables read in and assessed for 
# stagnation condition are 500, 1000 (taken as SLP) geostrophic winds, and 
# precipitation. 

# Read arguments passed via the command line if there are any. Otherwise use the
# interactive/development chioces list in else statement. 

# TODO: Should I consider adding the creation of other masks that reqiore daily 
# TODO: look into one file so that they are all made at the same time? 

import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	scenario    = sys.argv[1]
	wind1000Lim = sys.argv[2]   # m/s  
	wind500Lim  = sys.argv[3]   # m/s
	precLim     = sys.argv[4]   # inches/day
else:
	print 'Using default arguments.'
	# These are the defualt definitions of stagnation defined:
	# http://www.arl.noaa.gov/documents/reports/atlas.pdf
	scenario    = '2000Base'
	wind1000Lim = 8.   # m/s  
	wind500Lim  = 13.  # m/s
	precLim     = 0.01 # inches/day


import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import numpy.ma as ma
from datetime import date
import matplotlib.ticker as tkr
import datetime
import cesm_nc_manager as cnm
import time as timer

startTime = timer.time()

# Load geostophic wind components
ugFile = cnm.makeAQNCFile('ug_P', scenario, 'daily')
ug_nc  = Dataset(ugFile, 'r')
ug     = ug_nc.variables['ug']

vgFile = cnm.makeAQNCFile('vg_P', scenario, 'daily')
vg_nc  = Dataset(vgFile, 'r')
vg     = vg_nc.variables['vg']

# Get the wind dimensions from vg file
plevel  = vg_nc.variables['plevel'][:]  
lat     = vg_nc.variables['lat'][:]
lon     = vg_nc.variables['lon'][:]
time	= vg_nc.variables['time'][:]
nTime   = len(time)
nLat    = len(lat)
nLon    = len(lon) 

# Open the precipitation file
precFile = cnm.makeAQNCFile('PRECT', scenario, 'daily')
prec_nc  = Dataset(precFile, 'r')
prec     = prec_nc.variables['PRECT'] # [m/s]

# prec needs to be converted from m/s to incher/day
inchPerM = 39.3701      # [inch/m]
secondsPerDay = 86400.  # [s/day]
mPersToInPerDay = inchPerM * secondsPerDay # [(inch s) / (m day)]

# Convert m/s to inches/day
prec = prec[:] * mPersToInPerDay

# Find dates that sea-level and 500 hPa geo winds are less than 8 m/s 
# 13 m/s respectivly
wind1000Index = np.where(plevel==1000.)[0][0]
wind500Index  = np.where(plevel==500.)[0][0]

# Create numpy arrays to store mask of stagnation events
stagnationMask = np.zeros(prec.shape, dtype=int)

for i in range(nTime):  

	print i 

	# Handle 1000 mb geo winds and crazy high values that are near eq.
	ug1000 = ug[i, wind1000Index, :, :]
	vg1000 = vg[i, wind1000Index, :, :]
	ug1000 = ma.masked_where(ug1000 > 200., ug1000)
	vg1000 = ma.masked_where(vg1000 > 200., vg1000)

	wind1000Mag = np.sqrt(ug1000**2 + vg1000**2)
	mask1000    = np.array(wind1000Mag < wind1000Lim, dtype=bool)

	# Handle 500 mb geo winds and crazy high values that are near eq.
	ug500 = ug[i, wind500Index, :, :]
	vg500 = vg[i, wind500Index, :, :]
	ug500 = ma.masked_where(ug500 > 200., ug500)
	vg500 = ma.masked_where(vg500 > 200., vg500)

	# Calculate wind speed magnitude
	wind500Mag = np.sqrt(ug500**2 + vg500**2)
	mask500    = np.array(wind500Mag < wind500Lim, dtype=bool) 

	# Convert m/s to inches per day
	maskRain     = prec[i,:,:] < precLim

	# Create a mask of where all three are true
	m = np.array(mask1000 & mask500 & maskRain, dtype=int)
	print np.unique(m)
	stagnationMask[i, :, :] = m

###############################################################################
# Sanity check the output before writing the mask to an nc file 
###############################################################################
stagPrec  = ma.masked_where(stagnationMask == 0, prec, False)

print 'The max value of stagPrec is: ' + str(stagPrec.max())
print 'The unique values in stagnationMask are: ' + str(np.unique(stagnationMask))

# Make sure min values of these masked arrays is above the desired threshold
if stagPrec.max() >= precLim:
	print 'The maximum value of precip on stangation days exceeds threshold!'
	raise ValueError("This means creating the mask has failed!!!!!")


###############################################################################
# Write the stagnation mask as daily netCDF data
###############################################################################

saveName = cnm.makeAQNCFile('stagnation_mask', scenario, 'daily')
outPutSavename = saveName

ncFile = Dataset(outPutSavename, 'w', format='NETCDF4')
ncFile.description = 'Mask indicating stangation days'
ncFile.location = 'Global'
ncFile.createDimension('time', nTime )
ncFile.createDimension('lat', nLat )
ncFile.createDimension('lon', nLon )

# Create variables on the dimension they live on 
maskVar = ncFile.createVariable('stagnationMask', 'i', ('time','lat','lon'))
maskVar.units = '1 indicates stagnition condition is present. 0 means it is not.'

time_var = ncFile.createVariable('time', 'i4', ('time',))
time_var.units = 'days from origin'

latitude = ncFile.createVariable('lat', 'f4', ('lat',))
latitude.units = 'degrees north'
 
longitude = ncFile.createVariable('lon', 'f4', ('lon',))
longitude.units = 'degrees east'

# Write the actual data to these dimensions
maskVar[:]       = stagnationMask
latitude[:]      = lat
longitude[:]     = lon
time_var[:]      = time

ncFile.close()

writingComplete = timer.time()
dt = (writingComplete - startTime) / 60. 
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to run the entire script.'
print '----------------------------------------------------------------------'













	
