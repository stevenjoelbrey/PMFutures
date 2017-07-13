#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to generate daily (and 6 hourly) fields of Air
# qualility variables that are only available in hourly time steps.

# TODO: Make this work for emissions also. 

# NOTE: Performance testing indicates that this script should take about 
# NOTE: ~100 minutes to average a decades worth of hourly data to daily. 

import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	hourlyVAR =  str(sys.argv[1]) # e.g. 'PSL'
	scenario  =  str(sys.argv[2]) # e.g. '2000Base'

else: 
	# Development environment. Set variables by hand here. 
	hourlyVAR =  'PSL'
	scenario  =  '2000Base'


import cesm_nc_manager as cnm
import os
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import numpy.ma as ma
from datetime import date
import matplotlib.ticker as tkr
import datetime
import time as timer

# Start the timer on this code
timeStart = timer.time()

# Load hourly PSL data
HourlyFile  = cnm.makeAQNCFile(hourlyVAR, scenario, 'hourly')
Hourly_nc   = Dataset(HourlyFile, 'r')    # This is an 11GB file. Soo, holy smokes
time_nc     = Hourly_nc.variables['time'] # days since YYYY-01-01 00:00:00 
lon         = Hourly_nc.variables['lon']
lat         = Hourly_nc.variables['lat']

# Divide by hours to get days yo
nDays = len(time_nc)/24 
nLat  = len(lat)
nLon  = len(lon)

# Create array to store daily averaged data
dailyVAR = np.zeros((nDays, nLat, nLon))

# Load the date file that defines the dates. 
DateFile  = cnm.makeAQNCFile('date', scenario, 'daily') 
date_nc   = Dataset(DateFile, 'r')
d         = date_nc.variables['date'][:]
t         = cnm.dateNumToDate(d)

# Close the date file. 
date_nc.close()

# Change time_nc to an int, so each int clearly represents a calendar 
# day
dayInt = np.array(time_nc[:], dtype=int)
uniqueDayInts = np.unique(dayInt)


print 'Working on the large loop averaging hourly values for each day'
for i in uniqueDayInts: 
#for i in range(0,10):

	# We have to make a mask for each day of these hourly data.
	# Since the units are time are days from origin, all integers
	# are unique days at some hour (expressed by decimal). Se we 
	# need to mask each 24 hours. [0,1), [1,2), ...
	indexMask = np.where(dayInt == i)[0]

	if len(indexMask) != 24:
		# This is exspected for the last day. Make sure this is last day
		if i != uniqueDayInts[-1]:
			raise ValueError('Day index does not have 24 hours. Code Broken')
	else:
		# This day index looks good. Take the daily mean fo 24 hour values. 
		dayVAR            = Hourly_nc.variables[hourlyVAR][indexMask, :, :]
		dayMeanVAR        = np.mean(dayVAR, axis=0)
		dailyVAR[i, :, :] = dayMeanVAR


meansCompleteTime = timer.time()
dt = (meansCompleteTime - timeStart) / 60. 
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to create daily averages.'
print '----------------------------------------------------------------------'

# Make sure the time array matches the daily mean array in length and
# date. This dimension needs to be described correctly. 
t_subset = d[0:nDays]
if len(t_subset) != dailyVAR.shape[0]:
	raise ValueError('Labels for dates and dates averaged to no match!')

###############################################################################
# Write the new netcdf data with the exact same formatting as the
# data read here. 
# Create the Save name and make sure that this file does not exist
###############################################################################
outPutSavename = cnm.makeAQNCFile(hourlyVAR, scenario, 'daily')

# Will have to change where this is placed becuase I do not own the output
# This is temporary untill I can get permission to write in the scratch

# Directory names have different charactor lengths between base and RCP 
if scenario == '2000Base':
	firstChar = 69
else:
	firstChar = 70
	
outPutSavenameTail = outPutSavename[firstChar:]
outPutDir = '/home/sbrey/projects/PMFutures/Python_output/' + scenario + '/'
outPutSavename =  outPutDir + outPutSavenameTail

fileCount = -1
while os.path.isfile(outPutSavename):
	fileCount += 1
	print outPutSavename
	print 'File name already used. outPutSavename will be altered.'
	outPutSavename =  outPutDir + str(fileCount) + '_' + outPutSavenameTail

# This is the most important check. Do not want to overwrite orignal output. 
if outPutSavename == HourlyFile:
	raise ValueError('You are going to overwrite orignal data. Killing program.')

print '----------------------------------------------------------------------'
print 'outPutSavename used:'
print outPutSavename
print '----------------------------------------------------------------------'

###############################################################################
# Write the daily averaged netCDF data
###############################################################################
ncFile = Dataset(outPutSavename, 'w', format='NETCDF4')
ncFile.description = 'Daily average of hourly data'
ncFile.location = 'Global'
ncFile.createDimension('time',  len(t_subset) )
ncFile.createDimension('lat', nLat )
ncFile.createDimension('lon', nLon )

# Create variables on the dimension they live on 
dayVAR_ = ncFile.createVariable(hourlyVAR, 'f4', ('time','lat','lon'))
dayVAR_.units = Hourly_nc.variables[hourlyVAR].units

time_var = ncFile.createVariable('time', 'i4', ('time',))
time_var.units = 'daily dates'

latitude = ncFile.createVariable('lat', 'f4', ('lat',))
latitude.units = 'degrees north'
 
longitude = ncFile.createVariable('lon', 'f4', ('lon',))
longitude.units = 'degrees east'

# Write the actual data to these dimensions
dayVAR_[:]       = dailyVAR
latitude[:]      = lat
longitude[:]     = lon
time_var[:]      = t_subset

ncFile.close()
Hourly_nc.close()

writingComplete = timer.time()
dt = (writingComplete - meansCompleteTime) / 60. 
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to write the netCDF data of daily averages.'
print '----------------------------------------------------------------------'

dt = (timer.time() - timeStart) / 60.
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to run entire script.'
print '----------------------------------------------------------------------'


# ncFile = '../Python_output/2000Base/2_cesm122_fmozsoa_f09f09_2000_fires_00.PSL.daily.200001-201012.nc'




