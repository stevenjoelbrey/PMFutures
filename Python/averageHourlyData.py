#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to generate daily (and 6 hourly) fields of Air
# qualility variables that are only available in hourly time steps.
# TODO: Make generic for any hourly variable. Have that variable read in in 
# TODO: as a passed argument from a shel script since this will be run as
# TODO: a qsub job. 

# TODO: Make this work for emissions also. 

hourlyVAR = 'PSL'
scenario  = '2000Base'
year      = 2000

# NOTE: Performance testing indicates that this script should take about 
# NOTE: ~100 minutes to average a decades worth of hourly data to daily. 

import cesm_nc_manager as cnm
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

print 'Working on the large loop averaging hourly values for each day'
for i in np.arange(0, 10): # nDays

	# We have to make a mask for each day of these hourly data.
	# Since the units are time are days from origin, all integers
	# are unique days at some hour (expressed by decimal). Se we 
	# need to mask each 24 hours. [0,1), [1,2), ...
	mask  = (time_nc[:] >= i) & (time_nc[:] < (i+1) ) 
	index = np.where(mask ==  True)[0]

	if len(index) != 24:
		raise ValueError('Day index does not have 24 hours. Code Broken')

	# Now that we have the index of a day, we need to take the 
	# average value from the data for those 24 hours 
	dayVAR            = Hourly_nc.variables[hourlyVAR][index, :, :]
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

###############################################################################
# Write the new netcdf data with the exact same formatting as the
# data read here. 
# Create the Save name and make sure that this file does not exist
###############################################################################
outPutSavename = cnm.makeAQNCFile(hourlyVAR, scenario, 'daily')

# Will have to change where this is placed becuase I do not own the output
# This is temporary untill I can get permission to write in the scratch
outPutSavenameTail = outPutSavename[69:]
outPutDir = '/home/sbrey/projects/PMFutures/Python_output/' + '2000Base/'
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

###############################################################################
# Write the daily averaged netCDF data
###############################################################################
ncFile = Dataset(outPutSavename, 'w', format='NETCDF4')
ncFile.description = 'Daily average of hourly data'
ncFile.location = 'Global'
ncFile.createDimension('ncols', nLon )
ncFile.createDimension('nrows', nLat )
ncFile.createDimension('time',  dailyVAR.shape[0] )

# Create variables on the dimension they live on 
dayVAR_ = ncFile.createVariable(hourlyVAR, 'f4', ('nrows','ncols','time'))
dayVAR_.units = Hourly_nc.variables[hourlyVAR].units

latitude = ncFile.createVariable('lat', 'f4', ('nrows',))
latitude.units = 'degrees north'
 
longitude = ncFile.createVariable('lon', 'f4', ('ncols',))
longitude.units = 'degrees east'

time_var = ncFile.createVariable('time', 'i4', ('time',))
time_var.units = 'daily dates'

# Write the actual data to these dimensions
dayVAR_[:,:,:]   = dailyVAR
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


 




