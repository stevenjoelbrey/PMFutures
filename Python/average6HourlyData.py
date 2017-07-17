#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to generate daily met fields from hourly nc files. 
# The 6-houly data to average live in /barnes-scratch/sbrey/era_interim_nc_6_hourly

# TODO: Merge daily variable yearly files into a single 1997-2016 file!

dataDir = "/barnes-scratch/sbrey/era_interim_nc_6_hourly"
outputDir = "/barnes-scratch/sbrey/era_interim_nc_daily"

import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	hourlyVAR =  str(sys.argv[1]) # e.g. 'Z'

else: 
	# Development environment. Set variables by hand here. 
	hourlyVAR =  'Z'
	year      = '2012'

import cesm_nc_manager as cnm
import os
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr
import datetime
import time as timer
import os.path

# Start the timer on this code
timeStart = timer.time()

# Load 6-hourly  data
HourlyFile  = os.path.join(dataDir, hourlyVAR + "_" + year + ".nc")	
nc          = Dataset(HourlyFile, 'r')
VAR         = nc.variables['z']   
time        = nc.variables['time'] 
time_hour   = np.array(time[:], dtype=int)
lon         = nc.variables['longitude']
lat         = nc.variables['latitude']
level       = nc.variables['level']

##########################################
# Handle date from hourly time dimensions
##########################################
t0 = datetime.datetime(year=1900, month=1, day=1, hour=0, minute=0, second=0) 
date0 = date(year=1900, month=1, day=1)
t = []
dates = []
for i in range(len(time_hour)):
	dt = timedelta(hours=time_hour[i])
	date_new = date0 + dt
	t_new  = t0 + dt 
	t.append(t_new)
	dates.append(date_new)
t = np.array(t)
dates = np.array(dates)

# Create structure to save daily data and do the averaging
unique_dates = np.unique(dates)
nDays = len(unique_dates)
nLon  = len(lon)
nLat  = len(lat)
nLevel= len(level)


# Create array to store daily averaged data
dailyVAR = np.zeros((nDays, nLevel, nLat, nLon))
VAR_array = VAR[:]


print 'Working on the large loop averaging 6-hourly values for each day'
for i in range(nDays):

	print str(i/float(nDays)*100.) + " % complete"

	# find unique day to work on
	indexMask = np.where(dates == unique_dates[i])[0]
	VAR_array_subset = VAR_array[indexMask, :, :, :]
	day_time_mean = np.mean(VAR_array_subset, 0)

	dailyVAR[i, :, : , :] = day_time_mean

	

meansCompleteTime = timer.time()
dt = (meansCompleteTime - timeStart) / 60. 
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to create daily averages.'
print '----------------------------------------------------------------------'



###############################################################################
# Write the new netcdf data with the exact same formatting as the
# data read here. 
# Create the Save name and make sure that this file does not exist
###############################################################################

outputFile = os.path.join(outputDir, hourlyVAR + "_" + year + ".nc")

# See if the file already exists
# os.path.isfile(outPutSavename):

print '----------------------------------------------------------------------'
print 'outputFile used:'
print outputFile
print '----------------------------------------------------------------------'

###############################################################################
# Write the daily averaged netCDF data
###############################################################################
ncFile = Dataset(outputFile, 'w', format='NETCDF4')
ncFile.description = 'Daily average of 6 hourly data'
ncFile.location = 'Global'
ncFile.createDimension('time',  nDays )
ncFile.createDimension('level',  nLevel )
ncFile.createDimension('latitude', nLat )
ncFile.createDimension('longitude', nLon )

# Create variables on the dimension they live on 
dailyVAR_ = ncFile.createVariable(hourlyVAR, 'f4',('time','level','latitude','longitude'))
dailyVAR_.units = VAR.units

time_ = ncFile.createVariable('time', 'i4', ('time',))
time_.units = time.units

level_ = ncFile.createVariable('level', 'i4', ('level',))
level_.units = level.units

latitude_ = ncFile.createVariable('lat', 'f4', ('latitude',))
latitude_.units = lat.units
 
longitude_ = ncFile.createVariable('lon', 'f4', ('longitude',))
longitude_.units = lon.units

# Write the actual data to these dimensions
dailyVAR_[:]     = dailyVAR
latitude_[:]     = lat[:]
longitude_[:]    = lon[:]
time_[:]         = time[0::4] # Take every 4th element, starting at 0th

ncFile.close()


dt = (timer.time() - timeStart) / 60.
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to run entire script.'
print '----------------------------------------------------------------------'



