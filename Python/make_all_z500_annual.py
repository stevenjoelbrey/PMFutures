#!/usr/bin/env python2

###############################################################################
# ------------------------- Description ---------------------------------------
###############################################################################
# This script will be used to generate daily met fields from 6-hourly nc files.
# The 6-houly data to average live in /barnes-scratch/sbrey/era_interim_nc_6_hourly



###############################################################################
# ---------------------- Set analysis variables--------------------------------
###############################################################################
domain = "_NA_" # "_NA_" subsets domain to North America lat lon. Any other
				# string returns the whole world




# Import the required modules
import sys
import os
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as pl
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr
import datetime
import time as timer
import os.path
import cesm_nc_manager as cnm


# Handle local file structure pathways
drive = cnm.getDrive()
dataDir = drive + "era_interim_nc_6_hourly/"
outputDir = drive + "era_interim_nc_daily/"

# The 6-hourly data to load
ncFile = dataDir  + 'z_all.nc'

nc = Dataset(ncFile, 'r')
allTime = nc.variables['time'][:]
VAR         = nc.variables['z']
time        = nc.variables['time']
time_hour   = np.array(time[:], dtype=int)
lon         = nc.variables['longitude']
lat         = nc.variables['latitude']
lonUnits    = lon.units
latUnits    = lat.units


# get unique year sequence
t, months, years = cnm.get_era_interim_time(time)
uniqueYears = np.unique(years)

# want just boring old date though, for looping
dates = []
for i in range(len(t)):
	d = datetime.datetime(t[i].year, t[i].month, t[i].day)
	dates.append(d)
dates = np.array(dates)

unique_dates = np.unique(dates)
nDays = len(unique_dates)
nLat = len(lat)
nLon = len(lon)

###############################################################################
print 'Working on the large loop averaging 6-hourly values for each day'
###############################################################################
timeStart = timer.time()

dailyVAR  = np.zeros((nDays, nLat, nLon), dtype=float)
ii = 0
for i in range(nDays): # nDays
	ii = ii + 1.
	print str(ii/float(nDays)*100.) + " % complete"

	# find unique day to work on
	indexMask = np.where(dates == unique_dates[i])[0]

	VAR_array_subset = VAR[indexMask, :, :]
	day_time_mean = np.mean(VAR_array_subset, 0)
	dailyVAR[i, :, : ] = day_time_mean


meansCompleteTime = timer.time()
dt = (meansCompleteTime - timeStart) / 60.
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to create daily averages.'
print '----------------------------------------------------------------------'

###############################################################################
# Handle spatially subsetting the data to be written.
###############################################################################
if domain == "_NA_":

	# These are the same _NA_ "North America" subset region used by
	# merge_yearly_nc.py

	minLat     = 30.
	maxLat     = 50.
	minLon     = 234.
	maxLon     = 259.

	dailyVAR, lat, lon = cnm.mask2dims(dailyVAR, lon[:], lat[:], tDim=0,
										xmin=minLon, xmax=maxLon,
										ymin=minLat, ymax=maxLat)

	outputFile = os.path.join(outputDir, 'z_all_NA_daily.nc')

else:
	outputFile = os.path.join(outputDir, 'z_all_daily.nc')


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
ncFile.createDimension('latitude', len(lat) )
ncFile.createDimension('longitude', len(lon) )

# Create variables on the dimension they live on
dailyVAR_ = ncFile.createVariable('z500',\
            'f4',('time','latitude','longitude'))

# Assign the same units as the loaded file to the main variable
dailyVAR_.units = VAR.units

# Create time variable
time_ = ncFile.createVariable('time', 'i4', ('time',))
time_.units = time.units
time_.calendar = time.calendar

# create lat variable
latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
latitude_.units = latUnits

# create longitude variable
longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
longitude_.units = lonUnits

# Write the actual data to these dimensions
dailyVAR_[:]     = dailyVAR
latitude_[:]     = lat[:]
longitude_[:]    = lon[:]


# NOTE: In general, every 4th element, starting at 0th, since there
# NOTE: are 4 sets of 6 hourly data for any given date.
tstep = len(time) / nDays
time_[:] = time[0::tstep]
# The difference in each time_[:] element in hours must be 24 or
# this is not working.
if np.unique(np.diff(time_[:])) != 24:
	ValueError('Difference in hours between days not all 24! Broken.')

# The data is written, close the ncFile and move on to the next year!
ncFile.close()


dt = (timer.time() - timeStart) / 60.
print '----------------------------------------------------------------------'
print 'It took ' + str(dt) + ' minutes to run entire script.'
print '----------------------------------------------------------------------'


