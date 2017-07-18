#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to generate daily met fields from hourly nc files. 
# The 6-houly data to average live in /barnes-scratch/sbrey/era_interim_nc_6_hourly

# TODO: Merge daily variable yearly files into a single 1997-2016 file!



###############################################################################
# ---------------------- Set analysis variables--------------------------------
###############################################################################
import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	hourlyVAR =  str(sys.argv[1]) # e.g. 'z'
	startYear =  int(sys.argv[2])
	endYear   =  int(sys.argv[3])

else: 
	# Development environment. Set variables by manually here. 
	hourlyVAR =  't'
	startYear = 1997
	endYear   = 1997


dataDir = "/barnes-scratch/sbrey/era_interim_nc_6_hourly"
outputDir = "/barnes-scratch/sbrey/era_interim_nc_daily"
	

print '----------------------------------------------------------------------'
print 'Working on variable: ' + hourlyVAR + ' for years: ' + str(startYear) +\
      '-'+ str(endYear)
print '----------------------------------------------------------------------'


# Import the required modules
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

# Loop over the selected years 
years  = np.arange(startYear, endYear+1)
ii     = 0 # for counting total iterations	
for year in years:

	year = str(year)
	print '---------------------------------------------------------------'
	print 'Working on : ' + year
	print '---------------------------------------------------------------'	
	

	# Load 6-hourly data
	HourlyFile  = os.path.join(dataDir, hourlyVAR + "_" + year + ".nc")	
	nc          = Dataset(HourlyFile, 'r')
	VAR         = nc.variables[hourlyVAR]   
	time        = nc.variables['time'] 
	time_hour   = np.array(time[:], dtype=int)
	lon         = nc.variables['longitude']
	lat         = nc.variables['latitude']

	# Some variables live on (t, level, lat, lon) grids, others (t, lat, lon)
	# Find out which one using dimension keys
	# e.g. [u'longitude', u'latitude', u'time'] for 'sp'
	#      [u'longitude', u'latitude', u'level', u'time'] for 'z'
	dims = nc.dimensions.keys()
	if len(dims) == 4:
		level = nc.variables['level']

	#######################################################################
	# Handle date from hourly time dimensions
	#######################################################################
	if time.units == 'hours since 1900-01-01 00:00:0.0':
		# For time datetime array
		t0 = datetime.datetime(year=1900, month=1, day=1,\
                                       hour=0, minute=0, second=0) 
		# For simply getting dates 
		date0 = date(year=1900, month=1, day=1)

	else:
		raise ValueError('Unknown time origin! Code will not work.')		

	# Create arrays to store datetime objects
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


	# Create array to store daily averaged data, based on dimensions
	if len(dims) == 4:
		nLevel= len(level)
		dailyVAR  = np.zeros((nDays, nLevel, nLat, nLon))
	else:
		dailyVAR  = np.zeros((nDays, nLat, nLon))

	# Get all values numpy array into workspace
	VAR_array = VAR[:]

	print 'Working on the large loop averaging 6-hourly values for each day'
	for i in range(nDays):
		ii = ii + 1.
		print str(ii/(float(nDays)*len(years))*100.) + " % complete"
		
		# find unique day to work on
		indexMask = np.where(dates == unique_dates[i])[0]

		if len(dims) == 4:

			VAR_array_subset = VAR_array[indexMask, :, :, :]
			day_time_mean = np.mean(VAR_array_subset, 0)
			dailyVAR[i, :, : , :] = day_time_mean
		else:

			VAR_array_subset = VAR_array[indexMask, :, :]
			day_time_mean = np.mean(VAR_array_subset, 0)
			dailyVAR[i, :, : ] = day_time_mean	

	
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
	ncFile.createDimension('latitude', nLat )
	ncFile.createDimension('longitude', nLon )
	
	# Create variables on the dimension they live on 
	if len(dims) == 4:

		ncFile.createDimension('level',  nLevel )

		dailyVAR_ = ncFile.createVariable(hourlyVAR,\
		            'f4',('time','level','latitude','longitude'),
			    fill_value=VAR._FillValue)

		# While here create the level dimesion
		level_ = ncFile.createVariable('level', 'i4', ('level',))
		level_.units = level.units

	else:
		dailyVAR_ = ncFile.createVariable(hourlyVAR,\
		            'f4',('time','latitude','longitude'),
                            fill_value=VAR._FillValue)

	# Assign the same units as the loaded file to the main variable
	dailyVAR_.units = VAR.units
	dailyVAR_.scale_factor = VAR.scale_factor
	dailyVAR_.add_offset = VAR.add_offset
	dailyVAR_.missing_value = VAR.missing_value

	# Create time variable	
	time_ = ncFile.createVariable('time', 'i4', ('time',))
	time_.units = time.units
	time_.calendar = time.calendar

	# create lat variable
	latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
	latitude_.units = lat.units
	
	# create longitude variable
	longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
	longitude_.units = lon.units

	# Write the actual data to these dimensions
	dailyVAR_[:]     = dailyVAR
	latitude_[:]     = lat[:]
	longitude_[:]    = lon[:]
	if len(dims) == 4:
		level_[:] = level[:]

	# NOTE: In general, every 4th element, starting at 0th, since there
	# NOTE: are 4 sets of 6 hourly data for any given date. However, tp
	# NOTE: (total precip) only has two chunks of 6 hourly data per day.
	# NOTE: so this needs to be handled seperately.
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



