#!/usr/bin/env python2

###############################################################################
# ------------------------- Description ---------------------------------------
###############################################################################
# This script will be used to generate daily met fields from hourly nc files.
# The 6-houly data to average live in /barnes-scratch/sbrey/era_interim_nc_6_hourly

# Follows ----------------------------------------
# 	- get_ERA_Interim_data.py
# Precedes
#	- merge_yearly_nc.py

# TODO: Handle fg10 (wind gust) daily value creation.

###############################################################################
# ---------------------- Set analysis variables--------------------------------
###############################################################################
import sys
import cesm_nc_manager as cnm
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	hourlyVAR =  str(sys.argv[1]) # e.g. 'z'
	startYear =  int(sys.argv[2])
	endYear   =  int(sys.argv[3])

else:
	# Development environment. Set variables by manually here.
	hourlyVAR =  'rh2m'
	startYear = 2003
	endYear   = 2016



drive = cnm.getDrive()

dataDir   = drive + "era_interim_nc_6_hourly"
outputDir = drive + "era_interim_nc_daily"


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
if startYear == 'all':
	years = ['all']
else:
	years  = np.arange(startYear, endYear+1)

ii     = 0 # for counting total iterations
for year in years:

	year = str(year)
	print '---------------------------------------------------------------------'
	print 'Working on : ' + year
	print '---------------------------------------------------------------------'


	# Load 6-hourly data
	HourlyFile  = os.path.join(dataDir, hourlyVAR + "_" + year + ".nc")
	print HourlyFile
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
	yearList = []
	monthList = []
	hourList = []
	for i in range(len(time_hour)):
		dt = timedelta(hours=time_hour[i])
		date_new = date0 + dt
		t_new  = t0 + dt
		year_new = t_new.year
		t.append(t_new)
		dates.append(date_new)
		yearList.append(year_new)
		monthList.append(t_new.month)
		hourList.append(t_new.hour)

	t = np.array(t)
	dates = np.array(dates)
	dateYears = np.array(yearList)
	dateMonths = np.array(monthList)
	dateHours = np.array(hourList)


	# NOTE: Accumulation parameters (total precip (tp) and e) represent
	# NOTE: accumulated values from intitialization time. For these data those
	# NOTE: times are 00:00:00 and 12:00:00. I downloaded the data in 12 hour steps.
	# NOTE: So for tehse parameters, each time in the data represents a total for the
	# NOTE: previous 12 hours. This is why time series start at 12:00:00 for
	# NOTE: these fields and 00:00:00 for analysis fields.




	# Get all values numpy array into workspace
	print '---------------------------------------------------'
	print 'Loading the large variable array into the workspace'
	print '---------------------------------------------------'
	VAR_array = VAR[:]

	print 'Working on the large loop averaging 6-hourly values for each day'

	if (hourlyVAR != 'tp') & (hourlyVAR != 'e'): # Try checking first time hour?
		# these are the analysis variables that always require averages for a
		# given calendar date.

		print '---------------------------------------------------------------------'
		print 'Working with an analysis parameter whos first hour is 0. '
		print '---------------------------------------------------------------------'
		if dateHours[0] != 0:
			raise ValueError('The first hour of analysis field was not 0Z.')

		# Create structure to save daily data and do the averaging
		unique_dates = np.unique(dates)
		nDays = len(unique_dates) # might have to do a - 1 here now. Or search for feb 29th and set length based on that.
		nLon  = len(lon)
		nLat  = len(lat)


		# Create array to store daily averaged data, based on dimensions
		if len(dims) == 4:
			nLevel= len(level)
			dailyVAR  = np.zeros((nDays, nLevel, nLat, nLon))
		else:
			dailyVAR  = np.zeros((nDays, nLat, nLon))


		for i in range(nDays):

			# find unique day to work on
			indexMask = np.where(dates == unique_dates[i])[0]

			if len(dims) == 4:

				VAR_array_subset = VAR_array[indexMask, :, :, :]
				day_time_mean = np.mean(VAR_array_subset, 0)
				dailyVAR[i, :, : , :] = day_time_mean

			else:

				# Non-precip variables of this size need an average.
				VAR_array_subset = VAR_array[indexMask, :, :]
				day_time_mean    = np.mean(VAR_array_subset, 0)
				dailyVAR[i, :, : ] = day_time_mean


	elif (dateHours[0] == 12) & (dateHours[-1]) == 0 & (dateYears[-1] > int(year)):

		print '---------------------------------------------------------------------'
		print 'Working with an accumulation parameter with start time hour == 12. '
		print '---------------------------------------------------------------------'

		# These strange time conditions are all true when we are working with
		# tp and e accumulation forecast fields.

		# Precip units of m per 12 hour window. Requires a sum NOT an average.
		# Need matching dates noon and next dates midnight to get a days total.
		# e.g. total precip for Jan 1 is sum of tp at 01-01-year 12:00:00 AND
		# 01-02-year 00:00:00.

		# the last date in the time array will be the next year, since midnight or
		# 0Z.

		# In order for the code to work for these variables the same as the
		# analysis fields, we are going to subtract 12 hours from each time
		# element.
		t = t - timedelta(hours=12)

		nTime = len(t)
		if nTime % 2 != 0:
			raise ValueError("There is something wrong. Somehow there is a date without two 23 hour chuncks. ")

		nDays = len(t)/2
		nLon  = len(lon)
		nLat  = len(lat)

		# To make a mask of unique dates, we need to make timetime.datetime objects
		# a more simple datetime.date object.
		dates = []
		for i in range(nDays*2):
			dates.append(t[i].date())
		dates = np.array(dates)
		unique_dates = np.unique(dates)

		# Now that these strange time contrains have been met, we know we can
		# sum the values of every other. Create an array to store daily data in.

		dailyVAR = np.zeros((nDays, nLat, nLon))
		for j in range(nDays):

			# The hours that match our date.
			indexMask = np.where(dates == unique_dates[j])[0]

			if (dateHours[indexMask[0]] == 12) & (dateHours[indexMask[1]] == 0):

				# Subset the dataframe to include the two twelve hour slices we
				# want.
				timeSlice  = VAR[indexMask, :, :]

				# This statement makes sure we are really getting a time slice with
				# a dimension of 2, e.g. 2 12 hour segments.
				if timeSlice.shape[0] == 2:

					dailySum = np.sum(timeSlice, axis=0)
					dailyVAR[j,:,:] = dailySum

				else:

					raise ValueError("The time size was not two deep in time dim.")

				# if the sum of the dailyVAR array for this date is still zero,
				# no data was assigned.
				if np.sum(dailyVAR[j, :,:]) == 0:

					raise ValueError("No data was assigned to day index j = " + str(j))


	meansCompleteTime = timer.time()
	dt = (meansCompleteTime - timeStart) / 60.
	print '---------------------------------------------------------------------'
	print 'It took ' + str(dt) + ' minutes to create daily averages.'
	print '---------------------------------------------------------------------'


	# Check to see if the total amount of precip was conserved.
	if hourlyVAR == 'tp':
		originalSum = np.sum(VAR, axis=0)
		dailySum   = np.sum(dailyVAR, axis=0)
		dtp = np.abs(originalSum - dailySum)
		# ideally dtp is all zero. With float rounding issues it could be slightly
		# different. This matrix needs to be examined.
		maxDiff = np.max(dtp)
		print '------------------------------------------------------------------'
		print 'Maximum annual difference in rainfall is: ' + str(maxDiff)
		print '------------------------------------------------------------------'
		if maxDiff > 1e-10:
			raise ValueError("Total rainfall depth in meters not conserved within tolerance")


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

		dailyVAR_ = ncFile.createVariable(hourlyVAR,'f4',
									('time','level','latitude','longitude'),
									fill_value=VAR._FillValue)

		# While here create the level dimesion
		level_ = ncFile.createVariable('level', 'i4', ('level',))
		level_.units = level.units

	else:

		dailyVAR_ = ncFile.createVariable(hourlyVAR,
									'f4',('time','latitude','longitude'),
									fill_value=VAR._FillValue)

	# Assign the same units as the loaded file to the main variable
	# NOTE: Does this change values? Since netCDF4 package "accounts" for this?
	dailyVAR_.units = VAR.units
	dailyVAR_.scale_factor = VAR.scale_factor
	dailyVAR_.add_offset = VAR.add_offset
	dailyVAR_.missing_value = VAR.missing_value
	# NOTE: tests on the loaded values of nc files I write show that it does
	# NOTEL not actually matter if these off_set and scale factors are used.

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
	# NOTE: are 4 analysis snapshots space by 6 hours for any given date.
	# NOTE: However, tp (total precip) only has two chunks of 12 hourly data per
	# NOTE: day so this needs to be handled seperately. Because tp and e fields
	# NOTE: time were adjusted by minus 12 hours, all daily mean or sum fields
	# NOTE: have a time stamp of the 0Z for the date of the data.
	if (hourlyVAR == 'tp') | (hourlyVAR == 'e'):
		time = time[:] - 12

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