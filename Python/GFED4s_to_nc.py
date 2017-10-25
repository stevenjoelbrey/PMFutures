#!/usr/bin/env python2

###############################################################################
# ------------------------- Description ---------------------------------------
###############################################################################
# The purpose of this script is read the hdf5 GFED4s data and save out the same
# data as daily gridded data in nc file format for a single species. These are
# saved as yearly files.
#
# This script reads in data downloaded from the web where no changes have yet
# been made.
# When the startYear and endYear argument are different the different years
# data are merged and saved to an .nc file that has the year range in the file
# name.
# Functions
# 	getYearlyEmissions()
#	getMonthlyBurnArea()

# TODO: Function that just gets the monthly emissions?

# Daily emissions estimates made possible by
# http://onlinelibrary.wiley.com/doi/10.1029/2011JD016245/abstract

# Follows ----------------------------------------
# 	- NA
# Precedes
#	- any script that reads in GFED4s in .nc format.

# Datasource: http://www.geo.vu.nl/~gwerf/GFED/GFED4/
# Data README: http://www.geo.vu.nl/~gwerf/GFED/GFED4/Readme.pdf

# Clear all  variables before running this script.
#import sys
#sys.modules[__name__].__dict__.clear()
import sys
import h5py # if this creates an error please make sure you have the h5py library
import os
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import datetime
from datetime import date
from datetime import timedelta
from datetime import datetime
import cesm_nc_manager as cnm


# TODO: estimate daily burn area and save this out!
# TODO: include 'basis_regions' in nc output?
# TODO: Save all the two dimensional attributes as thier own NETCDF file
startYear = 2003
endYear   = 2003  # If different than startYear, they will be appended.
species   = 'DM'  # 'C' , 'DM' # (These have daily fraction est.)
getDaily  = False # execute code to create daily nc
getMonthly= True # execute code to create monthly nc

# Figure out what machine this code is running on. Set file paths.
drive = cnm.getDrive()
dataDir   = drive + 'GFED4s/'


# Months to loop over
months = ['01', '02', '03', '04', '05','06',\
          '07', '08', '09', '10', '11', '12']


def getDailyEmissions(dataDir, year, months, species):
	"""This function gets all the data for a species for a given year and returns
	   time, lat, lon arrays that describe daily species emission data on.

		return: time, latitude, longitude, yearData
	"""

	yearFile = 'GFED4.1s_' + str(year) + '.hdf5'
	f = h5py.File(dataDir + yearFile, 'r')

	# Get dimensions
	latitude  = f['lat'][:]
	longitude = f['lon'][:]
	nLat = latitude.shape[0]
	nLon = longitude.shape[1]

	# Get grid cell area [m**2]
	grid_cell_area_m2 = f['ancill/grid_cell_area'][:]

	# Create an array with the correct lat and lon dimension to append data
	# NOTE: will trim 366th day if no assignment is made
	yearData = np.zeros((366, latitude.shape[0], latitude.shape[1]))
	yearData[365,:,:] = -1

	# Create an array to append datetime.date objects to
	date0 = date(year=year, month=1, day=1) # reference date in Jan 1 of year
	time  = []

	jDay = 0 # Be careful about leap yaers?
	for m in months:

		print 'Getting ' + str(year) + ' ' + m + ' month daily data'

		# Write species emissions path
		speciesDir = '/emissions/' + m + '/' + species + '/'

		# Get the months emission array
		month_emission = f[speciesDir][:]

		# How many days?
		days         = f['/emissions/' + m + '/daily_fraction/']
		nDaysInMonth = len(days.keys())


		# because keys() does not put them in order
		dayNumber    = np.arange(1,nDaysInMonth+1)
		month_daily_frac  = np.zeros((nDaysInMonth, nLat, nLon))

		# loop through the days the monthly emissions are distributed over
		# keep track of daily_fraction

		for i in range(nDaysInMonth):

			# Advance the JDay Count (after adding dt to date0, since origin jan1)
			dt = timedelta(days=jDay)
			time.append(date0+dt)
			jDay = jDay + 1

			# Get fraction of monthly emissions that occured on THIS day
			dayString = 'day_' + str(dayNumber[i])
			dayFraction = days[dayString][:]
			month_daily_frac[i,:,:] = dayFraction

			# apply fraction to monthly data
			daily_data = month_emission * dayFraction * grid_cell_area_m2

			# Append the daily data to 'yearData' array
			yearData[jDay-1, :, :] = daily_data # -1 for python 0 based index

		# At the end of looping through each months days data, make sure the
		# daily fraction at each location adds up to 1 or 0.
		month_daily_frac_sum = np.sum(month_daily_frac, axis=0)
		# At locations not equal to zero, how different are values from 1?
		notZero = month_daily_frac_sum != 0.
		notZeroValues = month_daily_frac_sum[notZero]
	   	diff = np.abs(notZeroValues - 1.)
	   	test = diff > 1e-4
	   	if np.sum(test) > 0:
	   		print 'These is a monthly fraction sum equal to: ' + str(np.max(diff))
	   		raise ValueError('Monthly Fraction array non 0 or 1 at a location.')


	# Check for leap year, if 366 day of year is still all -1 get rid of it
	if np.unique(yearData[365,:,:])[0] == -1:
		yearData = yearData[0:365,:,:]

	# Make this a much more useful array
	time = np.array(time)

	return time, latitude, longitude, yearData, grid_cell_area_m2


################################################################################
# Get monthly emissions too.
################################################################################
def getMonthlyEmissions(dataDir, year, months, species):
	"""This function gets all the monthly burn area.

		return: time, latitude, longitude, yearData
	"""

	yearFile = 'GFED4.1s_' + str(year) + '.hdf5'
	f = h5py.File(dataDir + yearFile, 'r')

	# Get dimensions
	latitude  = f['lat'][:]
	longitude = f['lon'][:]

	# Get grid cell area [m**2]
	grid_cell_area_m2 = f['ancill/grid_cell_area'][:]

	# Create an array with the correct lat and lon dimension to append data
	# NOTE: will trim 366th day if no assignment is made
	dims = (12, latitude.shape[0], latitude.shape[1])
	emissions = np.zeros(dims) # to store emissions
	AGRI      = np.zeros(dims) # to store fraction from this type of source
	BORF      = np.zeros(dims) # ...
	DEFO      = np.zeros(dims) # ..
	PEAT      = np.zeros(dims) # .
	SAVA      = np.zeros(dims)
	TEMF      = np.zeros(dims)

	# To store burn area fraction
	area_fraction = np.zeros(dims)


	# Create a list where time string year-month can be stored
	timeString = []

	monthCount = -1
	for m in months:

		timeString.append(str(year) + m)

		monthCount = monthCount + 1
		print 'Getting ' + str(year) + ' ' + m + ' monthly '+species+' and source data'

		# Write species emissions path
		EPath = '/emissions/' + m + '/' + species
		emissions[monthCount, :, :] = f[EPath][:]

		# Get the source fraction for these emission data
		sourceBase = '/emissions/' + m + '/partitioning/' + species + '_'
		AGRI[monthCount, :, :] = f[sourceBase + 'AGRI']
		BORF[monthCount, :, :] = f[sourceBase + 'BORF']
		DEFO[monthCount, :, :] = f[sourceBase + 'DEFO']
		PEAT[monthCount, :, :] = f[sourceBase + 'PEAT']
		SAVA[monthCount, :, :] = f[sourceBase + 'SAVA']
		TEMF[monthCount, :, :] = f[sourceBase + 'TEMF']

		# Get burn area fraction
		areaPath = '/burned_area/' + m + '/burned_fraction'
		area_fraction[monthCount,:,:] = f[areaPath]

	# Make the return of these many variables easier with a dictionary
	yearData = {}
	yearData['emissions'] = emissions
	yearData['AGRI'] = AGRI
	yearData['BORF'] = BORF
	yearData['DEFO'] = DEFO
	yearData['PEAT'] = PEAT
	yearData['SAVA'] = SAVA
	yearData['TEMF'] = TEMF
	yearData['area_fraction'] = area_fraction



	timeString = np.array(timeString) # will make easier to append and handle later

	return timeString, latitude, longitude, yearData, grid_cell_area_m2


################################################################################
# Append the yearly data matricies by executing yearly data function
################################################################################
if getDaily:

	years = np.arange(startYear, endYear+1)
	for year in years:

		print 'Appending: ' + str(year)

		if year == years[0]:
			timeBase, lat, lon, dataBase, a = getDailyEmissions(dataDir, year, months, species)
		else:
			time, lat, lon, yearData, a = getDailyEmissions(dataDir, year, months, species)

			# append the new data to the existing base
			dataBase = np.append(dataBase, yearData, axis=0)
			timeBase = np.append(timeBase, time)

	# go back to the nice names
	yearlyData = dataBase
	time = timeBase

	# Create origin object that matches ecmwf era-interm
	t0 = datetime(year=1900, month=1, day=1, hour=0, minute=0, second=0)
	secondsPerHour = 60**2
	hoursFromOrigin = []
	for i in range(len(time)):
		# assume midnight on each date
		time_datetime = datetime.combine(time[i], datetime.min.time())
		# create a difftime object so we can extract an absolute diff in seconds
		diff = time_datetime - t0
		# convert difference in seconds to difference in hours
		hoursFromOrigin.append(diff.total_seconds()/secondsPerHour)

	# Make into nice array
	hoursFromOrigin = np.array(hoursFromOrigin)

	# make sure time step is ALWAYS 1 day, or something when wrong somewhere
	if len(np.unique(np.diff(time))) > 1.:
		raise ValueError('There is a missing timestep in the datamerge.')

	################################################################################
	# Write the NETCDF data. Make sure to include all relevant units for a given
	# species!
	################################################################################
	print 'Working on writing the output as netCDF data'
	nLat = lat.shape[0]
	nLon = lon.shape[1]
	nTime = len(time)

	# When the start year is the same as the end year, only assign one year for
	# file name
	if startYear == endYear:
		outputFile = dataDir + 'GFED4.1s_' + species + '_' +\
				str(startYear) + '.nc'
	else:
		outputFile = dataDir + 'GFED4.1s_' + species + '_' +\
				str(startYear) + '_' + str(endYear) + '.nc'

	ncFile = Dataset(outputFile, 'w', format='NETCDF4')
	ncFile.description = 'Data downloaded converted from hdf5 format'
	ncFile.location = 'Global'
	ncFile.createDimension('time',  nTime )
	ncFile.createDimension('latitude', nLat )
	ncFile.createDimension('longitude', nLon )

	VAR_ = ncFile.createVariable(species,\
		            'f4',('time','latitude','longitude'))

	grid_area_ = ncFile.createVariable("grid_area", 'f4', ('latitude', 'longitude'))
	grid_area_.units = 'm**2'

	if species == 'C':
		VAR_.units = 'g ' + species + ' per grid cell per day'
	if species == 'DM':
		VAR_.units = 'kg ' + species + ' per grid cell per day'

	# Create time variable
	time_ = ncFile.createVariable('time', 'i4', ('time',))
	time_.units = 'hours since ' + str(t0)

	# create lat variable
	latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
	latitude_.units = 'degrees north'

	# create longitude variable
	longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
	longitude_.units = 'degrees east'

	# Write the actual data to these dimensions
	VAR_[:]       = yearlyData[:]
	grid_area_[:] = a
	latitude_[:]  = lat[:,0]
	longitude_[:] = lon[0,:]
	time_[:]      = hoursFromOrigin[:]

	ncFile.close()




################################################################################
# Get all years mothly emissions and write the nc
################################################################################
if getMonthly:

	years = np.arange(startYear, endYear+1)
	for year in years:

		print 'Appending: ' + str(year)

		if year == years[0]:

			timeBase, lat, lon, yearData, a = getMonthlyEmissions(dataDir, year, months, species)

			area_fraction_base = yearData['area_fraction']
			emissions_base     = yearData['emissions']
			PEAT_base          = yearData['PEAT']
			TEMF_base          = yearData['TEMF']
			AGRI_base          = yearData['AGRI']
			BORF_base          = yearData['BORF']
			DEFO_base          = yearData['DEFO']
			SAVA_base          = yearData['SAVA']

		else:

			time, lat, lon, yearData, a = getMonthlyEmissions(dataDir, year, months, species)

			# append the new data to the existing base
			area_fraction_base = np.append(area_fraction_base, yearData['area_fraction'])
			emissions_base     = np.append(emissions_base, yearData['emissions'])
			PEAT_base          = np.append(PEAT_base, yearData['PEAT'])
			TEMF_base          = np.append(TEMF_base, yearData['TEMF'])
			AGRI_base          = np.append(AGRI_base, yearData['AGRI'])
			BORF_base          = np.append(BORF_base, yearData['BORF'])
			DEFO_base          = np.append(DEFO_base, yearData['DEFO'])
			SAVA_base          = np.append(SAVA_base, yearData['SAVA'])


			# Append time too, simple 1D append
			timeBase = np.append(timeBase, time)

	# Time needs to be type int in order to be stored in nc data as an array
	timeBase = np.array(timeBase, 'int')
	nLat = lat.shape[0]
	nLon = lon.shape[1]
	################################################################################
	# Write nc burn area
	################################################################################

	# When the start year is the same as the end year, only assign one year for file name
	if startYear == endYear:
		outputFile = dataDir + 'GFED4.1s_monthly_'+species+'_' +\
				str(startYear) + '.nc'
	else:
		outputFile = dataDir + 'GFED4.1s_monthly_'+species+'_' +\
				str(startYear) + '_' + str(endYear) + '.nc'

	ncFile = Dataset(outputFile, 'w', format='NETCDF4')
	ncFile.description = 'Data downloaded converted from hdf5 format'
	ncFile.location = 'Global'
	ncFile.createDimension('time',  len(timeBase) )
	ncFile.createDimension('latitude', nLat )
	ncFile.createDimension('longitude', nLon )

	# Burn area
	burn_area_fraction_ = ncFile.createVariable('burn_area_fraction',\
		            'f4',('time','latitude','longitude'))
	burn_area_fraction_.units = 'fraction of grid cell that burned'

	# Emissions
	emissions_ = ncFile.createVariable(species,\
		            'f4',('time','latitude','longitude'))
	if species == 'DM':
		emissions_.units = 'kg DM m**-2 month**-1'
	else:
		emissions_.units = 'g C m**-2 month**-1'

	# The source partition fractions
	PEAT_base_ = ncFile.createVariable('PEAT_fraction','f4',('time','latitude','longitude'))
	TEMF_base_ = ncFile.createVariable('TEMF_fraction','f4',('time','latitude','longitude'))
	AGRI_base_ = ncFile.createVariable('AGRI_fraction','f4',('time','latitude','longitude'))
	BORF_base_ = ncFile.createVariable('BORF_fraction','f4',('time','latitude','longitude'))
	DEFO_base_ = ncFile.createVariable('DEFO_fraction','f4',('time','latitude','longitude'))
	SAVA_base_ = ncFile.createVariable('SAVA_fraction','f4',('time','latitude','longitude'))

	PEAT_base_.units = 'fraction of emissions'
	TEMF_base_.units = 'fraction of emissions'
	AGRI_base_.units = 'fraction of emissions'
	BORF_base_.units = 'fraction of emissions'
	DEFO_base_.units = 'fraction of emissions'
	SAVA_base_.units = 'fraction of emissions'

	# Area
	grid_area_ = ncFile.createVariable("grid_area", 'f4', ('latitude', 'longitude'))
	grid_area_.units = 'm**2'

	# Create time variables
	time_ = ncFile.createVariable('time', 'f4', ('time',))
	time_.units = 'YYYYMM for monthly data'

	# create lat variable
	latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
	latitude_.units = 'degrees north'

	# create longitude variable
	longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
	longitude_.units = 'degrees east'

	# Write the actual data to these dimensions
	burn_area_fraction_[:]  = area_fraction_base[:]
	emissions_[:] = emissions_base[:]
	PEAT_base_[:] = PEAT_base
	TEMF_base_[:] = TEMF_base
	AGRI_base_[:] = AGRI_base
	BORF_base_[:] = BORF_base
	DEFO_base_[:] = DEFO_base
	SAVA_base_[:] = SAVA_base

	grid_area_[:] = a[:]
	latitude_[:]  = lat[:,0]
	longitude_[:] = lon[0,:]
	time_[:]      = timeBase[:]

	ncFile.close()