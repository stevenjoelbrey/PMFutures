#!/usr/bin/env python2


###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# The purpose of this script is to regrid GFED data to ecmwf era-interim grid
# (or any lat lon passed) and save as a single daily nc file. 

# TODO: include 'basis_regions' in nc output?

dataDir   = '/barnes-scratch/sbrey/GFED4s/'
startYear = 2003
endYear   = 2004
species   = 'C' # 'C' , 'DM', 'small_fire_fraction' (These have daily fraction est.)

# Months to loop over
months = ['01', '02', '03', '04', '05','06',\
          '07', '08', '09', '10', '11', '12']

import h5py # if this creates an error please make sure you have the h5py library
import os
import numpy as np
import sys
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import datetime
from datetime import date
from datetime import timedelta
from datetime import datetime


def getYearlyData(dataDir, year, months, species):
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

	jDay = 0
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
		month_daily  = np.zeros((nDaysInMonth, nLat, nLon))

		for i in range(nDaysInMonth):
	
			# Advance the JDay Count (after adding dt to date0, since origin jan1)
			dt = timedelta(days=jDay)
			time.append(date0+dt)
			jDay = jDay + 1

			# Get fraction of monthly emissions that occured on THIS day
			dayString = 'day_' + str(dayNumber[i])
			dayFraction = days[dayString][:]	

			# apply fraction to monthly data
			daily_data = month_emission * dayFraction * grid_cell_area_m2

			# Append the daily data to 'yearData' array
			yearData[jDay-1, :, :] = daily_data # -1 for python 0 based index

	# Check for leap year, if 366 day of year is still all -1 get rid of it
	if np.unique(yearData[365,:,:])[0] == -1:
		yearData = yearData[0:365,:,:]

	# Make this a much more useful array 
	time = np.array(time)

	return time, latitude, longitude, yearData, grid_cell_area_m2

####################################################################
# Append the yearly data matricies by executing yearly data function 	
####################################################################
years = np.arange(startYear, endYear+1)
for year in years:
	
	print 'Appending: ' + str(year)

	if year == years[0]:
		timeBase, lat, lon, dataBase, a = getYearlyData(dataDir, year, months, species)		
	else:	
		time, lat, lon, yearData, a = getYearlyData(dataDir, year, months, species)
		
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

######################################################################################
# Write the NETCDF data. Make sure to include all relevant units for a given species! 
######################################################################################
nLat = lat.shape[0]
nLon = lon.shape[1]
nTime = len(time)

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
VAR_.units = 'g ' + species + ' per grid cell per day'

#area_ = ncFile.createVariable('grid_cell_area_m2', 'f4', ('latitude','longitude',))
#area_.units = 'meters square in grid box'

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
latitude_[:]  = lat[:,0]
longitude_[:] = lon[0,:]
time_[:]      = hoursFromOrigin[:]


ncFile.close()	

