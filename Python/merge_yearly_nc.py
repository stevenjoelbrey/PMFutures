#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to merge yearly daily era-interim nc files 


# import required modules
import os
import numpy as np
import sys
from netCDF4 import Dataset 
import matplotlib.pyplot as plt

import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	ncVAR     =  str(sys.argv[1]) # e.g. 'z'
	startYear =  int(sys.argv[2])
	endYear   =  int(sys.argv[3])

else: 
	# Development environment. Set variables by manually here. 
	ncVAR     =  'tp'
	startYear = 1997
	endYear   = 1998


dataDir = '/barnes-scratch/sbrey/era_interim_nc_daily/'
outDir  = '/barnes-scratch/sbrey/era_interim_nc_daily_merged/'


#def mergeSurfaceData(var, startYear=1997, endYear=2016):
#""""if x1 is (time1, lat, lon) and x2 is (time2, lat, lon) we want to create
#    y such that y is (time1:time2, lat, lon)
#"""

years = np.arange(startYear, endYear+1)

for year in years:


	loadFile = dataDir + ncVAR + '_' + str(year) + '.nc'

	# Use the first years data to get dimension information 
	# and array to be appended too.
	if year == years[0]:

		print 'Creating base data and dimensions with year = ' + str(year)

		ncBase = Dataset(loadFile, 'r')
		varBase = ncBase.variables[ncVAR]
		# TODO: change to 'latitude' and 'longitude' once daily values are
		# TODO: rerun 
		latitude = ncBase.variables['latitude']
		longitude = ncBase.variables['longitude']
		tBase = ncBase.variables['time']

		dims = ncBase.dimensions.keys()
		if len(dims) == 4:
			level = ncBase.variables['level']

	else:

		print 'appending: ' + ncVAR + ' ' + str(year)	
		
		nc  = Dataset(loadFile, 'r')
		var = nc.variables[ncVAR]
		t   = nc.variables['time']

		# Append this nc file to the array loaded by ncBase, join
		# via the time axis, which is always 0 with era-interim data
		varBase = np.append(varBase, var, axis=0)
		tBase = np.append(tBase, t)


if (len(tBase)) != (varBase.shape[0]):
	ValueError('Time dimensions of data array and time array do not match!')

if np.unique(np.diff(tBase)) != 24:
	ValueError('Some time jump not equal to 24 hours present!')

print 'Final merged var size: ' + str(varBase.shape)
print 'Final merged time array: ' + str(len(tBase))
	
#######################################################################	
# Write the merged file
#######################################################################	
outputFile = outDir + ncVAR + '_' +\
             str(startYear) + '_' + str(endYear) + '.nc'

print '----------------------------------------------------------------------'
print 'Writing the output file. outputFile used:'
print outputFile
print '----------------------------------------------------------------------'

# Get lengths of dimensions that are always present 
nDays = len(tBase)
nLat = len(latitude)
nLon = len(longitude)


ncFile = Dataset(outputFile, 'w', format='NETCDF4')
ncFile.description = 'Daily data created by average6HourlyData.py and merged with merge_year_nc.py'
ncFile.location = 'Global'
ncFile.createDimension('time',  nDays )
ncFile.createDimension('latitude', nLat )
ncFile.createDimension('longitude', nLon )

# Create variables on the dimension they live on 
# TODO: preserve fill_value
if len(dims) == 4:

	ncFile.createDimension('level',  len(level) )
	mergedVAR_ = ncFile.createVariable(ncVAR,\
	            'f4',('time','level','latitude','longitude'),
                     fill_value=ncBase.variables[ncVAR]._FillValue)

	# While here create the level dimesion
	level_ = ncFile.createVariable('level', 'i4', ('level',))
	level_.units = level.units

else:
	mergedVAR_ = ncFile.createVariable(ncVAR,\
	            'f4',('time','latitude','longitude'),
                     fill_value=ncBase.variables[ncVAR]._FillValue)

# Assign the same units as the loaded file to the main variable
mergedVAR_.units = ncBase.variables[ncVAR].units
mergedVAR_.scale_factor = ncBase.variables[ncVAR].scale_factor
mergedVAR_.add_offset = ncBase.variables[ncVAR].add_offset
mergedVAR_.missing_value = ncBase.variables[ncVAR].missing_value

# Create time variable	
time_ = ncFile.createVariable('time', 'i4', ('time',))
time_.units = ncBase.variables['time'].units
time_.calendar = ncBase.variables['time'].calendar

# create lat variable
latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
latitude_.units = latitude.units

# create longitude variable
longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
longitude_.units = longitude.units

# Write the actual data to these dimensions
mergedVAR_[:] = varBase[:]
latitude_[:]  = latitude[:]
longitude_[:] = longitude[:]
time_[:] = tBase[:]
if len(dims) == 4:
	level_[:] = level[:]

ncFile.close()	





