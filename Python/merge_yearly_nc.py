#!/usr/bin/env python2

###############################################################################
# ------------------------- Description ---------------------------------------
###############################################################################
# This script will be used to merge yearly daily era-interim nc files or GFED4s
# fire emissions files.

# Follows ----------------------------------------
# 	- average6HourlyData.py | process_FINNY.py
# Precedes ----------------------------------------
#	- any function that reads in multi-year nc data of the global domain or any
#     domain subset.

# import required modules
import os
import numpy as np
import sys
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import cesm_nc_manager as cnm



# Figure out what machine this code is running on. Set file paths.
drive = cnm.getDrive()
dataDirBase = drive + "era_interim_nc_daily_merged/"
dirRoot = drive


###############################################################################
# ------------------------- Handle Args ---------------------------------------
###############################################################################
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

# TODO: Add met grid as a variable forth paths, for now all ecmwf

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	ncVARType =  str(sys.argv[1])
	ncVAR     =  str(sys.argv[2]) # e.g. 'z'
	region    =  str(sys.argv[3])
	startYear =  int(sys.argv[4])
	endYear   =  int(sys.argv[5])


else:
	# Development environment. Set variables by manually here.
	ncVARType = 'era_interim'   # 'era_interim' | 'GFED4s' | 'FINN'
	ncVAR     = 'v10'           # 'tp' | 'C' | 'CO2' | 'SPDH'
	region    = '_'             #  "_" = global | any region in cnm.getRegionBounds()
	startYear = 2003            # HMS only spans 2006 - 2015
	endYear   = 2016


# Set directory paths based on the type of data to be merged
if ncVARType == 'era_interim':
	dataDir = dirRoot + 'era_interim_nc_daily/'
	outDir  = dirRoot + 'era_interim_nc_daily_merged/'

elif ncVARType == 'GFED4s':
	dataDir = dirRoot + 'GFED4s/'
	outDir  = dirRoot + 'GFED4s/'

elif ncVARType == 'FINN':
	dataDir = dirRoot + 'FINN/'
	outDir  = dirRoot + 'FINN/'

elif ncVARType == 'HMS':
	dataDir = dirRoot + 'HMS/'
	outDir  = dirRoot + 'HMS/'

###############################################################################
# ----------------------- Begin main script -----------------------------------
###############################################################################

years = np.arange(startYear, endYear+1)

for year in years:

	if ncVARType == 'era_interim':
		loadFile = dataDir + ncVAR + '_' + str(year) + '.nc'
	elif ncVARType == 'GFED4s':
		loadFile = dataDir + 'GFED4.1s_ecmwf_' + ncVAR + '_' + str(year) + '.nc'
	elif ncVARType == 'FINN':
		loadFile = dataDir + 'FINN_ecmwf_' + ncVAR + '_' + str(year) + '.nc'
	elif ncVARType == 'HMS':
		loadFile = dataDir + 'HYSPLITPoints_ecmwf_' + ncVAR + '_' + str(year) + '.nc'

	# Use the first years data to get dimension information
	# and array to be appended too.
	if year == years[0]:

		print 'Creating base data and dimensions with year = ' + str(year)

		ncBase = Dataset(loadFile, 'r')
		varBase = ncBase.variables[ncVAR]

		# Get dimension information one time
		latitude = ncBase.variables['latitude']
		longitude = ncBase.variables['longitude']
		latitudeUnits = latitude.units
		longitudeUnits = longitude.units

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
# Handle making the spatial subset
#######################################################################
if region  != '_':

	# Get the chosen region bounds
	minLat, maxLat, minLon, maxLon, resolution = cnm.getRegionBounds(region)

	# Subset the data based on the bounds of this region.
	varBase, latitude, longitude = cnm.mask2dims(varBase, longitude[:], latitude[:], 0,\
							                     xmin=minLon, xmax=maxLon,\
							                     ymin=minLat, ymax=maxLat)

	print 'Data subset to region ' + region



#######################################################################
# Write the merged file
#######################################################################
if ncVARType == 'era_interim':
	outputFile = outDir + ncVAR + region + \
		     str(startYear) + '_' + str(endYear) + '.nc'

elif ncVARType == 'GFED4s':
	outputFile = dataDir + 'GFED4.1s_ecmwf_' + ncVAR + region + \
			str(startYear)+ '_' + str(endYear) + '.nc'

elif ncVARType == 'FINN':
	outputFile = dataDir + 'FINN_ecmwf_' + ncVAR + region + \
			str(startYear)+ '_' + str(endYear) + '.nc'

elif ncVARType == 'HMS':
	outputFile = dataDir + 'HMS_ecmwf_' + ncVAR + region + \
			str(startYear)+ '_' + str(endYear) + '.nc'

else:
	raise ValueError('ncVARType unknown. Please use known ncVARType.')


print '----------------------------------------------------------------------'
print 'Writing the output file. outputFile used:'
print outputFile
print '----------------------------------------------------------------------'

# Get lengths of dimensions that are always present
nDays = len(tBase)
nLat = len(latitude)
nLon = len(longitude)


ncFile = Dataset(outputFile, 'w', format='NETCDF4')
if ncVARType == 'era_interim':
	ncFile.description = 'Daily data created by average6HourlyData.py and merged with merge_year_nc.py'
else:
	ncFile.description = 'Merged daily' + ncVARType + 'emissions from yearly files after regridding'


ncFile.location = 'Global'
ncFile.createDimension('time',  nDays )
ncFile.createDimension('latitude', nLat )
ncFile.createDimension('longitude', nLon )

# Create variables on the dimension they live on
if len(dims) == 4:

	ncFile.createDimension('level',  len(level) )
	mergedVAR_ = ncFile.createVariable(ncVAR,\
	            'f4',('time','level','latitude','longitude'),
                     fill_value=ncBase.variables[ncVAR]._FillValue)

	# While here create the level dimesion
	level_ = ncFile.createVariable('level', 'i4', ('level',))
	level_.units = level.units


elif (len(dims)==3) & (ncVARType == 'era_interim'):

	mergedVAR_ = ncFile.createVariable(ncVAR,\
	            'f4',('time','latitude','longitude'))

else:
	mergedVAR_ = ncFile.createVariable(ncVAR,\
	            'f4',('time','latitude','longitude'))

# Assign the same units as the loaded file to the main variable
mergedVAR_.units = ncBase.variables[ncVAR].units

# Create time variable
time_ = ncFile.createVariable('time', 'i4', ('time',))
time_.units = ncBase.variables['time'].units
if ncVARType == 'era_interim':
	time_.calendar = ncBase.variables['time'].calendar

# create lat variable
latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
latitude_.units = latitudeUnits

# create longitude variable
longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
longitude_.units = longitudeUnits

# Write the actual data to these dimensions
mergedVAR_[:] = varBase[:]
latitude_[:]  = latitude[:]
longitude_[:] = longitude[:]
time_[:] = tBase[:]
if len(dims) == 4:
	level_[:] = level[:]

ncFile.close()





