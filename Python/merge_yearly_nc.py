#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to merge yearly daily era-interim nc files or GFED4s
# fire emissions files.


# import required modules
import os
import numpy as np
import sys
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import cesm_nc_manager as cnm

dirRoot = '/barnes-scratch/sbrey/'


###############################################################################
# ------------------------- Handle Args --------------------------------------- 
###############################################################################
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:
	print 'Using arguments passed via command line.'
	ncVARType =  str(sys.argv[1])
	ncVAR     =  str(sys.argv[2]) # e.g. 'z'
	subsetNA  =  str(sys.argv[3])
	startYear =  int(sys.argv[4])
	endYear   =  int(sys.argv[5])

else: 
	# Development environment. Set variables by manually here. 
	ncVARType = 'era_interim'       # 'era_interim' | 'GFED4s'
	ncVAR     = 'tp'            # 'tp' | 'C'
	subsetNA  = '_'         # "_NA_" | "_"
	startYear = 2003
	endYear   = 2016


# Set directory paths based on the type of data to be merged
if ncVARType == 'era_interim':
	dataDir = dirRoot + 'era_interim_nc_daily/'
	outDir  = dirRoot + 'era_interim_nc_daily_merged/'

elif ncVARType == 'GFED4s':
	dataDir = dirRoot + 'GFED4s/'
	outDir  = dirRoot + 'GFED4s/'

###############################################################################
# ------------------------Begin main script -----------------------------------
###############################################################################

years = np.arange(startYear, endYear+1)

for year in years:

	if ncVARType == 'era_interim':
		loadFile = dataDir + ncVAR + '_' + str(year) + '.nc'
	else: 
		loadFile = dataDir + 'GFED4.1s_METGrid_' + ncVAR + '_' + str(year) + '.nc'

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
if subsetNA  == '_NA_':
	ymin   = 10   
	ymax   = 90.  
	xmin   = -130. + 180.
	xmax   = -60.  + 180.
	varBase, latitude, longitude = cnm.mask2dims(varBase, longitude[:], latitude[:], 0,\
							xmin, xmax, ymin, ymax)

	print 'Final merged var NA size: ' + str(varBase.shape)



#######################################################################	
# Write the merged file
#######################################################################	
if ncVARType == 'era_interim':
	outputFile = outDir + ncVAR + subsetNA + \
		     str(startYear) + '_' + str(endYear) + '.nc'

elif ncVARType == 'GFED4s':
	outputFile = dataDir + 'GFED4.1s_METGrid_' + ncVAR + subsetNA + \
			str(startYear)+ '_' + str(endYear) + '.nc'	

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
elif ncVARType == 'GFED4s':
	ncFile.description = 'Merged daily GFED4s emissions from yearly files after regridding'	

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

	mergedVAR_.scale_factor = ncBase.variables[ncVAR].scale_factor
	mergedVAR_.add_offset = ncBase.variables[ncVAR].add_offset
	mergedVAR_.missing_value = ncBase.variables[ncVAR].missing_value

elif (len(dims)==3) & (ncVARType == 'era_interim'):

	mergedVAR_ = ncFile.createVariable(ncVAR,\
	            'f4',('time','latitude','longitude'),
                     fill_value=ncBase.variables[ncVAR]._FillValue)

	mergedVAR_.scale_factor = ncBase.variables[ncVAR].scale_factor
	mergedVAR_.add_offset = ncBase.variables[ncVAR].add_offset
	mergedVAR_.missing_value = ncBase.variables[ncVAR].missing_value

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





