#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script is used to identify high wind days, precipitation days, and
# high temperature days in era_interim data.

import sys
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

print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

if len(sys.argv) != 1:

	print 'Using arguments passed via command line.'
	TThresh     = sys.argv[1] # K
	WindThresh  = sys.argv[2] # m/s
	PRECTThresh = sys.argv[3] # inches/day

else:

	print 'Using defualt arguments. None passed via command line.'
	TThresh     = 95 # 298 # K, 90th percentile of all 20-60 North days daily temeperate
	WindThresh  = 95 # 9  #m/s = 20 mph 1hr, NWS fire weather warning
	PRECTThresh = 0.01  # 0.01  # percentile
    
dataDirBase= os.path.join("/barnes-scratch","sbrey")


startTime = timer.time()

# Make temperature, wind, and precip file connections
U = cnm.getSelf(findHighValueDays, scenario, "U")
V = cnm.getSelf(findHighValueDays, scenario, "V")
T = cnm.getSelf(findHighValueDays, scenario, "T")
PRECT = cnm.getSelf(findHighValueDays, scenario, "PRECT")

# Use T to get dimension information
TFile  = cnm.makeAQNCFile(findHighValueDays, 'T', scenario, 'daily')
Tnc    = Dataset(TFile, 'r')
t      = Tnc.variables['time'][:]
lat    = Tnc.variables['lat'][:]
lon    = Tnc.variables['lon'][:]
lev    = Tnc.variables['lev'][:]
Tnc.close()

# Get the surface index from lev. This is index used to make masks.
# NOTE: U, V, and T all on lev NOT ilev
srf_index = np.where(lev == lev.max())[0][0]

####################################################
# Get all variables massive arrays into useful forms
####################################################

# convert precip to inches per day 
PRECT = cnm.metersPerSecToMetersPerDay(PRECT)

print '-----------------------------------------------------------------'
print 'working on getting full size surface arrays'
print '-----------------------------------------------------------------'
T = T[:, srf_index, :, :]
print 'done with T'

U = U[:, srf_index, :, :]
print 'Done with U'

V = V[:, srf_index, :, :]
print 'Done with V'

print 'Done getting wind, temperature, and precip variables into environment'

# Create surface wind mag
windMag = np.sqrt(U**2 + V**2)


# Get a more friendly time dimension that is useful for precentile
# maximum values
Time     = cnm.getSelf(dataDirBase, scenario, "date")
dateTime = cnm.dateNumToDate(Time)
nTime    = len(Time)
nLon     = len(lon)
nLat     = len(lat)

# We need to identify the threshold value and mask for each variable. 
# Using a single value is the simple method. Using a percentile value
# for a given grid cell is much more complicated and still in development
if threshType == 'usePercentile':

	print 'using percentile thresholds for masking'

	TMask, TLimVals = cnm.findHighValueDays(T,\
										    dateTime,\
			                                TThresh)

	TMask = np.array(TMask, dtype=int)

	WindMask, WindLimVals = cnm.findHighValueDays(windMag,\
						                          dateTime,\
			                                      WindThresh)
	WindMask = np.array(WindMask, dtype=int)

	PRECTMask, PrecLimVals = cnm.findHighValueDays(PRECT,\
						       					   dateTime,\
					                               PRECTThresh)

	PRECTMask = np.array(PRECTMask, dtype=int)

elif threshType == 'useValue':

	print 'using hard value thresholds for masking'

	# Predefine the mask arrays
	TMask     = np.zeros((len(t), len(lat), len(lon)), dtype=int)
	WindMask  = TMask
	PRECTMask = TMask

	print 'Creating surface temperature Mask'
	# Mask high temperature days for this date
	TMask = np.array(T >= TThresh, dtype=int)

	print 'Creating srf windMag Mask'
	WindMask = np.array(windMag >= WindThresh, dtype=int)

	print 'Create the low precip day Mask'
	PRECTMask = np.array(PRECT <= PRECTThresh, dtype=int)


######################################################################
# Write all three masks with detailed description of chosen thresholds 
######################################################################
masksDict = {}
masksDict['T']     = TMask
masksDict['Wind']  = WindMask
masksDict['PRECT'] = PRECTMask

for var in masksDict.keys():

	if var!='PRECT':
		condition = 'high'+var

	else:
		condition = 'low'+var

	saveString = condition +'Mask_'+str(threshDict[var])

	# connect this descriptive name to the directory of the scenario and var
	saveName = cnm.makeAQNCFile(dataDirBase, saveString, scenario, 'daily')

	ncFile = Dataset(saveName, 'w', format='NETCDF4')
	ncFile.description = 'Mask indicating ' + condition + ' condition.'
	ncFile.location    = 'Global'
	ncFile.createDimension('time', nTime )
	ncFile.createDimension('lat', nLat )
	ncFile.createDimension('lon', nLon )

	# Create variables on the dimension they live on 
	maskVar = ncFile.createVariable(saveString, 'i', ('time','lat','lon'))
	maskVar.units = '1=True, 0=False. ' + 'Limit used: ' + str(threshDict[var])

	time_var = ncFile.createVariable('time', 'i4', ('time',))
	time_var.units = 'days from origin'

	latitude = ncFile.createVariable('lat', 'f4', ('lat',))
	latitude.units = 'degrees north'
	 
	longitude = ncFile.createVariable('lon', 'f4', ('lon',))
	longitude.units = 'degrees east'

	# Write the actual data to these dimensions
	maskVar[:]       = masksDict[var]
	latitude[:]      = lat
	longitude[:]     = lon
	time_var[:]      = t

	ncFile.close()




