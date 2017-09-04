#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script is used to grid daily HMS points to GFED4s grid. This will be 
# regridded to an ecmwf grid soon after but needs to be on the GFEDF4s grid
# first for sanity checks. 


import sys # for reading command line arguments
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, cm, shiftgrid
import numpy as np
import scipy.interpolate
import time as timer
import pandas as pd
import cesm_nc_manager as cnm
import datetime

# Regridding degrees tolerance
TOLERANCE = 0.25/2.

print 'Loading all data, setting the fields...'

# Load GFED4s grid to assign data to.
ncFile = '/barnes-scratch/sbrey/GFED4s/GFED4.1s_C_2003_2016.nc'
nc     = Dataset(ncFile, 'r')

glat = nc.variables['latitude'][:]
glon = nc.variables['longitude'][:]
grid_area = nc.variables['grid_area'][:]
# time = nc.variables['time']
# 
# # Make time useful (nc file connection required)
# t, month, year = cnm.get_era_interim_time(time)

nc.close()


######################################################################################
# Load HMS data
######################################################################################
f = '/barnes-scratch/sbrey/HMS/hysplitPoints_land_all.csv'
# want Dur column preserved as a string because of leading zeros
df = pd.read_csv(f, dtype={'Dur':'str'}) 
nRow = df.shape[0]
hlat = df.Lat
hlon = df.Lon	

######################################################################################
# Process time field (make hours from origin) and DUR (make hour)
######################################################################################
# NOTE: make a nice date object for masking, also handle dur. 
# NOTE: DUR is format HHMM, we need to make this into a numeric hour quantity. 
# NOTE: Some of these values are '****', these will be ignored. 
HMSTime = [None] * nRow
HMSHour = np.zeros(nRow) 

# Will use same origin as GFED data, unit in hours, calling universal for project
t0_universal = datetime.datetime(year=1900, month=1, day=1, hour=0, minute=0, second=0)

# Make a place to store formatted DUR values
DUR = np.zeros(nRow)

# Loop over each row processing the time and 
for i in range(nRow):
	
	YYYYMMDD = str(df.YearMmDd[i])
	year_  = int(YYYYMMDD[0:4]) 
	month_ = int(YYYYMMDD[4:6])
	day_   = int(YYYYMMDD[6:8])		
	
	t_ = datetime.datetime(year=year_, month=month_, day=day_,\
						   hour=0, minute=0, second=0)
	HMSTime[i] = t_					  
						  
	# Create and store desired origin version of time array. Find out how many hows
	# this date is from origin date.           
	dt = t_ - t0_universal 
	dt_hours = dt.total_seconds() / 60**2	  
	HMSHour[i] = dt_hours		   

	thisDur = df.Dur[i]
	# If it is not nan, get the value, otherwise leave as zero
	# len needs to be 4 for following code to work
	if (type(thisDur) != float) & (thisDur != "****"):
		if (len(thisDur) == 4):
			# Make HHMM into a numeric hourly quantity
			HH = thisDur[0:2]
			MM = thisDur[2:4]
			hourFrac = float(MM)/60.
			DurNew = float(HH) + hourFrac
		
			DUR[i] = DurNew
	
HMSTime = np.array(HMSTime)					   
df['DUR'] = DUR	
# These are the unique date we are going to loop over and save out for HMS data 	
uniqueDates = np.unique(HMSTime)

# Skip the first date that occurs in 2003. Preserve all other dates
uniqueDates_new = uniqueDates[1::]
uniqueDates = uniqueDates_new

# This is our new dimension for dates that will be used as the nc dimension
nDates = len(uniqueDates)

# We need to save out the unique hour starts for the unique dates to define the
# time dimension for the nc data 
uniqueHours = np.zeros(nDates)
for i in range(nDates):
	dt = uniqueDates[i] - t0_universal 
	dt_hours = dt.total_seconds() / 60**2	  
	uniqueHours[i] = dt_hours	

######################################################################################
# Loop through each day. Loop through each point. Assign duration hours a home
# based on what ecmwf grid cell center they are closest too. 
######################################################################################
nLat = len(glat)
nLon = len(glon)

# Place to store the new gridded HMS data 
data = np.zeros((nDates, nLat, nLon))

# Loop through the HYSPLIT points one at a time 
for i in range(nRow): 

	dur = df.DUR[i]

	# Subset the data for this date
	dateMask = uniqueDates == HMSTime[i]
	
	# Figure out if this is a non-zero dur value to assign, otherwise, leave data
	# array as zero, unchanged, since adding zero does nothing 
	if dur != 0:
		
		fire_lon = df.Lon[i]
		fire_lat = df.Lat[i]
		
		dx = np.abs(fire_lon - glon)
		dy = np.abs(fire_lat - glat)
	
		# Index of best match
		xi = np.argmin(dx)
		yi = np.argmin(dy)
		
		# Make sure the gridded index actually matches the HMS point location 
		dx_test = glon[xi] - fire_lon
		dy_test = glat[yi] - fire_lat

		# Assign the HYSPLIT point duration to the correct location and date
		if (dx_test <= TOLERANCE) & (dy_test <= TOLERANCE):
			data[dateMask, yi, xi] = data[dateMask, yi, xi] + dur
		else:
			print 'dx mismatch = ' + str(dx_test)
			print 'dy mismatch = ' + str(dy_test)
			raise ValueError('HYSPLIT Point location assignment not within tolerance.')

	percentComplete = np.round(float(i)/nRow*100., 5)
	print str(percentComplete) + '% complete gridding'



##########################################################################################
# Save with names matching GFED as closely as possible, e.g., latitude not lat for dim 
# name 
##########################################################################################
print 'Working on writing converted data to large nc file...'
	
fout = '/barnes-scratch/sbrey/HMS/HYSPLITPoints_SPDH.nc'
	
ncFile = Dataset(fout, 'w', format='NETCDF4')
ncFile.description = 'Gridded HYSPLIT point duration data. See Brey et al. 2017 for details'
ncFile.location = 'Global'
ncFile.createDimension('time',  nDates )
ncFile.createDimension('latitude', nLat )
ncFile.createDimension('longitude', nLon )

VARDIMS = ('time','latitude','longitude')

SPDH_ = ncFile.createVariable('SPDH','f4', VARDIMS)
SPDH_.units = 'smoke production duration hours'	            
	
grid_area_ = ncFile.createVariable("grid_area", 'f4', ('latitude', 'longitude')) 	            
grid_area_.units = 'm**2'
		
# Create time variables	
time_ = ncFile.createVariable('time', 'f4', ('time',))
time_.units = 'hours since 1900-01-01 00:00:00'

# create lat variable
latitude_ = ncFile.createVariable('latitude', 'f4', ('latitude',))
latitude_.units = 'degrees north'

# create longitude variable
longitude_ = ncFile.createVariable('longitude', 'f4', ('longitude',))
longitude_.units = 'degrees east'

# Write the actual data to these dimensions
SPDH_[:]      = data
grid_area_[:] = grid_area
latitude_[:]  = glat
longitude_[:] = glon
time_[:]      = uniqueHours


ncFile.close()