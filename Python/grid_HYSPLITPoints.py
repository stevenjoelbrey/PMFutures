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

# Handle strange duration quantities and put them back into df as correct numeric
# value 




######################################################################################
# Loop through each day. Loop through each point. Assign duration hours a home
# based on what ecmwf grid cell center they are closest too. 
######################################################################################
# nDates = len(uniqueDates)
# nLat = len(glat)
# nLon = len(glon)
# 
# data = np.zeros((nDates, nLat, nLon))
# 
# 
# for i in range(nTime):
# 
# 	dateMask = HMSTime == uniqueDates[i]
# 	
# 	_subset = df[dateMask]
# 




# COPY NC WRITING FORMAT USED FOR GFED PROCESSING






