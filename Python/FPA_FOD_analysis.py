#!/usr/bin/env python2

# The purpose of this script is to explore the relationships between summer met,
# burn area and FPA-FOD and GFED4s burn area. 

import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import datetime
import matplotlib.ticker as tkr
import cesm_nc_manager as cnm
import pandas as pd
from scipy.stats import pearsonr
import scipy.stats as stats
import statsmodels

region = "_west_"
minLat, maxLat, minLon, maxLon, resolution  = cnm.getRegionBounds(region)


# 
# dataFile = "/home/sbrey/projects/PMFutures/Dataframes/fireOccurrence_t.csv"
# df = pd.read_csv(dataFile)
# 
# 
# ###############################################################################
# # TODO: Make this faster... This is stupid...
# ###############################################################################
# # 
# # t = ['none']*len(startDate)
# # year = ['none']*len(startDate)
# # month = ['none']*len(startDate)
# # day = ['none']*len(startDate)
# # for i in range(len(startTime)):
# # 	DATE = datetime.datetime.strptime(startDate[i], "%Y-%m-%d")
# # 	t[i] = DATE
# # 	year[i] = DATE.year
# # 	month[i] = DATE.month
# # 	day[i] = DATE.day
# # 	
# # t = np.array(t)
# # year = np.array(year)
# # month = np.array(month)
# # day = np.array(day)
# # 
# # df['datetime'] = t
# # df['year'] = year
# # df['month'] = month
# # df['day'] = day
# # df.to_csv('/home/sbrey/projects/PMFutures/Dataframes/fireOccurrence_t.csv')
# 
# # subset by dates of interest
# dateMask = (df.year >= 2003) & (df.year <= 2016)
# df_subset = df[dateMask]
# 
# # TODO: convert latitude longitude
# lon = df_subset.LONGITUDE
# lonAdjust = lon + 360.
# 
# # TODO: subset by space of interest
# lonMask = (lonAdjust >= minLon) & (lonAdjust <= maxLon)
# latMask = (df_subset.LATITUDE >= minLat) & (df_subset.LATITUDE <= maxLat)
# spatialMask =  lonMask & latMask
# 
# df_subset =  df_subset[spatialMask]
# 
# ###############################################################################
# # Make summer summary dataframes
# ###############################################################################
# summerBlank = np.zeros(shape=(nYears, len(colnames)))
# FOD_summer_df   = pd.DataFrame(data=summerBlank, columns=colnames)


#!/usr/bin/env python2

# plot_interannual_variability.py 

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to explore the monthly and seasonal relationships
# between emissions (TODO: Burn area) and meteorology. 

# currently, this script contains much of the same functionality and load 
# statements as plotObservationParameterSpace.py and 
# plot_GFED4s_emission_summary.py. These will be consolidated at a later time
# but I am separating here to make organization more clear. 


################################################################################
#------------- How well correlated are the emissions of different regions ------
################################################################################


# Get region lat lon range and basemap plotting resolution based on the 
# chosen region
# TODO: Make a SE region! 
minLat, maxLat, minLon, maxLon, resolution  = cnm.getRegionBounds(region)

# Figure out what machine this code is running on
pwd = os.getcwd()
mac = '/Users/sbrey/Google Drive/sharedProjects/PMFutures/Python'
if pwd == mac:
	drive = "/Volumes/Brey_external/"
else:
	drive = "/barnes-scratch/sbrey/"

# Set directory paths
metDataDirBase = drive + "era_interim_nc_daily_merged/"
figureDir = '../Figures/GFED_era_interm_analysis/' # always relativ to py
	
# Get emissions, use this to get dimensions
ncFile  = drive + "GFED4s/GFED4.1s_METGrid_C_2003_2016.nc"
nc = Dataset(ncFile, 'r')
latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
time = nc.variables['time'][:]
C = nc.variables['C'][:]
nc.close()

# Make time into datetime arrays
time, month, year = cnm.get_era_interim_time(time)

def regionDailyTimeSeries(C, region):
	""""
	Returns a summed time series the length of axis 0 for C subet by the passed
	region name
	"""
	minLat, maxLat, minLon, maxLon, resolution  = cnm.getRegionBounds(region)
	C_PNW, ynew, xnew = cnm.mask2dims(C, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
	C_PNW_daily = np.sum(C_PNW, axis=(1,2))
	return C_PNW_daily
	
C_PNW_daily	= regionDailyTimeSeries(C, "_PNW_")
C_Rocky_daily	= regionDailyTimeSeries(C, "_CentralRockies_")
C_CAL_daily	= regionDailyTimeSeries(C, "_CAL_")




figure()
plt.plot(time, C_PNW_daily, label="PNW", linewidth=2)
plt.plot(time, C_CAL_daily, label="CAL", linewidth=2)
plt.plot(time, C_Rocky_daily, label="Central Rockies", linewidth=2)

plt.legend()

