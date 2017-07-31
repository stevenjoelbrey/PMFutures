#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to generate daily met fields from hourly nc files. 
# The 6-houly data to average live in /barnes-scratch/sbrey/era_interim_nc_6_hourly



###############################################################################
# ---------------------- Set analysis variables--------------------------------
###############################################################################
import sys
 
	# Development environment. Set variables by manually here. 
hourlyVAR =  'z'


dataDir = "/barnes-scratch/sbrey/era_interim_nc_6_hourly"
outputDir = "/barnes-scratch/sbrey/era_interim_nc_daily"
	

# Import the required modules
import os
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset 
import matplotlib.pyplot as pl
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr
import datetime
import time as timer
import os.path
import cesm_nc_manager as cnm

ncFile = dataDir + '/' + 'z_all.nc'

nc = Dataset(ncFile, 'r')
allTime = nc.variables['time'][:]




