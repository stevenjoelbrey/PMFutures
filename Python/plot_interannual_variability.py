#!/usr/bin/env python2

# plot_interannual_variability.py 

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to explore the monthly and seasonal relationships
# between emissions and meteorology. 

# currently, this script contains much of the same functionality and load 
# statements as plotObservationParameterSpace.py and 
# plot_GFED4s_emission_summary.py. These will be consolidated at a later time
# but I am separating here to make organization more clear. 


################################################################################
#------------- Arguments to Subset model emissions in space and time -----------
################################################################################

cutoffPercentile = 80.
startMonth = 6
endMonth   = 9
region     = "_west_" # "_west_"| "_PNW_" | "_CAL_" | "_NorthRockies_" 


# Load resources
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

# Get region lat lon range	
minLat, maxLat, minLon, maxLon  = cnm.getRegionBounds(region)

# Figure out what machine this code is running on
pwd = os.getcwd()
mac = '/Users/sbrey/Google Drive/sharedProjects/PMFutures/Python'
if pwd == mac:
	dataDirBase = "/Volumes/Brey_external/"
else:
	dataDirBase = "/barnes-scratch/sbrey/"