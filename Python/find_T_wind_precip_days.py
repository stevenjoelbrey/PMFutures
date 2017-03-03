#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script is used to identify high wind days, precipitation days, and
# high temperature days. 

import sys
print 'Number of arguments:', len(sys.argv), 'arguments.'
print 'Argument List:', str(sys.argv)

# TODO: Make it possible to pass absolute thresholds. 

if len(sys.argv) != 1:

	threshType = sys.argv[1]
	scenario   = sys.argv[2]
	print 'Using arguments passed via command line.'

	if threshType == 'usePercentile':

		TThresh    = sys.argv[3] # percentile
		WindThresh = sys.argv[4] # percentile
		PRECThresh = sys.argv[5] # percentile

	if threshType == 'useValue':

		TThresh    = sys.argv[3] # K
		WindThresh = sys.argv[4] # m/s
		PRECThresh = sys.argv[5] # inches/day

else:

	print 'Using defualt arguments. None passed via command line.'
	threshType = 'usePercentile'
	TThresh    = 95 # percentile
	WindThresh = 95 # percentile
	PRECThresh = 95 # percentile

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
import cesm_nc_manager as cnm
import time as timer

startTime = timer.time()

# Make temperature, wind, and precip file connections
U = cnm.getSelf(scenario, "U")
V = cnm.getSelf(scenario, "V")
T = cnm.getSelf(scenario, "T")
Prec = cnm.getSelf(scenario, "T")

# COPPIED FROM cesm_nc_manager.py, review to make sure it makes sense and use
# it as needed to create these three masks. 
def findHighValueDays(M, t, percentile):
    """This function accepts a numeric array (t, lat, lon) and returns
    an equal size array of 1s and 0s. The 1s indicate dates where M exceeded
    monthly percentile threshold for that location and 0s dates when they did
    not.
        Parameters:
            M:           Met variuable array to be passed. format (t, lat, lon)
            t:           datetime.date array that defines time diemsion of T
            percentile:  value between 0 and 100 that determines value used for
                         a given month and grid cells threshold. 
                         
        return: 
            mask:       True and False values indicating if temperature for day 
                        and location exceeds percentile threshold value for 
                        that month. True equals high T days. 
            threshold:  An array (month, lat, lon) that contains the percentile
                        threshold value for each month and location used to 
                        make 'mask'.    
    """
    # set the size of output based on the size of passed temperature array
    nDays = len(t)
    nLat = M.shape[1]
    nLon = M.shape[2]    
    
    # make arrray that is month of t    
    mon = np.zeros(nDays,dtype=int)
    for i in range(nDays):
        mon[i] = t[i].month
        
    # monthly threshold array
    threshold = np.zeros((12, nLat, nLon))
    mask      = np.zeros(M.shape, dtype=bool) 
    
    # Loop through each month of the year and each grid cell to find
    # that months locations percentile thresh value. Create a mask.
    for i in range(12):
        month = i + 1
        monthMask = np.where(mon == month)[0]       
        
        for LAT in range(nLat):
            for LON in range(nLon):
                M_cell = M[monthMask, LAT, LON]
                threshValue = np.percentile(M_cell, percentile)
                threshold[i, LAT, LON] = threshValue 
                mask[monthMask, LAT, LON] = M_cell >= threshValue
                
    return mask, threshold









