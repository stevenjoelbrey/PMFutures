#!/usr/bin/env python2

# met_event_finder_m.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This module is the home of functions used to identify meteorogy events 
# relevent for air quality and climate change impact assesment in the EPA
# "Planning for an uncertain future" project (dubbed PMFutures). 

# TODO: Deal with output that has been placed on pressure coordinates? 


import cesm_nc_manager as cnm
import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr


def detectStringInList(string, some_list):
    matching = [s for s in some_list if string in s]   
    return matching

def getGroundAirQaulityData(variable, ncFile, scenario):
    """This function loads nc data of given variable and scenario into
    the working environment using numpy arrays.
        Parameters:
            variable:   The AQ variable to be loaded (T, U, V, Precip?)
            ncFile:     The path of the ncFile to be loaded
            scenario:   The RCP scenario of the data to be loaded
            
        return: M, MUnits, MLongName, t, lat, lon     
    """
        
    nc = Dataset(ncFile, 'r')
    
    # Need to determine if variable on ilev or lev dimension
    ncDims = nc.variables[variable].dimensions
    test   = detectStringInList('lev', ncDims)
    if test == ['lev']:
        levelDim = 'lev'
    elif test == ['ilev']:
        levelDim = 'ilev'
    elif test == ['Pressure']: 
        levelDim = 'Pressure'
    else:
        print 'unknown dimension type for variable '
    
    # Find the index of the levelDim that is closest to surface
    print 'levelDim being used is : ' + levelDim    
    print ncDims    
    
    lev = nc.variables[levelDim][:]
    surface_index = np.where(lev == lev.max())[0][0]
    
    print 'loading large file.....'
    # Get the data and the units 
    M         = nc.variables[variable][:, surface_index, :, :] 
    MUnits    = nc.variables[variable].units
    MLongName = nc.variables[variable].long_name
    print 'done loading large file'
    
    # get the spatial dims
    lat = nc.variables['lat'][:]
    lon = nc.variables['lon'][:]
    
    nc.close()
    
    # Get the time dimension 
    ncTimeFile = cnm.makeAQNCFile('date', scenario)
    ncTime = Dataset(ncTimeFile, 'r')
    date = ncTime['date'][:] # Gets date string
    ncTime.close()
    
    t = cnm.dateNumToDate(date) # converts date string to date object 
    
    return M, MUnits, MLongName, t, lat, lon 
    

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
    
def getSurfaceWindScaler(scenario):
    """This function uses getGroundAirQaulityData() to load u and v winds and
    make a scaler wind for threshold analysis threshold. 
        Parameters:
            scenario: The year and RCP scenario for the wind data 
    """
    
    uFile = cnm.makeAQNCFile("U", scenario)
    u, MUnits, MLongName, t, lat, lon = getGroundAirQaulityData("U", \
                                                                uFile,\
                                                                scenario)
    
    vFile = cnm.makeAQNCFile("V", scenario)
    v, MUnits, MLongName, t, lat, lon = getGroundAirQaulityData("V", \
                                                                vFile,\
                                                                scenario)
    
    # Now make a scaler wind from the vector coomponents
    scaler = np.sqrt(v**2 + u**2)

    return scaler, MUnits, 'Scaler Wind', t, lat, lon                                                            
                                                                

###############################################################################
# ---------------- Environment for testing functions -------------------------- 
###############################################################################
variable = 'T'
scenario = '2000Base'


#ncFile = cnm.makeAQNCFile(variable, scenario)
#M, MUnits, MLongName, t, lat, lon = getGroundAirQaulityData(variable, ncFile, scenario)

#wind, wUnits, wName, t, lat, lon = getSurfaceWindScaler(scenario)
