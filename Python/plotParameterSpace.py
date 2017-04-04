#!/usr/bin/env python2

# plotEmissionSummary.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to make scatterplots that allow us to explore the
# relationship between different parameters and the emissions. 


# Load resources
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

EScenario = "RCP85" # NOTE: RCP85 2010 correponds to 2000Base
EYear     = "2010"  # NOTE: For emissions this gives 2000-2010 decade
AScenario = "2000Base" 

# Get Particulate species emissions
BC, BCunits, BCLongName, t, ELat, ELon = cnm.getEmissionVariableData("BC",EScenario,EYear)
OC, OCunits, OCLongName, t, ELat, ELon = cnm.getEmissionVariableData("OC",EScenario,EYear)

# Combine emissions
secondsPerDay = 24. * 60**2
E = (OC + BC) * secondsPerDay # units of E now kg/m2/day

################################################################################
#------------------ Subset model emissions in space and time ------------------
################################################################################
# For report use western U.S. only

startMonth = 6
endMonth   = 9
minLat     = 30.    # 10.
maxLat     = 50.    # 90.
minLon     = 234.   # 190.
maxLon     = 259.   # 320


E, t, ELat, ELon  =  cnm.subsetModelEmissions(E, t, ELat, ELon,\
                                              startMonth, endMonth,\
                                              minLat, maxLat, minLon, maxLon)

# Set up a relevant map to use later
# Set up the projection etc. Things that can be outside time loop
m = Basemap(projection='merc',llcrnrlat=minLat,urcrnrlat=maxLat,\
            llcrnrlon=minLon,urcrnrlon=maxLon,lat_ts=1,resolution='c')

lons, lats = np.meshgrid(ELon, ELat)
x, y = m(lons, lats)

# Load grid cell area once final domain of emissions analysis is set
# TODO: Look into run time warning
area, landMask = cnm.loadEmissionGridAttributes(ELat, ELon)





# Get the met parameters we want to explore the parameter space of 

TFile = cnm.makeAQNCFile('T', AScenario, 'daily')
T, TUnits, TLongName, t, lat, lon = cnm.getGroundAirQaulityData("T", TFile, AScenario)

VFile = cnm.makeAQNCFile('V', AScenario, 'daily')
V, VUnits, VLongName, vt, lat, lon = cnm.getGroundAirQaulityData("V", VFile, AScenario)

UFile = cnm.makeAQNCFile('U', AScenario, 'daily')
U, UUnits, ULongName, ut, lat, lon = cnm.getGroundAirQaulityData("U", UFile, AScenario)

windSpeed = np.sqrt(V**2 * U**2)

RHFile = cnm.makeAQNCFile('RELHUM', AScenario, 'daily')
RH, RHUnits, RHLongName, RHt, lat, lon = cnm.getGroundAirQaulityData("U", RHFile, AScenario)


# Subset the parameters to the same limited time and space domain as the emissions
TSubset, tSubset, latSubset, lonSubset  =  cnm.subsetModelEmissions(T, t, lat, lon,\
                                              startMonth, endMonth,\
                                              minLat, maxLat, minLon, maxLon)

windSpeedSubset, tSubset, latSubset, lonSubset  =  cnm.subsetModelEmissions(windSpeed, t, lat, lon,\
                                              startMonth, endMonth,\
                                              minLat, maxLat, minLon, maxLon)


# TODO: Plot all parameters in space with abline for event thresholds
# Make a version of all these plots that shows the points that make up 99% of total 
# emissions
fig = plt.figure(figsize=(10,10))
c = plt.scatter(windSpeedSubset, TSubset, c = E, marker=".", lw = 0)
cb = plt.colorbar(c, size="5%", pad="0%")
plt.show()

# Plot row of all against wind speed
 









