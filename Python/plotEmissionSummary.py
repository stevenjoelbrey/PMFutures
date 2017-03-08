#!/usr/bin/env python2

# plotEmissionSummary.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to summarize and plot emissions from CESM. 
# A Multipanel plot will sumarize the present day, 2050, 2100 emissions for a 
# selected variable. Later on, this script will be able to make the standard 
# plots for desired time and season subsets. Much later, this code will be
# incorperated into a shiny app that will allow easy exploration of the model 
# output. 
# TODO: Make this into a class. Copy the stucture of Sam Atwoods code. 
# TODO: Update function documentation to be the accepted format. 

# As needed for the April 1 EPA report. 

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

################################################################################
# Load particulate species emissions 
################################################################################
EmissionScenario = "RCP85"
EmissionYear     = "2010"

# Get Particulate species emissions
BC, BCunits, BCLongName, t, ELat, ELon = cnm.getEmissionVariableData("BC","RCP85","2010")
OC, OCunits, OCLongName, t, ELat, ELon = cnm.getEmissionVariableData("OC","RCP85","2010")

# Combine emissions
secondsPerDay = 24. * 60**2
E = OC + BC

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

# Load area onces final domain of emissions analysis is set
area = cnm.loadEmissionGridAttributes(ELat, ELon)

################################################################################
#----------------------- Load desired met event masks --------------------------
################################################################################

highWindMask = cnm.getSelf('2000Base', 'highWindMask_9')
highTMask    = cnm.getSelf('2000Base', 'highTMask_298')
lowPrecMask  = cnm.getSelf('2000Base', 'lowPRECTMask_0.01')
#stagMask     = cnm.getSelf('2000Base', 'stagnation_mask') # needs to be rerun 

ncFile = cnm.makeAQNCFile('stagnation_mask','2000Base', 'daily')
nc     = Dataset(ncFile, 'r')
stagMask = nc.variables['stagnationMask'][:]
AQLat  = nc.variables['lat'][:]
AQLon  = nc.variables['lon'][:]
nc.close()

# Get mask dates so we can match with emissions
dateNum      = cnm.getSelf('2000Base', 'date') 
mTime        = cnm.dateNumToDate(dateNum)

################################################################################
# Subset each of the masks to match the emissions area
################################################################################
highWindMask, mt, mLat, mLon  =  cnm.subsetModelEmissions(
                                              highWindMask, mTime, AQLat, AQLon,\
                                              startMonth, endMonth,\
                                              minLat, maxLat, minLon, maxLon)

highTMask, mt, mLat, mLon  =  cnm.subsetModelEmissions(
                                              highTMask, mTime, AQLat, AQLon,\
                                              startMonth, endMonth,\
                                              minLat, maxLat, minLon, maxLon)

lowPrecMask, mt, mLat, mLon  =  cnm.subsetModelEmissions(
                                              lowPrecMask, mTime, AQLat, AQLon,\
                                              startMonth, endMonth,\
                                              minLat, maxLat, minLon, maxLon)

stagMask, mt, mLat, mLon  =  cnm.subsetModelEmissions(
                                              stagMask, mTime, AQLat, AQLon,\
                                              startMonth, endMonth,\
                                              minLat, maxLat, minLon, maxLon)

################################################################################
# Count the number of events at each grid cell for each type 
################################################################################
nHighWindGrid = np.sum(highWindMask, axis=0)
nHighTGrid    = np.sum(highTMask,    axis=0)
nLowPrecGrid  = np.sum(lowPrecMask,  axis=0)
nStagMaskGrid = np.sum(stagMask,     axis=0)

nHighWind = np.sum(highWindMask)
nHighT    = np.sum(highTMask)
nLowPrec  = np.sum(lowPrecMask)
nStagMask = np.sum(stagMask)


# Make masked arrays for each emission day type
highWindE = ma.masked_where(highWindMask == 1, E)
HighTE    = ma.masked_where(highTMask == 1, E)
LowPrecE  = ma.masked_where(lowPrecMask == 1, E)
stagE     = ma.masked_where(stagMask == 1, E)
# NOTE: This implies that E is all summer E

# We also need to mask where emissions are zero 
highWindE = ma.masked_where(highWindE == 0, highWindE)
HighTE    = ma.masked_where(HighTE == 0, HighTE)
LowPrecE  = ma.masked_where(LowPrecE == 0, LowPrecE)
stagE     = ma.masked_where(stagE == 0, stagE)
E         = ma.masked_where(E == 0, E)

# Flatten all of these values, trailing _ indicates flattened to 1D
highWindE_ = ma.compressed(highWindE)
HighTE_    = ma.compressed(HighTE)
LowPrecE_  = ma.compressed(LowPrecE)
stagE_     = ma.compressed(stagE)
E_         = ma.compressed(E)


boxLabel = ['All days', 'High Wind', 'High T', 'No Precip', 'Stagnation']

# Make the kg/s/m2 raw plot
data = [E_, highWindE_, HighTE_, LowPrecE_, stagE_]

fig = plt.figure()
ax = plt.subplot(1, 1, 1)
ax.boxplot(data, sym="", labels=boxLabel, whis=[10, 90])
ax.set_ylabel(BCunits)
ax.set_yscale('log')
plt.show()

# Now sum over the days dimensions
highWindE_days = np.sum(highWindE, 0)
HighTE_days    = np.sum(HighTE, 0)
LowPrecE_days  = np.sum(LowPrecE, 0)
stagE_days     = np.sum(stagE, 0)
E_days         = np.sum(E, 0)

# Flatten and make a boxplot again
highWindE_days_ = ma.compressed(highWindE_days)
HighTE_days_    = ma.compressed(HighTE_days)
LowPrecE_days_  = ma.compressed(LowPrecE_days)
stagE_days_     = ma.compressed(stagE_days)
E_days_         = ma.compressed(E_days)

data = [ E_days_, highWindE_days_, HighTE_days_, LowPrecE_days_, stagE_days_]

fig = plt.figure()
ax = plt.subplot(1, 1, 1)
ax.boxplot(data,  sym="", labels=boxLabel, whis=[10, 90])
ax.set_yscale('log')
plt.show()

# TODO Make a map of the total emissions by day type! 


