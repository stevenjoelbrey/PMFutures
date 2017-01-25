#!/usr/bin/env python2

# createLinearModel.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to create a linear model of the relationship between
# meteorlogy parameters and wildfire emissions. 

import cesm_nc_manager as cnm
import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr
import pandas as pd

# The current goals of this script are to:
# 1) Get North American Emissions and Air Quality data onto a single grid so 
#    can be easily compared
# 2) Make a scatterplot of T vs. emissions for a given month. Find a way to 
#    visualize and report the summary of the observed relationship. 


###############################################################################
# Load Emissions grid for a selected variable
# Emissions variables all have leading 'E'
###############################################################################
EVar      = "BC"
EScenario = "RCP85" # NOTE: RCP85 all that exists for base period 2000s
EYear     = "2010"

E, EUnits, ELongName, Et, Elat, Elon = cnm.getEmissionVariableData(EVar,\
                                                                   EScenario,\
                                                                   EYear)
# Take care of propert area weighting of emissions and chance from per second
# average for the day to total for day
area              = cnm.loadEmissionGridAttributes(Elat, Elon) # m^2 of each cell
kgPerDay, kgTotal = cnm.makeTotalEmissions(E, area, Et)                                                                   
                                                                   
###############################################################################
# Load Emissions grid for a selected variable
# Meteorology varibles all have leading 'M'                                                                   
###############################################################################
MVar      = 'T'                                                                   
MScenario = '2000Base'

                                              
ncFile = cnm.makeAQNCFile(MVar, MScenario)                                                                   
M, MUnits, MLongName, Mt, Mlat, Mlon = cnm.getGroundAirQaulityData(MVar,\
                                                                   ncFile,\
                                                                   MScenario) 

# Check the sizes of the returned data arrays. Make sure dimensions and dates
# match before any analysis is considered. 

# Regridding information link:
# http://stackoverflow.com/questions/25544110/regridding-regular-netcdf-data

# Figure out a way to make the two grids comparable
# Emissions are on a smaller grid than air quality parameters which are global.
# Subset air quality grid by the four corners of the emissions grid

minLon = Elon[0]
maxLon = Elon[-1]
maxLat = Elat[-1]
minLat = Elat[0]             

startMonth = 1
endMonth   = 12

MSubset, MtSubset, MlatSubset, MlonSubset = cnm.subsetModelEmissions(M, Mt, Mlat, Mlon, startMonth, endMonth, 
                                                                     minLat, maxLat, minLon, maxLon)

# Deal with the Jan 1 day in 2011 that Met data has for some reason 
if ( (MtSubset[-1] != Et[-1]) | (MtSubset[0] != Et[0]) ):
    
    # See where the first and last dates fall    
    Mfirst_t = np.where(MtSubset == Et[0])[0][0]
    Mlast_t  = np.where(MtSubset == Et[-1])[0][0] 

    # Subset the time dimension and array directly 
    MtSubset = MtSubset[Mfirst_t:(Mlast_t+1)]
    MSubset  = MSubset[Mfirst_t:(Mlast_t+1),:,:]
    
# Make sure times and shapes align before allowing code to proceed beyond this
# point
if (MSubset.shape == E.shape):
    
    M = MSubset        ; del MSubset
    Mlat = MlatSubset  ; del MlatSubset
    Mlon = MlonSubset  ; del MlonSubset
    
    # Make sure all times are the same also 
    if np.unique(Et- MtSubset)[0] == timedelta(0):
        Mt = MtSubset ; del MtSubset
        t = Mt
    else:
        raise ValueError('Emission t(ime) and Air Quality t(ime) do not match') 
else:
    raise ValueError('Emission grid and Air Quality grids do not match')        
        

# I want an array that labels the month of each time value
nDays = len(t)
# make arrray that is month of t    
mon = np.zeros(nDays,dtype=int)
for i in range(nDays):
	mon[i] = t[i].month

###############################################################################
# Explore relationship between the two variables                                                                   
###############################################################################

# Make the Raw plot
if (False):
	# TODO: Color each point by month 
	plt.figure(figsize=(12,9))

	# Remove the plot frame lines because they are useless
	ax = plt.subplot(111)
	ax.spines["top"].set_visible(False)
	ax.spines["right"].set_visible(False)

	# Also remove tick marks from these sides
	ax.get_xaxis().tick_bottom()
	ax.get_yaxis().tick_left()

	plt.scatter(M.ravel(), kgPerDay.ravel(), marker='.', color='k')
	plt.title("Raw Scatterplot of " + EVar + ' vs. ' + MVar, fontsize=30)
	plt.xlabel(MVar + ' [' + MUnits + ']', fontsize=20)
	plt.ylabel(EVar + ' [kg/day]', fontsize=20)

	plt.savefig('testRaw.png')


###############################################################################
# Now plot daily totals for the chosedomain                                                               
###############################################################################
# to consider 
# http://stackoverflow.com/questions/17450313/summing-over-months-with-pandas

dailyDomainMMean   = np.mean(M, axis=(1,2))
dailyDomainESum    = np.sum(kgPerDay, axis=(1,2))

plt.figure(figsize=(12,9))

# Remove the plot frame lines because they are useless
ax = plt.subplot(111)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
#ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
#ax.yaxis.get_major_formatter().set_powerlimits((0, 1))
ax.yaxis.set_major_formatter(mpl.ticker.ScalarFormatter(useMathText=True, useOffset=False))

# Also remove tick marks from these sides
ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()


plt.scatter(KtoF(dailyDomainMMean), dailyDomainESum, marker='.', c=mon, s=40, lw=0)
cbar = plt.colorbar(ticks = [1,2,3,4,4,5,6,7,8,9,10,11,12])
cbar.set_label('Month', rotation=270, fontsize=20)

plt.title("Scatterplot of daily summed " + EVar + ' vs. monthly mean ' +\
          MVar, fontsize=24,  y=1.04)
plt.xlabel(MVar + ' [' + MUnits + ']', fontsize=20)
plt.ylabel(EVar + ' [kg/day]', fontsize=20)
plt.ylim([0, dailyDomainESum.max()])
plt.savefig('testDaily.png')

###############################################################################
# Also explore this mon shaded time series to explore why there are some days 
# in unexspected months that have really high emissions. Could it be Mexico? 
# Subset the domain being plotted and see if the signal remains.                
###############################################################################
plt.scatter(t, dailyDomainESum, c=mon)

###############################################################################
# Now plot monthy totals for the chosedomain, first make a handy pandas dataframe                                                                
###############################################################################
d = {'time': t, 'dailyDomainMMean': dailyDomainMMean}
        
    
                                                 
