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

# NOTE: In this scrip E will generally stand for "emissions", m for mask,
# NOTE: and _ will represent a flattened multidimensional array. 

# TODO: Double check time handling functions for emissions and airquality to 
# TODO: make sure the times actually overlap. 

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
# Select and load particulate species emissions and AirQualityData Masks
################################################################################
EScenario = "RCP85" # NOTE: RCP85 2010 correponds to 2000Base
EYear     = "2010"

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

# Load area once final domain of emissions analysis is set
area = cnm.loadEmissionGridAttributes(ELat, ELon)

################################################################################
#----------------------- Load desired met event masks --------------------------
################################################################################

highWindMask = cnm.getSelf(AScenario, 'highWindMask_9')
highTMask    = cnm.getSelf(AScenario, 'highTMask_298')
lowPrecMask  = cnm.getSelf(AScenario, 'lowPRECTMask_0.01')
#stagMask     = cnm.getSelf('2000Base', 'stagnation_mask') # needs to be rerun for 2000Base

ncFile = cnm.makeAQNCFile('stagnationMask',AScenario, 'daily')
nc     = Dataset(ncFile, 'r')
stagMask = nc.variables['stagnationMask'][:]
AQLat  = nc.variables['lat'][:]
AQLon  = nc.variables['lon'][:]
nc.close()

# Get mask dates so we can match with emissions
dateNum      = cnm.getSelf(AScenario, 'date') 
mTime        = cnm.dateNumToDate(dateNum)

################################################################################
# Subset each of the masks to match the emissions area and seasonality
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

# Count total occurance
nHighWind = np.sum(highWindMask)
nHighT    = np.sum(highTMask)
nLowPrec  = np.sum(lowPrecMask)
nStagMask = np.sum(stagMask)

################################################################################
# Map the occurance of each type 
################################################################################

fig = plt.figure(figsize=(10,10))

ax = fig.add_subplot(2,2,1)

m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, ma.masked_where(nHighWindGrid == 0, nHighWindGrid) )
plt.colorbar(c)
plt.title('Days Wind > 9 m/s')

ax = fig.add_subplot(2,2,2)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, ma.masked_where(nHighTGrid == 0, nHighTGrid) )
plt.colorbar(c)
plt.title('Days T > 25 C')

ax = fig.add_subplot(2,2,3)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, ma.masked_where(nLowPrecGrid == 0, nLowPrecGrid) )
plt.colorbar(c)
plt.title('Days with precip < 0.01 inches')

ax = fig.add_subplot(2,2,4)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, ma.masked_where(nStagMaskGrid == 0, nStagMaskGrid) )
plt.colorbar(c)
plt.title('Stagnation Days')

plt.savefig('../Figures/metEventMaskCounts_' + AScenario + '.png')


################################################################################
# Make masked arrays for each emission day type. DO NOT show where mask == 0
################################################################################
highWindE = ma.masked_where(highWindMask == 0, E)
HighTE    = ma.masked_where(highTMask == 0, E)
LowPrecE  = ma.masked_where(lowPrecMask == 0, E)
stagE     = ma.masked_where(stagMask == 0, E)
# NOTE: This implies that E is all summer E

# We also need to mask where emissions are zero, because zero is non-interesting 
highWindE = ma.masked_where(highWindE == 0, highWindE)
HighTE    = ma.masked_where(HighTE == 0, HighTE)
LowPrecE  = ma.masked_where(LowPrecE == 0, LowPrecE)
stagE     = ma.masked_where(stagE == 0, stagE)
E         = ma.masked_where(E == 0, E)

################################################################################
# TODO: Map total emissions for each met event type 
################################################################################
# cnm.makeTotalEmissions


# Flatten all of these values, trailing _ indicates flattened to 1D
highWindE_ = ma.compressed(highWindE)
HighTE_    = ma.compressed(HighTE)
LowPrecE_  = ma.compressed(LowPrecE)
stagE_     = ma.compressed(stagE)
E_         = ma.compressed(E)

medianValues = [np.percentile(E_, 50), 
				np.percentile(highWindE_, 50),
				np.percentile(HighTE_, 50),
				np.percentile(LowPrecE_, 50),
				np.percentile(stagE_ , 50)]

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

# Calculate total emissions for each grid cell by multiplying kg/m^2 by area
highWindEkg = highWindE_days * area
HighTEkg    = HighTE_days    * area
LowPrecEkg  = LowPrecE_days  * area
stagEkg     = stagE_days     * area
Ekg         = E_days         * area

################################################################################
# Make Box plot 
################################################################################
boxLabel = ['All Summer days', 'Wind > 9 m/s', 'T > 25C', 'Precip < 0.01', 'Stagnation']

# Make the kg/day/m2 raw plot
boxData = [E_, highWindE_, HighTE_, LowPrecE_, stagE_]

fig = plt.figure(figsize=(7,12))

ax = fig.add_subplot(3, 1, 1)
bpRaw = ax.boxplot(boxData, sym="",labels=['','','','',''], whis=[10, 90])
ax.set_ylabel("BC + OC kg/m$^{2}$/day")
plt.xticks(rotation=45)
ax.set_yscale('log')
plt.title('Black Carbon and Organic Carbon Emissions')
fig.subplots_adjust(bottom=0.3) # Increase bottom margin

# Change the units to kg/m
dayBoxData = [ E_days_, highWindE_days_, HighTE_days_, LowPrecE_days_, stagE_days_]

ax = fig.add_subplot(3, 1, 2)
bpDay = ax.boxplot(dayBoxData,  sym="", labels=['','','','',''], whis=[10, 90])
plt.xticks(rotation=45)
ax.set_yscale('log')
ax.set_ylabel('kg/m$^{2}$')
fig.subplots_adjust(bottom=0.3)


kgPerCell = [Ekg, highWindEkg, HighTEkg, LowPrecEkg, stagEkg]

ax = fig.add_subplot(3, 1, 3)
bpkg = ax.boxplot(kgPerCell,  sym="", labels=boxLabel, whis=[10, 90])
plt.xticks(rotation=45)
ax.set_yscale('log')
ax.set_ylabel('kg')
fig.subplots_adjust(bottom=0.3)

plt.savefig('../Figures/metEmissionBoxplot_' + AScenario + '.pdf')


# TODO Make a map of the total emissions by day type! 
highWindTotal = np.sum(highWindEkg)
HightTTotal   = np.sum(HighTEkg)
LowPrecTotal  = np.sum(LowPrecEkg)
stagETotal    = np.sum(stagEkg)


plt.figure(figsize=(6,4))

index = np.arange(4)
rects1 = plt.bar(left=index ,height= [highWindTotal, HightTTotal, LowPrecTotal, stagETotal], 
                 tick_label=['High Wind', 'Hight T', 'Low Precip', 'Stagnation'])

plt.savefig('../Figures/metEmissionTotals_' + AScenario + '.pdf')
















