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

# Load grid cell area once final domain of emissions analysis is set
# TODO: Look into run time warning
area, landMask = cnm.loadEmissionGridAttributes(ELat, ELon)


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

# cyclone & blocking days coming soon (I hope)



# Now we only care about the values over land. So mask out areas over water. 
# This is important for parameters like "highWindMask", which mainly shows up
# over the ocean, which is not particularly helpful. 
# TODO: Figure out how to do this without using a for loop
highWindMask = ma.core.MaskedArray(highWindMask)
highTMask    = ma.core.MaskedArray(highTMask)
lowPrecMask  = ma.core.MaskedArray(lowPrecMask)
stagMask     = ma.core.MaskedArray(stagMask)

for i in range(len(mt)):
	highWindMask[i,:,:]= ma.masked_where(landMask == 0, highWindMask[i,:,:])
	highTMask[i,:,:]   = ma.masked_where(landMask == 0, highTMask[i,:,:])
	lowPrecMask[i,:,:] = ma.masked_where(landMask == 0, lowPrecMask[i,:,:])
	stagMask[i,:,:]    = ma.masked_where(landMask == 0, stagMask[i,:,:])

################################################################################
# Count the number of events at each grid cell for each type 
################################################################################
year = []
for i in range(len(mt)):
	year.append(mt[i].year)

uniqueYears = np.unique(year)
nYears = len(uniqueYears)
daysPerSummer = np.sum(mt < date(uniqueYears[1],1,1))
n = daysPerSummer * nYears * 1.

# Sum across time dimension TODO: and make % June-September
nHighWindGrid = np.sum(highWindMask, axis=0) / n
nHighTGrid    = np.sum(highTMask,    axis=0) / n
nLowPrecGrid  = np.sum(lowPrecMask,  axis=0) / n
nStagMaskGrid = np.sum(stagMask,     axis=0) / n

# Count total occurance
nHighWind = np.sum(highWindMask)
nHighT    = np.sum(highTMask)
nLowPrec  = np.sum(lowPrecMask)
nStagMask = np.sum(stagMask)

################################################################################
# Map the occurance of each type of met event 
################################################################################

fig = plt.figure(figsize=(8,10))

ax = fig.add_subplot(2,2,1)

m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighWindGrid, cmap="Reds", vmin=0., vmax=.1 )
cbar = m.colorbar(c, location='bottom', pad="1%", extend='max', ticks=[0, 0.05, .1])
cbar.set_label('proportion of days')
plt.title('Days Wind > 9 m/s')

ax = fig.add_subplot(2,2,2)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighTGrid, cmap="Purples", vmin=0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days T > 25 C')

ax = fig.add_subplot(2,2,3)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowPrecGrid, cmap="Purples", vmin=0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days with precip < 0.01 inches')

ax = fig.add_subplot(2,2,4)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nStagMaskGrid, cmap="Purples", vmin = 0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Stagnation Days')

plt.savefig('../Figures/metEventMaskCounts_' + AScenario + '.png')
plt.close()

################################################################################
# Now we want to make the exact same plot, only we want to mask ot the locations
# where emissions are low. We want to look at met events where emissions for 
# fires are high. Once those locations are known. Show exact same plot again
# only where emissions are high 
################################################################################
ETotal = np.sum(E, axis=0) * area # kg/grid

#ENew = np.zeros(E.shape)
#ENew = ma.core.MaskedArray(ENew)
#for i in range(len(mt)):
#	ENew[i, :, :] = (E[i,:,:] * area)
#
#ENewTotal = np.sum(ENew, axis=0)


fig = plt.figure(figsize=(8,8))
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, ETotal, cmap="OrRd")
cbar = m.colorbar(c, location='bottom',pad="1%")
cbar.set_label('BC + OC kg')
plt.title('Total Emissions, June-Sept 2000-2010', fontsize=24)
plt.savefig('../Figures/decadeEmissions_'+ AScenario + '.png')

# Pick a value that boxes must exceed if they are going to be considering
# TODO: Consider whether it would be more meaningful to look at the
# TODO: boxes that account for 50% of the emissions, will be fewer. 
cutoff    = np.percentile(ETotal, 75)
highEMask = ETotal > cutoff
ETotalSum = np.sum(ETotal)

# Estimate proportion of total emissions these account for after masking out
# the low emission amounts 
highETotal = ma.masked_where(highEMask==False, ETotal)
highEmittersSum = np.sum(highETotal)
percentTotalERetained = highEmittersSum / ETotalSum * 100. 
print 'Percent of domain emissions retained with chosen mask: ' + str(percentTotalERetained)


fig = plt.figure(figsize=(8,8))
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, highETotal, cmap="OrRd")
cbar = m.colorbar(c, location='bottom',pad="1%")
cbar.set_label('BC + OC kg')
plt.title('Total Emissions, June-Sept 2000-2010', fontsize=24)
plt.savefig('../Figures/decadeHighEmitters_'+ AScenario + '.png')

# Mask the met event totals by this high emission mask 
nHighWindGrid = ma.masked_where(highEMask==False, nHighWindGrid)
nHighTGrid    = ma.masked_where(highEMask==False, nHighTGrid)
nLowPrecGrid  = ma.masked_where(highEMask==False, nLowPrecGrid)
nStagMaskGrid = ma.masked_where(highEMask==False, nStagMaskGrid)

# Now do it with the emissions in time, we do not ever want to look at the 
# locations where the values are this low 
for i in range(len(mt)):
	E[i,:,:] = ma.masked_where(highEMask==False, E[i,:,:])

# Plot the met even counter again now that they are masked 
fig = plt.figure(figsize=(8,10))

ax = fig.add_subplot(2,2,1)

m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighWindGrid, cmap="Reds", vmin=0., vmax=.1 )
cbar = m.colorbar(c, location='bottom', pad="1%", extend='max', ticks=[0, 0.05, .1])
cbar.set_label('proportion of days')
plt.title('Days Wind > 9 m/s')

ax = fig.add_subplot(2,2,2)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nHighTGrid, cmap="Purples", vmin=0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days T > 25 C')

ax = fig.add_subplot(2,2,3)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nLowPrecGrid, cmap="Purples", vmin=0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Days with precip < 0.01 inches')

ax = fig.add_subplot(2,2,4)
m.drawcoastlines(linewidth=1)
m.drawstates(linewidth=1)
m.drawcountries(linewidth=1)
c = m.pcolor(x, y, nStagMaskGrid, cmap="Purples", vmin = 0, vmax=1 )
cbar = m.colorbar(c, location='bottom',pad="1%", ticks=[0, 0.5, 1])
cbar.set_label('proportion of days')
plt.title('Stagnation Days')

plt.savefig('../Figures/metEventMaskCountsHighE_' + AScenario + '.png')

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

# Flatten all of these values, trailing _ indicates flattened to 1D, removes 
# masked values 
highWindE_ = ma.compressed(highWindE)
HighTE_    = ma.compressed(HighTE)
LowPrecE_  = ma.compressed(LowPrecE)
stagE_     = ma.compressed(stagE)
E_         = ma.compressed(E)

# First make a multihist of these met events. 
bins = np.linspace(0, E_.max(), 100)
x = np.cumsum(np.diff(bins))

# Count the amount in each emission bin for each type of met event
E_heights         = plt.hist(E_, bins, label='E | No meteorology event', alpha=0.2)[0]
HighTE_heights    = plt.hist(HighTE_, bins, label='E | high T', alpha=0.2)[0]
highWindE_heights = plt.hist(highWindE_, bins, label = 'E | high wind', alpha=0.2)[0]
LowPrecE_heights  = plt.hist(LowPrecE_, bins, label = 'E | low precip', alpha=0.2)[0]
stagE_heights     = plt.hist(stagE_, bins, label = 'E | stagnation', alpha=0.2)[0]


# Create data structure to plot nice side by side by histograms
fig = plt.figure(figsize=(10,10))
lw=2
ax= fig.add_subplot(1,1,1)
plt.subfigure(1,1)
plt.plot(x, E_heights, label='No meteorology event', linewidth=lw)
plt.plot(x, HighTE_heights, label='E | high T', linewidth=lw)
plt.plot(x, highWindE_heights, label='E | high wind', linewidth=lw)
plt.plot(x, LowPrecE_heights, label='E | low precip', linewidth=lw)
plt.plot(x, stagE_heights, label='E | stagnation', linewidth=lw)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.xaxis.set_ticks_position('bottom')
ax.yaxis.set_ticks_position('left')
plt.xlim([0, 0.00004])
plt.xlabel('Kg m$^{-2}$', fontsize=24)
plt.ylabel('Total grid boxe days with given emission', fontsize=24)

plt.legend(loc='upperright')

plt.savefig('../Figures/lineHistogram_'+AScenario+'.pdf')





# Now normalize the occurance of these events and plot the lines
# of the heights of the histogram. _n stands for normalized.
def normalizeE(heights):
	heights_n = heights / np.sum(heights)
	return heights_n

E_heights_n         = normalizeE(E_heights)
HighTE_heights_n    = normalizeE(HighTE_heights)
highWindE_heights_n = normalizeE(highWindE_heights)
LowPrecE_heights_n  = normalizeE(LowPrecE_heights)
stagE_heights_n     = normalizeE(stagE_heights) 




fig = plt.figure(figsize=(8,6))

plt.plot(x, E_heights_n, label='E | No meteorology event', linewidth=3)
plt.plot(x, HighTE_heights_n, label='E | high T', linewidth=3)
plt.plot(x, highWindE_heights_n, label = 'E | high wind', linewidth=3)
plt.plot(x, LowPrecE_heights_n, label = 'E | low precip', linewidth=3)
plt.plot(x, stagE_heights_n, label = 'E | stagnation', linewidth=3)
#plt.yscale('log', nonposy='clip')
plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
plt.xlim([0, x[30]])

plt.xlabel('Emissions (kg m$^{-2}$ day$^{-1}$)', fontsize=18)
plt.ylabel('Proportion of event days', fontsize=18)
plt.legend(loc='upper right')
plt.legend(frameon=False)

plt.savefig('../Figures/EmissionProbabilityGivenMet.pdf')

plt.close()


FACTOR = 1e14

# Make kernal density estimates for these arrays of values. 
mu = np.mean(E_)
sigma = np.std(E_)
Xn =  (E_ - mu) / sigma

# Try normalizing a second one with this

Xn = E_ * FACTOR

from sklearn.neighbors.kde import KernelDensity
X = Xn[:, np.newaxis]
kde = KernelDensity(kernel='tophat', bandwidth=0.1).fit(X)

# Make a range of values to predict density on
xRange = np.linspace(Xn.min(), Xn.max(), 100)[:, np.newaxis]
log_dens = kde.score_samples(xRange)

dens = np.exp(log_dens)

plt.plot(xRange, dens)


xRangeBack = xRange * sigma + mu

plt.plot(xRangeBack, dens)
plt.xlabel('No event emissions (kg/m$^{2}$/day)')
plt.ylabel('Density Estimate')


# Try a poor mans estimate
#height = plt.hist([highWindE_, HighTE_], 
#				  bins=np.linspace(E_.min(), E_.max(),1e6)))


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
















