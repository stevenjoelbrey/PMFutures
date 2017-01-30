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

# TODO: Dynamic map object range extension based on the bounds of analysis
# TODO: Save the figure in a sensible place with a sensible name. 
# TODO: Integrate into a shiny app that I can run on ozone. 


################################################################################
## ------------------------- USER INPUT ---------------------------------------- 
################################################################################
## TODO: Make sure these can work as agruments passed by the command line. 

variable   = 'BC'    # Anything listed in outputFromYellowstone/FireEmissions
startMonth = 1
endMonth   = 12
minLat     = 37    # 10.
maxLat     = 40    # 90.
minLon     = -109   # 190.
maxLon     = -102   # 320

saveName = variable + "_" + str(startMonth) + "-" + str(endMonth) + "_" + \
           str(minLat) + "-" + str(maxLat) + "N_" + str(minLon) + "-" +\
           str(maxLon) + ".pdf"
# Remove any "." in the saveName 
#saveName = str_re	

################################################################################
## ------------------------- exe envire ---------------------------------------- 
################################################################################
if minLon < 0:
	minLon = minLon + 360.
	maxLon = maxLon + 360.

# Create loop[ that gives us output we want for plotting and saves
# everything into python dictionaries.
kgTotals     = {}; kgTotals_max     = []
tSeries      = {}; tSeries_max      = []
sMonthTotals = {}; sMonthTotals_max = []
ts           = {}
for scenario in ['RCP45', 'RCP85']:
    for year in ['2050', '2100']:

        # Key names
        NAME = scenario+year        
        
        # Load specified variable
        bb, bbUnits, bbLongName, t, lat, lon = cnm.getVariableData(variable, scenario, year)
#
#        # Subset variable
#        bb, t, lat, lon = cnm.subsetModelEmissions(bb, t, lat, lon, 
#                                               startMonth, endMonth, 
#                                               minLat, maxLat,
#                                               minLon, maxLon)
#        ts[year] = t
#                                              
#        # Get Associated area mask
#        area = cnm.loadEmissionGridAttributes(lat, lon)
#
#        # Summarize the data for plotting 
#        kgPerDay, kgTotal = cnm.makeTotalEmissions(bb, area, t)
#        kgTotals[NAME]    = kgTotal
#        kgTotals_max.append(kgTotal.max())
#        
#        s                  = np.sum(kgPerDay, (1,2)) 
#        tSeries[NAME]      = s 
#        tSeries_max.append(s.max())
#        
#        sMonthTotals[NAME] = cnm.monthlyTotals(s, t)
#        sMonthTotals_max.append(sMonthTotals[NAME].max())
#        
## Reduce max lists to single values
#kgTotals_max     = max(kgTotals_max)
#tSeries_max      = max(tSeries_max)
#sMonthTotals_max = max(sMonthTotals_max)
#
################################################################################
## north america plot setup area 
## Set up the projection etc. Things that can be outside time loop
## TODO: Make a map projection that adjust to tghe limts passed above.
## TODO: Will probably have to go with a square projection. 
################################################################################
#
#
##m = Basemap(width=9000000, height=6000000,
##            rsphere=(6378137.00, 6356752.3142),\
##            resolution='l',area_thresh=10000.,projection='lcc',\
##            lat_1=45.,lat_2=55,lat_0=50,lon_0=-105.)
#
#m = Basemap(projection='cyl', llcrnrlat=minLat, urcrnrlat=maxLat,\
#            llcrnrlon=minLon, urcrnrlon=maxLon, lat_ts=20, resolution='c')
#
#lons, lats = np.meshgrid(lon, lat) # get lat/lons of ny by nx evenly space grid.
#x, y = m(lons, lats) # compute map proj coordinates.
#
#fig = plt.figure(figsize=(12,12))
#
#ax1 = fig.add_subplot(421)
#cnm.makePcolorFig(ax1, m, x, y, z=kgTotals['RCP452050'], 
#              titleText='RCP45 2050', maxVal=kgTotals_max)
#
#ax2 = fig.add_subplot(422)
#cnm.makePcolorFig(ax2, m, x, y, z=kgTotals['RCP452100'], 
#              titleText='RCP45 2100', maxVal=kgTotals_max)
#
#ax3 = fig.add_subplot(423)
#cnm.makePcolorFig(ax3, m, x, y, z=kgTotals['RCP852050'], 
#              titleText='RCP85 2050', maxVal=kgTotals_max)
#
#ax4 = fig.add_subplot(424)
#cnm.makePcolorFig(ax4, m, x, y, z=kgTotals['RCP852100'], 
#              titleText='RCP85 2100', maxVal=kgTotals_max)
#
#ax5 = fig.add_subplot(425)
#cnm.timeSeries(ax5, tSeries['RCP452050'], tSeries['RCP852050'], 
#           s1Label='RCP45', s2Label='RCP85', t=ts['2050'], maxValue=tSeries_max, titleText='2050')
#
#ax6 = fig.add_subplot(426)
#cnm.timeSeries(ax6, tSeries['RCP452100'], tSeries['RCP852100'], 
#           s1Label='RCP45', s2Label='RCP85', t=ts['2100'], maxValue=tSeries_max, titleText='2100')
#
#ax7 = fig.add_subplot(427)
#cnm.makeHist(ax7, sMonthTotals['RCP452050'], sMonthTotals['RCP852050'],
#         s1Label='RCP45', s2Label='RCP85', maxValue=sMonthTotals_max,titleText='2050')
#
#ax8 = fig.add_subplot(428)
#cnm.makeHist(ax8, sMonthTotals['RCP452100'], sMonthTotals['RCP852100'],
#         s1Label='RCP45', s2Label='RCP85', maxValue=sMonthTotals_max, titleText='2100')
#
#
## Handle the amount of space between plots 
#plt.subplots_adjust(wspace=0.5, hspace=0.5)
#
#
## Make overarching title 
#plt.suptitle(bbLongName + ' Summary Figure')
#
#plt.show()
#plt.savefig('../Figures/'+saveName)






