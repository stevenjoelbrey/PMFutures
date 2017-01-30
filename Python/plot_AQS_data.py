#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 12 11:21:09 2016

@author: sbrey
"""

###############################################################################
# -------------------------Description  
###############################################################################
# This script will be used to calculate and visualize the difference in Air
# quality as putput from the CESM122 


###############################################################################
# User arguments to determine what is done by this script 
###############################################################################
# PS is surface pressure 

baseScenario = '2000Base'
deltaScenario = '2050RCP85'
dataDirBase = '/fischer-scratch/sbrey/outputFromYellowstone/AirQualityData/'
fileBase = 'cesm122_fmozsoa_f09f09_2000_fires_00.'
fileTail = '.daily.200001-201012.nc'


###############################################################################

import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
#import pickle
import numpy.ma as ma
from datetime import date
import matplotlib.ticker as tkr


os.chdir('/home/sbrey/projects/PMFutures/Python')




###############################################################################
# Load the model output 
# TODO: Load model dimensions in a dynamic way 
# TOD" how to read a subset of a variable 
# TODEO: convert time (days since origin ) to math friendly no leapyear
###############################################################################

# Get the date information and make python specific    
fileList = os.listdir(dataDirBase+baseScenario)
# TODO: Learn more about how to detect specific string patterns in Lists
baseScenarioDateF  = dataDirBase + baseScenario + '/' + fileBase + 'date' + fileTail
baseDate = Dataset(baseScenarioDateF).variables['date'][:]
# TODO: Convert to nice python time object, ideally without a for loop

variable = 'HEIGHT' 
baseFile  = dataDirBase + baseScenario + '/' + fileBase + variable + fileTail
ncFile = Dataset(baseFile, 'r')
#deltaFile = dataDirBase + deltaScenario + '/' + fileBase + variable + fileTail
ncFile.variables # in unicode i.g. leading u 


lat = ncFile.variables[u'lat'][:]
lon = ncFile.variables[u'lon'][:]

# latitude lower and upper index
lati = np.argmin( np.abs(lat - 40.) ) 
loni = np.argmin( np.abs( lon - 100.) ) 

# atmosphere_hybrid_sigma_pressure_coordinate
ilev = ncFile.variables[u'ilev'][:]

# Determine what pressure level you want to look for blocking events 
surface_index=np.argmin(np.abs(500 - ilev))


var = ncFile.variables[variable][0, surface_index, :, :]


ncFile.close()





    
###############################################################################
# north america plot setup area 
###############################################################################

# Set up the projection etc. Things that can be outside time loop
m = Basemap(width=8000000,height=6500000,
            rsphere=(6378137.00, 6356752.3142),\
            resolution='l',area_thresh=10000.,projection='lcc',\
            lat_1=45.,lat_2=55,lat_0=50,lon_0=-105.)

parallels = np.arange(0.,90,15.)
meridians = np.arange(180.,360.,15.)

lons, lats = np.meshgrid(lon, lat) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.


###############################################################################
# Specify and create subset of data for plotting monthly gridded values 
###############################################################################


# Create the figure 
fig = plt.figure(figsize=(6, 6))

m.drawcoastlines()
m.drawstates()
m.drawcountries()

m.drawparallels(parallels,labels=[1,0,0,0], fontsize=10)
m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)

cs = m.pcolor(x, y, var, cmap='Reds')

con = m.contour(x,y, var, 8,colors='blue', linewidths=2)


# add colorbar.
cbar = m.colorbar(cs, location='bottom',pad="5%")
cbar.set_label(variable, fontsize=18)
# add title
plt.title('test', fontsize=20)
#plt.savefig(figureName, format='png')
        
plt.show(fig)
        
################################################################################
## Create June-Sept total smoke hour for each year for each dataset
################################################################################    
#years = np.arange(2007,2015,1)
#nlat = len(lat)
#nlon = len(lon)
#nyears = len(years)
#HMSSummerTotal = np.zeros([nlat, nlon, nyears])
#GFEDSummerTotal = np.zeros([nlat, nlon, nyears])
#
#i = 0
#for year in years:
#    yearMask = gridYear == year
#    monthMask = (gridMonth >= 6) & (gridMonth <= 9)
#    timeMask = yearMask & monthMask
#    # Sum emissions over the time dimension
#    HMSSummerTotal[:,:,i] = np.sum(HMS_emissions[:,:,timeMask], axis=2)
#    GFEDSummerTotal[:,:,i] = np.sum(GFED_emissions[:,:,timeMask], axis=2)
#    i+=1
#
## Set percentile based limit for colorbars , mask out zero values 
#HMS_max = np.percentile(HMSSummerTotal[(HMSSummerTotal > 0)], 95) 
#GFED_max = np.percentile(GFEDSummerTotal[(GFEDSummerTotal > 0)], 95) 
#
## TODO: Consider plotting the percentile value to make these comparable. 
#formatter = tkr.ScalarFormatter(useMathText=True)
#formatter.set_scientific(True)
#formatter.set_powerlimits((-2, 2))
#
## Plot the May-Sept totals for each individual summer
#for i in range(nyears):    
#    
#    # Set up figure 
#    fig = plt.figure(figsize=(12, 6))
#    figureName = figureSaveDir + 'HMS+' + gfed + '_emissions_summer_' \
#                  + str(years[i]) + '.png'    
#    HMSTitle = 'HMS June-September Emissions ' + str(years[i])    
#    GFEDTitle = gfed + ' June-Sept Emissions ' + str(years[i]) 
#    
#    # Set up first panel
#    ax = fig.add_subplot(121)
#
#    m.drawcoastlines()
#    m.drawstates()
#    m.drawcountries()
#        
#    m.drawparallels(parallels,labels=[1,0,0,0], fontsize=10)
#    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#    
#    data =  HMSSummerTotal[:,:,i]
#    data = ma.masked_where(data < 1, data)
#    cs = m.pcolor(x, y, data, cmap='Reds', vmin=data.min(), vmax=np.percentile(data,99))
#    
#    # add colorbar.
#    cbar = m.colorbar(cs, location='bottom',pad="5%") 
#    cbar.set_label('Observed emissions (hrs)', fontsize=20)
#    cbar.formatter.set_powerlimits((0, 0))
#    cbar.update_ticks()
#    # TODO: make colorbar label scientific notation 
#    
#    # add title
#    plt.title(HMSTitle, fontsize=20)
#    
#    ############################
#    # Plot the GFEDSummerTotal 
#    ############################
#    ax = fig.add_subplot(122)
#    
#    m.drawcoastlines()
#    m.drawstates()
#    m.drawcountries()
#        
#    m.drawparallels(parallels,labels=[1,0,0,0], fontsize=10)
#    m.drawmeridians(meridians,labels=[0,0,0,1],fontsize=10)
#    
#    data =  GFEDSummerTotal[:,:,i]
#    data = ma.masked_where(data <= 0, data)
#    cs = m.pcolor(x, y, data, cmap='Reds', vmin=data.min(), vmax=np.percentile(data,99))
#    
#    # add colorbar.
#    cbar = m.colorbar(cs, location='bottom',pad="5%")
#    cbar.set_label(gfedUnits, fontsize=20)
#    # TODO: make colorbar label scientific notation 
#    plt.title(GFEDTitle, fontsize=20)
#
#    plt.savefig(figureName, format='png')        
#    plt.close(fig)   

    
################################################################################
## Make a normalized [0-1] summed may-sept 2007-2014 side by side comparison
################################################################################
#HMSAllSummersTotal  = np.sum(HMSSummerTotal, axis=2)
#GFEDAllSummersTotal = np.sum(GFEDSummerTotal, axis=2)
#
##HMSAllSummersTotal_norm = HMSAllSummersTotal / HMSAllSummersTotal.max()
##GFEDAllSummersTotal_norm = GFEDAllSummersTotal / GFEDAllSummersTotal.max()
#
#HMSAllSummersTotal_noZero = ma.masked_where(HMSAllSummersTotal == 0,\
#                                            HMSAllSummersTotal)
#HMSLimit = np.percentile(HMSAllSummersTotal_noZero, 99)
#
#GFEDAllSummersTotal_noZero = ma.masked_where(GFEDAllSummersTotal == 0,\
#                                             GFEDAllSummersTotal)
#GFEDLimit = np.percentile(GFEDAllSummersTotal_noZero, 99)
#
## Make a figure showing both datasets side by side 
#fig = plt.figure(figsize=(12, 6))
#figureName = figureSaveDir + 'HMS+' + gfed +\
#             '_summer_emission_totals_2007-2014.png'    
# 
## Set up first panel
#ax = fig.add_subplot(121)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#    
#
#cs = m.pcolor(x, y, HMSAllSummersTotal_noZero, cmap='Reds',\
#              vmin=0, vmax=HMSLimit)
#cbar = m.colorbar(cs, location='bottom',pad="5%")#, ticks=np.linspace(0,1,5)) 
#cbar.set_label('HMS-Hours', fontsize=17)
#cbar.ax.tick_params(labelsize=14)
#
## add title
##plt.title('HMS', fontsize=20)

#############################
## Plot the GFEDSummerTotal 
#############################
#ax = fig.add_subplot(122)
#
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#    
#cs = m.pcolor(x, y, GFEDAllSummersTotal_noZero, cmap='Reds',\
#              vmin=0,vmax=GFEDLimit)
#cbar = m.colorbar(cs, location='bottom',pad="5%")#, ticks=np.linspace(0,1,5))
#cbar.set_label(gfed + ' ' + gfedUnits, fontsize=17)
#cbar.ax.tick_params(labelsize=14)
#
## TODO: make colorbar label scientific notation 
##plt.title(gfed, fontsize=20)
#
#
#plt.suptitle('June-September 2007-2014 Total Emissions', fontsize=20)
#fig.tight_layout()
#plt.savefig(figureName, format='png', bbox_inches='tight')        
#plt.close(fig)   
#
#
################################################################################
## How many HMS-hours per g/C/m**2
## TODO: MAKE SURE THAT ZERO VALUES ARE BEING PLOTTED 
################################################################################
#HMSHourPerGC = HMSAllSummersTotal_noZero / GFEDAllSummersTotal_noZero
#meanHMSHourPerGC = HMSHourPerGC.mean()
#
#print '------------------------------------------------------------'
#print 'The mean value of HMSHourPerGC is: ' + str(meanHMSHourPerGC)
#print '------------------------------------------------------------'
#
## NOTE: Locations where HMSHourPerGC == nan means both emissions are zero.
## NOTE: Locations where HMSHourPerGC == inf is float divided by zero, i.e. HMS
## NOTE: had emissions where GFED has none. For the sake of a colorbar, these 
## NOTE: both must be excluded from this analysis. 
##HMSHourPerGC_rm_nan = ma.masked_where(np.isnan(HMSHourPerGC), HMSHourPerGC)
##HMSHourPerGC_rm_nan_inf = ma.masked_where(np.isinf(HMSHourPerGC_rm_nan),\
##                                          HMSHourPerGC_rm_nan)
#limit99 = np.percentile(HMSHourPerGC, 99)
#infMask = np.isinf(HMSAllSummersTotal / GFEDAllSummersTotal)
#
## Show the inf as the max value on the colorbar
#HMSHourPerGC[infMask] = limit99
#
#
#fig = plt.figure(figsize=(12, 12))
#figureName = figureSaveDir + 'HMSHourPer' + gfed +\
#             '_summer_emission_totals_2007-2014.png'    
# 
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#    
#cs = m.pcolor(x, y, HMSHourPerGC, cmap='jet', vmin=0, vmax=limit99)
#cbar = m.colorbar(cs, location='bottom',pad="3%")#, ticks=np.linspace(0,1,5)) 
#cbar.set_label('HMS-Hours Per ' + gfedUnits, fontsize=23)
#cbar.ax.tick_params(labelsize=24)
#
#plt.suptitle('June-September 2007-2014 ' + 'HMS-Hours Per ' + gfedUnits,\
#             fontsize=34)
#fig.tight_layout()
#plt.savefig(figureName, format='png', bbox_inches='tight')        
#plt.close(fig)  
#
################################################################################
## Correlation analysis against GFED. These correlations will be for the
## 0)  WHOLE DOMAIN TIME SERIES COMPARISON of summertime totals inter-annual 
##     variability.
## 1) Totals in regions roughly equal to those used in analysis based on 
##    interannual variability
## 2) The interannual and monthly correlation at EACH grid cell.    
## 3) Depending on how these look, we may take the correlation to new levels, 
##    like subsetting by land cover etc. 
##
################################################################################
#
################################################################################
## Create total monthly time series plot + correlation for a given domain,
## start with the entire dataset. 
################################################################################    
#
#nMonths = GFED_emissions.shape[2]
#GFED_monthly_total = np.zeros(nMonths)
#HMS_monthly_total  = np.zeros(nMonths)
#for i in range(nMonths):
#    GFED_monthly_total[i] = np.sum(GFED_emissions[:,:,i])
#    HMS_monthly_total[i] = np.sum(HMS_emissions[:,:,i])
#
#
#fig, ax1 = plt.subplots(figsize=(10,6))
#
#figureName = figureSaveDir + 'HMS+' + gfed +'_monthly_emissions' + \
#                 '.pdf'    
#
## Plot normalized values since the quantity for HMS is arbitrary 
#GFED_monthly_total_norm = GFED_monthly_total / GFED_monthly_total.max()
#HMS_monthly_total_norm  = HMS_monthly_total / HMS_monthly_total.max()
#
#ax2 = ax1.twinx()
#ax1.plot(t, GFED_monthly_total_norm, 'g-', linewidth=4)
#ax1.tick_params(axis='x', labelsize=20)
#ax1.tick_params(axis='y', labelsize=20)
#
#ax2.plot(t, HMS_monthly_total_norm, 'b-', linewidth=4)
#ax2.tick_params(axis='y', labelsize=20)
#ax2.tick_params(axis='x', labelsize=20)
#
#
#ax1.set_xlabel('Time (month)', fontsize=20)
#ax1.set_ylabel(gfed + ' Total ' + gfedUnits, color='g', fontsize=19)
#ax2.set_ylabel('HMS-Hours', color='b', fontsize=19)
#
## Calculate the r-squared and use it in the title
#rsquared = np.corrcoef(GFED_monthly_total, HMS_monthly_total)[0,1]
#rsquared = np.round(rsquared, decimals=5)
#
#titleText = "Normalized Domain Summed Monthly Totals, r=" + str(np.round(rsquared,3))
#plt.title(titleText, fontsize=20)
#
#plt.savefig(figureName, format='pdf')
#plt.close(fig)   
#

#
#fig = plt.figure(figsize=(6,6) )
#plt.plot(GFED_monthly_total, HMS_monthly_total, 'ro')
#plt.xlabel("GFED4s total monthly Carbon emissions")
#plt.ylabel("HMS monthly burn hours observed")
#plt.title('All North America')
#plt.show()

###############################################################################
# Make the same style correlation only use summer totals instead of monthly 
# NOT IN PAPER
###############################################################################    
#GFED_summer_total = np.zeros(nyears)
#HMS_summer_total = np.zeros(nyears)
#t_year = []
#
#for i in range(nyears):
#    GFED_summer_total[i] = np.sum(GFEDSummerTotal[:,:,i])
#    HMS_summer_total[i] = np.sum(HMSSummerTotal[:,:,i])
#    t_year.append( date(years[i], 7, 15) )
#    
#t_year = np.array(t_year)
#    
#fig, ax1 = plt.subplots(figsize=(10,6))
#figureName = figureSaveDir + 'HMS+' + gfed + '_summer_emissions' + \
#                 '.pdf'    
#
## Plot normalized values since the quantity for HMS is arbitrary 
#GFED_summer_total_norm = GFED_summer_total / GFED_summer_total.max()
#HMS_summer_total_norm = HMS_summer_total / HMS_summer_total.max()
#
#ax2 = ax1.twinx()
#ax1.plot(t_year, GFED_summer_total_norm, 'g-', linewidth=4)
#ax2.plot(t_year, HMS_summer_total_norm, 'b-', linewidth=4)
#
#ax1.set_xlabel('Year', fontsize=15)
#ax1.set_ylabel(gfed + ' ' + gfedUnits, color='g', fontsize=19)
#ax2.set_ylabel('HMS Smoke Hrs', color='b', fontsize=19)
#
## Calculate the r-squared and use it in the title
#rsquared = np.corrcoef(GFED_summer_total, HMS_summer_total)[0,1]
#rsquared = np.round(rsquared, decimals=2)
#
#titleText = "Total June-Sept emissions (normalized), r=" + str(rsquared)
#plt.title(titleText, fontsize=22)
#
#plt.savefig(figureName, format='pdf')
#plt.close(fig)   
#
#
################################################################################
## Point by point monthly correlation map
###############################################################################    

#monthCorMat = np.zeros([nlat, nlon])
#for ii in range(nlat):
#    for jj in range(nlon):
#        monthCorMat[ii, jj] = np.corrcoef(HMS_emissions[ii,jj,:],\
#                                          GFED_emissions[ii,jj,:])[0,1]
#        
## Mask out the nan, no real meaning in plotting. The physical meaning of nan
## is zero variability in the value of the time series for one or both of the 
## grid cells. 
#monthCorMat = ma.masked_where(np.isnan(monthCorMat), monthCorMat)
#
## Plot this correlation matrix 
#fig = plt.figure(figsize=(12, 12))
#figureName = figureSaveDir + 'HMS+' + gfed + '_monthly_corr.png'
# 
#m.drawcoastlines()
#m.drawstates()
#m.drawcountries()
#
#
#cs = m.pcolor(x, y, monthCorMat, cmap='winter')
#cbar = m.colorbar(cs, location='bottom', pad="3%")#, ticks=ticks) 
#cbar.set_label('Pearson correlation coefficient', fontsize=22)
#cbar.ax.tick_params(labelsize=20)
#
## add title
#plt.title(gfed + ' and HMS-hours Monthly Correlation', fontsize=23)
#
#plt.savefig(figureName, bbox_inches='tight')
#plt.close(fig)


###############################################################################
# TODO: Figure out which data source has data for the SE and which one does not
###############################################################################   

