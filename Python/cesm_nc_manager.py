#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will evantually be a class to be used to help analyze CESM air
# quality data. 

# NOTE on NC variables vertical dimensions: Z3 is the height of the middle
# model layers above sea level. The model uses hybrid vertical coordinates
# which are terrain-following, sigma coordinates in most of the troposphere
# and then gradually become pure pressure coordinates somewhere around 100 hPa

#os.chdir('/home/sbrey/projects/PMFutures/Python')

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
from datetime import date

###############################################################################
# ------------------------- USER INPUT ---------------------------------------- 
###############################################################################

scenario = '2000Base'
AQ_or_fire = 'AirQuality' # TODO: Decide if this is needed 
NCVariable = 'HEIGHT'
analysisDay = date(2005, 3, 1) # will be something that is eventually looped over 
analysisLevel = 100 # mb

# TODO: Make these functions into a class

def dateNumToDate(dateNum):
	"""This function takes dates of the form YYYYMMDD and turns them into python
		date objects.

	Keyword arguments:
	dateNum -- A single date string in the format YYYYMMDD or an array of datestring containing YYYYMMDD

	"""
	n = len(dateNum)
	Dates = []
	for i in range(n):    
		DS = str(dateNum[i])
		yearInt = int(DS[0:4])
		monthInt = int(DS[4:6])
		dayInt = int(DS[6:8])
		Dates.append(date(yearInt, monthInt, dayInt))
	Dates = np.array(Dates)
	return(Dates)


def makeAQNCFile(NCVariable, scenario):
	"""Function for getting path of desired AirQuality variable nc file path."""

	scenDec = scenario[0:4]
	scenarioEmissions = scenario[4:len(scenario)]	
	scenDate = {}
	scenDate['2000tail']   = '200001-201012'
	scenDate['2000center'] = '2000'
	scenDate['2050tail']   = '204001-205012'
	scenDate['2050center'] = '2050'
	scenDate['2100tail']   = '209001-209912'
	scenDate['2100center'] = '2100'	

	rcp = {}
	rcp['2000Base']  = ''
	rcp['2050RCP45'] = 'rcp45'
	rcp['2050RCP85'] = 'rcp85'
	rcp['2100RCP45'] = 'rcp45'
	rcp['2100RCP85'] = 'rcp85'
	
	# Data will always live in the same place on Yellowstone, static. 
	dataDirBase = '/fischer-scratch/sbrey/outputFromYellowstone/'
	fileHead = 'cesm122_fmozsoa_f09f09_' 
	fileMid  =  scenDate[scenDec+'center']+'_'+rcp[scenario]+'fires_00.'+ NCVariable
	fileTail = '.daily.' + scenDate[scenDec+'tail'] + '.nc'
	ncFile   = dataDirBase + "AirQualityData/" + scenario + '/' + fileHead + fileMid + fileTail

	return ncFile
 
def makeEmissionNCFile(NCVariable, scenario, year):
	"""function for getting path of desired variable nc file path"""

	dateSpan = {}
	dateSpan['2010']   = '20000101-20101231'
	dateSpan['2050']   = '20400101-20501231'
	dateSpan['2100']   = '20900101-21001231'
	
	# Data will always live in the same place on Yellowstone, static. 
	dataDirBase = '/fischer-scratch/sbrey/outputFromYellowstone/FireEmissions/'
	fileHead = 'cesm130_clm5_firemodule_' 
	fileMid  =  scenario + '_02.' + NCVariable
	fileTail = '.daily.NA.' + dateSpan[year] + '.nc'
	ncFile   = dataDirBase + fileHead + fileMid + fileTail

	return ncFile 


def getNCVarDims(ncFile, NCVariable):
	""""Function takes NCVariable and RCP scenario and returns the dimensions of
        that variable in a dictionary."""	
	
	nc     = Dataset(ncFile, 'r')
	ncDims = nc.variables[NCVariable].dimensions
	
	# Assign each dimension variable to a datastructure (list?) with keys 
	dimDict = {}	
	for dim in ncDims:
		dimDict[dim] = nc.variables[dim][:]

	nc.close()

	# Datefile is the same as ncFile only NCVariable needs to be replaced with 
	# 'date'
	nc = Dataset(ncFile.replace(NCVariable, 'date'), 'r')
	numDate = nc.variables[u'date'][:]
	Dates = dateNumToDate(numDate)
	nc.close()
	# Append nicely formatted dates of these data 
	dimDict['Dates'] = Dates
	
	return dimDict
	
 


def getNCData(NCVariable, scenario, analysisDay, analysisLevel):
	"""Get the data for a selected level of NCVariable. This function will 
       automatically determine if the dataset has multiple levels. If there
       is only one dimension in the vertical analysisLevel arugment is ignored.
       TODO: Make this work for a passed time mask instead of a single date.

	  Keyword arguments:
	  NCVarialbe -- String , the variable to be loaded. 
	  scemario -- String, RCP scenario, helps select data 
	  analysisDay -- String, YYYYMMDD for which to get data
	  analysisLevel -- float, the pressure surface of interest for data defined on
                     pressure surfaces.
	
	"""


	# Get model dimension information from other functions 
	dimDict = getNCVarDims(NCVariable, scenario)
	
	# Find date of interest index in the data 
	dt=np.abs(analysisDay - dimDict['Dates']) 
	dayIndex = np.argmin(dt)
	
	# Load desired layer 
	ncFile = makeNCFile(NCVariable, scenario)	
	nc = Dataset(ncFile, 'r')
	dimDict['units'] = nc.variables[NCVariable].units
	dimDict['longName'] = nc.variables[NCVariable].long_name
	
	print nc.variables

	# Find the ilev to pull from the data, if relevent
	# TODO: Make dynamic based on what dimensions exist and in that order
    # subcalls 

		
	if 'ilev' in dimDict.keys():
		diff = np.abs(dimDict['ilev'] - analysisLevel)
		levIndex = np.argmin(diff)
		data = nc.variables[NCVariable][dayIndex, levIndex, :, :]
	elif 'lev' in dimDict.keys():
		diff = np.abs(dimDict['lev'] - analysisLevel)
		levIndex = np.argmin(diff)
		data = nc.variables[NCVariable][dayIndex, levIndex, :, :]
	else:
		data = nc.variables[NCVariable][dayIndex, :, :]

	nc.close()

	print 'the level chosen is: ' + str(levIndex)
	
	# Now get the land mask
	fireBase = '/fischer-scratch/sbrey/outputFromYellowstone/FireEmissions/'
	nc = Dataset(fireBase + 'cesm130_clm5_firemodule_landmask_f09x125.nc')
	landMask = nc.variables[u'landmask'][:]
	nc.close

	dimDict[NCVariable] = data
	dimDict['landMask'] = landMask

	return dimDict


def getUSGSTopo():
	"""This function returns USGS ncobject"""
	
	baseDir = '/fischer-scratch/sbrey/outputFromYellowstone/'
	ncFile= baseDir + 'USGS-gtopo30_0.9x1.25_remap_c051027.nc'
		
	nc = Dataset(ncFile,'r')
	srf_geo = {}
	srf_geo['PHIS'] = nc.variables[u'PHIS'][:]
	srf_geo['units'] = nc.variables[u'PHIS'].units
	srf_geo['longName'] = nc.variables[u'PHIS'].long_name
	nc.close()
	
	srf_geo['topo'] = srf_geo['PHIS'] / 9.81 # (m**2/s**2) / (m/s**2)
	
	return srf_geo


def getCESMPHIS():
	"""This function extracts and returns surface geopotential 
	   m**2/s**2."""

	baseDir = '/fischer-scratch/sbrey/outputFromYellowstone/'
	ncFile= baseDir + 'cesm122_fmozsoa_f09f09_2000_fires_00.PHIS.monthly.200001-201012.nc'
		
	nc = Dataset(ncFile,'r')
	srf_PHIS = {}
	# All times are the same, does not change, so just grab first
	srf_PHIS['PHIS'] = nc.variables[u'PHIS'][0,:,:] 
	srf_PHIS['Z_srf'] = srf_PHIS['PHIS'] / 9.81 # (m**2/s**2) / (m/s**2) = [m]
	srf_PHIS['units'] = nc.variables[u'PHIS'].units
	srf_PHIS['longName'] = nc.variables[u'PHIS'].long_name
	nc.close()

	return srf_PHIS

def getZ3():
	"""This function gets a month geopotential heights"""
	
	baseDir = '/fischer-scratch/sbrey/outputFromYellowstone/'
	ncFile= baseDir + 'cesm122_fmozsoa_f09f09_2000_fires_00.Z3.monthly.200001-201012.nc'
		
	nc = Dataset(ncFile,'r')
	Z = {}
	Z['Z3'] = nc.variables[u'Z3'][0,:,:,:] 
	Z['units'] = nc.variables[u'Z3'].units
	Z['longName'] = nc.variables[u'Z3'].long_name
	nc.close()

	return Z	

def contourPlot(Z, titleText, units):
	"""Plots a Single date model layer on global projection centered over NA"""
	# TODO: PLACE ALL ARGUMENTS INTO DICTIONARY

	# set up orthographic map projection with
	# perspective of satellite looking down at 50N, 100W.
	# use low resolution coastlines.

	fig = plt.figure(figsize=(6, 6))

	m = Basemap(projection='ortho',lat_0=45,lon_0=-100,resolution='l')
	# draw coastlines, country boundaries, fill continents.
	m.drawcoastlines(linewidth=0.25)
	m.drawcountries(linewidth=0.25)
	m.fillcontinents(color='coral',lake_color='aqua')
	# draw the edge of the map projection region (the projection limb)
	m.drawmapboundary(fill_color='aqua')
	# draw lat/lon grid lines every 30 degrees.
	m.drawmeridians(np.arange(0,360,30))
	m.drawparallels(np.arange(-90,90,30))

	# make up some data on a regular lat/lon grid.
	lons, lats = np.meshgrid(data['lon'], data['lat']) 
	x, y = m(lons, lats) # compute map proj coordinates.

	# contour data over the map.
	cs = m.contourf(x, y, Z, 16,linewidths=3)
	csFill = m.pcolor(x, y, Z, visible=False)
	cbar = m.colorbar(csFill, location='bottom', pad="1%")
	cbar.set_label(units, fontsize=17)
	#titleText= dateLab+' '+scenario + ' ' + str(analysisLevel) + 'hPa ' +NCVariable
	plt.title(titleText, fontsize=15)
	#plt.draw()
	plt.show(block=False)


###############################################################################
# Area for testing by calling functions, will eventually be a different script 
# of its own 
###############################################################################
#data = getNCData(NCVariable,scenario,analysisDay, analysisLevel)

###############################################################################
# HOW TO CONVERT TO Z500 from sigma pressure coordinate system
###############################################################################
# In any case for "lev" >100 hPa you should EXPECT to see the signature of
# topography on Z3 and HEIGHT contoured at any particular level.  Picking a
# level in the vertical that corresponds to around 500 hPa in the reference
# pressure profile will not be the same as plotting Z500.  To actually obtain
# contours of Z500 you should interpolate the daily 3D HEIGHT fields in pressure
# to 500 hPa using a 3D field of edge pressures constructed from a formula like

#P(x,y,L)= a(L)*100000. + b(L)*PS(x,y)

#Here a(L) and b(L) are "hyai" and "hybi" fields on most CAM history files, and PS is the surface pressure.  Please check the formula above to see where the a's and b's should actually go with respect to 100000 (Pa) and PS(x,y)
	

#contourPlot(data[NCVariable], 'geopotential height above surface at interfaces', '')
#contourPlot(Z, 'Ground Level Geopotential Height [m]', 'meters')
#contourPlot(data[NCVariable] - Z, 'Z - srf_Z', 'meters')

#Z_sea = getZ3()
#contourPlot(Z_sea['Z3'][18,:,:] , 'Geopotential Height of 500mb (above sea level, Z3)', 'meters')

#if Name == '__main__':
	
