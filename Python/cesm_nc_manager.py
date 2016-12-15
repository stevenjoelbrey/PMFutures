###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will evantually be a class to be used to help analyze CESM air
# quality data. 

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
# ------------------------- USE INPUT ----------------------------------------- 
###############################################################################
scenario = '2000Base'
AQ_or_fire = 'AirQuality' # TODO: Decide if this is needed 
NCVariable = 'HEIGHT'
analysisDay = date(2005, 2, 1) # will be something that is eventually looped over 
analysisLevel = 500 # mb

# TODO: Make these functions into a class

def dateNumToDate(dateNum):
	"""This function takes dates of the form YYYYMMDD and turns them into python
		date objects."""
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


def makeNCFile(NCVariable, scenario):
	"""function for getting path of desired variable nc file path\n
       TODO: Make work for emissions data also, dynamically"""

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
		
	dataDirBase = '/fischer-scratch/sbrey/outputFromYellowstone/AirQualityData/'
	fileHead = 'cesm122_fmozsoa_f09f09_' 
	fileMid  =  scenDate[scenDec+'center']+'_'+rcp[scenario]+'fires_00.'+ NCVariable
	fileTail = '.daily.' + scenDate[scenDec+'tail'] + '.nc'
	ncFile   = dataDirBase + scenario + '/' + fileHead + fileMid + fileTail
	
	return ncFile


def getNCVarDims(NCVariable, scenario):
	""""Function takes NCVariable and RCP scenario and returns the dimensions of
        that variable in a dictionary."""	
	
	ncFile = makeNCFile(NCVariable, scenario)	
	nc = Dataset(ncFile, 'r')
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
       TODO: Make this work for a passed time mask instead of a single date."""


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
	
	# Find the ilev to pull from the data, if relevent
	# TODO: Make dynamic based on what dimensions exist and in that order subcalls 
	print dimDict.keys()
		
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
	
	# Now get the land mask
	fireBase = '/fischer-scratch/sbrey/outputFromYellowstone/FireEmissions/'
	nc = Dataset(fireBase + 'cesm130_clm5_firemodule_landmask_f09x125.nc')
	landMask = nc.variables[u'landmask'][:]
	nc.close

	return data, landMask, dimDict


def contourPlot(data, dimDict, scenario, NCVariable, dateLab, analysisLevel):
	"""Plots a Single date model layer on global projection centered over NA"""
	if(type(dateLab) != str):
		dateLab = str(dateLab)

	# set up orthographic map projection with
	# perspective of satellite looking down at 50N, 100W.
	# use low resolution coastlines.

	fig = plt.figure(figsize=(10, 10))

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
	lons, lats = np.meshgrid(dimDict['lon'], dimDict['lat']) 
	x, y = m(lons, lats) # compute map proj coordinates.

	# contour data over the map.
	cs = m.contour(x, y, data, 15,linewidths=3)
	cbar = m.colorbar(cs, location='bottom',pad="1%")
	titleText= dateLab+' '+scenario + ' ' + str(analysisLevel) + 'hPa ' +NCVariable
	plt.title(titleText, fontsize=30)
	#plt.draw()
	plt.show(block=False)


###############################################################################
# Area for testing by calling functions, will eventually be a different script 
# of its own 
###############################################################################
data,landMask,dimDict=getNCData(NCVariable,scenario,analysisDay, analysisLevel)


contourPlot(data, dimDict, scenario, NCVariable, dateLab=analysisDay,analysisLevel=analysisLevel)
	
	








