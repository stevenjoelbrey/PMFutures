#!/usr/bin/env python2

# plot_interannual_variability.py 

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to explore the monthly and seasonal relationships
# between emissions and meteorology. 

# currently, this script contains much of the same functionality and load 
# statements as plotObservationParameterSpace.py and 
# plot_GFED4s_emission_summary.py. These will be consolidated at a later time
# but I am separating here to make organization more clear. 


################################################################################
#------------- Arguments to Subset model emissions in space and time -----------
################################################################################

cutoffPercentile = 80.
startMonth = 6
endMonth   = 9
region     = "_CAL_" # "_west_"| "_PNW_" | "_CAL_" | "_CentralRockies_" 


# Load resources
import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import datetime
import matplotlib.ticker as tkr
import cesm_nc_manager as cnm
import pandas as pd
from scipy.stats import pearsonr
import scipy.stats as stats
import statsmodels

# Get region lat lon range	
minLat, maxLat, minLon, maxLon, resolution  = cnm.getRegionBounds(region)

# Figure out what machine this code is running on
# Figure out what machine this code is running on
pwd = os.getcwd()
mac = '/Users/sbrey/Google Drive/sharedProjects/PMFutures/Python'
if pwd == mac:
	drive = "/Volumes/Brey_external/"
else:
	drive = "/barnes-scratch/sbrey/"


metDataDirBase = drive + "era_interim_nc_daily_merged/"
figureDir = '../Figures/GFED_era_interm_analysis/'
	
# Get emissions, use this to get dimensions
ncFile  = drive + "GFED4s/GFED4.1s_METGrid_C_NA_2003_2016.nc"
nc = Dataset(ncFile, 'r')
latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
time = nc.variables['time'][:]
C = nc.variables['C'][:]
nc.close()

# Make time into datetime arrays
time, month, year = cnm.get_era_interim_time(time)

# Spatially subset the data 
C, ynew, xnew = cnm.mask2dims(C, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)

################################################################################
# Show emissions time series for the domain
################################################################################
C_daily_total = np.sum(C,axis=(1,2))
C_cumulative = np.cumsum(C_daily_total)

fig = plt.figure(figsize=(12,8))
ax = plt.subplot(111)
plt.plot(time, C_daily_total)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
plt.xlabel("date", fontsize=26)
plt.ylabel("grams carbon emitted", fontsize=26)
plt.title("GDED4.1s Daily Emissions", fontsize=29)
plt.savefig(figureDir + "daily_timeSeries" +region+".png")
plt.close()

# TODO: Also plot different region contributions lines! That would be dope!

################################################################################
# show emissions monthly histogram
################################################################################
uniqueMonths = np.unique(month)
monthTotal = np.zeros(12)
for i in range(12):
	monthMask = month == uniqueMonths[i]
	monthTotal[i] = np.sum(C_daily_total[monthMask])

fig = plt.figure(figsize=(12,8))
ax = plt.subplot(111)
plt.bar(uniqueMonths, monthTotal)
plt.xticks(uniqueMonths, uniqueMonths)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
plt.xlabel("Month", fontsize=26)
plt.ylabel("grams carbon emitted", fontsize=26)
plt.title("GDED4.1s Seasonality", fontsize=29)
plt.savefig(figureDir + 'emissions_seasonality' + region + ".png")
plt.close()


################################################################################
# Show total summer emissions
################################################################################
uniqueYears = np.unique(year)
nYears = len(uniqueYears)
C_summer = np.zeros(nYears)
C_june = np.zeros(nYears)
C_july = np.zeros(nYears)
C_aug = np.zeros(nYears)
C_sept = np.zeros(nYears)

for i in range(nYears):
	summerMask = (month >= 6.) & (month <= 9.) & (year == uniqueYears[i])
	juneMask = (month == 6.) & (year == uniqueYears[i])
	julyMask = (month == 7.) & (year == uniqueYears[i])
	augMask  = (month == 8.) & (year == uniqueYears[i])
	septMask  = (month == 9.) & (year == uniqueYears[i])
	
	C_summer[i] = np.sum(C_daily_total[summerMask])
	C_june[i] = np.sum(C_daily_total[juneMask])
	C_july[i] = np.sum(C_daily_total[julyMask])
	C_aug[i] = np.sum(C_daily_total[augMask])
	C_sept[i] = np.sum(C_daily_total[septMask])

# First plot with all one color 
fig = plt.figure(figsize=(12,8))
ax = plt.subplot(111)
p1 = plt.bar(uniqueYears, C_summer)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
plt.xlabel("date", fontsize=26)
plt.ylabel("grams carbon emitted", fontsize=26)
plt.title("GDED4.1s June-Sept Emissions", fontsize=29)
plt.savefig(figureDir + "summer_interannual_variability"+region+".png")
plt.close()

# Now stack the monthly contributions 
fig = plt.figure(figsize=(12,8))
ax = plt.subplot(111)
p1 = plt.bar(uniqueYears, C_june, color="blue")
p2 = plt.bar(uniqueYears, C_july, bottom=C_june, color="green")
p3 = plt.bar(uniqueYears, C_aug, bottom=(C_june+C_july), color="grey")
p4 = plt.bar(uniqueYears, C_sept, bottom=(C_june+C_july+C_aug), color="lightpink")
plt.legend( (p1[0], p2[0], p3[0], p4[0]), ("June", "July", "Aug", "Sept"), 
			frameon=False, loc="best", fontsize=27)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
plt.xlabel("date", fontsize=26)
plt.ylabel("grams carbon emitted", fontsize=26)
plt.title("GDED4.1s Summer Emissions", fontsize=29)
plt.savefig(figureDir + "summer_interannual_variability_months"+region+".png")
plt.close()


################################################################################
# Transitioning to looking at monthly 
################################################################################

# TODO: MASK LOCATIONS THAT DO NOT HAVE LOTS OF EMISSIONS

################################################################################
# Bring in meteorology and event masks
################################################################################
# make a handy little function for getting met 
def ncGetSelf(Dir, varname):
	nc = Dataset(Dir + varname +"_NA_2003_2016.nc", 'r')
	VAR = nc.variables[varname][:]
	nc.close()
	return VAR

t2m = ncGetSelf(metDataDirBase, 't2m')[:]
tp = ncGetSelf(metDataDirBase, 'tp')[:] 
u10 = ncGetSelf(metDataDirBase, 'u10')[:] 
v10 = ncGetSelf(metDataDirBase, 'v10')[:]
windSpeed = np.sqrt(u10**2 + v10**2) 
z   = ncGetSelf(metDataDirBase, 'z')[:,1,:,:] / 9.81
RH   = ncGetSelf(metDataDirBase, 'RH')[:] # created by make_era_interim_met_masks.py
d2m = ncGetSelf(metDataDirBase, 'd2m')[:]
# TODO: consider also using dew point strait up?

# spatially subset these data based on specified range
t2m, ynew, xnew = cnm.mask2dims(t2m, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
tp, ynew, xnew = cnm.mask2dims(tp, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
windSpeed, ynew, xnew = cnm.mask2dims(windSpeed, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
z, ynew, xnew = cnm.mask2dims(z, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
RH, ynew, xnew = cnm.mask2dims(RH, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)
d2m, ynew, xnew = cnm.mask2dims(d2m, longitude, latitude, 0, minLon, maxLon, minLat, maxLat)

# Make units America
d2m = cnm.KtoF(d2m)
t2m = cnm.KtoF(t2m)

# These are the only lat lon values moving forward
latitude = ynew
longitude = xnew 

# Make Monthly totals and or means for meteorology 
# TODO: Use pandas

# rows of dataframe
nMonths = nYears * 12.
# columns of dataframe
colnames = ['year', 'month', 'E', 't2m', 'tp', 'windSpeed', 'z', 'RH', 'd2m']
monthBlank = np.zeros(shape=(int(nMonths), len(colnames))) 
month_df = pd.DataFrame(data=monthBlank, columns=colnames)

summerBlank = np.zeros(shape=(nYears, len(colnames)))
summer_df = pd.DataFrame(data=summerBlank, columns=colnames)

# TODO: Consider all kinds of different lags, like spring weather and
# TODO: summer emissions 

df_i = -1 # this is a month counter
for y in range(nYears):

	# identify this years summer days to make summer average
	summerMask = (year == uniqueYears[y]) & (month >= 6) & (month <= 9)
	summer_df.year[y] = uniqueYears[y]
	#summer_df.month[y] = '6-9'
	# variables requiring summing first
	summer_df["E"][y] = np.sum(C[summerMask, :,:])
	summer_df["tp"][y] = np.sum(tp[summerMask, :,:])
	# then variables that require means
	summer_df["t2m"][y] = np.mean(t2m[summerMask, :,:])
	summer_df["windSpeed"][y] = np.mean(windSpeed[summerMask, :,:])
	summer_df["z"][y] = np.mean(z[summerMask, :,:])
	summer_df["RH"][y] = np.mean(RH[summerMask, :,:])
	summer_df["d2m"][y] = np.mean(d2m[summerMask, :,:])

	for m in range(12):
		df_i = df_i + 1
		# Mask out the days of this unique year month combo 
		monthMask = (uniqueMonths[m] == month) & (year == uniqueYears[y])
		month_df.year[df_i] = uniqueYears[y]
		month_df.month[df_i]= uniqueMonths[m]
		
		# Assign met values for this month mask (~30 days) and all of 
		# the spatial domain 
		
		# variables requiring summing first
		month_df["E"][df_i] = np.sum(C[monthMask, :,:])
		month_df["tp"][df_i] = np.sum(tp[monthMask, :,:])
		# then variables that require means
		month_df["t2m"][df_i] = np.mean(t2m[monthMask, :,:])
		month_df["windSpeed"][df_i] = np.mean(windSpeed[monthMask, :,:])
		month_df["z"][df_i] = np.mean(z[monthMask, :,:])
		month_df["RH"][df_i] = np.mean(RH[monthMask, :,:])
		month_df["d2m"][df_i] = np.mean(d2m[monthMask, :,:])
		
		
		
################################################################################		
# plot monthly and yearly values against each-other in a 7x7 figure
# where each frame is a scatterplot and each 		
################################################################################
units = {}
units['E'] = 'emission (g carbon)'
units['tp'] = 'precip (inches)'
units['t2m'] = 'temperature (F)'
units['windSpeed'] = 'windspeed (m s$^{-1}$)'
units['z'] = '500mb height (m)'
units['RH'] = "RH%"
units['d2m'] = "dew point (f)"

################################################################################		
# TODO: Make a lagged correlation plot! Maybe easier to do with R? 
################################################################################


# TODO: make it possible for emissions to lag weather in this plot 
def plotParameterSpace(df, units, region, figureDir, string="monthly"): 

	colnames = df.columns
	fig = plt.figure(figsize=(30, 30))

	figSaveName = figureDir + string + 'MeanCorrelations' + region + '.png'
				
	frameN = 0
	frameRow = 0
	for c1 in colnames[2:9]:
		
		frameRow = frameRow + 1
		frameColumn = 0
		
		for c2 in colnames[2:9]:
		    
		    frameN = frameN + 1
		    frameColumn = frameColumn + 1
		    
		    if frameColumn >= frameRow:
				# Get the data to plot against each-other
				yData = df[c1]
				xData = df[c2] # ~predictor
		
				# Get the linear fit 
				# https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.linregress.html
				# Just to the pierceson cor r for now
				r, p = pearsonr(xData, yData)
				r = round(r,3)
				p = round(p,5)
				
				
				ax = fig.add_subplot(7,7, frameN, frame_on=False)	
		
				ax.scatter(xData, yData, 
						   edgecolors='none',\
						   color="k"
						   )
						   
				# show linear fit if r is decent 
				if (np.abs(r) >= 0.5):
					lm = stats.linregress(xData, yData)		
					yhat = lm.slope * xData + lm.intercept   
					plt.plot(xData, yhat, linewidth=3)
						   
				plt.xlabel(units[c2], fontsize=18, weight='bold')
				plt.ylabel(units[c1], fontsize=18, weight='bold')
				titleString = 'r=' + str(r) + ', p=' + str(p)
				plt.title(titleString, 
							fontsize=19, weight='bold')
		
				ax.spines['top'].set_visible(False)
				ax.spines['right'].set_visible(False)
				ax.yaxis.set_ticks_position('left')
				ax.xaxis.set_ticks_position('bottom')
				ax.tick_params(axis='y', labelsize=18)
				ax.tick_params(axis='x', labelsize=18)
			
			# TODO: plot linear model and R**2
			# TODO: when not plotting emissions as x or y, colorcode by emissions? 
		
	fig.tight_layout()		
	plt.savefig(figSaveName)
	plt.close()
	
plotParameterSpace(month_df, units, region, figureDir, "monthly")	
plotParameterSpace(summer_df, units, region, figureDir, "summer")	

# TODO: make it possible for emissions to lag weather in this plot 
def plotEmissionVsMet(df, units, region, figureDir, string="monthly", plotType="scatter"): 

	colnames = df.columns
	fig = plt.figure(figsize=(22, 4))

	figSaveName = figureDir + string + 'EmissionsVsMet_' + plotType + region + '.png'
				
	frameN = 0
	for c2 in colnames[3:9]:
		
		frameN = frameN + 1
		
		# Get the data to plot against each-other
		yData = df["E"]
		xData = df[c2] # ~predictor

		# Just to the pierceson cor r for now
		r, p = pearsonr(xData, yData)
		r = round(r,3)
		p = round(p,5)
		
		ax = fig.add_subplot(1,6, frameN, frame_on=False)			

		if plotType=="scatter":
		
			
			ax.scatter(xData, yData, 
					   edgecolors='none',\
					   color="k"
					   )
				   
			# show linear fit if r is decent 
			if (np.abs(r) >= 0.5):
				lm = stats.linregress(xData, yData)		
				yhat = lm.slope * xData + lm.intercept   
				plt.plot(xData, yhat, linewidth=3)
				   
			plt.xlabel(units[c2], fontsize=16, weight='bold')
			plt.ylabel(units["E"], fontsize=16, weight='bold')
			titleString = 'r=' + str(r) 
			plt.title(titleString, 
						fontsize=15, weight='bold')
			
		elif plotType=="ccf": 		
		
			ccf = plt.xcorr(df.E, df[c2], maxlags=12, 
							usevlines=True, color="blue") 
			plt.plot(ccf[0], ccf[1], color="blue")
			plt.xlabel('lag (months)', fontsize=16, weight='bold')
			plt.ylabel('correlation', fontsize=16, weight='bold')
			plt.title('cor(E, ' + c2 + ')', weight='bold', fontsize=20)				
	
		# universal style 
		ax.spines['top'].set_visible(False)
		ax.spines['right'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		ax.tick_params(axis='y', labelsize=16)
		ax.tick_params(axis='x', labelsize=16)
			
					
	fig.tight_layout()		
	plt.savefig(figSaveName)
	plt.close()

plotEmissionVsMet(month_df, units, region, figureDir, "monthly", "scatter")
plotEmissionVsMet(summer_df, units, region, figureDir, "summer", "scatter")

# Make the lagged cross correlation version of the figure
plotEmissionVsMet(month_df, units, region, figureDir, "monthly", "ccf")


################################################################################
# Now see if the inter-annual variability in the occurrence of difference events
# predicts inter-annual variability in summer emissions. 
################################################################################
# https://stackoverflow.com/questions/40578701/acf-confidence-intervals-in-r-vs-python-why-are-they-different
#statsmodels.graphics.tsaplots.plot_acf()



#matplotlib.pyplot.xkcd(scale=1, length=100, randomness=2)


