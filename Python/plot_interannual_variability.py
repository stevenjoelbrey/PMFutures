#!/usr/bin/env python2

# plot_interannual_variability.py 

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to explore the monthly and seasonal relationships
# between emissions (TODO: Burn area) and meteorology. 

# currently, this script contains much of the same functionality and load 
# statements as plotObservationParameterSpace.py and 
# plot_GFED4s_emission_summary.py. These will be consolidated at a later time
# but I am separating here to make organization more clear. 

# TODO: There is something up with precip. 

################################################################################
#------------- Arguments to Subset model emissions in space and time -----------
################################################################################
region = "_west_" # "_west_"| "_PNW_" | "_CAL_" | "_CentralRockies_" 


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

# Get region lat lon range and basemap plotting resolution based on the 
# chosen region
# TODO: Make a SE region! 
minLat, maxLat, minLon, maxLon, resolution  = cnm.getRegionBounds(region)

# Figure out what machine this code is running on
pwd = os.getcwd()
mac = '/Users/sbrey/Google Drive/sharedProjects/PMFutures/Python'
if pwd == mac:
	drive = "/Volumes/Brey_external/"
else:
	drive = "/barnes-scratch/sbrey/"

# Set directory paths
metDataDirBase = drive + "era_interim_nc_daily_merged/"
figureDir = '../Figures/GFED_era_interm_analysis/' # always relativ to py
	
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
C_daily_total = np.sum(C,axis=(1,2)) # sums daily value for all lon lat
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
# Show total summer emissions by month
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
plt.title("GDED4.1s June-Sept Emissions, 2003-2016", fontsize=27)
plt.savefig(figureDir + "summer_interannual_variability"+region+".png")
plt.close()

# Now stack the monthly contributions 
fig = plt.figure(figsize=(12,8))
ax = plt.subplot(111)
p1 = plt.bar(uniqueYears, C_june, color="blue")
p2 = plt.bar(uniqueYears, C_july, bottom=C_june, color="green")
p3 = plt.bar(uniqueYears, C_aug, bottom=(C_june + C_july), color="grey")
p4 = plt.bar(uniqueYears, C_sept, bottom=(C_june + C_july + C_aug), color="brown")
plt.legend( (p1[0], p2[0], p3[0], p4[0]), ("June", "July", "Aug", "Sept"), 
			frameon=False, loc="best", fontsize=27)
# Style this plot like a pro
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.tick_params(axis='y', labelsize=20)
ax.tick_params(axis='x', labelsize=20)
plt.xlabel("date", fontsize=26)
plt.ylabel("grams carbon emitted", fontsize=26)
plt.title("GDED4.1s June-Sept Emissions, 2003-2016", fontsize=27)
plt.savefig(figureDir + "summer_interannual_variability_months"+region+".png")
plt.close()


################################################################################
# Transitioning to looking at monthly 
################################################################################

# TODO: MASK LOCATIONS THAT DO NOT HAVE LOTS OF EMISSIONS. 
# UPDATE: No probably not. 

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

# Make temperature units America
d2m = cnm.KtoF(d2m)
t2m = cnm.KtoF(t2m)
tp = tp * 39.3701 # 39.37 inches per meter 

# These are the only lat lon values moving forward
latitude = ynew
longitude = xnew 

################################################################################
# Make Monthly totals and or means for meteorology. This dataframe will be used
# for analysis for the rest of the script. 
################################################################################

# rows of dataframe
nMonths = nYears * 12.

# columns of dataframe, TODO: soil moisture
colnames = ['year', 'month', 'E', 't2m', 'tp', 'windSpeed', 'z', 'RH', 'd2m']
# Create monthly dataframe
monthBlank = np.zeros(shape=(int(nMonths), len(colnames))) 
month_df = pd.DataFrame(data=monthBlank, columns=colnames)
# Create summer (6-9) dataframe
summerBlank = np.zeros(shape=(nYears, len(colnames)))
summer_df   = pd.DataFrame(data=summerBlank, columns=colnames)

df_i = -1 # this is a month counter
for y in range(nYears):

	# identify this years summer days to make summer average
	summerMask = (year == uniqueYears[y]) & (month >= 6) & (month <= 9)
	summer_df.year[y] = uniqueYears[y]
	
	# variables requiring summing first
	summer_df["E"][y] = np.sum(C[summerMask, :,:])
	summer_df["tp"][y] = np.sum(tp[summerMask, :,:])
	# then variables that require means
	summer_df["t2m"][y] = np.mean(t2m[summerMask, :,:])
	summer_df["windSpeed"][y] = np.mean(windSpeed[summerMask, :,:])
	summer_df["z"][y] = np.mean(z[summerMask, :,:])
	summer_df["RH"][y] = np.mean(RH[summerMask, :,:])
	summer_df["d2m"][y] = np.mean(d2m[summerMask, :,:])

	# for each year, loop through each month
	for m in range(12):
		df_i = df_i + 1
		# Mask out the days of this unique year month combo (e.g. jan 2003)
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


# TODO: make it possible for emissions to lag weather in this plot 
def plotParameterSpace(df, units, region, figureDir, string): 

	colnames = df.columns
	col_i = np.where(colnames == "E")[0][0] # Start with E column
	col_f = len(colnames) # end with the last column
	
	fig = plt.figure(figsize=(30, 30))
	figSaveName = figureDir + string + 'MeanCorrelations' + region + '.png'
				
	frameN = 0
	frameRow = 0
	for c1 in colnames[col_i:col_f]:
		
		frameRow = frameRow + 1
		frameColumn = 0
		
		for c2 in colnames[col_i:col_f]:
		    
		    frameN = frameN + 1
		    frameColumn = frameColumn + 1
		    
		    if frameColumn >= frameRow:
				# Get the data to plot against each-other
				yData = df[c1]
				xData = df[c2] # ~predictor
		
				# Get the linear fit 
				# https://docs.scipy.org/doc/scipy-0.15.1/reference/generated/scipy.stats.linregress.html
				# Just to the pierceson cor r for now
				r = np.corrcoef(xData, yData)[0][1]
				r = round(r,3)
				#p = round(p,5)
				
				ax = fig.add_subplot(7,7, frameN, frame_on=False)	
				ax.scatter(xData, yData, 
						   edgecolors='none',\
						   color="k",
						   s=50)
						   
				# show linear fit if r is decent 
				if (np.abs(r) >= 0.5):
					lm = stats.linregress(xData, yData)		
					yhat = lm.slope * xData + lm.intercept   
					plt.plot(xData, yhat, linewidth=3)
						   
				plt.xlabel(units[c2], fontsize=20, weight='bold')
				plt.ylabel(units[c1], fontsize=20, weight='bold')
				titleString = 'r=' + str(r)
				plt.title(titleString, fontsize=20, weight='bold')
		
				ax.spines['top'].set_visible(False)
				ax.spines['right'].set_visible(False)
				ax.yaxis.set_ticks_position('left')
				ax.xaxis.set_ticks_position('bottom')
				ax.tick_params(axis='y', labelsize=20)
				ax.tick_params(axis='x', labelsize=20)
			
			# TODO: plot linear model and R**2
			# TODO: when not plotting emissions as x or y, colorcode by emissions? 
		
	fig.tight_layout()		
	plt.savefig(figSaveName)
	plt.close()
	
plotParameterSpace(month_df, units, region, figureDir, "monthly")	
plotParameterSpace(summer_df, units, region, figureDir, "summer")	



		
		
# TODO: make it possible for emissions to lag weather in this plot 
def plotEmissionVsMet(df, units, region, figureDir, string, plotType):
	"""This basically just gives the top row of plotParameterSpace()"""
	
	colnames = df.columns
	col_i = np.where(colnames == "E")[0][0] + 1 # Start with first column after E
	col_f = len(colnames) # end with the last column
    
	fig = plt.figure(figsize=(22, 4))
	figSaveName = figureDir + string + 'EmissionsVsMet_' + plotType + region + '.png'
				
	frameN = 0
	for c2 in colnames[col_i:col_f]:
		
		frameN = frameN + 1
		
		# Get the data to plot against each-other
		yData = df["E"] # This is always the y ~predicated
		xData = df[c2] # ~predictor

		# Just to the pierceson cor r for now
		r = np.corrcoef(xData, yData)[0][1]
		r = round(r,3)
		#p = round(p,5)
		
		ax = fig.add_subplot(1,6, frameN, frame_on=False)			

		if plotType=="scatter":
			# Sometimes we want to show a scatterplot of these data 
			ax.scatter(xData, yData, 
					   edgecolors='none',
					   color="k")
				   
			# show linear fit if r is decent 
			if (np.abs(r) >= 0.5):
				lm = stats.linregress(xData, yData)		
				yhat = lm.slope * xData + lm.intercept   
				plt.plot(xData, yhat, linewidth=3)
				   
			plt.xlabel(units[c2], fontsize=20, weight='bold')
			plt.ylabel(units["E"], fontsize=20, weight='bold')
			titleString = 'r=' + str(r) 
			plt.title(titleString, fontsize=20, weight='bold')
			
		elif plotType=="ccf": 		
			# Sometimes we want to do a cross correlation to see what 
			# lag provides the highest correlation pearson coef
			ccf = plt.xcorr(df['E'], df[c2], maxlags=12, normed=True,
							usevlines=True, color="blue")
			lag = ccf[0]
			cor = ccf[1]
			j = np.where(lag==0)[0][0]
			correlation = cor[j]			
		
			plt.plot(lag, cor, color="blue")
			plt.ylim( (0, 0.6) )
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
# We saw that total precip ~4 months before summer is well correlated with 
# summertime emissions. Show this relationship
################################################################################
spring_tp = np.zeros(nYears)
for y in range(nYears):
	springMask = (year == uniqueYears[y]) & (month >= 4) & (month <= 9)
	spring_tp[y] = np.sum(tp[springMask,:,:])

# Add this to the dataframe 
#summer_df["spring_tp"] = spring_tp
figureName = figureDir + 'spring+summer_tp_vs_E'+region+'.png'

fig = plt.figure(figsize=(5,5))
ax = fig.add_subplot(1,1,1, frame_on=True)
plt.scatter(spring_tp, summer_df.E)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

lm = stats.linregress(spring_tp, summer_df.E)		
yhat = lm.slope * spring_tp + lm.intercept   
r = round(lm.rvalue,3)

plt.plot(spring_tp, yhat, linewidth=3)

plt.xlabel('April-Sept Precip (in)', fontsize=15, weight='bold')
plt.ylabel('Emissions (grams carbon)', fontsize=15, weight='bold')
#plt.legend('r=' + str(r), loc='best', frameon=False, fontsize=15)
plt.title(region[1:-1] + ', r=' + str(r), fontsize=16, weight='bold')

fig.tight_layout()	

plt.savefig(figureName)
plt.close()

################################################################################
# Make emissions to precip relationship grid box specific
# For every day with emissions, I want to know how long it has been since it 
# rained in that grid box
################################################################################
daysSinceRain = np.zeros(C.shape)
daysSinceRain[:] =-9999
nLat = len(latitude)
nLon = len(longitude)
nTime = len(time) 
for x in range(nLon):
	for y in range(nLat):
		
		tp_ = tp[:,y,x]
		daysSinceRain_  = 0. # reset the counter for each new grid point 
		
		for t in range(nTime):
		
			if tp_[t] > 0.001:
				# If in here it rained, reset daily counter
				daysSinceRain_  = 0.
				daysSinceRain[t,y,x] = 0.
				
			elif tp_[t] <= 0.001:
				# If here it has not rained, advance the daily counter 
				daysSinceRain_ = daysSinceRain_ + 1
				daysSinceRain[t,y,x] = daysSinceRain_
			else
				print 'how did you get here?'
		
summerMask = (month >= 6) & (month <= 9)
#plt.scatter(daysSinceRain[summerMask,:,:], C[summerMask,:,:])

################################################################################
# Now see if the inter-annual variability in the occurrence of difference events
# predicts inter-annual variability in summer emissions. 
# TODO: look at these masks monthly too. 
################################################################################
maskFile = drive + 'era_interim_nc_daily_merged/met_event_masks_NA_2003_2016.nc'
nc = Dataset(maskFile, 'r')

latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
metMaskDict = {}
metMaskDict['stagnation_mask'] = nc.variables['stagnation_mask'][:]
metMaskDict['high_T_mask']     = nc.variables['high_T_mask'][:]
metMaskDict['low_precip_mask'] = nc.variables['low_precip_mask'][:]
metMaskDict['high_wind_mask']  = nc.variables['high_wind_mask'][:]
metMaskDict['low_RH_mask']     = nc.variables['low_RH_mask'][:]
metMaskDict['blocking_mask']   = nc.variables['blocking_mask'][:]
# TODO: cyclone days coming soon
nc.close()

# Spatially subset these met masks 
for key in metMaskDict.keys():
	summer_df[key] = np.zeros(nYears) # Need to place these columns in here now
	metMaskDict[key], ynew, xnew, = cnm.mask2dims(metMaskDict[key], longitude, latitude, 
											      0, minLon, maxLon, minLat, maxLat)

# Add summer totals of these identified events to summer_df
newColumns = metMaskDict.keys()

for y in range(nYears):

	# identify this years summer days to make summer average
	summerMask = (year == uniqueYears[y]) & (month >= 6) & (month <= 9)
	
	# Assign the new variables
	for key in newColumns:
		summer_df[key][y] = np.sum(metMaskDict[key][summerMask, :,:])



# plot up these relationships and the summary statistics
nMax = C.shape[1] * C.shape[2] * np.sum(summerMask) * 1. # max any mask value for a year

fig = plt.figure(figsize=(15,10))
i = 0
for key in newColumns:
	i = i + 1
	ax = fig.add_subplot(2,3,i, frame_on=True)		
	plt.scatter(summer_df[key], summer_df.E, label=key)
	plt.xlabel(key + ' metric', fontsize=15, weight='bold')
	plt.ylabel('Emissions (grams carbon)', fontsize=15, weight='bold')
	
	lm = stats.linregress(summer_df[key], summer_df.E)		
	yhat = lm.slope * summer_df[key] + lm.intercept   
	r = round(lm.rvalue,3)
	ax.text(0.45, 0.9, 'r='+str(r),
        verticalalignment='bottom', horizontalalignment='right',
        transform=ax.transAxes,
        color='green', fontsize=19)
	
	if r > 0.5:
		plt.plot(summer_df[key], yhat, linewidth=3)
	
	plt.title(key, fontsize=20, weight='bold')
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)

fig.tight_layout()		
plt.savefig(figureDir + 'mask_predictors'+region+'.png')
plt.close()



#matplotlib.pyplot.xkcd(scale=1, length=100, randomness=2)
###############################################################################
# bring on the FPA-FOD data for the same analysis 
###############################################################################
dataFile = "/home/sbrey/projects/PMFutures/Dataframes/fireOccurrence_t.csv"
df = pd.read_csv(dataFile)

# subset by dates of interest
dateMask = (df.year >= 2003) & (df.year <= 2016)
df_subset = df[dateMask]

# TODO: convert latitude longitude
lon = df_subset.LONGITUDE
lonAdjust = lon + 360.

# TODO: subset by space of interest
lonMask = (lonAdjust >= minLon) & (lonAdjust <= maxLon)
latMask = (df_subset.LATITUDE >= minLat) & (df_subset.LATITUDE <= maxLat)
spatialMask =  lonMask & latMask

df_subset =  df_subset[spatialMask]

df = df_subset

###############################################################################
# Make summer summary dataframes
###############################################################################
acres_lighting = np.zeros(nYears)
acres_other    = np.zeros(nYears)

lightningMask = df.STAT_CAUSE_DESCR == 'Lightning'

# Count summer total BA in acres for each time
for y in range(nYears):
	summerMask = (df.year == uniqueYears[y]) & (df.month >= 6) & (df.month <= 9)
	acres_lighting[y] = np.sum(df.FIRE_SIZE[summerMask & lightningMask])
	acres_other[y] = np.sum(df.FIRE_SIZE[summerMask & (lightningMask==False)])
	
acres_total = acres_lighting + acres_other

# Add these totals to the ongoing summer dataframe
summer_df['FOD_acres_lightning'] = acres_lighting
summer_df['FOD_acres_other'] = acres_other
summer_df['FOD_acres_total'] = acres_total


# plt.scatter(summer_df.E[0:10], summer_df.FOD_acres_total[0:10], label="All acres")
# plt.scatter(summer_df.E[0:10], summer_df.FOD_acres_lightning[0:10], label="lightning caused")
# plt.scatter(summer_df.E[0:10], summer_df.FOD_acres_other[0:10], label="human caused")
# 
# plt.legend()
# plt.show()



metCols = [u't2m', u'tp', u'windSpeed', u'z', u'RH', u'd2m']

for metName in metCols:

	fig = plt.figure(figsize=(6,6))

	xData = summer_df[metName][0:10]
	y_lightning = summer_df.FOD_acres_lightning[0:10]
	y_human =  summer_df.FOD_acres_other[0:10]
	y_all = summer_df.FOD_acres_total[0:10]


	lm_all = stats.linregress(xData, y_all)
	yhat_all = lm_all.slope * xData + lm_all.intercept	
	r_all = round(lm_all.rvalue, 3) 

	lm_lightning = stats.linregress(xData, y_lightning)		
	yhat_lightning = lm_lightning.slope * xData + lm_lightning.intercept   
	r_lightning = round(lm_lightning.rvalue,3)

	lm_human = stats.linregress(xData, y_human)
	yhat_human = lm_human.slope * xData + lm_human.intercept   
	r_human = round(lm_human.rvalue,3)

	plt.scatter(xData, y_all, label="Total reported burn area, r$^{2}$="+ str(r_all))
	plt.scatter(xData, y_lightning, label="lightning ignition, r$^{2}$="+ str(r_lightning))
	plt.scatter(xData, y_human, label="human ignition, r$^{2}$="+ str(r_human))

	plt.legend()

	plt.ylabel('Acres burned', weight='bold')
	plt.xlabel(metName, weight='bold')

	plt.plot(xData, yhat_all, color="k")
	plt.plot(xData, yhat_lightning, color="k")
	plt.plot(xData, yhat_human, color="k")

	fig.tight_layout()		

	plt.savefig(figureDir + metName + '_FPA' + region + '.png')
	plt.close()




#summerBlank = np.zeros(shape=(nYears, len(colnames)))
#FOD_summer_df = pd.DataFrame(data=summerBlank, columns=colnames)








