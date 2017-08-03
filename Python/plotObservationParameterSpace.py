#!/usr/bin/env python2

# plotEmissionSummary.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to make scatterplots that allow us to explore the
# relationship between different parameters and the emissions. 


################################################################################
#------------------ Subset model emissions in space and time ------------------
################################################################################

startMonth = 6
endMonth   = 9
minLat     = # "_NA_" file default  
maxLat     = # 
minLon     = #
maxLon     = # 


# TODO: spatial subset the analysis by the locations that have the highest 
# TODO: cumulative emissions, since places with few emissions dont really
# TODO: matter anyways

# TODO: Make sure that we are actually only plotting the non-masked values. 
# TODO: use ma.compressed() to be sure. 

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

dataDirBase = "/barnes-scratch/sbrey/era_interim_nc_daily_merged/"


# Get emissions, use this to get dimensions
ncFile  = "/barnes-scratch/sbrey/GFED4s/GFED4.1s_METGrid_C_NA_2003_2016.nc"
nc = Dataset(ncFile, 'r')
latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
time = nc.variables['time'][:]
C = nc.variables['C'][:]
nc.close()

t, month, year = cnm.get_era_interim_time(time)

month_matrix = np.zeros(C.shape)
for i in range(len(t)):
	month_matrix[i,:,:] = month[i]


# Subset all by months of interest
monthMask = (month >= 6) & (month <= 9)

C = C[monthMask, :,:]
t = t[monthMask]
month = month[monthMask]
year = year[monthMask]
month_matrix = month_matrix[monthMask,:,:]

#####################################################
# function for calculating the monthly mean 
#####################################################
def makeMonthlyAverage(data, month, year, mathType='mean'):
	"""
	This function take data of the from (t, x, y) where t is a daily time dimension
        and calculates monthly mean or summed values.
	"""
	# TODO: DOCUMENT ME

	# build the structure to save everything 
	yearMonths = np.array( [year.astype(str)[o] + month.astype(str)[o] for o in range(0, len(year)) ] )	
	nUniqueMonth = len(np.unique(yearMonths))
	monthlyValue = np.zeros((nUniqueMonth, data.shape[1], data.shape[2]))
	i=-1
	for y in np.unique(year):
		for m in np.unique(month):
			i = i + 1
			mask = (y == year) & (m == month)	
			if mathType == 'mean':
				monthlyValue[i, :,:] = np.mean(data[mask, :, :], axis=0)
			elif mathType == 'sum':
				monthlyValue[i, :,:] = np.sum(data[mask, :, :], axis=0)

	return monthlyValue

			
# make a handy little function for getting met 
def ncGetSelf(Dir, varname):
	nc = Dataset(Dir + varname +"_NA_2003_2016.nc", 'r')
	VAR = nc.variables[varname][:]
	nc.close()
	return VAR

# TODO: consider only showing relationships for where 
# TODO: emissions are greater than zero 

t2m = ncGetSelf(dataDirBase, 't2m')[monthMask, :,:]
tp = ncGetSelf(dataDirBase, 'tp')[monthMask, :,:] 
u10 = ncGetSelf(dataDirBase, 'u10')[monthMask, :,:] 
v10 = ncGetSelf(dataDirBase, 'v10')[monthMask, :,:]
windSpeed = np.sqrt(u10**2 + v10**2) 
z   = ncGetSelf(dataDirBase, 'z')[monthMask,1,:,:] 
RH   = ncGetSelf(dataDirBase, 'RH')[monthMask, :,:] # created by make_era_interim_met_masks.py
d2m = ncGetSelf(dataDirBase, 'd2m')[monthMask, :,:]
# TODO: consider also using dew point strait up?

# Make units America
d2m = cnm.KtoF(d2m)
t2m = cnm.KtoF(t2m)


# TODO: Plot all parameters in space with abline for event thresholds
# Make a version of all these plots that shows the points that make up 99% of total 
# emissions

# place all meteorology into a dictionary for easy forloop handling 
v = {}
v['T'] = t2m
v['windSpeed'] = windSpeed
v['precip'] = tp
v['Td'] = d2m
v['Z500'] = z

# Include/make monthly mean versions of the data! 
vMonth = {}
vMonth['T'] = makeMonthlyAverage(t2m, month, year, mathType='mean')
vMonth['windSpeed'] = makeMonthlyAverage(windSpeed, month, year, mathType='mean')
vMonth['precip'] = makeMonthlyAverage(tp, month, year, mathType='sum')
vMonth['Td'] = makeMonthlyAverage(d2m, month, year, mathType='mean')
vMonth['Z500'] = makeMonthlyAverage(z, month, year, mathType='mean')
C_monthly = makeMonthlyAverage(C, month, year, mathType='sum')

# Mask all variables (v[key]) where C == 0. We care about these relationships 
# where there are emissions! 
for v_key in v.keys():
	v[v_key] = np.ma.masked_where(C == 0, v[v_key])
 	vMonth[v_key] = np.ma.masked_where(C_monthly == 0, vMonth[v_key])

# Now Mask C itself
C = np.ma.masked_where(C == 0, C)
C_monthly = np.ma.masked_where(C_monthly == 0, C_monthly)

u = {}
u['T'] = 'temperature [F]'
u['windSpeed'] = 'windspeed [m s$^{-1}$]'
u['precip'] = 'precip [in day$^{-1}$]'
u['Td'] = "dew point [F]"
u['Z500'] = "Z500 [m$^{2}$ s$^{-2}$]"

uMonth = {}
uMonth['T'] = 'temperature [F]'
uMonth['windSpeed'] = 'windspeed [m s$^{-1}$]'
uMonth['precip'] = 'precip [in month$^{-1}$]'
uMonth['Td'] = "dew point [F]"
uMonth['Z500'] = "Z500 [m$^{2}$ s$^{-2}$]"

# Range of emissions is the same for all variables in the space so set a common
# colorbar range. 
fig = plt.figure(figsize=(30,22))

figSaveName = '../Figures/GFED_era_interim_'+str(startMonth)+'_'+ str(endMonth)+\
		'_daily_parameter_space.png'

nVar = len(v.keys())
frameN = 0

vmin = C.min()
vmax = C.max()


for i in range(nVar):

	# Get x-axis plotting parameter 
	x = v.keys()[i]

	for j in range(nVar):
		
		# Get y-axis plotting parameter 
		y = v.keys()[j]

		frameN += 1 
		
		# Recal these have been masked where emissions (C) is zero
		xData = v[x]
		yData = v[y]

		ax = fig.add_subplot(5,5,frameN, frame_on=False)

		c = ax.scatter(xData, yData, c = C, marker=".",\
				s=5, edgecolors='none',\
				vmin=vmin, vmax=vmax,\
				alpha=1,\
				norm=matplotlib.colors.LogNorm()
				#cmap="Reds"
	                       )
	
		ax.spines['top'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		ax.tick_params(axis='y', labelsize=15)
		ax.tick_params(axis='x', labelsize=15)


		#cb = plt.colorbar(c, cmap='jet', extend='both')
		#cb.set_label('Emissions [kg/day]')

		plt.xlim([v[x].min(), v[x].max()])
		plt.ylim([v[y].min(), v[y].max()])

		plt.xlabel(u[x], fontsize=20)
		plt.ylabel(u[y], fontsize=20)

		print 'Figure frameN: ' + str(frameN) + ' is complete'

fig.tight_layout()
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(c, cax=cbar_ax, extend="both")
cbar.set_label('Emissions [g C day$^{-1}$]', fontsize=30)
cbar.ax.tick_params(labelsize=24)
print 'Writing the figure to memory'
plt.savefig(figSaveName, dpi=600)
 
##################################################################
# Same figure done monthly
# TODO: make this space a function that is passed a type of data
##################################################################
fig = plt.figure(figsize=(30,22))

figSaveName = '../Figures/GFED_era_interim_'+str(startMonth)+'_'+ str(endMonth)+\
		'_monthly_parameter_space.png'

nVar = len(v.keys())
frameN = 0

vmin = C_monthly.min()
vmax = C_monthly.max()


for i in range(nVar):

	# Get x-axis plotting parameter 
	x = v.keys()[i]

	for j in range(nVar):
		
		# Get y-axis plotting parameter 
		y = v.keys()[j]

		frameN += 1 
		
		# Recal these have been masked where emissions (C) is zero
		xData = vMonth[x]
		yData = vMonth[y]

		ax = fig.add_subplot(5,5,frameN, frame_on=False)

		c = ax.scatter(xData, yData, c = C_monthly, marker=".",\
				s=6, edgecolors='none',\
				vmin=vmin, vmax=vmax,\
				alpha=1,\
				norm=matplotlib.colors.LogNorm()
				#cmap="Reds"
	                       )
	
		ax.spines['top'].set_visible(False)
		ax.yaxis.set_ticks_position('left')
		ax.xaxis.set_ticks_position('bottom')
		ax.tick_params(axis='y', labelsize=15)
		ax.tick_params(axis='x', labelsize=15)


		#cb = plt.colorbar(c, cmap='jet', extend='both')
		#cb.set_label('Emissions [kg/day]')

		plt.xlim([vMonth[x].min(), vMonth[x].max()])
		plt.ylim([vMonth[y].min(), vMonth[y].max()])

		plt.xlabel(uMonth[x], fontsize=20)
		plt.ylabel(uMonth[y], fontsize=20)

		print 'Figure frameN: ' + str(frameN) + ' is complete'

fig.tight_layout()
fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(c, cax=cbar_ax, extend="both")
cbar.set_label('Emissions [g C month$^{-1}$]', fontsize=30)
cbar.ax.tick_params(labelsize=24)
print 'Writing the figure to memory'
plt.savefig(figSaveName, dpi=600)



##################################################################
# TODO: next make plots of C vs met parameters, colorcoded by???
##################################################################



for v_key in v.keys():

	print v_key

	fig = plt.figure(figsize=(13,10))
	savename = '../Figures/' + 'C_vs_'+ v_key + '_' +\
			str(startMonth) + '_' +str(endMonth)+'_daily.png'

	ax = plt.subplot(111)
	xData=v[v_key]
	c = ax.scatter(xData, C, c=month_matrix,
		       edgecolors='none'
		       )

	# Make axis pro
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.tick_params(axis='y', labelsize=20)
	ax.tick_params(axis='x', labelsize=20)


	plt.xlabel(u[v_key], fontsize=26)
	plt.ylabel('Emissions [g C day$^{-1}$]', fontsize=26)

	cbar = fig.colorbar(c)
	cbar.set_label('Month', fontsize=30)
	cbar.ax.tick_params(labelsize=24)

	plt.savefig(savename, dpi=600)

	#############################
	# CREATE FOR MONTHLY AS WELL
	#############################
	fig = plt.figure(figsize=(10,10))
	savename = '../Figures/' + 'C_vs_'+ v_key + '_' +\
			str(startMonth) + '_' +str(endMonth)+'_monthly.png'

	ax = plt.subplot(111)
	xData=vMonth[v_key]
	c = ax.scatter(xData, C_monthly,
		       edgecolors='none'
		       )

	# Make axis pro
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	ax.tick_params(axis='y', labelsize=20)
	ax.tick_params(axis='x', labelsize=20)


	plt.xlabel(uMonth[v_key], fontsize=26)
	plt.ylabel('Emissions [g C month$^{-1}$]', fontsize=26)

	plt.savefig(savename, dpi=600)




##################################################################
# TODO: Calculate total emissions for a given parameter (e.g. T)
##################################################################


