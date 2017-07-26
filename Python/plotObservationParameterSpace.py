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
# For report use western U.S. only

startMonth = 6
endMonth   = 9

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

dataDirBase = "/barnes-scratch/sbrey/era_interim_nc_daily_merged/"

# Get emissions, use this to get dimensions
ncFile  = "/barnes-scratch/sbrey/GFED4s/GFED4.1s_METGrid_C_NA_2003_2016.nc"
nc = Dataset(ncFile, 'r')
latitude = nc.variables['latitude'][:]
longitude = nc.variables['longitude'][:]
time = nc.variables['time'][:]
C = nc.variables['C'][:]
nc.close()

# Make hourly time array into a nice datetime object 
t0 = datetime.datetime(year=1900, month=1, day=1,\
                       hour=0, minute=0, second=0)

# Make a nice month and time array for masking 
nTime = len(time)
t = []
month = []
for i in range(nTime):
	dt = timedelta(hours=int(time[i]))
	t_new  = t0 + dt 
	t.append(t_new)
	month.append(t_new.month)
t = np.array(t)
month = np.array(month)


# Subset all by months of interest
monthMask = (month >= 6) & (month <= 9)

C = C[monthMask, :,:]
t = t[monthMask]
month = month[monthMask]


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
RH   = ncGetSelf(dataDirBase, 'RH')[monthMask, :,:] 



# TODO: Plot all parameters in space with abline for event thresholds
# Make a version of all these plots that shows the points that make up 99% of total 
# emissions


v = {}
v['T'] = t2m
v['windSpeed'] = windSpeed
v['precip'] = tp
v['RH'] = RH
v['Z500'] = z

u = {}
u['T'] = 'temperature [k]'
u['windSpeed'] = 'windspeed [m s$^{-1}$]'
u['precip'] = 'precip [in day$^{-1}$]'
u['RH'] = "RH %"
u['Z500'] = "geopotential [m$^{2}$ s$^{-2}$]"


# Range of emissions is the same for all variables in the space so set a common
# colorbar range. 

fig = plt.figure(figsize=(30,22))
nVar = len(v.keys())
frameN = 0

vmin = 0
vmax = C.max()

#E_high = ma.masked_where(C < Emin, E)

for i in range(nVar):

	x = v.keys()[i]

	for j in range(nVar):
	
		y = v.keys()[j]

		frameN += 1 

		#xData = ma.masked_where(E < np.percentile(E, 95), v[x])
		#yData = ma.masked_where(E < np.percentile(E, 95), v[y])

		xData = v[x]
		yData = v[y]

		ax = fig.add_subplot(5,5,frameN, frame_on=False)

		c = ax.scatter(xData, yData, c = C, marker=".",\
				s=4, edgecolors='none',\
				vmin=vmin,vmax=vmax,\
				alpha=0.7,\
				norm=matplotlib.colors.LogNorm()
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

fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(c, cax=cbar_ax)
cbar.set_label('Emissions [g C day$^{-1}$]', fontsize=30)

plt.savefig('../Figures/GFED_era_interim_6_9_parameter_space.png', dpi=500)
 









