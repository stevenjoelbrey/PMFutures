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
minLat     = 30.    # 10.
maxLat     = 50.    # 90.
minLon     = 234.   # 190.
maxLon     = 259.   # 320

# Load resources
import cesm_nc_manager as cnm
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
import matplotlib.ticker as tkr

EScenario = "RCP85" # NOTE: RCP85 2010 correponds to 2000Base
EYear     = "2100"  # NOTE: For emissions this gives 2000-2010 decade
AScenario = "2100RCP85" 

# Get Particulate species emissions
BC, BCunits, BCLongName, t, ELat, ELon = cnm.getEmissionVariableData("BC",EScenario,EYear)
OC, OCunits, OCLongName, t, ELat, ELon = cnm.getEmissionVariableData("OC",EScenario,EYear)

# Combine emissions
secondsPerDay = 24. * 60**2
E = (OC + BC) * secondsPerDay # units of E now kg/m2/day


# Subset missions
E, t, ELat, ELon  =  cnm.subsetModelEmissions(E, t, ELat, ELon,\
                                              startMonth, endMonth,\
                                              minLat, maxLat, minLon, maxLon)

# Load grid cell area once final domain of emissions analysis is set
area, landMask = cnm.loadEmissionGridAttributes(ELat, ELon)




################################################################################
# Get the met parameters we want to explore the parameter space of 
################################################################################

TFile = cnm.makeAQNCFile('T', AScenario, 'daily')
T, TUnits, TLongName, t, lat, lon = cnm.getGroundAirQaulityData("T", TFile, AScenario)

VFile = cnm.makeAQNCFile('V', AScenario, 'daily')
V, VUnits, VLongName, vt, lat, lon = cnm.getGroundAirQaulityData("V", VFile, AScenario)

UFile = cnm.makeAQNCFile('U', AScenario, 'daily')
U, UUnits, ULongName, ut, lat, lon = cnm.getGroundAirQaulityData("U", UFile, AScenario)

windSpeed = np.sqrt(V**2 * U**2)

RHFile = cnm.makeAQNCFile('RELHUM', AScenario, 'daily')
RH, RHUnits, RHLongName, RHt, lat, lon = cnm.getGroundAirQaulityData("RELHUM", RHFile, AScenario)


ZFile = cnm.makeAQNCFile('Z3_P', AScenario, 'daily')
Z_nc  = Dataset(ZFile, 'r')
Z     = Z_nc.variables['Z3']
plevel= Z_nc.variables['plevel'][:] 
Z500Index = np.where(plevel == 500)[0][0]
Z500  = Z[:, Z500Index, :, :]
Z_nc.close()

# precip in meters per second
precip = cnm.getSelf(AScenario, 'PRECT')[:]
precip = cnm.metersPerSecToMetersPerDay(precip)

################################################################################
# Make subsetting these easier by looping through a dictionary
################################################################################
v = {}
v['T'] = T
v['windSpeed'] = windSpeed
v['precip'] = precip
v['RH'] = RH
v['Z500'] = Z500

u = {}
u['T'] = 'temperature [k]'
u['windSpeed'] = 'windspeed [m s$^{-1}$]'
u['precip'] = 'precip [in day$^{-1}$]'
u['RH'] = "RH %"
u['Z500'] = "geopotential [m$^{2}$ s$^{-2}$]"


# Subset the parameters to the same limited time and space domain as the emissions
for k in v.keys():

	v[k], tSubset, latSubset, lonSubset  =  cnm.subsetModelEmissions(v[k], t,\
													  lat, lon,\
		                                              startMonth, endMonth,\
		                                              minLat, maxLat, minLon, maxLon)

# TODO: Plot all parameters in space with abline for event thresholds
# Make a version of all these plots that shows the points that make up 99% of total 
# emissions

# Range of emissions is the same for all variables in the space so set a common
# colorbar range. 

fig = plt.figure(figsize=(30,22))
nVar = len(v.keys())
#nVar = 2
frameN = 0

#vmin = 10**-9 * 1.
Emin = np.percentile(E, 95)


E_high = ma.masked_where(E < Emin, E)

for i in range(nVar):

	x = v.keys()[i]

	for j in range(nVar):
	
		y = v.keys()[j]

		frameN += 1 

		xData = ma.masked_where(E < np.percentile(E, 95), v[x])
		yData = ma.masked_where(E < np.percentile(E, 95), v[y])

		ax = fig.add_subplot(5,5,frameN, frame_on=False)

		c = ax.scatter(xData, yData, c = E_high,\
	                   marker=".",\
                       s=4,\
	                   edgecolors='none',\
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

fig.subplots_adjust(right=0.80)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cbar = fig.colorbar(c,cax=cbar_ax, cmap='jet', extend='both')
cbar.set_label('Emissions [kg m$^{-2}$ day$^{-1}$]', fontsize=30)

plt.savefig('../Figures/' + AScenario + '_parameter_space.png', dpi=500)
 









