#!/usr/bin/env python2

# met_event_finder_m.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This module is the home of functions used to identify meteorogy events 
# relevent for air quality and climate change impact assesment in the EPA
# "Planning for an uncertain future" project (dubbed PMFutures). 

# TODO: Deal with output that has been placed on pressure coordinates? 


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



###############################################################################
# ---------------- Environment for testing functions -------------------------- 
###############################################################################
variable = 'T'
scenario = '2000Base'


ncFile = cnm.makeAQNCFile(variable, scenario)
M, MUnits, MLongName, t, lat, lon = cnm.getGroundAirQaulityData(variable, ncFile, scenario)

wind, wUnits, wName, t, lat, lon = cnm.getSurfaceWindScaler(scenario)
