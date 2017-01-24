#!/usr/bin/env python2

# createLinearModel.py

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# This script will be used to create a linear model of the relationship between
# meteorlogy parameters and wildfire emissions. 

import cesm_nc_manager as cnm
import met_event_finder_m as mef
import plotEmissionSummary as pme
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


# Make a scatterplot of T vs. emissions for a given month

# make these managers into one. update naming conventions and make direct
# comparisons easier. 

E, EUnits, ELongName, Et, Elat, Elon = pme.getEmissionVariableData("BC", \
                                                                   "RCP85",\
                                                                   "2010")
ncFile = cnm.makeAQNCFile("T", "2000Base")                                                                   
M, MUnits, MLongName, Mt, Mlat, Mlon = mef.getGroundAirQaulityData("T",\
                                                                ncFile,\
                                                                "2000Base") 

# Figure out a way to make the two firds comparable                                                                  