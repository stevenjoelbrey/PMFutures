#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# Blocking days will be found by identified using Barnes et al. 2012. This 
# searches for a persisten reversal in the 500 hPa gradient. 


import os
import numpy as np
import sys
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import numpy.ma as ma
from datetime import date
import matplotlib.ticker as tkr
from datetime import date












