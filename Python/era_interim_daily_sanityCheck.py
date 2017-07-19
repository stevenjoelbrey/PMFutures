#!/usr/bin/env python2

###############################################################################
# ------------------------- Description --------------------------------------- 
###############################################################################
# Sanity check on merged daily nc reanalysis data.
# I want to plot met fields together to make sure things are not all murphed up.

import os
import numpy as np
from mpl_toolkits.basemap import Basemap, cm
from netCDF4 import Dataset 
import matplotlib.pyplot as plt
import numpy.ma as ma
from datetime import date
from datetime import timedelta
import matplotlib.ticker as tkr
import datetime
import time as timer
import os.path



